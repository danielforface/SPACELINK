"""Pinned spectral reference and gated click receiver models."""
import pathlib,importlib.util,math,json,hashlib,urllib.request,itertools
import numpy as np
spec=importlib.util.spec_from_file_location("o7",pathlib.Path(__file__).resolve().parents[1]/"SL-E007"/"optical.py")
o7=importlib.util.module_from_spec(spec);spec.loader.exec_module(o7)
URL="https://raw.githubusercontent.com/pytroll/pyspectral/b10cc644148f78443aad75e9e12646b775e921c3/pyspectral/data/e490_00a.dat"
BLOB="8bd99d4ecb37b72125cc41bc3db907f782408b64"
def parse(data):
    if hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()!=BLOB:raise ValueError("spectrum hash")
    x,y=np.array([list(map(float,l.split())) for l in data.decode("ascii").splitlines() if l.strip() and not l.startswith("#")]).T
    if np.any(np.diff(x)<=0) or np.any(y<0) or not np.all(np.isfinite(x+y)):raise ValueError("spectrum")
    total=float(np.trapz(y,x))
    if abs(total/1366.1-1)>.01:raise ValueError("integral")
    return x*1000,y/1000*1361/total,total
def retrieve():
    with urllib.request.urlopen(URL,timeout=15) as r:data=r.read(2_000_001)
    if len(data)>2_000_000:raise ValueError("size")
    x,y,total=parse(data)
    return x,y,total,hashlib.sha256(data).hexdigest()
def band(x,y,lo,hi):
    if not x[0]<=lo<hi<=x[-1]:raise ValueError("band")
    z=np.r_[lo,x[(x>lo)&(x<hi)],hi]
    return float(np.trapz(np.interp(z,x,y),z))
def budget(x,y,mode,width,bgf,center=1064,dark=100):
    irradiance=band(x,y,center-width/2,center+width/2)
    s,b=o7.budget(mode,width,1,bgf,1)
    factor=irradiance/(o7.spectral()*width*1e-9)
    return (s*factor if mode=="mirror" else s),(b-100)*factor+dark
def probs(s,b,t):return -math.expm1(-b*t),-math.expm1(-(b+s)*t)
def binom(n,k,q):
    return math.comb(n,k)*q**k*(1-q)**(n-k) if 0<=k<=n else 0.
def ppm(M,q0,q1):
    if q0==q1:return 0.
    if q0==0:return q1*math.log2(M)
    if q0==1:return 0. # q1>=q0 in this channel
    I=0.
    for k in range(1,M):
        mix=q1*binom(M-1,k-1,q0)+(1-q1)*binom(M-1,k,q0)
        den=k*q1*(1-q0)+(M-k)*(1-q1)*q0
        if not den or not mix:continue
        p=k/M;d=k*(M-k)*(q1-q0)/(M*den)
        I+=mix*(p*o7.phi(d/p)+(1-p)*o7.phi(-d/(1-p)))/math.log(2)
    return I
def brute(M,q0,q1):
    I=0.
    for y in itertools.product((0,1),repeat=M):
        ps=[]
        for j in range(M):
            ps.append(math.prod((q1 if j==k else q0) if v else (1-(q1 if j==k else q0)) for k,v in enumerate(y)))
        p=sum(ps)/M
        if p:I+=sum(v/M*math.log2(v/p) for v in ps if v)
    return I
def rate(s,b,t,guard,M):
    q0,q1=probs(s,b,t)
    return dict(gate_s=t,mean_fraction=t/(M*(t+guard)),ook=o7.mutual(1/M,q0,q1)/(t+guard),ppm=ppm(M,q0,q1)/(M*(t+guard)))
def grid(x,y):
    for mode in ("mirror","solar_tem00","pv_laser_assumed"):
        for width in (1,10):
            for bg in (1e-10,1e-7,1e-4):
                s,b=budget(x,y,mode,width,bg)
                for guard in (20e-9,1e-6):
                    for M in (2,4,8,16):
                        rr=[rate(s,b,t,guard,M) for t in (1e-8,1e-7,1e-6,1e-5,.001,.01,.1)]
                        yield dict(mode=mode,band_nm=width,bg=bg,guard_s=guard,M=M,signal_s=s,background_s=b,ook_best=max(rr,key=lambda z:z["ook"]),ppm_best=max(rr,key=lambda z:z["ppm"]),slow=max(rr[-3:],key=lambda z:z["ppm"]))
def multi(x,y,bgf,dark):
    t=.001;guard=1e-6;M=2
    s,b=budget(x,y,"mirror",10,bgf,dark=dark)
    rs=[budget(x,y,"mirror",2.5,bgf,c,dark) for c in (1060.25,1062.75,1065.25,1067.75)]
    return dict(bg=bgf,dark_per_channel=dark,whole_signal=s,split_signal=sum(a for a,b in rs),whole_rate=rate(s,b,t,guard,M)["ook"],split_rate=sum(rate(a,b,t,guard,M)["ook"] for a,b in rs))
def main():
    x,y,total,sha=retrieve()
    print(json.dumps(dict(e009="SIMULATED reference spectral model",samples=len(x),total_raw=total,sha256=sha,band_1nm=band(x,y,1063.5,1064.5),band_10nm=band(x,y,1059,1069),planck_1nm=o7.spectral()*1e-9,multiplex=[multi(x,y,g,d) for g in (1e-10,1e-7,1e-4) for d in (0,100)])))
    for row in grid(x,y):print(json.dumps(dict(e009_row=row),separators=(",",":")))
if __name__=="__main__":main()
