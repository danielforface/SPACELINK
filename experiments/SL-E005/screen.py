"""Conditional state and ambient-carrier screening."""
import importlib.util, pathlib, math, json, hashlib, urllib.request
import numpy as np
spec=importlib.util.spec_from_file_location("link",pathlib.Path(__file__).resolve().parents[1]/"SL-E003"/"link.py")
e3=importlib.util.module_from_spec(spec);spec.loader.exec_module(e3)
P=e3.P
def mass_delta(S,p=P):
    T=abs(S)**2; d=np.zeros(3)
    for i in range(3):
        for j in range(3):
            if i!=j:
                v=p[i]*T[j,i]; d[i]-=v;d[j]+=v
    return d
def kl_term(x):
    if abs(x)<1e-4:return x*x*(.5+x*(-1/6+x*(1/12-x/20)))
    return (1+x)*math.log1p(x)-x if x>-1 else 1.
def mass_metrics(layers):
    d=mass_delta(e3.oracle(layers)); B=e3.old.bound(layers); tv=float(sum(abs(d))/2)
    chi=float(sum(d*d/P))/math.log(2)
    kl=sum(float(p)*kl_term(float(v/p)) for p,v in zip(P,d))/math.log(2)
    cap=B*B*sum(1/x for x in sorted(P)[:2])/math.log(2)
    return dict(tv=tv,B=B,kl_bits=kl,chi_bits=chi,column_ceiling_bits=cap)
MODELS=[("BS05-OP","bs2005opflux.dat","9435d546fc7079b1d82e6922ccbf0405c5a609dd",13,10),("B16-GS98","nudistr_b16_gs98.dat","f73c47cf6f2d77086634a5c180b50039e10805e7",12,7),("B16-AGSS09","nudistr_b16_agss09.dat","d9bd29f3374c63e8ea898733a55fb7aa566a2c96",12,7)]
REF="1783f9f976c2f445ccf686c25ba003baeee8f357"
def parse(data,sha,ncol,idx):
    if hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()!=sha:raise ValueError("blob hash mismatch")
    rows=[]
    for line in data.decode("ascii").splitlines():
        t=line.split()
        if len(t)!=ncol:continue
        try: x=list(map(float,t))
        except ValueError:continue
        rows.append((x[0],x[idx]))
    if len(rows)<100 or not all(math.isfinite(r) and math.isfinite(w) and 0<r<=1 and w>=0 for r,w in rows) or not all(a[0]<b[0] for a,b in zip(rows,rows[1:])):raise ValueError("invalid rows")
    norm=math.fsum(w for r,w in rows)
    if abs(norm-1)>.02:raise ValueError("unexpected normalization "+str(norm))
    return rows,norm
def shell_coverage(rows,norm,a,zmod=1,zrx=1.524):
    if zrx<=zmod:return 0.
    return math.fsum(w*e3.shell_fraction(a/((1-zmod/zrx)*r*e3.RS)) for r,w in rows)/norm
def get_models():
    result=[]
    for name,file,sha,ncol,idx in MODELS:
        url=f"https://raw.githubusercontent.com/michelelucente/PEANUTS/{REF}/Data/{file}"
        try:
            with urllib.request.urlopen(url,timeout=15) as res:data=res.read(2_000_001)
            if len(data)>2_000_000:raise ValueError("size limit")
            rows,norm=parse(data,sha,ncol,idx)
            result.append(dict(name=name,url=url,blob=sha,sha256=hashlib.sha256(data).hexdigest(),rows=len(rows),norm=norm,coverage=[dict(radius_m=a,f=shell_coverage(rows,norm,a),disk=e3.coverage(a,.1)) for a in (1,1000,100000,1e7)]))
        except Exception as ex:result.append(dict(name=name,url=url,blocked=str(ex)))
    return result
def laser(P,D,dt,dr,wavelength):
    return P*(math.pi*dt*dr/(4*wavelength*D))**2
def mirror(P,D,dr,theta=.00465):
    return P*(dr/(2*D*math.tan(theta)))**2
def photons(P,lam):return P*lam/(6.62607015e-34*299792458)
def main():
    rows=[]
    for L in (1,100,1000,100000):
        m=mass_metrics([(L,19.3,74/183.84)]);m["length_m"]=L;m["illustrative_tagged_ceiling_bits_day"]=m["column_ceiling_bits"]*48.3/1.524**2*e3.coverage(1000,.1);rows.append(m)
    n=6e10/2.99792458e10; V=2*7.6e-14*n/6.02214076e23
    D=225e9;mr=mirror(1361,D,5)
    print(json.dumps(dict(evidence="DERIVED/SIMULATED screening",mass=rows,models=get_models(),collective=dict(n_cm3=n,V_eV=V,ratio_vacuum=V/(7.5e-5/(2*862000)),phase_1au_constant=V*149597870700/1.973269804e-7),bragg=dict(wavelength_m=1.239841984e-6/862000,theta_rad=math.asin((1.239841984e-6/862000)/(2*.2e-9))),optics=dict(mirror_W=mr,mirror_equiv_photons_s=photons(mr,550e-9),laser=[dict(efficiency=eta,tx_W=1361*eta,rx_W=laser(1361*eta,D,.22,5,1550e-9),photons_s=photons(laser(1361*eta,D,.22,5,1550e-9),1550e-9)) for eta in (.001,.01,.03,.1)],dsoc_4w_rx_W=laser(4,D,.22,5,1550e-9))),indent=2))
if __name__=="__main__":main()
