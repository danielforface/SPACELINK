"""Adiabatic solar screen and continuous-spectrum phase-average bound."""
import pathlib,importlib.util,math,hashlib,urllib.request,json
import numpy as np
spec=importlib.util.spec_from_file_location("e5",pathlib.Path(__file__).resolve().parents[1]/"SL-E005"/"screen.py")
e5=importlib.util.module_from_spec(spec);spec.loader.exec_module(e5)
e3=e5.e3
W=e3.W; U=e3.U
def download(file,sha):
    url=f"https://raw.githubusercontent.com/michelelucente/PEANUTS/{e5.REF}/Data/{file}"
    with urllib.request.urlopen(url,timeout=15) as r: data=r.read(2_000_001)
    if len(data)>2_000_000 or hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()!=sha:raise ValueError("source integrity")
    return data
def mixture(logne,E=862000.):
    V=7.6e-14*10**logne
    h=np.diag(np.array([0,7.5e-5,2.5e-3])/(2*E))+V*np.outer(U,U)
    _,vec=np.linalg.eigh(h)
    return abs(vec.T@U)**2
def reduced(logne,E=862000.):
    c=.693-.307;s=2*math.sqrt(.693*.307)
    a=2*E*7.6e-14*10**logne*.978/7.5e-5
    cm=(c-a)/math.hypot(c-a,s)
    return np.array([.978*(1+cm)/2,.978*(1-cm)/2,.022])
def solar(file,sha,ncol,idx):
    data=download(file,sha);rows,norm=e5.parse(data,sha,ncol,idx)
    logs=[]
    for line in data.decode("ascii").splitlines():
        t=line.split()
        if len(t)!=ncol:continue
        try:x=list(map(float,t))
        except ValueError:continue
        logs.append(x[2])
    if len(logs)!=len(rows):raise ValueError("density mismatch")
    weights=np.array([w for r,w in rows])/norm
    f=np.array([e3.shell_fraction(1000/((1-1/1.524)*r*e3.RS)) for r,w in rows])
    m=np.array([mixture(d) for d in logs]);red=np.array([reduced(d) for d in logs])
    p=weights@m;ap=(weights*f)@m/(weights@f);p2=weights@red
    if abs(p.sum()-1)>1e-10 or abs(ap.sum()-1)>1e-10 or abs((p-p2)@W)>.003:raise ValueError("mixture validation")
    out=[]
    for L in (1,1000):
        S=e3.oracle([(L,19.3,74/183.84)]);d=e5.mass_delta(S,ap)
        out.append(dict(length_m=L,delta_pe=float(d@W),old_delta_pe=float(e5.mass_delta(S,e5.P)@W),**trace_metrics(S,ap)))
    return dict(file=file,sha256=hashlib.sha256(data).hexdigest(),p=p.tolist(),conditional_p=ap.tolist(),Pee=float(p@W),conditional_Pee=float(ap@W),reduced_Pee=float(p2@W),coverage=float(weights@f),modulator=out)
def trace_metrics(S,p):
    # Use commutator difference via subtracting identity on diagonal; no
    # cancellation issue for off-diagonals, stable transfer for diagonal.
    D=(S*p)@S.conj().T
    d=e5.mass_delta(S,p)
    np.fill_diagonal(D,d)
    return dict(coherent_trace_distance=float(sum(abs(np.linalg.eigvalsh(D)))/2),dephased_tv=float(sum(abs(d))/2))
def line(data):
    x=[]
    for l in data.decode("ascii").splitlines():
        if not l.strip() or l.startswith("#"):continue
        E,g=map(float,l.split(","));x.append((E*1e6,g))
    E,g=np.array(x).T
    if len(E)<10 or np.any(np.diff(E)<=0) or np.any(g<0) or not np.all(np.isfinite(g)):raise ValueError("line schema")
    n=float(np.trapz(g,E))
    return E,g/n
def envelope(E,g,K):
    if K==0:return 1.
    numerator=E[-1]**2*(g[0]+g[-1]+sum(abs(np.diff(g))))+2*E[-1]
    return min(1.,float(numerator/abs(K)))
def phase_data():
    data=download("be7ground_shape.csv","19d97e87e8b03678d2e4629598f7d7fd7794eae0")
    E,g=line(data);out=[]
    for L in (1000.,1e6,1e9,.524*149597870700):
        out.append(dict(L_m=L,**{name:envelope(E,g,dm*L/(2*1.973269804e-7)) for name,dm in (("bound21",7.5e-5),("bound31",2.5e-3))}))
    return dict(sha256=hashlib.sha256(data).hexdigest(),samples=len(E),energy_min_eV=float(E[0]),energy_max_eV=float(E[-1]),bounds=out)
def main():
    models=[]
    for name,file,sha,ncol,idx in e5.MODELS:
        try:models.append(dict(name=name,**solar(file,sha,ncol,idx)))
        except Exception as e:models.append(dict(name=name,error=str(e)))
    try:phase=phase_data()
    except Exception as e:phase=dict(error=str(e))
    print(json.dumps(dict(evidence="SIMULATED/DERIVED",models=models,phase=phase),indent=2))
    if any("error" in x for x in models) or "error" in phase:raise SystemExit("source or scientific validation failed")
if __name__=="__main__":main()
