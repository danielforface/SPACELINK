"""Optimistic aligned count-only screen; not an actual solar link."""
import importlib.util, pathlib, math, json
import numpy as np
spec=importlib.util.spec_from_file_location("e002",pathlib.Path(__file__).resolve().parents[1]/"SL-E002"/"audit.py")
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
RS=6.957e8
P=np.array([.68,.30,.02])
U=np.sqrt([.693*.978,.307*.978,.022])
W=U**2
def oracle(layers,E=862000.):
    s=np.eye(3,dtype=complex)
    for length,rho,ye in layers:
        h=np.diag(np.array([0,7.5e-5,2.5e-3])/(2*E))+7.6e-14*rho*ye*np.outer(U,U)
        e,v=np.linalg.eigh(h)
        s=(v*np.exp(-1j*e*length/1.973269804e-7))@v.T@s
    return s
def compare():
    absmax=relmax=0.; n=0
    cases=[(ls,862000.) for ls in old.sample_profiles()]
    cases += [([(L,19.3,74/183.84)],E) for L in (1,100,1000,100000) for E in (384000.,862000.,1442000.)]
    for ls,E in cases:
        a=abs(np.asarray(old.propagate(ls,E)))**2;b=abs(oracle(ls,E))**2
        absmax=max(absmax,float(np.max(abs(a-b))))
        mask=(b>1e-10)&(~np.eye(3,dtype=bool))
        if mask.any():relmax=max(relmax,float(np.max(abs(a[mask]/b[mask]-1))))
        n+=1
    return dict(cases=n,max_absolute_probability_error=absmax,max_relative_offdiagonal_error=relmax)
def coverage(a,rfrac,zmod=1.,zrx=1.524):
    if a<0 or rfrac<=0 or zmod<=0 or zrx<=0: raise ValueError("geometry")
    if zrx<=zmod:return 0.
    return min(1.,(a/(rfrac*RS*(1-zmod/zrx)))**2)
def shell_fraction(x):
    if x<0:raise ValueError("negative ratio")
    if x>=1:return 1.
    return x*x/(1+math.sqrt(1-x*x))
def xs(E=.862,T0=0.,electron=True):
    m=.51099895; r=.2317; l=r+(.5 if electron else -.5)
    t=2*E*E/(m+2*E)
    if T0<0 or T0>=t:raise ValueError("threshold")
    def primitive(T):return l*l*T+r*r*(T-T*T/E+T**3/(3*E*E))-l*r*m*T*T/(2*E*E)
    return primitive(t)-primitive(T0)
def xs_quad(E=.862,T0=0.,electron=True):
    m=.51099895;r=.2317;l=r+(.5 if electron else -.5);t=2*E*E/(m+2*E)
    x=np.linspace(T0,t,10001)
    y=l*l+r*r*(1-x/E)**2-l*r*m*x/E**2
    return float(np.trapz(y,x))
def info_bound(signal,background,relative):
    if signal<=0 or background<0:raise ValueError("rates")
    return (signal*relative)**2/((signal+background)*math.log(2))
def link(length,a,rfrac=.1,rate=48.3,bg=0.,full=False):
    layers=[(length,19.3,74/183.84)]
    s=oracle(layers); dp=float(sum(P[i]*(W[j]-W[i])*abs(s[j,i])**2 for i in range(3) for j in range(3) if i!=j))
    b=old.bound(layers)*float(np.ptp(W))
    ratio=xs(electron=False)/xs(); response=(1-ratio)/(ratio+(1-ratio)*float(P@W))
    f=1. if full else coverage(a,rfrac)
    sig=rate/1.524**2
    return dict(length_m=length,radius_m=a,source_radius_fraction=rfrac,coverage=f,delta_pe=dp,relative_count_change=f*response*dp,relative_upper=f*response*b,info_upper_bits_day=info_bound(sig,bg,f*response*b),background_day=bg,reference_rate_earth_day=rate)
def main():
    print(json.dumps(dict(evidence="SIMULATED/DERIVED count-only model",numpy=np.__version__,oracle=compare(),sigma_e_over_x=xs()/xs(electron=False),rows=[link(L,a,r) for L in (1,1000) for a in (1,1000,100000,10000000) for r in (.05,.1,.2)],full=[link(L,0,full=True) for L in (1,1000)],sensitivity=[link(1000,1000,rate=r,bg=b) for r in (46.1,48.3,50.5) for b in (0,100)]),indent=2))
if __name__=="__main__":main()
