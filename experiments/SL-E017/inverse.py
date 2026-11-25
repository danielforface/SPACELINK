"""Exploratory inverse-discovery screens, not a physical link."""
import pathlib,importlib.util,math,json
import numpy as np
s=importlib.util.spec_from_file_location("e6",pathlib.Path(__file__).resolve().parents[1]/"SL-E006"/"solar.py")
e6=importlib.util.module_from_spec(s);s.loader.exec_module(e6)
HB=1.973269804e-7;E0=862000.;DM=7.5e-5;C=.386;S=2*math.sqrt(.693*.307)
P=np.array([.603837916,.373943596]);P=P/P.sum();RHO=np.diag(P);V=7.6e-14*.978*19.3*74/183.84
D0=DM/(2*E0);VS=D0/C;LS=math.pi*HB/(VS*S);Q=math.pi*HB/(2*D0)
X=np.array([[0,1],[1,0]],complex);Z=np.diag([1.,-1.])
def slab(E,v,L):
    h=(v*C-DM/(2*E))*Z/2+v*S*X/2
    om=math.hypot(v*C-DM/(2*E),v*S)/2
    if om==0:return np.eye(2)
    t=om*L/HB
    return math.cos(t)*np.eye(2)-1j*math.sin(t)*h/om
def imprint(E,v=V):
    u=slab(E,v,1)
    return u@RHO@u.conj().T-RHO
def maxswap(E,v):
    return (v*S)**2/((v*C-DM/(2*E))**2+(v*S)**2)
def projected(d):return float(C*d[0,0].real+S*d[0,1].real)
def echo(E,M,offset,v=V):
    d=imprint(E,v);t=abs(M[0,1])**2
    diag=(1-2*t)*d[0,0].real
    off=M[0,1]*M[1,0].conjugate()*d[1,0]*np.exp(1j*DM*offset/(2*E*HB))
    return C*diag+S*off.real,S*off.real
def direct_phase(E,M,offset):
    d=imprint(E)
    total=np.zeros((2,2),complex)
    for theta in 2*np.pi*np.arange(16)/16:
        A=np.diag(np.exp(1j*np.array([1,-1])*theta/2))
        B=np.diag(np.exp(1j*np.array([1,-1])*(theta+DM*offset/(2*E*HB))/2))
        U=B@M@A
        total+=U@d@U.conj().T/16
    return projected(total)
def spectrum():
    E,g=e6.line(e6.download("be7ground_shape.csv","19d97e87e8b03678d2e4629598f7d7fd7794eae0"))
    E1=np.linspace(E[0],E[-1],math.ceil(E[-1]-E[0])+1)
    E1=np.unique(np.r_[E1,E]);g1=np.interp(E1,E,g);g1/=np.trapz(g1,E1)
    return E1,g1
def run():
    E,g=spectrum()
    dep=float(np.trapz(g*np.array([C*imprint(e)[0,0].real for e in E]),E))
    rows=[]
    for label in ("ideal","slab"):
        for off in (0,Q,1e5,1e6):
            vals=np.array([echo(e,X if label=="ideal" else slab(e,VS,LS),off) for e in E])
            total,coh=np.trapz(vals*g[:,None],E,axis=0)
            rows.append(dict(mode=label,offset_m=off,delta_Pee=float(total),coherent=float(coh),ratio_to_dephased=abs(float(coh)/dep)))
    z=650*149597870700.;R=2.5;be=math.sqrt(4*6.67430e-11*1.98847e30*z/299792458.**2)
    root=math.hypot(R,2*be);bm=(root-R)/2;area=math.pi*R*root;gain=area/(math.pi*R*R)
    lens=[]
    for a in (1,1000,1e5):
        patch=min(math.pi*a*a,area*math.asin(a/bm)/math.pi);f=patch/area
        lens.append(dict(a_m=a,patch_bound_m2=patch,f=f,tag_gain=patch/(math.pi*R*R),aggregate_coefficient_gain=gain*f*f))
    n=6e10/2.99792458e10;pot=7.6e-14*n/6.02214076e23;phase=pot/6.582119569e-16;signal=phase*4.8968563e-9*2e-7
    print(json.dumps(dict(e017="DERIVED/SIMULATED",spectral_nodes=len(E),dephased=dep,echo=rows,swap=dict(Veff_eV=VS,rhoYe_g_cm3=VS/(.978*7.6e-14),length_m=LS,earth_max=maxswap(E0,.978*7.6e-14*13*.47),tungsten_max=maxswap(E0,V),tuned_probability=abs(slab(E0,VS,LS)[0,1])**2),lens=dict(b_m=be,total_gain=gain,source_angle_rad=R/z,rows=lens),reciprocal=dict(n_cm3=n,V_eV=pot,phase_1s=phase,modulated_phase_1s=signal,ideal_independent_probes_1sigma=1/signal**2),rare_event_kl_bound_bits=(1e-3*(2e-7)**2)/math.log(2)),separators=(",",":")))
if __name__=="__main__":run()
