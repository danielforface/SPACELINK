"""Conditional amplitude bounds and resource scales."""
import pathlib,importlib.util,math,json
import numpy as np
spec=importlib.util.spec_from_file_location("solar",pathlib.Path(__file__).resolve().parents[1]/"SL-E006"/"solar.py")
e6=importlib.util.module_from_spec(spec);spec.loader.exec_module(e6)
P=np.array([.6038379159606956,.37394359550989803,.02221848852940717])
HB=1.973269804e-7;E0=862000.;DIST=.524*149597870700;V=7.6e-14*19.3*74/183.84
def delta(E,L):
    S=e6.e3.oracle([(L,19.3,74/183.84)],E)
    D=(S*P)@S.conj().T
    np.fill_diagonal(D,e6.e5.mass_delta(S,P))
    return D
def amplitudes(L,Emin):
    phi=abs(V)*L/HB
    return 2*max(P)*phi,2*max(P)*phi*(2.5e-3*L/(2*HB*Emin**2))
def numerator(E,g):
    return E[-1]**2*(g[0]+g[-1]+sum(abs(np.diff(g))))+2*E[-1]
def entry_bound(E,g,L,dist,dm):
    A,M=amplitudes(L,E[0]);K=dm*dist/(2*HB)
    return min(A,(A*numerator(E,g)+E[-1]**2*M)/abs(K)) if K else A
def scales(dist,dm):
    K=dm*dist/(2*HB)
    w=2*E0**2/(math.hypot(K,E0)+K) if K else 2*E0
    dx=dist*dm/(2*E0**2)
    return dict(distance_m=dist,dm2=dm,energy_full_width_eV=w,position_one_rad_m=2*HB*E0/dm,packet_separation_m=dx,packet_sigma_half_overlap_m=dx/math.sqrt(8*math.log(2)),gaussian_sigmaE_half_contrast_eV=math.sqrt(2*math.log(2))*E0**2/K if K else None)
def run(E,g,L):
    A,M=amplitudes(L,E[0]);D=np.array([delta(e,L) for e in E]);mask=~np.eye(3,dtype=bool)
    peak=float(np.max(abs(D[:,mask])));der=float(np.max(abs(np.diff(D[:,mask],axis=0)/np.diff(E)[:,None])))
    if peak>A*(1+1e-9) or der>M*(1+1e-9):raise ValueError("amplitude check")
    return dict(material_m=L,A_bound=float(A),derivative_bound_per_eV=float(M),sample_max=peak,sampled_derivative=der,entry_bounds=[dict(pair=p,bound=float(entry_bound(E,g,L,DIST,dm))) for p,dm in (("21",7.5e-5),("31",2.5e-3),("32",2.425e-3))])
def main():
    E,g=e6.line(e6.download("be7ground_shape.csv","19d97e87e8b03678d2e4629598f7d7fd7794eae0"))
    print(json.dumps(dict(e008="DERIVED/SIMULATED",rows=[run(E,g,L) for L in (1,1000)],scales=[scales(d,m) for d in (1e6,1e9,DIST) for m in (7.5e-5,2.5e-3,2.425e-3)]),separators=(",",":")))
if __name__=="__main__":main()
