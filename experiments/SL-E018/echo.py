"""E018: conditional three-state echo, geometry, and delay-only screens."""
import pathlib,importlib.util,math,json
import numpy as np
sp=importlib.util.spec_from_file_location("e17",pathlib.Path(__file__).resolve().parents[1]/"SL-E017"/"inverse.py")
e17=importlib.util.module_from_spec(sp);sp.loader.exec_module(e17)
e6=e17.e6;e5=e6.e5;e3=e6.e3
HB=e17.HB;E0=e17.E0;DM=np.array([0.,7.5e-5,2.5e-3]);U=e3.U;W=e3.W
P=np.array([.603837916,.373943596,.022218488]);P/=P.sum()
RHO=19.3;YE=74/183.84;V=7.6e-14*RHO*YE
AU=149597870700.;C_LIGHT=299792458.;DIST=.524*AU;Q=e17.Q
D0=DM[1]/(2*E0);VA=.5*.978*V
PER=2*math.pi*HB/math.hypot(D0-VA*e17.C,VA*e17.S)
X=np.array([[0,1,0],[1,0,0],[0,0,1]],complex)
def layer(E,v,L):
    h=np.diag(DM/(2*E))+v*np.outer(U,U)
    ev,vec=np.linalg.eigh(h)
    return (vec*np.exp(-1j*ev*L/HB))@vec.T
def cell(E):
    return layer(E,0,PER/2)@layer(E,V,PER/2)
def choose():
    m=cell(E0);a=np.eye(3,dtype=complex);best=(0.,0)
    for n in range(1,601):
        a=m@a
        best=max(best,(abs(a[1,0])**2,n))
    return best[1]
def imprint(E):
    a=layer(E,V,1)
    d=(a*P)@a.conj().T
    np.fill_diagonal(d,e5.mass_delta(a,P))
    return d
def echo(E,M,d,offset=Q):
    diag=(abs(M)**2)@d.diagonal().real
    delta=DM[:,None]-DM[None,:]
    off=M*M.T.conj()*d.T*np.exp(-1j*delta*offset/(2*E*HB))
    np.fill_diagonal(off,0)
    co=float(np.real(U@off@U))
    return float(W@diag)+co,co
def direct(E,M,d,offset=Q):
    idx=np.array([0,3,100])
    extra=np.exp(-1j*DM*offset/(2*E*HB))
    a=0.
    for t in 2*np.pi*np.arange(401)/401:
        v=np.exp(-1j*idx*t)
        T=(v*extra)[:,None]*M*v[None,:]
        a+=float(np.real(U@T@d@T.conj().T@U))/401
    return a
def modes(E,n):
    return {"none":np.eye(3),"ideal":X,
            "dense":layer(E,e17.VS/.978,e17.LS),
            "periodic":np.linalg.matrix_power(cell(E),n),
            "uniform":layer(E,0,PER*n/2)@layer(E,V,PER*n/2)}
def twostate(E,n):
    m=e17.slab(E,0,PER/2)@e17.slab(E,.978*V,PER/2)
    return float(abs(np.linalg.matrix_power(m,n)[1,0])**2)
def integrate(E,g,n):
    out={}
    for e in E:
        d=imprint(e)
        for label,m in modes(e,n).items():
            total,co=echo(e,m,d)
            out.setdefault(label,[]).append([abs(m[1,0])**2,total,co,float(np.max(abs(m.conj().T@m-np.eye(3))))])
    norm=np.trapz(g,E)
    return {k:dict(zip(["swap21","delta_pe","coherent","mean_unitarity_error"],np.trapz(np.array(v)*g[:,None],E,axis=0)/norm)) for k,v in out.items()}
def coverage(rows,norm,D,a=1000):
    return math.fsum(w*e3.shell_fraction(a*1.524*AU/(D*r*e3.RS)) for r,w in rows)/norm
def delay(D,dm,E=E0):
    return D*dm/(2*E*E*C_LIGHT)
def tv(lam,dt):
    return -math.expm1(-2*lam*dt)
def budget(contrast,f):
    r=e3.xs(electron=False)/e3.xs()
    pe=.5264113600;pc=float(P@W)
    # Model baseline response: whole-flux average; tag subset gets its own response.
    den=r+(1-r)*pe;dc=r+(1-r)*pc
    rate=48.3/1.524**2
    rel=f*(1-r)*contrast/den
    tag_rate=rate*f*dc/den
    tag_rel=(1-r)*contrast/dc
    return dict(relative=rel,aggregate_bits_day=e3.info_bound(rate,0,rel),tagged_bits_day=e3.info_bound(tag_rate,0,tag_rel),tag_rate_day=tag_rate)
def run():
    n=choose();E,g=e6.line(e6.download("be7ground_shape.csv","19d97e87e8b03678d2e4629598f7d7fd7794eae0"))
    fine=np.sort(np.r_[E,(E[:-1]+E[1:])/2])
    coarse=integrate(E,g,n);refined=integrate(fine,np.interp(fine,E,g),n)
    # Independent half-energy-step refinement detects rapid residual phase dependence.
    finest=np.linspace(E[0],E[-1],math.ceil((E[-1]-E[0])/5)+1)
    finest=np.unique(np.r_[finest,E])
    final=integrate(finest,np.interp(finest,E,g),n)
    _,file,sha,ncol,idx=e5.MODELS[0]
    rows,norm=e5.parse(e6.download(file,sha),sha,ncol,idx)
    geo=[dict(distance_m=d,f=coverage(rows,norm,d)) for d in (1.514*AU,DIST,1e9,1e6,1e3)]
    f=geo[1]["f"]
    result=dict(e018="DERIVED/SIMULATED",period_m=PER,cycles=n,device_length_m=n*PER,matter_length_m=n*PER/2,lower_bound_m=math.pi*HB/V,
        material_mass_kg={str(a):RHO*1000*math.pi*a*a*n*PER/2 for a in (1,1000)},
        center={k:dict(swap21=float(abs(m[1,0])**2),delta_pe=echo(E0,m,imprint(E0))[0],oracle_error=abs(echo(E0,m,imprint(E0))[0]-direct(E0,m,imprint(E0)))) for k,m in modes(E0,n).items()},
        two_state_center=twostate(E0,n),coarse=coarse,refined=refined,final=final,nodes=[len(E),len(fine),len(finest)],geometry=geo,
        link={k:budget(v["delta_pe"],f) for k,v in final.items()},
        time=[dict(dm2=float(dm),delay_s=delay(DIST,dm),incident_tv=tv(6e14,delay(DIST,dm)),all_detected_tv=tv(48.3/1.524**2/86400,delay(DIST,dm)),covered_detected_tv=tv(48.3/1.524**2/86400*f,delay(DIST,dm)),delay_over_1ns=delay(DIST,dm)/1e-9,delay_over_1ps=delay(DIST,dm)/1e-12) for dm in DM[1:]])
    print(json.dumps(result,separators=(",",":")))
if __name__=="__main__":run()
