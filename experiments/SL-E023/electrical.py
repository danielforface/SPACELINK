"""E023: electrical material records and final-state occupation screen."""
import math,json,hashlib,urllib.request
import numpy as np
NA=6.02214076e23; LN2=math.log(2); KB=8.617333262e-5; HB=1.973269804e-5
GA=69.723; AS=74.921595; RHO=5.32; FRAC=GA/(GA+AS)
RATE_GA=8.382666805346743e-6 # per kg/day, Mars
LIFE=11.43/LN2; EPS=4.749291240080495e-18
R=RATE_GA*FRAC*RHO/1000 # per cm3/day
M=71*931.49410242e6; D=.524*149597870700.; F=4.8968563e-9
Q=862000.; ME=510998.95; ALPHA=1/137.035999084
def material():
    vol=100000*1000/RHO
    return dict(ga_mass_fraction=FRAC,captures_cm3_day=R,live_density_cm3=R*LIFE,
        spacing_cm=(R*LIFE)**(-1/3),volume_100t_cm3=vol,captures_100t_day=R*vol,
        live_100t=R*LIFE*vol,background_ratios=[dict(concentration_cm3=b,ratio=R*LIFE/b) for b in (1e10,1e13,1e15)],
        pairs=[dict(E_MeV=e,pairs=(e-.232443)*1e6/4.2) for e in (.42,.862,1.442)])
def cells(side_cm):
    v=side_cm**3; r=R*v; n=100000*1000/RHO/v
    return dict(side_cm=side_cm,cells=n,yearly_captures_per_cell=r*365,
        mean_false_hazard_per_cell_day=.1*r,yearly_false_probability=-math.expm1(-.1*r*365),
        aggregate_true_year=n*r*365)
def coulomb(r_nm,t):
    e=1.439964/(12.9*r_nm)
    return dict(r_nm=r_nm,T_K=t,energy_eV=e,over_kBT=e/(KB*t))
def latch(background_ratio):
    s=material()['captures_100t_day'];b=s*background_ratio
    return dict(background_ratio=background_ratio,baseline_day=s+b,chi_bits_day=(s*EPS)**2/((s+b)*LN2))
def spectrum():
    url='https://raw.githubusercontent.com/michelelucente/PEANUTS/1783f9f976c2f445ccf686c25ba003baeee8f357/Data/be7ground_shape.csv'
    with urllib.request.urlopen(url,timeout=20) as r:b=r.read(2000001)
    if hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()!='19d97e87e8b03678d2e4629598f7d7fd7794eae0':raise ValueError('hash')
    E,g=np.loadtxt(b.decode().splitlines(),delimiter=',').T
    E*=1e6;g/=np.trapz(g,E)
    return E,g,hashlib.sha256(b).hexdigest()
def weight(K):
    ee=ME+K;p=np.sqrt(K*(K+2*ME))
    beta=p/ee; a=2*math.pi*ALPHA*4
    z=np.divide(a,beta,out=np.full_like(beta,np.inf),where=beta>0)
    return a*(ee/ME)**2/(-np.expm1(-z))
def interval(E,g,lo,hi,step,mode='trap',mu=None,T=0):
    lo=max(lo,Q,E[0]);hi=min(hi,E[-1])
    if hi<=lo:return 0.
    x=np.unique(np.r_[np.linspace(lo,hi,math.ceil((hi-lo)/step)+1),E[(E>lo)&(E<hi)]])
    def y(x):
        K=x-Q; v=np.interp(x,E,g)*weight(K)
        if mu is not None:
            z=np.clip((K-mu)/(KB*T),-700,700)
            v*=1/(1+np.exp(z))
        return v
    if mode=='mid':return float(np.sum(y((x[1:]+x[:-1])/2)*np.diff(x)))
    return float(np.trapz(y(x),x))
def blocking(E,g,step):
    den=interval(E,g,Q,E[-1],step)
    base=interval(E,g,Q,E[-1],step,mu=0,T=300)
    rows=[]
    for mu in (1,5,10,100,1000):
        a=interval(E,g,Q,Q+mu,step)/den
        oracle=interval(E,g,Q,Q+mu,.07,mode='mid')/interval(E,g,Q,E[-1],.07,mode='mid')
        b=(interval(E,g,Q,E[-1],step,mu=mu,T=300)-base)/(den-base)
        rows.append(dict(mu_eV=mu,blocked_fraction_T0=a,midpoint_reference=oracle,
          blocked_change_T300=b,free_e_density_cm3=(2*ME*mu)**1.5/(3*math.pi**2*HB**3)))
    return dict(weight_integral=den,base_blocked_T300=base/den,rows=rows)
def momentum(sigma):
    req=862000*1000/D
    return dict(sigma_eV_c=sigma,required_eV_c=req,angular_rms_rad=sigma/862000,
        patch_rad=1000/D,within_patch=-math.expm1(-req**2/(2*sigma*sigma)))
def run():
    E,g,sha=spectrum(); coarse=blocking(E,g,1);fine=blocking(E,g,.1)
    a=fine['rows'][2]['blocked_fraction_T0'];tau=.534/7.016*NA*100*1e-42
    print(json.dumps(dict(e023='DERIVED/SIMULATED',sha256=sha,material=material(),
        cells=[cells(s) for s in (1e-6,1e-5,1e-4,.01,.1,1)],
        coulomb=[coulomb(r,t) for r in (10,100,1000) for t in (4,77,300)],
        latch=[latch(b) for b in (0,1,100)],blocking_coarse=coarse,blocking=fine,
        toy_sender=dict(tau=tau,delta_absorption=tau*a,diluted_delta=tau*a*F),
        direction_thermal=[dict(T_K=t,**momentum(math.sqrt(M*KB*t))) for t in (1,4,77,300)],
        direction_zero_point=[dict(mode_meV=e,**momentum(math.sqrt(M*e*1e-3/2))) for e in (1,10,30)]),separators=(',',':')))
if __name__=='__main__':run()
