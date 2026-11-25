"""E022: recoil scars and harvested carrier; idealized screens only."""
import math,json
import numpy as np
Q=.232443; ME=.51099895; M=71*931.49410242
RATE=1.946938073809502; HALF=11.43; LAM=math.log(2)/HALF
D=.524*149597870700.
def momentum(k):return math.sqrt(max(0.,k*(k+2*ME)))
def recoil(q2):return q2/(math.sqrt(M*M+q2)+M)
def exact(E,u):
    lo,hi=0.,E-Q
    for i in range(70):
        k=(lo+hi)/2;p=momentum(k);t=recoil(E*E+p*p-2*E*p*u)
        if k+t>E-Q:hi=k
        else:lo=k
    return (E-Q-(lo+hi)/2)*1e6
def range_recoil(E):
    p=momentum(E-Q)
    return ((E-p)**2/(2*M)*1e6,(E+p)**2/(2*M)*1e6)
def fraction(E,threshold,a):
    k=E-Q;p=momentum(k);beta=p/(k+ME)
    u=max(-1.,min(1.,(E*E+p*p-2*M*threshold/1e6)/(2*E*p)))
    return (u+1)/2+a*beta*(u*u-1)/4
def tag(accept,tau,b):
    w=-tau*math.log1p(-accept)
    return dict(accept=accept,window_s=w,background_s=b,accidental=-math.expm1(-b*w))
def geom(area):
    x=area/(math.pi*D*D)
    return x/(2*math.sqrt(1+x)*(1+math.sqrt(1+x)))
def h(p):
    if p<=0 or p>=1:return 0.
    return (-p*math.log(p)-(1-p)*math.log1p(-p))/math.log(2)
def onoff(mu):
    p=-math.expm1(-mu)
    if p<1e-8:return p/2+p*p/(8*math.log(2))
    return h(p/2)-h(p)/2
def harvest(T,area):
    n=RATE*T;s=RATE*(-math.expm1(-LAM*T))/LAM;f=geom(area)
    return dict(days=T,area_m2=area,stored_ideal=n,stored_normal=s,geometry=f,arrivals_ideal=n*f,arrivals_normal=s*f,MI_ideal=onoff(n*f))
def gate(T):
    n0=RATE/LAM;excess=RATE*T
    return dict(off_days=T,initial_inventory=n0,inventory_after_off=n0+excess,missed_emissions=excess,excess_after_resume=excess,initial_rate_after_resume=RATE+LAM*excess,release90_days=math.log(10)/LAM,net_integrated_difference=0.,one_percent_initial=1.01*RATE,one_percent_final=RATE)
def run():
    rows=[]
    for E in (.3,.42,.862,1.442,5,10):
        lo,hi=range_recoil(E);ex=[exact(E,1),exact(E,-1)]
        rows.append(dict(E_MeV=E,leading_eV=[lo,hi],exact_eV=ex,endpoint_error_eV=max(abs(lo-ex[0]),abs(hi-ex[1])),thresholds=[dict(eV=t,isotropic=fraction(E,t,0),GT_toy=fraction(E,t,-1/3)) for t in (9,13,25,39)]))
    print(json.dumps(dict(e022='DERIVED/SIMULATED',recoil=rows,EC_recoil_K_eV=recoil(.2221**2)*1e6,EC_recoil_max_eV=recoil(Q*Q)*1e6,K_gap_keV=(Q-.2221)*1000,tag=[dict(kind=k,**tag(a,t,b)) for k,t in [('fast',114e-9),('slow',HALF*86400/math.log(2))] for a in (.9,.99) for b in (1/86400,1)],branches=[dict(f=f,expected_day=.8382666805346743*f) for f in (.01,.1,1)],harvest=[harvest(t,a) for t in (28,365) for a in (1,100)],gate=gate(28)),separators=(',',':')))
if __name__=='__main__':run()
