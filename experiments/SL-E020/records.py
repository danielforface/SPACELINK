"""E020: persistent effect-carrier record, marked counts, threshold, resonance."""
import math,json,hashlib,urllib.request
import numpy as np
NA=6.02214076e23; DAY=86400.; LN2=math.log(2); HALF=11.43
F=4.8968563e-9; DP=5.062866248111545e-10; PE=.5220174386
def rate(kg,au=1):
    return kg*1000/69.723*.39892*NA*65.4e-36*DAY/au**2
def effective(T,half=HALF):
    if math.isinf(half):return T
    lam=LN2/half
    return -math.expm1(-lam*T)/lam
def klterm(x):
    if abs(x)<1e-4:return x*x*(.5+x*(-1/6+x*(1/12-x/20)))
    if x==-1:return 1.
    return (1+x)*math.log1p(x)-x
def kl(mu,eps):return mu*klterm(eps)/LN2
def chi(mu,eps):return mu*eps**2/LN2
def h(p):
    if p<=0 or p>=1:return 0.
    return (-p*math.log(p)-(1-p)*math.log1p(-p))/LN2
def exclusive(mu):
    p=-math.expm1(-mu)
    # I(X; event present) with X uniform and P(event|0)=0
    return h(p/2)-h(p)/2
def binary_poisson_mi(a,b):
    if a==b:return 0.
    if min(a,b)<0 or max(a,b)>1000:raise ValueError('mean range')
    nmax=math.ceil(max(a,b)+16*math.sqrt(max(a,b)+1)+50)
    total=0.; tails=[0.,0.]
    for n in range(nmax+1):
        probs=[math.exp(-v+n*math.log(v)-math.lgamma(n+1)) if v else float(n==0) for v in (a,b)]
        p,q=probs;s=p+q
        if s:
            # use x=(q-p)/s, series avoids cancellation for small contrast
            x=(q-p)/s
            if abs(x)<1e-4: val=x*x*(.5+x*x*(1/12+x*x/30))
            elif abs(x)==1:val=LN2
            else:val=((1+x)*math.log1p(x)+(1-x)*math.log1p(-x))/2
            total+=s*val/(2*LN2)
        tails[0]+=p;tails[1]+=q
    if min(tails)<1-1e-10:raise ValueError('tail')
    return total
def marked(d,r=100.):
    return dict(d=d,pooled_kl=0.,marked_kl=kl(r/2,2*d)+kl(r/2,-2*d),exclusive_presence_mi=exclusive(r*d),exclusive_event_probability=-math.expm1(-r*d))
def spectrum():
    url='https://raw.githubusercontent.com/michelelucente/PEANUTS/1783f9f976c2f445ccf686c25ba003baeee8f357/Data/be7ground_shape.csv'
    with urllib.request.urlopen(url,timeout=20) as r:b=r.read(2000001)
    if hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()!='19d97e87e8b03678d2e4629598f7d7fd7794eae0':raise ValueError('hash')
    E,g=np.loadtxt(b.decode().splitlines(),delimiter=',').T
    E*=1e6;g/=np.trapz(g,E)
    return E,g,hashlib.sha256(b).hexdigest()
def tail(E,g,q):
    if q<=E[0]:return float(np.trapz(g,E))
    if q>=E[-1]:return 0.
    x=np.r_[q,E[E>q]]
    return float(np.trapz(np.interp(x,E,g),x))
def threshold(E,g):
    q=862000.;base=tail(E,g,q);mx=max(g)
    tau=.534/7.016*NA*100*1e-42
    out=[]
    for d in (-100,-10,-1,1,10,100):
        v=tail(E,g,q+d); diff=v-base
        x=np.linspace(min(q,q+d),max(q,q+d),1001)
        quad=-math.copysign(float(np.trapz(np.interp(x,E,g),x)),d)
        out.append(dict(shift_eV=d,weight=v,delta=diff,relative=diff/base,absolute_bound=float(mx*abs(d)),oracle_error=abs(diff-quad),absorption_difference_toy=tau*abs(diff)))
    return dict(base=base,tau=tau,rows=out)
def resonance(gmax):
    width=6.582119569e-16*LN2/(53.22*DAY)
    peak=4*math.pi*(1.973269804e-5/861815)**2
    area=peak*math.pi*width/2
    avg=area*gmax
    # integral in u=atan(2(E-Eres)/Gamma), exact Lorentzian measure
    u=np.linspace(-math.pi/2,math.pi/2,1001)
    integral=float(np.trapz(np.ones_like(u)*peak*width/2,u))
    return dict(width_eV=width,peak_cm2=peak,area_cm2_eV=area,average_upper_cm2=avg,peak_to_avg=peak/avg,captures_per_tonne_day_upper=1e6/7.016*NA*4.9e9*avg*DAY,integral_rel_error=abs(integral/area-1))
def run():
    E,g,sha=spectrum(); r=rate(100000,1.524);eps=DP/PE*F
    out=dict(e020='DERIVED/SIMULATED',spectrum_sha256=sha,rate=[dict(mass_kg=m,earth_day=rate(m),mars_day=rate(m,1.524)) for m in (1,1000,100000)],
      memories=[dict(days=T,effective_days=effective(T),retained_mars_100t=r*effective(T),stable_mars_100t=r*T,fraction=effective(T)/T,chi_bits=chi(r*effective(T),eps)) for T in (1,7,28,90,365)],
      single_atom=dict(one_year=-math.expm1(-65.4e-36*365*DAY),years_for_95=-math.log(.05)/(65.4e-36*365*DAY)),
      mars=dict(rate=r,eps=eps,full_coverage_eps=DP/PE,chi_bits_day=chi(r,eps),full_coverage_chi_bits_day=chi(r,DP/PE),expected_inventory_steady=r*HALF/LN2,day_exclusive_full_onoff_mi=exclusive(r)),
      two_slots=[dict(T_days=T,late=effective(T),early=effective(T)*math.exp(-LN2*T/HALF)) for T in (1,7,28)],
      marked=[marked(d) for d in (.1,1e-4,1e-8)],mi_checks=[dict(a=a,b=b,mi=binary_poisson_mi(a,b)) for a,b in ((0,.01),(0,1),(1,1),(100,100.01),(100,110))],
      threshold=threshold(E,g),resonance=resonance(float(max(g))),gmax_per_eV=float(max(g)))
    print(json.dumps(out,separators=(',',':')))
if __name__=='__main__':run()
