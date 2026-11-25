"""E021 radioactive temporal memory; idealized computational models only."""
import math,json
import numpy as np
HALF=11.43; LN2=math.log(2); LAM=LN2/HALF; RATE=.8382666805346743; EPS=4.749291240080495e-18
def survive(age,half=HALF): return math.exp(-LN2*age/half) if not math.isinf(half) else 1.
def pulse_age(n,age,half=HALF):
    p=survive(age,half); lam=LN2/half
    return dict(n=n,age=age,p=p,sd_days=math.sqrt((1-p)/(n*p))/lam,n_for_1day=(1-p)/(p*lam**2))
def decay_filter(period,half=HALF):
    lam=LN2/half; w=2*math.pi/period
    h=lam/complex(lam,w)
    return dict(period_days=period,amplitude=abs(h),information_ratio=abs(h)**2,phase_lag_deg=-math.degrees(np.angle(h)),time_lag_days=-np.angle(h)/w)
def quad_filter(period,n,half=HALF):
    x=np.linspace(0,40,n+1); w=2*math.pi/period/(LN2/half)
    return complex(np.trapz(np.exp(-x-1j*w*x),x))
def dual(A,age,h2=35.):
    lam2=0. if math.isinf(h2) else LN2/h2
    u=A*.5*survive(age);v=A*.5*survive(age,h2)
    J=u*v/(u+v)*(LAM-lam2)**2
    return dict(A=A,age=age,half2=None if math.isinf(h2) else h2,means=[u,v],sd_days=1/math.sqrt(J),births_for_1day=A/J)
def kernel(halves,ages):
    return np.exp(-LN2*np.array(ages)[None,:]/np.array(halves)[:,None])/len(halves)
def collision():
    K=kernel([HALF,35.],[0,HALF,2*HALF]);d=np.cross(K[0],K[1]);d/=max(abs(d))
    a=10*np.maximum(d,0);b=10*np.maximum(-d,0)
    return dict(a=a.tolist(),b=b.tolist(),mean_a=(K@a).tolist(),mean_b=(K@b).tolist(),error=float(max(abs(K@(a-b)))))
def matrix(m):
    halves=np.geomspace(2,64,m) if m>1 else np.array([HALF])
    edges=np.arange(9)*7.;lam=LN2/halves
    K=(np.exp(-lam[:,None]*edges[:-1])-np.exp(-lam[:,None]*edges[1:]))/lam[:,None]/7/m
    normalized=K/np.linalg.norm(K,axis=1,keepdims=True)
    s=np.linalg.svd(normalized,compute_uv=False)
    return dict(clocks=m,rank=int(np.linalg.matrix_rank(normalized)),nullity=8-int(np.linalg.matrix_rank(normalized)),singular_values=s.tolist(),condition_nonzero=float(s[0]/s[-1]))
def pairing(w,b):
    return dict(window_days=w,background_day=b,true_accept=-math.expm1(-LAM*w),accidental=-math.expm1(-b*w))
def monte():
    rng=np.random.default_rng(213);out=[]
    for n in (100,10000):
        k=rng.binomial(n,.5,100000);v=-np.log(k[k>0]/n)/LAM
        out.append(dict(n=n,samples=len(k),zero=int(sum(k==0)),mean_days=float(v.mean()),sd_days=float(v.std(ddof=1)),fisher_sd=pulse_age(n,HALF)['sd_days']))
    return out
def run():
    periods=[1,7,HALF,28,90,365];rows=[]
    for p in periods:
        r=decay_filter(p);h=LAM/complex(LAM,2*math.pi/p)
        r['quadrature_error_100k']=abs(quad_filter(p,100000)-h);r['quadrature_error_200k']=abs(quad_filter(p,200000)-h)
        r['amplitude_2024_half']=decay_filter(p,11.468)['amplitude'];rows.append(r)
    print(json.dumps(dict(e021='DERIVED/SIMULATED',rate_day=RATE,mean_lifetime_days=1/LAM,terminal_collision=[10*survive(0),20*survive(HALF)],age=[pulse_age(n,t) for n in (10,100,10000) for t in (1,HALF,2*HALF)],monte=monte(),filter=rows,dual=[dual(A,t,h) for A in (10,100,10000) for t in (1,HALF,28,90) for h in (35.,math.inf)],collision=collision(),matrix=[matrix(m) for m in (1,2,4,8)],pair=[pairing(w,b) for w in (1,HALF,2*HALF,5*HALF) for b in (.01,.1,1)],prompt_chi_bits_day=RATE*EPS**2/LN2,prompt_chi_bits_year=365*RATE*EPS**2/LN2),separators=(',',':')))
if __name__=='__main__':run()
