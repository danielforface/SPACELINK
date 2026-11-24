"""Mass-basis screening; no physical signal or fitted solar model."""
import math, cmath, random, json
HBARC=1.973269804e-7
VFACTOR=7.6e-14
def eye(n): return [[complex(i==j) for j in range(n)] for i in range(n)]
def mul(a,b): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def scale(a,s): return [[x*s for x in row] for row in a]
def add(a,b): return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]
def expm(a):
    n=len(a); norm=max(sum(abs(x) for x in row) for row in a)
    k=max(0,math.ceil(math.log2(norm))) if norm else 0
    a=scale(a,2.**-k); t=eye(n); out=eye(n)
    for j in range(1,32):
        t=scale(mul(t,a),1/j); out=add(out,t)
    for _ in range(k): out=mul(out,out)
    return out
def electron(s12=.307,s13=.022):
    return [math.sqrt((1-s12)*(1-s13)),math.sqrt(s12*(1-s13)),math.sqrt(s13)]
def propagate(layers,E=862000.,s13=.022):
    u=electron(s13=s13); m=[0,7.5e-5,2.5e-3]; S=eye(3)
    for length,rho,ye in layers:
        v=VFACTOR*rho*ye
        h=[[v*u[i]*u[j]+(m[i]/(2*E) if i==j else 0) for j in range(3)] for i in range(3)]
        S=mul(expm(scale(h,-1j*length/HBARC)),S)
    return S
def leave(S,i): return sum(abs(S[j][i])**2 for j in range(3) if i!=j)
def bound(layers):
    phi=sum(VFACTOR*rho*ye*length/HBARC for length,rho,ye in layers)
    return math.sin(min(phi/2,math.pi/2))**2
def contrast(S,p=(.68,.30,.02),s13=.022):
    w=[v*v for v in electron(s13=s13)]
    return sum(p[i]*(w[j]-w[i])*abs(S[j][i])**2 for i in range(3) for j in range(3) if i!=j)
def analytic2(length,rho,ye,E=862000.):
    a,b,_=electron(s13=0); v=VFACTOR*rho*ye
    off=v*a*b; dz=(v*(a*a-b*b)-7.5e-5/(2*E))/2
    om=math.hypot(dz,off)
    return (off/om)**2*math.sin(om*length/HBARC)**2
def h2(x): return -x*math.log2(x)-(1-x)*math.log2(1-x) if 0<x<1 else 0.
def zchannel(lam=.81):
    q=-math.expm1(-lam); p=1/(q*(1+2**(h2(q)/q)))
    f=lambda p:h2(p*q)-p*h2(q)
    return dict(optimal_prior=p,capacity=f(p),grid_capacity=max(f(i/10000) for i in range(10001)),uncoded_ber=.5*math.exp(-lam),reported_rate=40/184/2.2)
def coefficient(): return math.sqrt(2)*1.1663787e-5*1e-18*6.02214076e23*(1.973269804e-5)**3
def sample_profiles():
    r=random.Random(20261124)
    for _ in range(100):
        yield [(10**r.uniform(0,5),r.uniform(0,20),r.uniform(.4,.6)) for _ in range(5)]
def diagnostics(S):
    return max(abs(sum(S[k][i].conjugate()*S[k][j] for k in range(3))-(i==j)) for i in range(3) for j in range(3))
def main():
    rows=[]
    for L in (1,100,1000,100000):
        ls=[(L,19.3,74/183.84)]; S=propagate(ls)
        rows.append(dict(length_m=L,bound=bound(ls),max_leave=max(leave(S,i) for i in range(3)),persistent_delta_pe=contrast(S),analytic2=analytic2(*ls[0]),numeric2=leave(propagate(ls,s13=0),0)))
    err=0.; excess=0.
    for ls in sample_profiles():
        S=propagate(ls); err=max(err,diagnostics(S)); excess=max(excess,max(leave(S,i)-bound(ls) for i in range(3)))
    print(json.dumps(dict(evidence="SIMULATED/DERIVED",coefficient=coefficient(),zchannel=zchannel(),rows=rows,random_profiles=100,max_unitarity_error=err,max_positive_bound_excess=excess),indent=2))
if __name__=="__main__": main()
