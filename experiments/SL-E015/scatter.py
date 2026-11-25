"""Heavy spin-zero nucleus neutral-current angular reference."""
import numpy as np, math,json
GF=1.1663787e-5;HB=197.3269804;CONV=.389379338e-27
def ff(x):
    x=np.asarray(x,dtype=float)
    out=np.empty_like(x);small=abs(x)<.05
    z=x[small]**2;out[small]=1-z/10+z*z/280-z*z*z/15120
    y=x[~small];out[~small]=3*(np.sin(y)-y*np.cos(y))/y**3
    return out
def point(E,A,Z):
    Q=A-Z-(1-4*.2317)*Z
    return GF**2*(E/1000)**2*Q**2/(4*math.pi)*CONV
def integ(E,A,lo,hi,n=128,point_only=False,backward=False):
    # h=1-cos(theta) for forward, h=1+cos(theta) for backward.
    x,w=np.polynomial.legendre.leggauss(n)
    h=lo+(hi-lo)*(x+1)/2
    q=E*np.sqrt(2*(2-h if backward else h))
    f=1 if point_only else ff(q*1.2*A**(1/3)/HB)**2
    return float((hi-lo)/2*np.dot(w,(h/2 if backward else 1-h/2)*f))
def run(E,A,Z):
    t=integ(E,A,0,2);t2=integ(E,A,0,2,256)
    if abs(t-t2)>1e-10:raise ValueError("quadrature")
    sigma=point(E,A,Z)
    tau=19.3*6.02214076e23/A*sigma*1e5 if A==184 else None
    rows=[]
    for a in (1e-8,1e-4,.01,.1):
        h=2*math.sin(a/2)**2
        f=integ(E,A,0,h);b=integ(E,A,0,h,backward=True)
        rows.append(dict(angle_rad=a,forward_fraction=f/t,backward_fraction=b/t,forward_point=h-h*h/4,backward_point=h*h/4,W_1km_forward_probability=tau*f if tau is not None else None))
    return dict(E_MeV=E,A=A,Z=Z,sigma_point_cm2=sigma,sigma_cm2=sigma*t,form_factor_integral=t,convergence=abs(t-t2),W_1km_tau=tau*t if tau is not None else None,cones=rows)
def main():
    for A,Z in ((40,18),(184,74),(208,82)):
        for E in (.862,10,50):print(json.dumps(dict(e015_row=run(E,A,Z)),separators=(",",":")))
if __name__=="__main__":main()
