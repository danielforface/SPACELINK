"""Ideal tagging and interaction scales, not a physical link."""
import importlib.util, pathlib, math, json, random
spec=importlib.util.spec_from_file_location("link",pathlib.Path(__file__).resolve().parents[1]/"SL-E003"/"link.py")
e3=importlib.util.module_from_spec(spec);spec.loader.exec_module(e3)
AU=149597870700.
def tagged(S,f,r,t,l,B=0):
    if S<=0 or B<0 or not all(0<=v<=1 for v in (f,t,l)):raise ValueError("invalid channel")
    vals=[]
    for a,b in [(t,l),(1-t,1-l)]:
        lam=S*(f*a+(1-f)*b)+B/2
        d=S*f*r*a
        vals.append(d*d/lam/math.log(2) if lam else 0.)
    return sum(vals)
def count(S,f,r,B=0):return (S*f*r)**2/(S+B)/math.log(2)
def row(L,a):
    x=e3.link(L,a);f=x["coverage"];r=e3.link(L,0,full=True)["relative_upper"];S=48.3/1.524**2
    return dict(length_m=L,radius_m=a,f=f,angle_rad=a/(.524*AU),count=count(S,f,r),perfect=tagged(S,f,r,1,0),partial=[dict(false_positive=l,bound=tagged(S,f,r,.9,l)) for l in (1e-3,1e-6,1e-9,0)],perfect_bg100=tagged(S,f,r,1,0,100))
def cevns(E=.862,Z=74,A=184,rho=19.3,L=1):
    q=A-Z-(1-4*.2317)*Z
    sigma=1.1663787e-5**2*(E/1000)**2*q*q/(4*math.pi)*.389379338e-27
    tau=rho*6.02214076e23/A*sigma*L*100
    return dict(sigma_cm2=sigma,tau=tau,scatter_probability=-math.expm1(-tau))
def spin(mu,BL):
    x=mu*5.7883818e-5*BL/1.973269804e-7
    return math.sin(x)**2
def gravity(M,b):return 4*6.67430e-11*M/(b*299792458.**2)
def main():
    print(json.dumps(dict(evidence="DERIVED screening only",rows=[row(L,a) for L in (1,1000) for a in (1,1000,100000)],cevns=[dict(length_m=L,**cevns(L=L)) for L in (1,1000)],spin=[dict(BL_Tm=bl,sm=spin(3.2e-19*.05,bl),hypothetical=spin(6.4e-12,bl)) for bl in (1,1e4)],matter_index_contribution=7.6e-14*19.3*(74/183.84)/862000,gravity_rad=dict(earth_limb=gravity(5.9722e24,6.371e6),sun_limb=gravity(1.98847e30,6.957e8),local=gravity(1000,1))),indent=2))
if __name__=="__main__":main()
