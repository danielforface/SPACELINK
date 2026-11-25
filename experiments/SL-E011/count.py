"""Ideal Poisson count channel with bounded numerical support."""
import pathlib,importlib.util,math,json,random
import numpy as np
s=importlib.util.spec_from_file_location("e9",pathlib.Path(__file__).resolve().parents[1]/"SL-E009"/"optical.py")
e9=importlib.util.module_from_spec(s);s.loader.exec_module(e9)
def support(mu):
    r=14*math.sqrt(mu+1)+20
    return range(max(0,int(mu-r)),math.ceil(mu+r)+1)
def pmf(k,mu):
    if mu==0:return float(k==0)
    return math.exp(k*math.log(mu)-mu-math.lgamma(k+1))
def term(a,b,p):
    m=(1-p)*a+p*b
    if not m:return 0.
    d=b-a
    return m*((1-p)*e9.o7.phi(-p*d/m)+p*e9.o7.phi((1-p)*d/m))/math.log(2)
def count_mi(mu0,mu1,p):
    if not 0<p<1 or min(mu0,mu1)<0:raise ValueError("channel")
    ks=sorted(set(support(mu0))|set(support(mu1)))
    a=[pmf(k,mu0) for k in ks];b=[pmf(k,mu1) for k in ks]
    mass0=math.fsum(a);mass1=math.fsum(b)
    I=math.fsum(term(x,y,p) for x,y in zip(a,b))
    ber=math.fsum(min((1-p)*x,p*y) for x,y in zip(a,b))
    return I,ber,max(abs(1-mass0),abs(1-mass1))
def reference(mu0,mu1,p,n=300):
    a=math.exp(-mu0);b=math.exp(-mu1);z=0.
    for k in range(n):
        # Normalize the two probabilities before mixing to avoid a
        # subnormal mixture rounding to zero in the far tail.
        scale=max(a,b)
        if scale:
            u,v=a/scale,b/scale
            m=(1-p)*u+p*v
            if u:z+=(1-p)*a*math.log2(u/m)
            if v:z+=p*b*math.log2(v/m)
        a*=mu0/(k+1);b*=mu1/(k+1)
    return z
def row(s,b,t,g,p):
    c,ber,err=count_mi(b*t,(b+s)*t,p);q0,q1=e9.probs(s,b,t)
    click=e9.o7.mutual(p,q0,q1)
    return dict(gate_s=t,guard_s=g,prior=p,signal_s=s,background_s=b,count_bits_s=c/(t+g),click_bits_s=click/(t+g),count_ber=ber,mass_error=err)
def main():
    x,y,_,sha=e9.retrieve()
    for mode in ("mirror","solar_tem00","pv_laser_assumed"):
        for width in (1,10):
            for bg in (1e-10,1e-7,1e-4):
                s,b=e9.budget(x,y,mode,width,bg)
                for t in ((1e-6,1e-5,.001,.01) if mode=="mirror" else (1e-8,1e-7,1e-6)):
                    for p in (.25,.5):
                        print(json.dumps(dict(e011_row=dict(mode=mode,width_nm=width,bg=bg,**row(s,b,t,1e-6 if mode=="mirror" else 20e-9,p))),separators=(",",":")))
if __name__=="__main__":main()
