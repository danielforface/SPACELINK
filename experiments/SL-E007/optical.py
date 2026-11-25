"""Illustrative optical budgets and binary click-channel information."""
import math,json
import numpy as np
H=6.62607015e-34;C=299792458.;K=1.380649e-23;T=5772.;LAM=1064e-9;D=225e9;DT=.22;DR=5.;TH=.00465
def spectral(lam=LAM):
    # Planck spectral shape normalized to total irradiance, W/m² per m.
    x=H*C/(lam*K*T)
    return 1361*15/math.pi**4*x**4/(lam*math.expm1(x))
def airy(x,n=128):
    if abs(x)<1e-4:return (1-x*x/8+x**4/192)**2
    t=np.linspace(0,math.pi,n+1)
    j=float(np.trapz(np.cos(t-x*np.sin(t)),t)/math.pi)
    return (2*j/x)**2
def budget(mode,band_nm,offset_urad,bg_fraction,fov_factor,distance=D):
    band=band_nm*1e-9;A=math.pi*(DR/2)**2
    E=spectral();L=E/(math.pi*TH**2)
    fov=1.22*LAM/DR*fov_factor
    bg=L*bg_fraction*A*math.pi*fov**2*band
    if mode=="mirror":
        P=E*band*(DR/(2*distance*math.tan(TH)))**2
    else:
        eta={"solar_tem00":.0206,"solar_cw_optimistic":.0396,"pv_laser_assumed":.09}[mode]
        P=1361*eta*(math.pi*DT*DR/(4*LAM*distance))**2*airy(math.pi*DT*offset_urad*1e-6/LAM)
    factor=.5*.6*LAM/(H*C)
    return P*factor,bg*factor+100
def h2(x):
    return -(x*math.log2(x)+(1-x)*math.log2(1-x)) if 0<x<1 else 0.
def phi(x):
    if abs(x)<1e-4:return x*x*(.5+x*(-1/6+x*(1/12-x/20)))
    if x<=-1:return 1.
    return (1+x)*math.log1p(x)-x
def mutual(p,q0,q1):
    if p==0 or p==1 or q0==q1:return 0.
    d=q1-q0
    # Center relative to the mixture, avoid subtracting nearly equal entropies.
    total=0.
    for m,delta in (((1-p)*q0+p*q1,d),((1-p)*(1-q0)+p*(1-q1),-d)):
        if m>0:total+=m*((1-p)*phi(-p*delta/m)+p*phi((1-p)*delta/m))
    return total/math.log(2)
def joint(p,q0,q1):
    o=(1-p)*q0+p*q1;v=0.
    for px,q in ((1-p,q0),(p,q1)):
        for a,b in ((q,o),(1-q,1-o)):
            if px*a>0:v+=px*a*math.log2(a/b)
    return v
def channel(s,b,t):
    q0=-math.expm1(-b*t);q1=-math.expm1(-(b+s)*t)
    I,p=max((mutual(i/100,q0,q1),i/100) for i in range(101))
    return dict(slot_s=t,bits_slot=I,prior_on=p,bits_s=I/t,ber_equal=(q0+1-q1)/2)
def grid():
    result=[]
    for mode in ("mirror","solar_tem00","solar_cw_optimistic","pv_laser_assumed"):
        for band in (.1,1,10):
            for offset in (0,1,3):
                for bg in (1e-10,1e-7,1e-4):
                    for fov in (1,5):
                        s,b=budget(mode,band,offset,bg,fov)
                        best=max((channel(s,b,t) for t in (1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3)),key=lambda x:x["bits_s"])
                        result.append(dict(mode=mode,band_nm=band,offset_urad=offset,bg_fraction=bg,fov_factor=fov,signal_s=s,background_s=b,**best))
    return result
def main():
    print(json.dumps(dict(optical_header="SIMULATED illustrative click channel",spectral_W_m2_nm=spectral()*1e-9)))
    for row in grid():print(json.dumps(dict(optical_row=row),separators=(",",":")))
if __name__=="__main__":main()
