"""Energy-dependent adiabatic screen and scalar grating analogy."""
import pathlib,importlib.util,math,json
import numpy as np
s=importlib.util.spec_from_file_location("e6",pathlib.Path(__file__).resolve().parents[1]/"SL-E006"/"solar.py")
e6=importlib.util.module_from_spec(s);s.loader.exec_module(e6)
e5=e6.e5;e3=e6.e3;HB=1.973269804e-7
def profile():
    name,file,sha,ncol,idx=e5.MODELS[0];data=e6.download(file,sha);rows,norm=e5.parse(data,sha,ncol,idx)
    den=[]
    for l in data.decode("ascii").splitlines():
        t=l.split()
        if len(t)!=ncol:continue
        try:z=list(map(float,t))
        except ValueError:continue
        den.append(z[2])
    r,w=np.array(rows).T;w=w/norm
    f=np.array([e3.shell_fraction(1000/((1-1/1.524)*a*e3.RS)) for a in r])
    return np.array(den),w,f
def mixture(E,den,w,f):
    H=np.diag(np.array([0,7.5e-5,2.5e-3])/(2*E))[None,:,:]+(7.6e-14*10**den)[:,None,None]*np.outer(e6.U,e6.U)
    _,v=np.linalg.eigh(H)
    p=abs(np.einsum("rji,j->ri",v,e6.U))**2
    full=w@p;cond=(w*f)@p/(w@f)
    red=np.array([e6.reduced(d,E) for d in den])
    err=max(abs((full-w@red)@e6.W),abs((cond-(w*f)@red/(w@f))@e6.W))
    return full,cond,float(err)
def run(E,g,den,w,f):
    pairs=[mixture(e,den,w,f) for e in E];full,cond,errors=map(np.array,zip(*pairs))
    f0=float(w@f);frozen=mixture(862000.,den,w,f)[1]
    out=[]
    for L in (0,1,1000):
        ss=[e3.oracle([(L,19.3,74/183.84)],e) for e in E]
        dp=np.array([e5.mass_delta(s,p)@e6.W for s,p in zip(ss,cond)])
        dp0=np.array([e5.mass_delta(s,frozen)@e6.W for s in ss])
        for t in (0,.25,.5):
            ce=np.array([e3.xs(e/1e6,t,True) for e in E]);cx=np.array([e3.xs(e/1e6,t,False) for e in E])
            B=cx+(ce-cx)*(full@e6.W);Bcond=cx+(ce-cx)*(cond@e6.W)
            num=float(np.trapz(g*(ce-cx)*dp,E));num0=float(np.trapz(g*(ce-cx)*dp0,E))
            denom=float(np.trapz(g*B,E));denomc=float(np.trapz(g*Bcond,E))
            out.append(dict(length_m=L,threshold_MeV=t,conditional_relative=num/denomc,full_relative=f0*num/denom,frozen_relative=f0*num0/denom))
    return dict(coverage=f0,max_reduction_error=float(max(errors)),conditional_p_min=cond.min(axis=0).tolist(),conditional_p_max=cond.max(axis=0).tolist(),rows=out)
def diffraction(L,alpha=.001):
    k=862000/HB;theta=math.asin(math.pi/(k*.2e-9));kap=7.6e-14*19.3*74/183.84/(HB*math.sin(theta))
    a=kap*L;d=k*math.cos(theta)
    return dict(normal_length_m=L,kappa_per_m=kap,theta_rad=theta,kappa_L=a,resonant_scalar=math.tanh(a)**2,first_zero_angle_rad=math.pi/(d*L),weak_born_average_bound=min(a*a,math.pi*kap*kap*L/(2*alpha*d)) if a<=.01 else None)
def main():
    E,g=e6.line(e6.download("be7ground_shape.csv","19d97e87e8b03678d2e4629598f7d7fd7794eae0"));den,w,f=profile()
    a=run(E,g,den,w,f);fine=np.sort(np.r_[E,(E[1:]+E[:-1])/2]);b=run(fine,np.interp(fine,E,g),den,w,f)
    print(json.dumps(dict(e012="SIMULATED/DERIVED",coarse=a,refined=b,diffraction=[diffraction(L) for L in (.001,.01,.1,1,10)]),separators=(",",":")))
if __name__=="__main__":main()
