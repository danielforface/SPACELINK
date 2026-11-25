"""E019: natural Earth preparation and information invariant."""
import hashlib, urllib.request, pathlib, math, json, sys
import numpy as np
REF='1783f9f976c2f445ccf686c25ba003baeee8f357'
BLOBS={'Earth_Density.csv':'403bc19d86ec3d942af056bce9d9d09d8bfbcb57','be7ground_shape.csv':'19d97e87e8b03678d2e4629598f7d7fd7794eae0'}
HB=1.973269804e-7; DM=np.array([0.,7.5e-5,2.5e-3]); U=np.sqrt([.693*.978,.307*.978,.022]); W=U**2
P=np.array([.603837916,.373943596,.022218488]); P/=P.sum(); V=7.6e-14*19.3*74/183.84
R=6371000.; E0=862000.; Q=math.pi*HB*E0/DM[1]; F=4.8968563e-9
DLOG=np.log2(P[:,None]/P[None,:]) # initial i, final j

def source(name,local=None):
    if local is None:
        with urllib.request.urlopen(f'https://raw.githubusercontent.com/michelelucente/PEANUTS/{REF}/Data/{name}',timeout=20) as r: b=r.read(2000001)
    else: b=(pathlib.Path(local)/name).read_bytes()
    if len(b)>2000000 or hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()!=BLOBS[name]: raise ValueError('source hash')
    return np.loadtxt(b.decode().splitlines(),delimiter=','),hashlib.sha256(b).hexdigest()

def density(r,tab):
    i=np.minimum(np.searchsorted(tab[:,0],r),len(tab)-1)
    _,a,b,c=tab[i].T
    return a+b*r*r+c*r**4

def dagger(a): return a.conj().swapaxes(-1,-2)
def vac(E,L): return np.exp(-1j*DM*L/(2*E[:,None]*HB))
def layer(E,pot,L):
    h=np.eye(3)*DM/(2*E[:,None,None])+pot*np.outer(U,U)
    ev,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*ev[:,None,:]*L/HB))@dagger(v)

def earth(E,chord,step,tab):
    if chord<0 or chord>2*R: raise ValueError('chord')
    a=np.broadcast_to(np.eye(3,dtype=complex),(len(E),3,3)).copy()
    if chord==0:return a
    n=math.ceil(chord/step); dx=chord/n
    for k in range(n):
        x=(k+.5)*dx-chord/2
        r=math.sqrt(1-(chord/2/R)**2+(x/R)**2)
        a=layer(E,7.6e-14*float(density(r,tab)),dx)@a
    return a

def bound(L): return math.sin(min(math.pi/2,abs(V*L/HB)/2))**2*math.log2(max(P)/min(P))

def measure(E,A,L):
    # Rrel removes vacuum phase; diagonal readout unchanged by it.
    m=vac(E,L).conj()[:,:,None]*layer(E,V,L)
    rho=(A*P)@dagger(A)
    # stable diagonal-to-diagonal transfer and coherence cross terms
    t=abs(m)**2; pp=rho.diagonal(axis1=1,axis2=2).real
    d=np.zeros_like(pp)
    for i in range(3):
        for j in range(3):
            if i!=j:
                z=pp[:,i]*t[:,j,i]; d[:,i]-=z; d[:,j]+=z
    off=rho.copy(); off[:,np.arange(3),np.arange(3)]=0
    cross=(m@off@dagger(m)).diagonal(axis1=1,axis2=2).real
    rel=dagger(A)@m@A
    kl=np.sum(abs(rel)**2*(P[:,None]*DLOG).T,axis=(1,2))
    return (d+cross)@W,d@W,cross@W,kl,np.max(abs(np.sum(d+cross,axis=1)))

def average(E,g,v): return float(np.trapezoid(g*v,E)/np.trapezoid(g,E)) if hasattr(np,'trapezoid') else float(np.trapz(g*v,E)/np.trapz(g,E))

def grid(E,g,tab,step):
    rows=[]; unit=0.; conservation=0.; bound_ratio=0.
    for chord in (0,30000,100000,300000,1000000):
        # Three equally weighted angular proxy paths, no angular claim.
        paths=[earth(E,chord+x,step,tab) for x in ((-5000,0,5000) if chord else (0,))]
        unit=max(unit,max(float(np.max(abs(dagger(a)@a-np.eye(3)))) for a in paths))
        for gap in (0,Q,2*Q,3*Q):
            a=[vac(E,gap)[:,:,None]*m for m in paths]
            for L in (.1,1.,10.):
                z=[measure(E,m,L) for m in a]
                vals=np.mean([x[:4] for x in z],axis=0)
                pe,diag,coh,kl=[average(E,g,v) for v in vals]
                conservation=max(conservation,max(x[4] for x in z)); bound_ratio=max(bound_ratio,max(float(np.max(x[3]))/bound(L) for x in z))
                rows.append(dict(chord_m=chord,gap_m=gap,L_m=L,delta_pe=pe,diagonal=diag,coherence=coh,relative_entropy_bits=kl,bound_bits=bound(L)))
    return dict(rows=rows,unitarity_error=unit,conservation_error=conservation,max_bound_ratio=bound_ratio)

def two_state_oracle():
    # independent closed-form two-state uniform density vs eigensolution
    p=.5; E=E0; L=100000.; v=7.6e-14*1.67*.978
    c=.386; s=2*math.sqrt(.693*.307); z=(v*c-DM[1]/(2*E))/2; x=v*s/2; om=math.hypot(x,z)
    h=np.array([[z,x],[x,-z]])
    ev,a=np.linalg.eigh(h); m=(a*np.exp(-1j*ev*L/HB))@a.T
    exact=(x/om*math.sin(om*L/HB))**2
    return abs(float(abs(m[1,0])**2)-exact)

def sidebands(E,g):
    tv=float(g[0]+g[-1]+sum(abs(np.diff(g))))
    # density normalized to eV by caller
    phi=V/HB
    return [dict(f_Hz=f,shift_eV=4.135667696e-15*f,TV_bound=4.135667696e-15*f*tv*phi**2/4) for f in (1e6,1e8,1e9)]

def main(local=None):
    tab,sh1=source('Earth_Density.csv',local); a,sh2=source('be7ground_shape.csv',local)
    assert np.all(np.diff(tab[:,0])>0) and tab[-1,0]==1 and np.min(density(np.linspace(0,1,1001),tab))>0
    E=a[:,0]*1e6; g=a[:,1]; g/=average(E,np.ones_like(E),g)*(E[-1]-E[0])
    fine=np.sort(np.r_[E,(E[:-1]+E[1:])/2]); gf=np.interp(fine,E,g)
    coarse=grid(E,g,tab,5000); final=grid(fine,gf,tab,2500)
    err=max(abs(a['delta_pe']-b['delta_pe']) for a,b in zip(coarse['rows'],final['rows']))
    assert final['max_bound_ratio']<=1+1e-8 and final['unitarity_error']<1e-10 and final['conservation_error']<1e-12
    assert two_state_oracle()<1e-14
    best=max((r for r in final['rows'] if r['L_m']==1),key=lambda r:abs(r['delta_pe']))
    # zero-background Poisson KL chi-square reference, not receiver-independent bound
    ratio=1/4.76123463182996; rel=F*(1-ratio)*best['delta_pe']/(ratio+(1-ratio)*.52641136)
    rate=48.3/1.524**2
    print(json.dumps(dict(e019='DERIVED/SIMULATED',sha256=[sh1,sh2],numpy=np.__version__,coarse=coarse,final=final,spatial_spectral_max_abs_error=err,two_state_oracle_error=two_state_oracle(),best=best,aggregate_relative=rel,aggregate_bits_day=rate*rel**2/math.log(2),sidebands=sidebands(E,g)),separators=(',',':')))
if __name__=='__main__':main(sys.argv[1] if len(sys.argv)>1 else None)
