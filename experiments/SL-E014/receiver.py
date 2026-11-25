"""Held-out simulated pair receiver; no physical data."""
import numpy as np, json
def sample(seed,n,signal,drift,clip,repeat=1):
    rng=np.random.default_rng(seed)
    truth=rng.integers(0,2,n)
    bits=np.repeat(truth,repeat)
    j=np.arange(len(bits)*2);phase=rng.uniform(0,2*np.pi)
    bg=20*(1+drift*np.sin(2*np.pi*j/157+phase))
    on=np.stack([bits,1-bits],axis=1).reshape(-1)
    counts=rng.poisson(bg+signal*(.1+.8*on)).reshape(-1,2)
    return truth,np.minimum(counts,clip)
def calibrate(obs,bits,clip):
    best=min(range(clip+2),key=lambda t:np.mean((obs[:,0]>=t)!=bits))
    hist=np.ones((2,clip+1,clip+1))
    np.add.at(hist,(bits,obs[:,0],obs[:,1]),1)
    hist/=hist.sum(axis=(1,2),keepdims=True)
    return best,np.log(hist[1]/hist[0])
def decode(obs,model):
    threshold,ll=model
    d=obs[:,0].astype(float)-obs[:,1]
    return dict(fixed=(obs[:,0]>=threshold).astype(int),differential=(d>0).astype(int),empirical=(ll[obs[:,0],obs[:,1]]>0).astype(int))
def repeat_decode(obs,k=5):
    d=(obs[:,0].astype(float)-obs[:,1]).reshape(-1,k)
    hard=(d>0).sum(axis=1)>k//2
    soft=d.sum(axis=1)>0
    return dict(hard=hard,soft=soft)
def summarize(z):
    return dict(mean=float(np.mean(z)),block_sd=float(np.std(z,ddof=1)),min=float(min(z)),max=float(max(z)))
def run(signal,drift,clip):
    y,x=sample(7100,1800,signal,drift,clip);model=calibrate(x,y,clip)
    err={k:[] for k in ("fixed","differential","empirical","hard","soft")}
    for i in range(20):
        y,x=sample(8100+i,600,signal,drift,clip)
        pred=decode(x,model)
        for k,v in pred.items():err[k].append(float(np.mean(v!=y)))
        y,x=sample(9100+i,120,signal,drift,clip,5)
        for k,v in repeat_decode(x).items():err[k].append(float(np.mean(v!=y)))
    return dict(signal=signal,drift=drift,clip=clip,uncoded_bits=12000,coded_bits=2400,slots_each=24000,errors={k:summarize(v) for k,v in err.items()})
def main():
    for signal in (0,4,12,30):
        for drift in (0,.6):
            for clip in (31,63,255):
                r=run(signal,drift,clip)
                if signal==0 and any(abs(v["mean"]-.5)>.025 for v in r["errors"].values()):raise ValueError("null check")
                print(json.dumps(dict(e014_row=r),separators=(",",":")))
if __name__=="__main__":main()
