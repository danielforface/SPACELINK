"""Public source metadata probe; never stores raw data. Not a license grant."""
import urllib.request, hashlib, json, math
URL="https://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs2005opflux.dat"
def main():
    try:
        with urllib.request.urlopen(URL,timeout=12) as r:
            data=r.read(2_000_001)
        if len(data)>2_000_000:raise ValueError("size limit")
        rows=[]
        for line in data.decode("ascii").splitlines():
            p=line.split()
            if len(p)!=13:continue
            try:row=[float(x) for x in p]
            except ValueError:continue
            rows.append(row)
        if len(rows)<100 or not all(0<x[0]<=1 and x[10]>=0 and all(math.isfinite(v) for v in x) for x in rows) or not all(a[0]<b[0] for a,b in zip(rows,rows[1:])):raise ValueError("unexpected schema")
        print(json.dumps(dict(source_probe="retrieved_metadata_only",url=URL,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data),rows=len(rows),be7_sum=sum(x[10] for x in rows),license="unverified; not used in experiments")))
    except Exception as e:
        print(json.dumps(dict(source_probe="blocked",url=URL,error=type(e).__name__,detail=str(e))))
if __name__=="__main__":main()
