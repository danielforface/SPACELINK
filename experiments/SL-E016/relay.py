"""Ideal conditional two-hop link reference."""
import math,json
C=299792458.;K=1.380649e-23;LAM=1064e-9
def optical(d,dt,dr):
    return min(1.,(math.pi*dt*dr/(4*LAM*d))**2)
def relay(d,f,diam):
    p=4.;first=p*optical(d*f,.22,diam)
    second=optical(d*(1-f),diam,5)
    direct=p*optical(d,.22,5)
    return dict(distance_m=d,fraction=f,relay_diameter_m=diam,first_rx_W=first,passive_rx_W=.5*first*second,direct_rx_W=direct,passive_ratio=.5*first*second/direct,active_second_rx_W=p*second,active_local_output_W=p,one_way_latency_s=d/C)
def radio(d):
    B=1e8;T=100.;a=1.;rx=math.pi*2.5**2;lam=C/2.8e9
    factor=min(1.,a*rx/(lam*lam*d*d))
    p=[70e-22*B*a,1.,1361*.3*.3]
    return dict(distance_m=d,transmitter_W=p,rx_W=[v*factor for v in p],awgn_bits_s=[B*math.log1p(v*factor/(K*T*B))/math.log(2) for v in p])
def main():
    for d in (1e6,1e9,225e9):
        for f in (.01,.1,.5,.9,.99):
            for a in (.22,1,5):print(json.dumps(dict(e016_relay=relay(d,f,a))))
        print(json.dumps(dict(e016_radio=radio(d))))
if __name__=="__main__":main()
