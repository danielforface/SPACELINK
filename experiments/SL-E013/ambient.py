"""Conditional ambient RF and thermal reference calculations."""
import math,json
C=299792458.;K=1.380649e-23;SIG=5.670374419e-8;F=2.8e9;B=1e8;LAM=C/F;RX=math.pi*2.5**2
def transfer(A,D):return min(1.,A*RX/(LAM*LAM*D*D))
def radio(sfu,A,D,T):
    pin=sfu*1e-22*B*A;factor=transfer(A,D);pr=pin*factor
    def cap(p):return B*math.log1p(p/(K*T*B))/math.log(2)
    return dict(sfu=sfu,aperture_m2=A,distance_m=D,T_K=T,incident_W=pin,rx_W=pr,awgn_ref_bits_s=cap(pr),active_1W_ref_bits_s=cap(factor))
def thermal(eta,r,T,mc,mr,P=4):
    a=P/(eta*1361);q=r*(1361*a-P);rad=q/(.85*SIG*T**4)
    return dict(eta=eta,heat_fraction=r,T_K=T,collector_kg_m2=mc,radiator_kg_m2=mr,collector_m2=a,heat_W=q,radiator_m2=rad,area_mass_kg=mc*a+mr*rad)
def main():
    for s in (64,70,200,1000):
        for A in (1,100):
            for D in (1e6,1e9,225e9):
                for T in (20,100,300):
                    print(json.dumps(dict(e013_radio=radio(s,A,D,T)),separators=(",",":")))
    for r in (.1,.5,1):
        for T in (280,320,360):
            for mc in (1,5):
                for mr in (2,10):
                    rows=[thermal(e,r,T,mc,mr) for e in (.0206,.0396,.09)]
                    print(json.dumps(dict(e013_thermal=dict(rows=rows,pv_extra_mass_break_even_kg=[x["area_mass_kg"]-rows[-1]["area_mass_kg"] for x in rows[:-1]])),separators=(",",":")))
    print(json.dumps(dict(e013_asymmetric=[dict(direct=thermal(.0206,r,320,5,10),pv=thermal(.09,1,320,5,10)) for r in (.01,.1,.5,1)])))
if __name__=="__main__":main()
