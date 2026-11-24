"""Screening calculation, not measurement or full flavor simulation."""
import math
import json
HBARC = 1.973269804e-7  # eV m
V_FACTOR = 7.6e-14  # eV per g/cm3 times Ye
def phase(rho, ye, length_m):
    if rho < 0 or not 0 <= ye <= 1 or length_m < 0:
        raise ValueError("invalid material or length")
    return V_FACTOR * rho * ye * length_m / HBARC
def transition_bound(phi):
    return math.sin(min(abs(phi), math.pi)/2)**2
def integration_days(modulation, mass_t=100, au=1, z=5, rate=134):
    if not 0 < modulation <= 1 or min(mass_t,au,z,rate) <= 0:
        raise ValueError("invalid detector parameters")
    return (z/modulation)**2 / (rate * mass_t/100 / au**2)
def main():
    rows=[]
    for name,rho,ye in [("water",1,.555),("tungsten",19.3,74/183.84)]:
        for length in (1,100,1000,100000):
            p=phase(rho,ye,length)
            rows.append(dict(material=name,length_m=length,phase_rad=p,transition_upper_bound=transition_bound(p)))
    print(json.dumps(dict(evidence="DERIVED screening only",material_rows=rows,ideal_days_at_1pct=integration_days(.01)),indent=2))
if __name__ == "__main__":
    main()
