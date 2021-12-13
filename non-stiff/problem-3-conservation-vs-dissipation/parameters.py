import numpy as np

# Bulk modulus water, air
Kw = 2e9
Ka = 142e3


# Density water, air
rw = 999.7
ra = 1.225

ca = np.sqrt(Ka / ra)
cw = np.sqrt(Kw / rw)

print("Water density: ", rw, "kg/m^3")
print("Air density: ", ra, "kg/m^3")
print("Water speed of sound: ", cw, "m/s")
print("Air speed of sound: ", ca, "m/s")
print("Density contrast: ", rw/ra,  "water / air")
print("Wave speed contrast: ", cw/ca,  "water / air")
