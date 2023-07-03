import numpy as np

mu_0 = 1.0
epsilon_0 = 1.0

# From Clayton Paul's book: Analysis of Multiconductor Transmission Lines. P. 175.
# Assumes wide separation of wires. Distances in mm.
rS = 50 # shield radius.
d1 = 25 # distance of wire 1 to center.
rw1 = 2 # wire 1 radius.
d2 = 25 # distance of wire 2 to center.
rw2 = 2  # wire 2 radius.
theta12 = np.pi # angle between conductors.
mu = mu_0
epsilon = epsilon_0

l11 = mu / (2.0*np.pi) * np.log((rS**2 - d1**2)/(rw1*rS)) 
l22 = l11
l12 = mu / (2.0*np.pi) * np.log(d1/rS * \
    np.sqrt(
        ((d1*d2)**2 + rS**4 - 2 * d1 * d2* rS**2 *np.cos(theta12))/ \
        ((d1*d2)**2 + d2**4 - 2 * d1 * d2**3 * np.cos(theta12))
    )
)
l21 = l12

Lmat = np.array([[l11, l12], 
                 [l21, l22]])
Cmat = epsilon * mu * np.linalg.inv(Lmat)

print("Inductance matrix [natural units]: ")
print(Lmat)

print ("Capacitance matrix [natural units]: ")
print(Cmat)