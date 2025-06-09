#%%
import numpy as np
import scipy.constants as sc

Cmat0 = np.array([
    [22.494, -11.247],
    [-11.247, 16.581]
])*1e-12

Lmat = np.linalg.inv(Cmat0)*sc.epsilon_0*sc.mu_0

print ("Capacitance matrix: ")
print(Cmat0)

print("Inductance matrix: ")
print(Lmat)