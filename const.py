#------------------------------------------------------------------------------
# Declare some constants
#------------------------------------------------------------------------------
# Experimental/empirical
from math import pi  # 3.14159...
a     = 6.371e6      # earth radius, meters
Omega = 7.292115e-5  # rotation rate, radians/second
H     = 7.0e3        # scale height rule-of-thumb
g     = 9.80665      # standard gravity
G     = 6.67408e-11  # gravitational constant; m3 kg-1 s-2
R     = 8.31446      # Ideal gas constant; J K-1 mol-1
Na    = 6.02214e23   # Avo's number; mol-1
h     = 6.62607e-34  # Planck constant; J s
c     = 2.99792458e8 # Speed of light (exact because we define meters in terms of this); m s-1
Md    = 28.9645e-3   # dry air molar mass; kg mol-1
Mw    = 18.0153e-3   # water vapor molar mass; kg mol-1
cp    = 1.005e3      # specific heat at T=300K, surface press; J kg-1 K-1
p0    = 1013.25
# Derivative constants
kb    = R/Na  # Boltzmann constant; J K-1
Rd    = R/Md  # dry air gas constant; J K-1 kg-1
Rm    = R/Mw  # water vapor gas constant; J K-1 kg-1
kappa = Rd/cp # poisson constant
sigma = (2*(pi**5)*(kb**4))/(15*(h**3)*(c**2)) # Stefan-boltzmann constant (see wikipedia); W m-2 K-4
