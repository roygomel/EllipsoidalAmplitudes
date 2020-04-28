# EllipsoidalAmplitudes

Python module that calculates the expected semi amplitudes of the first three harmonics, based on the PHOEBE simulated light curves.

Use the following steps:

Download from GitHub "Aell_Harmonics" and save it in a local folder.

Open Python:

// First import sys:
1. import sys

// Second append the "Aell_Harmonics" path:
2. sys.path.insert(0, '/the/local/folder/path/Aell_Harmonics/')

// Third import the package "Harmonics":
3. import Harmonics as H

// Forth use the "Aell_PH" function to derive the semi-amplitudes of the first three harmonics:
4. H.Aell_PH(teff,logg,f,i,q,band)

Function description:
Aell_PH derives the ellipsoidal-semi-amplitudes of the first three harmonics.

Inputs:
teff - effective temperature of the primary star
logg - surface gravity of the primary star
f - Roche-lobe filling factor of the primary
i - orbital inclination
q - binary mass ratio
band - observing band

Outputs:
coeffA - semi-amplitudes of the first three harmonics (a1c, a2c, a3c).



