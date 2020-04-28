# EllipsoidalAmplitudes

Python module that calculates the expected semi amplitudes of the first three harmonics of the ellipsoidal effect.

In order to use this module, complete the following steps:

Download from GitHub "Aell_Harmonics" and save it in a local folder.

Open Python:

1. Import sys:
>> import sys

2. Append the "Aell_Harmonics" path:
>> sys.path.insert(0, '/the/local/folder/path/Aell_Harmonics/')

3. Import the package "Harmonics":
>> import Harmonics as H

4. Use the "Aell_PH" function to derive the semi-amplitudes of the first three harmonics:
>> H.Aell_PH(teff,logg,f,i,q,band)  
  
  
Function description:
Aell_PH derives the ellipsoidal-semi-amplitudes of the first three harmonics.

Inputs:  
teff - effective temperature of the primary star  
logg - surface gravity of the primary star  
f - Roche-lobe filling factor of the primary  
i - orbital inclination  
q - binary mass ratio  
band - observing band.  

Outputs:  
coeffA - semi-amplitudes of the first three harmonics (a1c, a2c, a3c).



