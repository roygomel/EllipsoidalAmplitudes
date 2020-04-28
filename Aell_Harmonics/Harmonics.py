# Python Packages:
import numpy as np
from scipy.interpolate import griddata
import os

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# These functions return the linear limb and gravity coefficients of a primary star 
# with logg and teff, observed in a specific band.
def claret2019(Tab,logg,teff,band,Type):
    
    """
    Inputs:
        Tab - LD/GD table that was pre loaded with:
    'logg','Teff','Z','xi','u','Filt','Met','Mod' columns (LD table)
    'logg','logTeff','Z','xi','y','Filt','Mod' columns (GD table)     
        logg,teff,band - parameters for which the value is interpolated (A sub table is taken for Z=0) 
        Type - Limb/Gravity Darkening Coefficient ('LD'/'GD')
    
    Output:
        Limb/Gravity Darkening Coefficient
    """
    
    teff = min(max(teff,3500),40000)
    logg = min(max(logg,0),5)
    if len(band) == 1:
        band = band + ' '
    
    if Type == 'LD':
        coeffString = 'u'
        teffString = 'Teff'
        teffVal = teff
    else:
        coeffString = 'y'
        teffString = 'logTeff'
        teffVal = np.log10(teff)
    
    claret = Tab[(Tab[:]['Z']==0.) & (Tab[:]['Filt']==band.encode('UTF-8'))]
    points = np.column_stack((claret[:][teffString],claret[:]['logg']))       
    vals = claret[:][coeffString]
    coeff= np.asscalar( griddata(points,vals,(teffVal,logg),method = 'linear') )

    return coeff

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Obtaining a1c using PHOEBE simulations:
def a1c_PH(Tab,teff,logg,q,f,i):
    
    """
    Inputs:
        Tab - a1c PHOEBE table that was pre loaded with:
    'q','f','sin2i','M','R' columns
        teff,logg,q,f,i - parameters for which the value is interpolated: 
        Effective temperature and surface gravity of the primary, binary mass ratio, 
        Roch-Lobe filling factor of the primary and the orbital inclination.

    Output:
        a1c - first harmonic semi-amplitude in one of the following bands: B,V,R,I.
    """

    sin2i = (np.sin(np.radians(i)))**2
    sin2i = min(max(sin2i,0.2),1)
    
    coeff = np.nan
    if all([q>=0.1 ,q<=10 ,f>=0.3 ,f<=0.9 ,sin2i>=0.2 ,sin2i<=1 ,teff>=4935 ,teff<=16700, logg>=4.0722, logg<=4.5963 ]):
        points = np.column_stack((Tab[:]['q'],(Tab[:]['F'])**4,Tab[:]['sin2i'],Tab[:]['Teff'],Tab[:]['logg']))       
        vals = Tab[:]['a1c']
        coeff= np.asscalar( griddata(points,vals,(q,f**4,sin2i,teff,logg),method = 'linear') )

    return coeff

# Obtaining the semi-amplitudes of the first three harmonics, using our analytical correction terms:
def Aell_Analytic(teff,logg,u1,tau1,q,i,f,band,Home_path):
    
    """
    Inputs:
        teff - effective temperature of the primary
        logg - surface gravity of the primary
        u1 - linear limb-darkening coefficient
        tau1 - linear gravity-darkening coefficient
        q - binary mass ratio
        i - orbital inclination
        f - Roch-Lobe filling factor of the primary
        band - observing band
        Home_path - Home-directory path

    Output:
        a1c, a2c, a3c - semi-amplitudes of the first three harmonics.
    """
    
    a1c, a2c, a3c = np.nan, np.nan, np.nan

    # Calculating the MN93 semi-amplitudes of the second and third harmonics:
    sini = np.sin(np.radians(i));
    a2cMN = -3*(15+u1)*(1+tau1)/(20*(3-u1))*(f*E(q))**3*q*sini**2
    Correction = -15*(1-u1)*(3+tau1)*(f*E(q))**5*q*(6*sini**2-7*sini**4)/(64*(3-u1))
    Lave = 1 + (15+u1)*(1+tau1)*(f*E(q))**3*(2+5*q)*(2-3*sini**2)/(60*(3-u1)) + 9*(1-u1)*(3+tau1)*(f*E(q))**5*q*(8-40*sini**2+35*sini**4)/(256*(3-u1))    
    
    a2cMN = (a2cMN + Correction)/Lave    
    a3cMN = -25*u1*(2+tau1)*(f*E(q))**4*q*sini**3/32/(3-u1)/Lave
    
    # Second-harmonic Correction-term parameters (for large filling factors):
    a2 = 1.09090812767358000
    b2 = 0.03791808778223120
    c2 = 0.00504447910111281
    d2 = 0.04464586115427450
    
    # Third-harmonic Correction-term parameters (for large filling factors):
    a3 = 0.2074876857271310;
    b3 = 0.0698076131977346;
    c3 = 2.0222797170372300;
    d3 = 0.3879775920855080;

    # analytical corrections:
    if band in {'B','V','R','I'}:
        FN = 'A1C_' + band + '.txt'
        a1cTab = np.loadtxt(os.path.join(Home_path,FN),
                       dtype ={'names':('q','F','sin2i','Teff','logg','a1c'),
                               'formats':('float','float','float','float','float','float')}, delimiter = ',',skiprows=1)   
        a1c = a1c_PH(a1cTab,teff,logg,q,f,i)
    
    a2c = (1 + (-f/(f - a2)) * (b2 + c2/(d2 + q))) * a2cMN
    
    if all([q >= 0.1, q<=10, f<=0.9]):
        a3c = (1 + (f**6 + a3*f**2 + b3*q*(sini)**2*f**6)/(c3*f + sini**4 + d3*f*np.log(q))) * a3cMN
    
    coeff = np.asarray([a1c,a2c,a3c])
    
    return(coeff)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# RL/a (=E) using Eggleton approximation
def E(q):
    E = 0.49 * q**(-2./3.) / (0.6*q**(-2./3.) + np.log(1.+q**(-1./3.)))
    return(E)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def Aell_PH(teff,logg,f,i,q,band):      

    """
    Inputs:
        teff - effective temperature of the primary
        logg - surface gravity of the primary
        f - Roche-lobe filling factor of the primary
        i - orbital inclination
        q - binary mass ratio
        band - observing band

    Output:
        coeffA - semi-amplitudes of the first three harmonics (a1c, a2c, a3c).
    """
    
    # Home Directory
    Home_path = str(os.path.dirname(os.path.abspath('Harmonics')))
    
    # Reading Claret11 tables into numpy array:
    LD_Claret = np.loadtxt(os.path.join(Home_path,'Claret_LD.tsv'),
                           dtype ={'names':('logg','Teff','Z','xi','u','Filt','Met','Mod'),
                                   'formats':('f4','f4','f4','f4','f4','S2','S1','S1')}, delimiter = ';',skiprows=1)
    
    GD_Claret = np.loadtxt(os.path.join(Home_path,'Claret_GD.tsv'),
                           dtype ={'names':('logg','logTeff','Z','xi','y','Filt','Mod'),
                                   'formats':('f4','f4','f4','f4','f4','S2','S1')}, delimiter = ';',skiprows=1)
    
    # Deriving limb and gravity darkening coefficients:
    U1 = claret2019(LD_Claret,logg,teff,band,'LD')
    TAU1 = claret2019(GD_Claret,logg,teff,band,'GD')
         
    # Deriving the semi amplitudes of the first three harmonics:
    coeffA = Aell_Analytic(teff,logg,U1,TAU1,q,i,f,band,Home_path)
    return(coeffA)    
        



