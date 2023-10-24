# ======================== Import Packages ==========================

import sys,os,pdb,glob
import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.table import Table, join, MaskedColumn
from astroquery.vizier import Vizier
import warnings
from astropy.logger import AstropyWarning
warnings.filterwarnings('ignore', category=AstropyWarning)


# ===================== Define Functions ===================

def get_data(catalog, join_key='Name', join_type='inner'):

    """
    PURPOSE:    Get data from literature with Vizier

    INPUT:      catalog = ctalog name on Vizier (str)
                join_key = column header to join tables, if multiple (str; optional)
                join_type = way to join tables, if multiple (str; optional)

    OUTPUT:     t = data table (AstroPy Table)

    """

    ### GET FULL CATALOG (ALL COLUMNS, ALL ROWS)
    viz = Vizier(catalog=catalog, columns=['**'])
    viz.ROW_LIMIT = -1
    tv = viz.get_catalogs(catalog)

    ### IF MULTIPLE TABLES, JOIN THEN
    for i, val in enumerate(tv.keys()):
        if i == 0:
            t = tv[val]
        else:
            tt = tv[val]
            if join_key in tt.columns:
                t = join(t, tt, join_type=join_type, keys=join_key)

    return t

    
def fix_units(b, l, f, d, t):
    
    """
    PURPOSE:    Change units of inputs to work with
                dust mass calculations

    INPUT:      b = dust opacity power-law index (unitless; float)
                l = wavelength of observations (microns; float)
                f = observed flux (mJy; float)
                d = istance to disk (pc; float)
                t = isk temperature (K; float)

    OUTPUT:     Inputs but now with correct units

    """

    b = float(b) 
    t = float(t) * u.K 
    d = float(d) * u.pc.to(u.cm) * u.cm
    l = float(l) * u.um.to(u.cm) * u.cm 
    f = float(f) * u.mJy

    return b, l, f, d, t


    
def get_planck(l, t):

    """
    PURPOSE:    Calculate Planck function

    INPUT:      l = wavelength of observations in cm (float)
                t = disk temperature in K (float)

    OUTPUT:     planck = Planck function (float)

    """

    ### CONVERT WAVELENGTH TO FREQUENCY
    nu = (const.c.cgs / l)

    ### CALCULATE PLANCK FUNCTION
    plank = 2 * const.h.cgs * nu**3 / const.c.cgs**2 / (np.exp(const.h.cgs *nu / (const.k_B.cgs * t)) -1)

    return plank


def get_kappa(l, b):

    """
    PURPOSE:    Calculate kappa parameter for dust mass calculation

    INPUT:      l = wavelength of observations in cm (float)
                b = dust opacity power-law index (float)

    OUTPUT:     kappa = kappa parameter for dust mass calculation (float)

    """

    ### CONVERT WAVELENGTH TO FREQUENCY
    nu = (const.c.cgs / l)

    ### ASSUMES NORMALIZATION OF 2.3 cm^2/g at 230 GHz
    nf = 2.3 / ( (230. / 1000.)**b )

    ### CALCULATE KAPPA PARAMETER
    kappa = (nf * (nu / (1e9 * 1000 * u.Hz))**b ) * ((u.cm**2)/(u.g))
        
    return kappa 


def get_mult(d, kappa, b_nu):

    """
    PURPOSE:    Calculate multiplication factor that produces dust mass
                when applied to source flux in mJy
                (matches Eq. 1 in Ansdell+2016 when d=150pc)

    INPUT:      d = distance to source in cm (float)
                kappa = kappa parameter from get_kappa() (float)
                b_nu = Planck function from get_planck() (float)

    OUTPUT:     mult = mltiplication factor that produces dust mass
                       when applied to source flux in mJy (float)

    """

    mult = (1e-26) * (d**2) / (kappa * b_nu) / const.M_sun.cgs
        
    return mult


def get_dustmass(b, l, f, d, t):
    
    """
    PURPOSE:    Calculate disk dust mass using prescription from
                Hildebrand 1983 (1983QJRAS..24..267H)

    INPUT:      b = dust opacity power-law index (unitless; float)
                l = wavelength of observations (microns; float)
                f = observed flux (mJy; float)
                d = distance to disk (pc; float)
                t = disk temperature (K; float)

    OUTPUT:     Disk dust mass in Earth masses (float)

    """

    b, l, f, d, t = fix_units(b, l, f, d, t)
    k_nu = get_kappa(l, b)
    b_nu = get_planck(l, t)
    mult = get_mult(d, k_nu, b_nu)
    mdust = mult * f * (const.M_sun.cgs / const.M_earth.cgs)

    return round(mdust.value, 5)


# ============================= Code ==================================

def main():

    #### LOAD IN LUPUS DATA
    T = get_data("J/ApJ/828/46")

    ### GET DUST MASSES
    DM, e_DM = [], []
    for i, val in enumerate(T['Name']):

        ### SET INPUTS TO DUST MASS CALCULATION
        Beta_Index = 1.                   # DUST OPACITY POWER-LAW INDEX (UNITLESS)
        Lambda_Obs = 890.                 # WAVELENGTH OF OBSERVATIONS (MICRONS)
        Temp_Disk = 20.                   # DISK TEMPERATURE (KELVIN)
        Dist_Disk = T['Dist'][i]           # DISTANCE TO DISK (PC)
        Flux_Disk = T['F890'][i]         # OBSERVED FLUX (mJy)
        e_Flux_Disk = T['e_F890'][i]     # ERROR ON OBSERVED FLUX (mJy)

        ### CALCULATE DUST MASS
        DM.append(get_dustmass(Beta_Index, Lambda_Obs, Flux_Disk, Dist_Disk, Temp_Disk))
        e_DM.append(get_dustmass(Beta_Index, Lambda_Obs, e_Flux_Disk, Dist_Disk, Temp_Disk))

    TD = Table()
    TD['Name'] = np.copy(T['Name'])
    TD.add_column(MaskedColumn(name='MDust', data=DM))
    TD.add_column(MaskedColumn(name='e_MDust', data=e_DM))
    TD.write('../output/dustmasses.txt', format='ascii.ipac')

if __name__ == "__main__":
    main()