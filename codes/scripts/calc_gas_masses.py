# ======================== Import Packages ==========================

import os, sys, pdb, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join, MaskedColumn
import matplotlib as mpl
from astropy import constants as const
from astropy import units as u
from scipy.stats import spearmanr, linregress
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from astroquery.vizier import Vizier
import warnings
from astropy.logger import AstropyWarning
warnings.filterwarnings('ignore', category=AstropyWarning)

# ========================== Define Functions ==========================

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


def get_model_grid_old(f):

    ### LOAD FILE FROM WILLIAMS & BEST 2014 (2014ApJ...788...59W)
    gasgrid = Table.read(f, format='ascii.csv')

    ### LOAD RELEVANT MODEL OUTPUTS
    gasgrid['M_gas'].name = 'Mgas'
    gasgrid['f_2-1_13co'].name = 'F13CO21'
    gasgrid['f_2-1_12co'].name = 'F12CO21'

    ### REMOVE MASKED VALUES (1 BAD END-OF-LINE IN GRID FROM JPW)
    gasgrid = gasgrid[np.where(~gasgrid['F13CO21'].mask)]
    gasgrid = gasgrid[np.where(~gasgrid['F12CO21'].mask)]

    return gasgrid['M_star', 'Mgas', 'gamma', 'R_c', 'T_m1', 'T_a1','q', 'incl', 'F13CO21', 'F12CO21']


def get_cal_err(f, e, e_cal=0.10):

    """
    PURPOSE:    Calculate multiplication factor for adding 
                calibration error to flux measurement error
                (assumed to be 10% unless otherwise specified)

    INPUT:      f = measured flux (float)
                e = measurement error (float)
                e_cal = calibration error fraction (float; optional)

    OUTPUT:     Multiplication factor for adding calibration error

    """

    mult = np.sqrt((f * e_cal)**2 + (e)**2)

    return mult


def get_model_idx(g13, g12, f13, e13, f12, e12, d13, d12):

    """
    PURPOSE:    Index model grid points for a given flux measurement 

    INPUT:      g13 = model grid points for 13CO flux
                g18_hi = model grid points for C18O flux (ISM abundance)
                g18_lo = model grid points for C18O flux (low abundance)
                f13 = 13CO flux measurement (float)
                e13 = 13CO flux measurement error(float)
                f18 = C18O flux measurement (float)
                e18 = C18O flux measurement error(float)
                d13 = detection flags for 13CO flux (masked array)
                d18 = detection flags for C18O flux (masked array)

    OUTPUT:     ifit = indexes of model grid
                inote = note indicating type of detection

    """

    ### IF BOTH LINES DETECTED DETECTED
    ### INDEX GRID WITHIN ERRORS (MEASUREMENT + FLUX CAL)
    if (d13 == True) and (d12 == True):
        i13 = ( abs(g13 - f13)  < get_cal_err(f13, e13) )
        i12 = ( abs(g12 - f12)  < get_cal_err(f12, e12) )
        inote = "GD"

    ### IF ONLY 12CO DETECTED
    ### INDEX GRID WITHIN UPPER LIMIT FOR 13CO
    elif (d12 == True) and (d13 == False):
        i12 = (abs(g12 - f12) < get_cal_err(f12, e12) )
        i13 = (g13 < f13)
        inote = "D12"

    ### IF BOTH LINES UNDETECTED
    ### INDEX GRID WITHIN UPPER LIMITS FOR BOTH
    elif (d12 == False) and (d13 == False):
        i13 = (g13 < f13)
        i12 = (g12 < f12)
        inote = "ND"

    ### STEP CODE IF UNKNOWN RESULT
    else:
        print("Unknown result")
        pdb.set_trace()

    ### KEEP ONLY THOSE IN BOTH LINES
    ### I.E., CREATING BOX AROUND MEASUREMENT OR UPPER LIMIT
    ifit = i13 & i12

    return ifit, inote


def get_model_gasmass(gm, ifit, inote):

    """
    PURPOSE:    Get gas mass from model grid 

    INPUT:      gm = all gas masses from model grid 
                ifit = model grid indexes for a given flux measurement
                inote = note indicating type of detection

    OUTPUT:     mgas_fit = gas mass based on model fit
                mgas_min = lower limti of gas mass based on model fit
                mgas_max = upper limit of gas mass based on model fit
                mgas_note = note indicating type of gas mass estimate
    """

    nfit = np.count_nonzero(ifit)
    if (nfit > 0):
        
        mgasfit = gm[ifit]

        ### BOTH LINES DETECTED
        if (d13==True) and (d12==True):
            mgas_fit = 10**np.mean(np.log10(mgasfit))
            mgas_min, mgas_max = np.min(mgasfit), np.max(mgasfit)
            mgas_note = "GF"

        ### ONLY 12CO DETECTED
        elif (d12==True) and (d13==False):
            mgas_fit = 10**np.mean(np.log10(mgasfit))
            mgas_min, mgas_max = np.min(mgasfit), np.max(mgasfit)
            mgas_note = "GL"

        ### BOTH LINES UNDETECTED
        elif (d13==False) and (d12==False):
            mgas_fit = -99.0
            mgas_min, mgas_max = np.min(mgasfit), np.max(mgasfit)
            mgas_note = "UL"

        ### STOP CODE IF UNKNOWN RESULT
        else:
            print("Unknown flux measurement result")
            pdb.set_trace()
                  
    else:

        ### STOP CODE IF NO MATCHES TO MODEL GRID
        print("No matches to model grid")
        pdb.set_trace()

    ### COMBINE NOTES
    mgas_note = (', ').join([inote, mgas_note])
                            
    return mgas_fit, mgas_min, mgas_max, mgas_note



def get_gasmass(g, f13, e13, d13, f12, e12, d12):

    """
    PURPOSE:    Calculate gas mass

    INPUT:      g = model grid 
                f13 = 13CO flux measurement (float)
                e13 = 13CO flux measurement error(float)
                d13 = detection flags for 13CO flux (masked array)
                f12 = 12CO flux measurement (float)
                e12 = 12CO flux measurement error(float)
                d12 = detection flags for 12CO flux (masked array)

    OUTPUT:     mg_fit = gas mass based on model fit
                mgas_min = lower limti of gas mass based on model fit
                mgas_max = upper limit of gas mass based on model fit
                mgas_note = note indicating type of gas mass estimate
    """

    ### INDEX MODEL GRID FOR THIS FLUX MEASUREMENT
    i_fit, i_note = get_model_idx(g['F13CO21'], g['F12CO21'], f13, e13, f12, e12, d13, d12)

    ### CALCULATE GAS MASS
    mgas_fit, mgas_lo, mgas_hi, mgas_note = get_model_gasmass(g['Mgas'], i_fit, i_note)
    
    return mgas_fit, mgas_lo, mgas_hi, mgas_note


# ========================== Code ==========================

### LOAD IN SIGMA ORIONIS DATA
T = get_data("J/AJ/153/240")
### T[['__HHM2007_','F12CO','e_F12CO','l_F12CO','F13CO','e_F13CO','l_F13CO']]


### LOAD GAS MODEL GRID FROM WILLIAMS & BEST 2014 (2014ApJ...788...59W)
# G = get_model_grid('../input/apj495435t3_mrt.txt')
G = get_model_grid_old('../data/gasgrid.csv')

### GET GAS MASSES
mg_f, mg_m, mg_l, mg_h, mg_n = [], [], [], [], []
for i, val in enumerate(T['__HHM2007_']):

    ### GET GAS FLUXES IN JY, SCALED TO 140 PC TO MATCH MODEL GRID
    f2l = (385. / 140.)**2 / 1000.0
    f13, e13 = f2l * T['F13CO'][i], f2l * T['e_F13CO'][i]
    f12, e12 = f2l * T['F12CO'][i], f2l * T['e_F12CO'][i]

    ### FLAG (NON-)DETECTIONS 
    d13 = T['l_F13CO'][i] != '<'
    d12 = T['l_F12CO'][i] != '<'

    ### CALCULATE GAS MASSES
    mgf, mgl, mgh, mgn = get_gasmass(G, f13, e13, d13, f12, e12, d12)

    ### SAVE GAS MASSES (M_JUP), LIMITS, AND FLAGS
    ### WEIRD ROUNDING / FORMATTING IS TO MATCH OUTPUT IN PAPER
    if mgn == 'ND, UL':
        mg_n.append('<')
        mg_f.append(float("{0:.2e}".format(-99.)))
        mg_l.append(float("{0:.2e}".format(-99.)))
        mg_h.append(round(float("{0:.2e}".format(mgh))*(const.M_sun.cgs/const.M_jup.cgs).value, 1))

    elif mgn == 'D12, GL':
        mg_n.append('<')
        mg_f.append(float("{0:.2e}".format(-99.)))
        mg_h.append(round(float("{0:.2e}".format(mgh))*(const.M_sun.cgs/const.M_jup.cgs).value, 1))
        mg_l.append(float("{0:.2e}".format(-99.)))

    elif mgn == 'GD, GF':
        mg_n.append('')
        mg_f.append(round(float("{0:.2e}".format(mgf))*(const.M_sun.cgs/const.M_jup.cgs).value, 1))
        mg_l.append(round(float("{0:.2e}".format(mgl))*(const.M_sun.cgs/const.M_jup.cgs).value, 1))
        mg_h.append(round(float("{0:.2e}".format(mgh))*(const.M_sun.cgs/const.M_jup.cgs).value, 1))

    else:
        print("Unknown gas mass result")
        pdb.set_trace()

TG = Table()
TG['__HHM2007_'] = np.copy(T['__HHM2007_'])
TG.add_column(MaskedColumn(name='l_Mgas', data=mg_n))
TG.add_column(MaskedColumn(name='Mgas', data=mg_f, mask=[x==-99.0 for x in mg_f]))
TG.add_column(MaskedColumn(name='b_Mgas', data=mg_l, mask=[x==-99.0 for x in mg_l]))
TG.add_column(MaskedColumn(name='B_Mgas', data=mg_h, mask=[x==-99.0 for x in mg_h]))
TG.write('../output-test/gasmasses.txt', format='ascii.ipac')

### NEED TO FIGURE OUT WHY 1327 & 936 ARE DIFFERENT FROM PAPER
# for i, val in enumerate(T['__HHM2007_']):
#     if str(T['B_Mgas'][i]) != str(TG['B_Mgas'][i]):
#         print(val, str(T['B_Mgas'][i]), TG['__HHM2007_'][i], str(TG['B_Mgas'][i]))
    # pdb.set_trace()
# T[['__HHM2007_','F12CO','e_F12CO','l_F12CO','F13CO','e_F13CO','l_F13CO']][T['__HHM2007_'] == 1327]









