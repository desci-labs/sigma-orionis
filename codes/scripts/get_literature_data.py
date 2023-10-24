# ======================== Import Packages ==========================

import sys, os, pdb, glob
import numpy as np
from astropy.table import Table, join
from astroquery.vizier import Vizier
import calc_dust_masses
import warnings
from astropy.logger import AstropyWarning
warnings.filterwarnings('ignore', category=AstropyWarning)


# ===================== Define Functions ===================

def get_data(catalog, join_key='Name'):

    """
    PURPOSE:    Get data from literature with Vizier

    INPUT:      catalog = ctalog name on Vizier (str)
                join_key = column header to join tables, if multiple (str; optional)

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
                t = join(t, tt, join_type='inner', keys=join_key)

    return t


def fix_table(t, region, col_name, col_spt, col_flux, col_eflux, col_mstar, col_emstar, dist, obs_wave, col_lflux=None, log_mstar=False):

    """
    PURPOSE:    Make tables for different regions consistent

    INPUT:      t = table from literature (AstroPy Table)
                region = name of region (str)
                col_X = column name for data values (str)
                dist = distance to source in pc (float arr)
                obs_wave = wavelength of observation in microns (float)
                col_lflux = column name if fluxes for non-detections are upper limits (str; optional)
                log_mstar = boolean flag if stellar mass is log (boolean; optional)

    OUTPUT:     t = fixed data table (AstroPy Table)

    """
   
    ### DO SOME EXTRA FIXING FOR SOME REGIONS
    t['ObsWave'] = obs_wave
    det_lim = 3.0
    if region == 'tau':

        ### USE 1.3mm CONTINUUM FLUXES THAT WERE MEASURED RATHER THAN ESTIMATED
        ind_1300 = np.where(t['Notes'] == 'e, m')
        t[col_flux][ind_1300] = t['F1.3'][ind_1300]
        t[col_lflux][ind_1300] = t['l_F1.3'][ind_1300]
        t['ObsWave'][ind_1300] = 1300.

        ### CONVERT TO mJy
        t[col_flux] *= 1e3
        t[col_eflux] *= 1e3

        ### FIX SPECTRAL TYPES
        t[col_spt] = [x.replace(')', '').replace('(', '').split('+/-')[0].split('-')[0] for x in t['SpT']]

        ### TAKE AVERAGE OF UPPER/LOWER STELLAR MASS UNCERTAINTIES
        t[col_emstar] = np.mean([abs(t['logM_3'] - t['b_logM_3']), abs(t['B_logM_3'] - t['logM_3'])], axis=0)

        ### FIX FLUX ERRORS
        t[col_eflux][t[col_name] == 'IRAS 04216+2603'] = 0.1 * t[col_flux][t[col_name] == 'IRAS 04216+2603'].data[0]

    if region == 'cha':

        ### ASSUME 20% STELLAR MASS ERRORS FOR THOSE WITH UNKNOWN ERRORS
        t.remove_column('Name')
        ind = np.where(t['b_logM_'] == 0.0)
        t['b_logM_'][ind] = t['logM_'][ind] - 0.434 * 0.2
        t['B_logM_'][ind] = t['logM_'][ind] - 0.434 * 0.2

        ### TAKE AVERAGE OF UPPER/LOWER STELLAR MASS UNCERTAINTIES
        t[col_emstar] = np.mean([abs(t['logM_'] - t['b_logM_']), abs(t['B_logM_'] - t['logM_'])], axis=0)

    if region == 'usc':

        ### TAKE AVERAGE OF UPPER/LOWER STELLAR MASS UNCERTAINTIES
        t[col_emstar] = np.mean(np.array([t['E_logM'], t['e_logM']]), axis=0)

    if region == 'sor':

        det_lim = 2.9

    ### FIX COLUMN NAMES
    t[col_name].name = 'Name'
    t[col_spt].name = 'SpT'
    t[col_mstar].name = 'Mstar'
    t[col_emstar].name = 'e_Mstar'
    t[col_flux].name = 'Flux'
    t[col_eflux].name = 'e_Flux'
    t['Dist'] = dist

    ### FLAG (NON-)DETECTIONS 
    ### CAN REPLACE WITH 3-SIGMA UPPER LIMITS; USE 3x ERROR IF TRUE ERRORS REPORTED
    if col_lflux is None:
        t['Det'] = np.repeat(0, len(t))
        t['Det'][t['Flux'] / t['e_Flux'] >= det_lim] = 1
        # t['Flux'][np.where(t['Det'] == 0)] = 3.0 * t['e_Flux'][np.where(t['Det'] == 0)]
    ### USE REPORTED FLUX IF UPPER LIMITS PROVIDED
    else:
        t['Det'] = np.repeat(1, len(t))
        t['Det'][np.where(t[col_lflux] == '<')] = 0

    ### FOR UPPER SCO, ONLY KEEP "PRIMORDIAL" DISKS 
    ### TO MATCH SAMPLES OF LUPUS & TAURUS
    if region == 'usc':
        t = t[np.where( (t['Type'] == 'Full') | (t['Type'] == 'Transitional') | (t['Type'] == 'Evolved'))] 
    # print(len(t['Det'][t['Det'] == 0]))

    ### ONLY KEEP STARS > 0.1 SOLAR MASS
    ### AND KEEP LUPUS SOURCES THAT HAVE UNKNOWN STELLAR MASSES
    if log_mstar:
        t['Mstar'] = 10**t['Mstar']
        t['e_Mstar'] = t['Mstar'] * t['e_Mstar'] / 0.434
    t = t[np.where((t['Mstar'] >= 0.1) | (t['Mstar'].mask == True))]
    # print(len(t['Det'][t['Det'] == 0]))

    ### CALCULATE DUST MASSES USING SAME METHOD AS LUPUS
    mdust, e_mdust = [], []
    for i, val in enumerate(t):
        mdust.append(calc_dust_masses.get_dustmass(1.0, t['ObsWave'][i], t['Flux'][i], t['Dist'][i], 20.))
        ### USE MEASURED ERRORS IF PROVIDED
        if col_lflux is None:
            e_mdust.append(calc_dust_masses.get_dustmass(1.0, t['ObsWave'][i], t['e_Flux'][i], t['Dist'][i], 20.))
        elif t[col_lflux][i] == '':
            e_mdust.append(calc_dust_masses.get_dustmass(1.0, t['ObsWave'][i], t['e_Flux'][i], t['Dist'][i], 20.))
        ### USE FLUX/3 IF ONLY UPPER LIMITS PROVIDED FOR NON-DETECTIONS
        else:
            e_mdust.append(calc_dust_masses.get_dustmass(1.0, t['ObsWave'][i], t['Flux'][i] / 3.0, t['Dist'][i], 20.))   
    t['MDust'] = mdust
    t['e_MDust'] = e_mdust

    return t['Name', 'SpT', 'Dist', 'Mstar', 'e_Mstar', 'Flux', 'ObsWave', 'Det', 'MDust', 'e_MDust']


# ========================== Code ==========================

### GET SIGMA ORI DATA
TS = get_data("J/AJ/153/240")
TS = fix_table(TS, 'sor', '__HHM2007_', 'SpT', 'F1.33', 'e_F1.33', 'Mass', 'e_Mass', np.repeat(385., len(TS)), 1330.)
TS.write('../output/data_sor.txt', format='ascii.ipac')

### GET LUPUS DATA
TL = get_data("J/ApJ/828/46")
TL = fix_table(TL, 'lup', 'Name', 'SpT', 'F890', 'e_F890', 'Mass', 'e_Mass', TL['Dist'], 890.)
TL.write('../output/data_lup.txt', format='ascii.ipac')

### GET TAURUS DATA
TT = get_data("J/ApJ/771/129/")
TT = fix_table(TT, 'tau', 'Name', 'SpT', 'F0.89', 'e_F0.89', 'logM_3', 'e_logM', np.repeat(140., len(TT)), 890., col_lflux='l_F0.89', log_mstar=True)
TT.write('../output/data_tau.txt', format='ascii.ipac')

### GET CHAM-I DATA
TC = get_data("J/ApJ/831/125")
TC = fix_table(TC, 'cha', '_2MASS', 'SpT', 'Fnu', 'e_Fnu', 'logM_', 'e_logM', np.repeat(160., len(TC)), 887., log_mstar=True)
TC.write('../output/data_cha.txt', format='ascii.ipac')

### GET UPPER SCO DATA
TU = get_data("J/ApJ/827/142")
TU = fix_table(TU, 'usc', '_2MASS', 'SpT', 'Snu', 'e_Snu', 'logM', 'e_logM', np.repeat(145., len(TU)), 880., log_mstar=True)
TU.write('../output/data_usc.txt', format='ascii.ipac')






