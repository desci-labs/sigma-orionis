# ======================== Import Packages ==========================

import sys, os, pdb, glob
import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import SkyCoord
from decimal import Decimal, ROUND_UP
from astropy.table import Table, join, MaskedColumn
from astroquery.vizier import Vizier
import warnings
from astropy.logger import AstropyWarning
warnings.filterwarnings('ignore', category=AstropyWarning)


# ===================== Define Functions ===========================

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

    
def calc_sep(ra1, de1, dist1, ra2, de2, dist2):

    """
    PURPOSE:    Calculate separation between two objects in 2D and 3D space

    INPUT:      ra, de = RA, DE of object in degrees (AstroPy coordinates)
                dist = distance of object in parsecs (AstroPy units)

    OUTPUT:     sep.arcsec = projected (2D) separation in arcseconds (float)
                dis.pc = 3D separation in parsecs (float)

    """

    ### 2D & 3D COORDINATES FOR OBJECT 2
    c2d2 = SkyCoord(ra=ra2, dec=de2, frame='icrs')
    c3d2 = SkyCoord(ra=ra2, dec=de2, distance=dist2, frame='icrs')

    ### SEPARATION IN ARCSEC FROM 2D COORDINATES OF OBJECT
    c2d1 = SkyCoord(ra=ra1, dec=de1, frame='icrs')
    sep = c2d2.separation(c2d1)

    ### SEPARATION IN PARSEC FROM 3D COORDINATES OF OBJECT
    c3d1 = SkyCoord(ra=ra1, dec=de1, frame='icrs', distance=dist1)
    dis = c3d2.separation_3d(c3d1)

    return sep.arcsec, dis.pc  


# ============================= Code ==================================

#### LOAD IN SIGMA ORIONIS DATA
T = get_data("J/AJ/153/240")

### GET SEPARATIONS
r_AS, r_PC = [], []
for i, val in enumerate(T['__HHM2007_']):

    ### GET ASTROPY COORDS OF THIS OBJECT
    coord = SkyCoord(str(T['RAJ2000'][i])+' '+str(T['DEJ2000'][i]), unit=(u.hourangle, u.deg))

    ### CALCULATE SEPARATION FROM SIGMA ORI SYSTEM
    ras, rpc = calc_sep(coord.ra, coord.dec, 385.*u.pc, 84.68658*u.degree, -2.60003*u.degree, 385.*u.pc)

    ### SAVE OUTPUT
    r_AS.append(str(Decimal(str(ras)).quantize(Decimal('.01'), rounding=ROUND_UP)))
    r_PC.append(str(Decimal(str(rpc)).quantize(Decimal('.01'), rounding=ROUND_UP)))

### SAVE TABLE
TD = Table()
TD['__HHM2007_'] = np.copy(T['__HHM2007_'])
TD.add_column(MaskedColumn(name='R_as', data=r_AS))
TD.add_column(MaskedColumn(name='R_pc', data=r_PC))
TD.write('../output/sep_OB.txt', format='ascii.ipac')
