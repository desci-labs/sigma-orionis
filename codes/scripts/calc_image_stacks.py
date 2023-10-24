# ======================== Import Packages ==========================

import sys, os, pdb, glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.vizier import Vizier
import warnings
from astropy.logger import AstropyWarning
warnings.filterwarnings('ignore', category=AstropyWarning)


# ===================== Define Functions ===================

def readfits(file):

    """
    PURPOSE:    Read in FITS file and header info

    INPUT:      Path to FITS file (str)

    OUTPUT:     img = image (float arr)
                xcen, ycen = image center coordinates in pixels (float)
                xpix, ypix = image pixel width in deg/pix units (float)
                bmaj, bmin, bpa = beam major axis, minor axis, position angle (float)
                xcen_ra, ycen_ra = image center coordinates in deg units (float)

    """

    ### READ IN FITS FILE
    hdulist = fits.open(file)
    data = hdulist[0].data[0, 0, :, :]
    head = hdulist[0].header
    hdulist.close()

    ### GET HEADER INFO
    xcen = head['CRPIX1']
    ycen = head['CRPIX2']
    xpix = head['CDELT1']
    ypix = head['CDELT2']
    xcen_ra = head['CRVAL1']
    xcen_de = head['CRVAL2']
    bmaj = head['bmaj']
    bmin = head['bmin']
    bpa  = head['bpa']

    return(data, xcen, ycen, xpix, ypix, bmaj, bmin, bpa, xcen_ra, xcen_de)


def write_fits(img, line, stype, xpix_all, ypix_all, bmaj_all, bmin_all, bpa_all):

    """
    PURPOSE:    Write FITS file and header info

    INPUT:      Image to write (array)
                Line that was stacked (str)
                Pixel width in deg/pix units for all images in stack (array)
                Beam major axis, minor axis, position angle for all images in stack (array)

    OUTPUT:     Stacked image (FITS file)

    """

    os.system('rm ../output/stack_nd_'+line+'_'+str(i)+'.fits')
    hdu = fits.PrimaryHDU()
    hdu.data = img

    hdu.header['CRPIX1'] = hdu.header['NAXIS1']/2
    hdu.header['CRPIX2'] = hdu.header['NAXIS2']/2
    hdu.header['bmaj'] = bmaj_all.mean()
    hdu.header['bmin'] = bmin_all.mean()
    hdu.header['bpa'] = bpa_all.mean()
    hdu.header['cdelt1'] = xpix_all.mean()
    hdu.header['cdelt2'] = ypix_all.mean()

    hdu.writeto('../output/stack_'+line+'_'+stype+'.fits')


def crop_img(file_img, hw_as, c_obj):

    """
    PURPOSE:    Crop & center an image

    INPUT:      file_img = path to FITS file (str)
                hw_as = half-width of desired cropped image size in arcsec (float)
                c_obj = coordinates of object (AstroPy SkyCoords)

    OUTPUT:     img = cropped & centered image (float arr)
                width_pix = half-width of cropped image in pixels (int)
                xpix_img, ypix_img = image pixel width in deg/pix units (float)
                Beam major axis, minor axis, position angle (float)
                bmaj, bmin, bpa = beam major axis, minor axis, position angle (float)

    """

    ### LOAD IMAGE AND GET CENTER COORDINATES
    img, xcen_img, ycen_img, xpix_img, ypix_img, bmaj_img, bmin_img, bpa_img, xcen_ra_img, ycen_de_img = readfits(file_img)
    c_img = SkyCoord(xcen_ra_img, ycen_de_img, frame='icrs', unit='deg')
    
    ### CENTER IMAGE ON OBJECT LOCATION 
    dra, ddec = c_img.spherical_offsets_to(c_obj)
    width_pix = int(round(hw_as / (ypix_img * 3600.0)))
    xctr = xcen_img + dra.value / xpix_img
    yctr = ycen_img + ddec.value / ypix_img

    ### CROP IMAGE
    img = img[int(round(yctr - width_pix)):int(round(yctr + width_pix)),
              int(round(xctr - width_pix)):int(round(xctr + width_pix))]
    
    return img, width_pix, xpix_img, ypix_img, bmaj_img, bmin_img, bpa_img


def stackme(t, line):

    """
    PURPOSE:    Stack image

    INPUT:      Table of sources to be stacked (AstroPy Table)
                Line name (str; must be 'cont', '13CO', C18O')

    OUTPUT:     Stacked image (array)
                Pixel width in deg/pix units for all images in stack (array)
                Beam major axis, minor axis, position angle for all images in stack (array)

    """

    xpix_all, ypix_all = np.empty(len(t)), np.empty(len(t))
    bmaj_all, bmin_all, bpa_all = np.empty(len(t)), np.empty(len(t)), np.empty(len(t))
    
    for i, val in enumerate(t['__HHM2007_']):

        ### GET IMAGE FILE NAME & CHECK FILE EXISTS
        if (line == 'C18O'): suffix = '_C18O_mom0.fits'
        if (line == '13CO'): suffix = '_13CO_mom0.fits'
        if (line == '12CO'): suffix = '_12CO_mom0.fits'
        if (line == 'cont'): suffix = '_cont.fits'
        file_img = '../data/FITS/S_' + str(val) + suffix
        if os.path.isfile(file_img) is False:
            print('missing FITS file for ' + str(val), line)
            pdb.set_trace()

        ### GET COORDINATES OF OBJECT FROM PAPER TABLE
        de_obj = str(t['DEJ2000'][i].split(' ')[0][0]) + str(t['DEJ2000'][i].split(' ')[0][1:]) + 'd' + str(t['DEJ2000'][i].split(' ')[1]) + 'm' + str(t['DEJ2000'][i].split(' ')[2]) + 's'
        ra_obj = str(t['RAJ2000'][i].split(' ')[0]) + 'h' + str(t['RAJ2000'][i].split(' ')[1]) + 'm' + str(t['RAJ2000'][i].split(' ')[2]) + 's'
        c_obj = SkyCoord(ra_obj, de_obj, frame='icrs')
    
        ### CROP IMAGE
        img_cont, width_pix, xpix_img, ypix_img, bmaj_img, bmin_img, bpa_img = crop_img(file_img, 8.0, c_obj)
        xpix_all[i], ypix_all[i] = xpix_img, ypix_img
        bmaj_all[i], bmin_all[i], bpa_all[i] = bmaj_img, bmin_img, bpa_img

        ### PUT INTO MJY UNITS
        img_cont = 1e3 * img_cont 

        ### ADD IMAGE TO STACK ARRAY
        if (i==0):
            img_all = np.zeros([2 * width_pix, 2 * width_pix, 1])
            temp = img_cont.reshape((2 * width_pix, 2 * width_pix, 1))
            img_all = temp
            
        else:
            temp = img_cont.reshape((2 * width_pix, 2 * width_pix, 1))
            img_all = np.append(img_all, temp, axis=2)
            

    ### COLLAPSE IMAGE STACK ARRAY
    stacked = np.sum(img_all, 2) / len(t)

    return stacked, xpix_all, ypix_all, bmaj_all, bmin_all, bpa_all


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


def meanerr_idl(val, err):

   weight = 1.0 / (err**2)
   wsum   = np.sum(weight)
   
   sumx   = np.sum(weight * val)
   xmean  = sumx / wsum
   sigmam = np.sqrt(1.0 / wsum)

   return xmean, sigmam


# ========================== Code ==========================

### GET SIGMA ORIONIS DATA
T = get_data("J/AJ/153/240")

### INDEX DUST NON-DETECTIONS
ind_dust_nd = T['F1.33'] / T['e_F1.33'] < 2.9

### INDEX dust detections, 13CO non-detections, C18O non-detections
ind_gas_nd = ((T['F1.33'] / T['e_F1.33'] >= 2.9)  & (T['l_F12CO'] == '<') )

### CALCULATE MEAN + STANDARD ERROR
xmc, smc = meanerr_idl(np.array(T[ind_dust_nd]['F1.33']), np.array(T[ind_dust_nd]['e_F1.33']))
print("Continuum Mean = {0:.3f}".format(xmc),'; ', "Standard Error = {0:.3f}".format(smc))

### STACK IMAGES
lines = ['cont', '12CO', '13CO', 'C18O']
stype = ['CND', 'GND']
for n, nval in enumerate(lines):
    for i, val in enumerate([T[ind_dust_nd], T[ind_gas_nd]]):
        stacked, xpix_all, ypix_all, bmaj_all, bmin_all, bpa_all = stackme(val, lines[n])
        write_fits(stacked, lines[n], stype[i], xpix_all, ypix_all, bmaj_all, bmin_all, bpa_all)
