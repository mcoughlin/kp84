### Tomas Ahumada
# Functions to reduce KPED, LCO data
# and to interact with the GROWTH Marshal
# version July 28 2019

import warnings
import glob
import matplotlib.pyplot as plt 
# %matplotlib inline
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

from astroML.crossmatch import crossmatch_angular
from astropy.coordinates import SkyCoord

from astropy.io import fits
from astropy import wcs

import urllib
import urllib.request

import astropy.units as u
import astropy.coordinates as coord

from astroquery.vizier import Vizier

def get_ztf_cand(url_report_page, username, password):

    '''
    Query the HTML ZTF report page to get the name & coord of the candidates
    params:
        url_report_page: can modify the date to get more specific results. tip: copy/paste url from ztf
        username, password : marshal GROWTH user and password
    retunrs:
        name_,coord: names and coordinates for the candidates

    '''

    # create a password manager
    password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()

    # Add the username and password.
    # If we knew the realm, we could use it instead of None.
    top_level_url = "http://skipper.caltech.edu:8080/"
    password_mgr.add_password(None, top_level_url, username, password)

    handler = urllib.request.HTTPBasicAuthHandler(password_mgr)

    # create "opener" (OpenerDirector instance)
    opener = urllib.request.build_opener(handler)

    # use the opener to fetch a URL
    opener.open(url_report_page)

    # Install the opener.
    # Now all calls to urllib.request.urlopen use our opener.
    urllib.request.install_opener(opener)

    with urllib.request.urlopen(url_report_page) as url:
        data = url.read().decode()


    print('Loaded :',len(data),'ZTF objects')

    df_list = pd.read_html(data,header=0)

    coord=[]
    name_=[]
    for i in range(len(df_list[1]['Name (age)'])):
        if pd.notna(df_list[1]['Name (age)'][i]):
            name_.append(df_list[1]['Name (age)'][i][:12])
            coord.append(df_list[1]['RA  Dec'][i]) 

    coord=np.array(coord)
    name_=np.array(name_)
    
    ra_transient,dec_transient=[],[]

    for i in range(len(coord)):
        try:
            c = SkyCoord(coord[i].split('+')[0],'+'+coord[i].split('+')[1], unit=(u.hourangle, u.deg))
        except:
            c = SkyCoord(coord[i].split('-')[0],'-'+coord[i].split('-')[1], unit=(u.hourangle, u.deg))
        ra_transient.append( c.ra.deg)
        dec_transient.append( c.dec.deg)
        
    return name_,ra_transient,dec_transient

def gaia_query(ra_deg, dec_deg, rad_deg, maxmag=25,
               maxsources=1):
    """
    Query Gaia DR1 @ VizieR using astroquery.vizier
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                maxmag: upper limit G magnitude (optional)
                maxsources: maximum number of sources
    returns: astropy.table object
    """
    vquery = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS',
                             'phot_g_mean_mag','phot_r_mean_mag',
                             'e_Gmag',
                             'Plx', 'e_Plx', 'BP-RP', 'e_BPmag', 'e_RPmag',
                             'Teff', 'Rad', 'Lum'],
                    column_filters={"phot_g_mean_mag":
                                    ("<%f" % maxmag),
                                   "phot_r_mean_mag":
                                    ("<%f" % maxmag)},
                    row_limit = maxsources)

    field = SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    try:
        source = vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="I/345/gaia2")
        return source[0]
    except:
        return []

def ps1_query(ra_deg, dec_deg, rad_deg, maxmag=25,
               maxsources=1):
    """
    Query Pan-STARRS @ VizieR using astroquery.vizier
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                maxmag: upper limit G magnitude (optional)
                maxsources: maximum number of sources
    returns: astropy.table object
    """
    vquery = Vizier(columns=['Source', 'RAJ2000', 'DEJ2000',
                             'gmag','rmag','imag','zmag','ymag'],
                    column_filters={"gmag":
                                    ("<%f" % maxmag),
                                   "imag":
                                    ("<%f" % maxmag)},
                    row_limit = maxsources)

    field = SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')

    try:
        source = vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="II/349/ps1")
        return source[0]
    except:
        return []

def panstarrs_query(ra_deg, dec_deg, rad_deg, ndet=4, 
                    maxsources=10000,minmag = 16.5,maxmag = 19.5,diffmag = 0.05,
                    server=('https://archive.stsci.edu/'+
                            'panstarrs/search.php')): 
    """
    Query Pan-STARRS DR1 @ MAST
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field 
                                          radius in degrees
                ndet: minimum number of detection (optional)
                maxsources: maximum number of sources
                server: servername
    returns: astropy.table object
            with ra,dec,rmag,gmag,imag,zmag,Imag
    """
    import requests 
    from astropy.io.votable import parse_single_table 
 
    
    r = requests.get(server, 
    params= {'RA': ra_deg, 'DEC': dec_deg, 
             'SR': rad_deg, 'max_records': maxsources, 
             'outputformat': 'VOTable', 
             'ndetections': ('>%d' % ndet)}) 
 
    # write query data into local file 
    outf = open('panstarrs.xml', 'w') 
    outf.write(r.text) 
    outf.close() 
 
    # parse local file into astropy.table object 
    data = parse_single_table('panstarrs.xml')
    table = data.to_table(use_names_over_ids=True) 
    
    mask = (table["nDetections"]>ndet) * (table["rMeanPSFMag"] > minmag) * (table["rMeanPSFMag"] < maxmag) * (table["gMeanPSFMag"] > minmag) * (table["gMeanPSFMag"] < maxmag) * (table["iMeanPSFMag"] > minmag) * (table["iMeanPSFMag"] < maxmag) *  (table["iMeanPSFMag"] - table["iMeanKronMag"] < diffmag)
    table = table[mask]
    newtable = np.zeros(len(table), dtype=[("ra", np.double), ("dec", np.double), ("rmag", np.float), ("gmag", np.float), ("imag", np.float),("Imag", np.float),("zmag", np.float)])
    newtable["ra"] = table["raMean"]
    newtable["dec"] = table["decMean"]
    newtable["gmag"] = table["gMeanPSFMag"]
    newtable["rmag"] = table["rMeanPSFMag"]
    newtable["imag"] = table["iMeanPSFMag"]
    newtable["Imag"] = table["iMeanPSFMag"]-0.386*(table["iMeanPSFMag"]-table["zMeanPSFMag"])-0.397
    newtable["zmag"] = table["zMeanPSFMag"]
    return newtable
 


def sdss_query(ra,dec,radius=5*u.arcmin,minmag=16.5,maxmag=20):
    from astroquery.sdss import SDSS
    from astropy import coordinates as coords
    import astropy.units as u

    pos = coords.SkyCoord(ra,dec,unit="deg", frame='icrs')
    table = SDSS.query_region(pos,fields=['type','ra','dec','u','g','r','i','z','err_u','err_g','err_r','err_i','err_z'],radius=5*u.arcmin)
    mask = (table["type"]==6) * (table["r"] > minmag) * (table["r"] < maxmag) * (table["g"] > minmag) * (table["g"] < maxmag)* (table["i"] > minmag) * (table["i"] < maxmag)* (table["u"] > minmag) * (table["u"] < maxmag)
    table = table[mask]
    
    newtable = np.zeros(len(table), dtype=[("ra", np.double), ("dec", np.double), ("rmag", np.float),("umag", np.float), ("gmag", np.float), ("imag", np.float),("zmag", np.float),("Imag", np.float),("Umag", np.float)])
    newtable["ra"] = table["ra"]
    newtable["dec"] = table["dec"]
    newtable["gmag"] = table["g"]
    newtable["rmag"] = table["r"]
    newtable["umag"] = table["u"]
    newtable["imag"] = table["i"]
    newtable["zmag"] = table["z"]
    newtable["Imag"] = table["i"]-0.386*(table["i"]-table["z"])-0.397
    B = table["g"]+0.313*(table["g"]-table["r"])+0.2271
    newtable["Umag"] = B+0.77*(table["u"]-table["g"])-0.88
    
    return newtable


def get_ztf_object(name,username,password):
    '''
    Query the HTML ZTF single object page to get the name & coord of a candidate
    params:
        name: ztf id of the object of interest
        username, password : marshal GROWTH user and password
    retunrs:
        name_,ra,dec: names and coordinates for the candidates

    '''
    url_single = 'http://skipper.caltech.edu:8080/cgi-bin/growth/view_source.cgi?name='+name

    password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()

    # Add the username and password.
    # If we knew the realm, we could use it instead of None.
    top_level_url = "http://skipper.caltech.edu:8080/"
    password_mgr.add_password(None, top_level_url, username, password)
    handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
    opener = urllib.request.build_opener(handler) # create "opener" (OpenerDirector instance)
    opener.open(url_single) # use the opener to fetch a URL
    urllib.request.install_opener(opener) # Now all calls to urllib.request.urlopen use our opener.

    with urllib.request.urlopen(url_single) as url:
            data = url.read().decode()

    print('GROWTH data loaded')

    df_list = pd.read_html(data,header=0)
    name = df_list[1].keys()[1]
    print(name)	
    ra_tran = df_list[1].keys()[6].split(' ')[-2]
    dec_tran = df_list[1].keys()[6].split(' ')[-1]

    coord=np.array(ra_tran+dec_tran)
    name_=name
    print('Coordinates ',name, coord)

    return name,np.array([ra_tran]),np.array([dec_tran])

def do_crossmatch(ra_im,dec_im,ra_std,dec_std):   
    '''
    Do the crossmatching
        params: ra,dec of data set 1 and 2
    returns:
        index of the match 
 
    '''
    # get imaging data
   
    try:
        imX = np.empty((len(ra_im), 2), dtype=np.float64)
    except TypeError:
        imX = np.empty((1, 2), dtype=np.float64)
    imX[:, 0] = ra_im
    imX[:, 1] = dec_im

    # get standard stars
    stX = np.empty((len(ra_std), 2), dtype=np.float64)
    stX[:, 0] = ra_std
    stX[:, 1] = dec_std

    # crossmatch catalogs
    max_radius = 1.5 / 3600  # 1 arcsec
    dist, ind_im = crossmatch_angular(imX, stX, max_radius)
    match_im = ~np.isinf(dist)
    return ind_im,match_im

def do_sextractor(files_entire,workdir_fits_red,workdir_sex_out,workdir_cat,sextractor_dir,control=False):
    import os
    # run sextractor in the images from KPED 
    sex_command_0 = 'sex '+workdir_fits_red+files_entire+' -c '+sextractor_dir+'0_first_catalog.sex -CATALOG_NAME '+workdir_cat+files_entire[:-5]+'.cat'
    os.system(sex_command_0)
    if control:
        print (sex_command_0)
    #psf_command_1 = 'psfex '+workdir_sex_out+files_entire[:-5]+'.ldac -c '+sextractor_dir+'1_make.psfex'
    #os.system(psf_command_1)
    #if control:
    #    print (psf_command_1)
    #sex_command_2 = 'sex '+workdir_fits_red+files_entire+' -c '+sextractor_dir+'2_genZP.sex -CATALOG_NAME '+workdir_cat+files_entire[:-5]+'.cat -PSF_NAME '+workdir_sex_out+files_entire[:-5]+'.psf'
    #os.system(sex_command_2)
    #if control:
    #    print (sex_command_2)        

def do_sextractor_KPED(files_entire,workdir_fits_red,workdir_sex_out,workdir_cat,sextractor_dir,control=False):
    import os
    # run sextractor in the images from KPED 
    sex_command_0 = 'sex '+workdir_fits_red+files_entire+' -c '+sextractor_dir+'0_first_catalog.sex -CATALOG_NAME '+workdir_sex_out+files_entire[:-4]+'.ldac'
    os.system(sex_command_0)
    if control:
        print (sex_command_0)
    psf_command_1 = 'psfex '+workdir_sex_out+files_entire[:-4]+'.ldac -c '+sextractor_dir+'1_make.psfex'
    os.system(psf_command_1)
    if control:
        print (psf_command_1)
    sex_command_2 = 'sex '+workdir_fits_red+files_entire+' -c '+sextractor_dir+'2_genZP.sex -CATALOG_NAME '+workdir_cat+files_entire[:-4]+'.cat -PSF_NAME '+workdir_sex_out+files_entire[:-4]+'.psf'
    os.system(sex_command_2)
    if control:
        print (sex_command_2)
