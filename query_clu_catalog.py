import numpy as np
import os, glob
import requests, json
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy import constants as const
from astropy.table import Table, Column, MaskedColumn
import threading
import time
from tqdm import tqdm
import warnings
from datetime import date
from astropy.time import Time
from astropy.io import ascii

today = date.today()
warnings.filterwarnings('ignore')

# Define global parameters
global TOKEN, BASEURL
with open('user_info.json') as usr:
    usr_data = json.load(usr)

GETTOKEN = usr_data['user']['FritzToken']
BASEURL = 'https://fritz.science/'

# Empty dict.
db = dict()
db['source'] = []

def api(method, endpoint, data=None):
    ''' Info : Basic API query, takes input the method (eg. GET, POST, etc.), the endpoint (i.e. API url)
               and additional data for filtering
        Returns : response in json format
        CAUTION! : If the query doesn't go through, try putting the 'data' input in 'data' or 'params'
                    argument in requests.request call
    '''
    headers = {'Authorization': f'token {GETTOKEN}'}

    response = requests.request(method, endpoint, json=data, headers=headers)

    return response.json()

def dist_mod_mag(app_mag, distance):
    """
    Calculate the absolute magnitude via the distance modulus.

    Input:
    -----
    app_mag (float): apparent magnitude of the source
    distance (float): distance of the source in Mpc

    Output:
    ------
    abs_mag (float): absolute magnitude
    """
    return (app_mag - 5*np.log10(distance)-25)

def redshift_to_distance(z):
    """
    Convert a given redshift to distance in Mpc, assuming H0=67.7 km/s/Mpc, and T_cmb=2.725 K
    Current cosmological constants used are similar to those on fritz.science

    Input:
    -----
    z (float): redshift of the source

    Output:
    ------
    distance (float): distance in Mpc
    """
    cosmo = FlatLambdaCDM(H0=67.66 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.30966, Neff=3.05, m_nu=[0., 0., 0.06]*u.eV, Ob0=0.049)
    return cosmo.luminosity_distance(z).value

def get_all_sources(program_id, saved_after=None):
    """
    Fetch all the obj_id and save dates of all sources for a given program_id.
    program_id (int): program_id of your program (i.e 43 for Caltech Census of the Local Univese)
    """
    all_sources = []
    num_per_page = 500
    page = 1
    total_matches = None
    retry_attempts = 0
    max_retry_attempts = 10
    
    if saved_after:
        date_query = Time(saved_after, format='isot')
        mainurl = BASEURL + f'api/sources?group_ids={program_id}&saveSummary=true&numPerPage={num_per_page}&savedAfter={date_query.isot}'
    else:
        mainurl = BASEURL + f'api/sources?group_ids={program_id}&saveSummary=true&numPerPage={num_per_page}'
    
    while retry_attempts <= max_retry_attempts:
        url = mainurl + f'&pageNumber={page}'
        response = api('GET',url)

        if response["status"] != "success":
            print(response)  # log as appropriate
            retry_attempts += 1
            time.sleep(1)
            continue
        if retry_attempts != 0:
            retry_attempts = 0

        all_sources.extend(response["data"]["sources"])
        total_matches = response['data']['total_matches']

        print(f'Fetched {len(all_sources)}/{total_matches} sources.')

        if len(all_sources) >= total_matches:
            break
            
        page += 1
    
    return all_sources

def get_source_api(ztfname, comments=False):
    ''' Info : Query a single source, takes input ZTF name
        comments (bool): If True it will also include the source comments
        Returns : all basic data of that source
    '''
    if comments:
        url = BASEURL + f'api/sources/{ztfname}?includeComments=true'
    else:
        url = BASEURL+f'api/sources/{ztfname}'

    response = api('GET', url)
    return response['data']

def get_user_id(last_name, group_id=43):
    url = BASEURL + 'api/groups/%s'%group_id
    response = api("GET", url)
    usr_id = []
    for user in response['data']['users']:
        if user['last_name']==last_name:
            usr_id.append(user['id'])
    if len(usr_id)>0:
        return (usr_id)
    else:
        return (False)

def get_source_photometry(ztfname, extrapolate_photometry=False):
    """ Fetch the photometry issued by ZTF for a given ZTF source.
    Params:
    -------
    ztfname (str): obj_id of the given ZTF source
    extrapolate_photometry (bool): Extrapolate the photometry from the .json file and convert it to astropy.Table format

    Output:
    ------
    astropy.Table containing mjd, filter, mag, mag_err (apparent magnitude), if extrapolate_photometry: True
    .json containing all magnitude information (including upper limits), if extrapolate_photometry: False
    """
    url = BASEURL + "api/sources/" + ztfname + "/photometry" # base-url to fetch data
    response = api("GET", url)
    
    return (response['data'])

def extrapolate_photometry(phot):
    """ Extrapolate .json photometry into astropy.Table format
    Params:
    ---
    phot (json): .json file with photometric information, obtained via get_source_photometry function
    
    Output:
    astropy.Table containing mjd, filter, mag, mag_err (apparent magnitude)
    """
    mjd = [alert['mjd'] for alert in phot if alert['mag']!=False]
    filt = [alert['filter'] for alert in phot if alert['mag']!=False]
    mag = [alert['mag'] for alert in phot if alert['mag']!=False]
    mag_err = [alert['magerr'] for alert in phot if alert['mag']!=False]

    if len(mag):
        return (Table([mjd, filt, mag, mag_err], names=("mjd", "filter", 'mag', 'mag_err')))
    else:
        return (None)

def CLU_luminous(photometry, source_z, nearby_bool, luminosity_threshold=-17.0):
    """
    This function will return a bool (True/False) according a to a luminosity threshold.
    The estimates of the luminosity take into account the Peak absolute magnitude (peak light) + error.

    Input
    -----
    photometry (astropy.Table): extrapolated photometry of the ZTF source
    soure_z (float): redshift of the ZTF source
    luminosity_threshold (float): luminosity threshold you want to test your peak absolute magnitude

    Output:
    luminous event (bool): True if it's more luminous than the luminosity threshold & False if fainter of the luminosity threshold
    """
    if photometry:
            if source_z:
                # Convert redshift to distance if source isn't super nearby
                if nearby_bool==True:
                    dist = 1. # 1 Mpc is slightly greater than the distance to M31, M33, etc; a source would need app. mag < 8 to pass -17 threshold here, and we'd definitely pick up any source that bright
                elif nearby_bool==False:
                    dist = redshift_to_distance(source_z)

                nan_index = photometry['mag'].data!=None # select index where no None's

                # Fetch absolute magnitude of the source
                Mabs = dist_mod_mag(photometry['mag'].data[nan_index], dist)

                # Lightcurve parameters
                filt = photometry['filter'].data[nan_index]
                Mapp_err = photometry['mag_err'].data[nan_index]
                #phase = photometry['mjd'].data[nan_index] - photometry['mjd'].data[nan_index][0]

                if len(filt):
                    # Select peak in lightcurve (first maximum)
                    peak = np.argmin(Mabs)

                    # Check if source is luminous by considering the absolute magnitude + error
                    is_luminous_error = (Mabs[peak] + Mapp_err[peak])<=luminosity_threshold

                    return (is_luminous_error, Mabs[peak])
                else:
                    return (False, None)
            else:
                return (None, None)
    else:
        return (None, None)

def get_source_phot_count(j):
    """Collects and counts the number of photometric points. Returns the total number in rgi
    Input:
    -----
    j (json): unextrapolated photometry of source
    """
    mt, mf = [], []
    for val in j:
        k0 = False
        for l in val['groups']:
            if l['name']=='Sitewide Group':
                k0 = True

        if k0==False:
            if val['mag']:
                mt.append(val['mag'])
                mf.append(val['filter'])

    mt = np.array(mt) # magnitude test
    mf = np.array(mf) # magnitude test filter

    return (len(np.where(mf=='ztfr')[0]), len(np.where(mf=='ztfg')[0]), len(np.where(mf=='ztfi')[0]))

def get_source_lastphot(j):
    """Retuns the latest photometric detection alongside the band it was detected at
    j (json): unextrapolated photometry of source"""
    mt, mf, mjde = [], [], []
    for val in j:
        k0 = False
        for l in val['groups']:
            if l['name']=='Sitewide Group':
                k0 = True

        if k0==False:
            if val['mag']:
                mt.append(val['mag'])
                mf.append(val['filter'])
                mjde.append(val['mjd'])

    mt = np.array(mt) # magnitude test
    mf = np.array(mf) # magnitude test filter
    mjde = np.array(mjde) # magnitude jde

    delta_time = Time(f"{today.year}-{today.month}-{today.day}T17:00").mjd - mjde # time since today in MJD note time in in UTC
    try:
        min_dt = np.argmin(delta_time)

        return (mt[min_dt], mf[min_dt], delta_time[min_dt])
    except:
        return (0, None, 0)

def get_source_spectra(ztfname):
    """Given ZTF name, load the .json containing **all** the spectra.
    Input
    ----
    ztfname (str): ZTF id

    Output
    ----
    Spectrum response in .json format
    """
    url = BASEURL + "api/sources/" + ztfname + "/spectra"
    response = api('GET',url)
    return (response)

def count_telescope(s, inst='DBSP'):
    """Count the number of unique spectra from i.e DBSP, SEDM, Keck
    Input:
    ------
    s (json): spectra of source from get_source_spectra

    Output:
    ------
    Number of unique of spectra for a given telescope ID
    """
    t0 = []
    for ss in s['data']['spectra']:
        if ss['instrument_name']==inst:
            t0.append(ss['created_at'].split("T")[0])
    t0 = np.array(t0)
    return (len(np.unique(t0)))

def count_fritz_spec(spec):
    """ Count the number of spectra uploaded to Fritz, returning the number of spectra found.

    Input
    -----
    spec (json): spectra from get_source_spectra
    """
    spec_count = 0
    if spec['status']=="error":
        return (spec_count)
    elif len(spec['data']['spectra'])==0:
        return (spec_count)
    else:
        spectrum_data = spec['data']['spectra'] # load all spectra

        for ss in enumerate(spec['data']['spectra']):
            indx, _ = ss[0], ss[1] # counting index, spectrum dict.
            spec_count = indx+1 # +1 since index starts at 0

        return (spec_count)

def fetch_likely_CLU_z(source):
    """Fetch the closest CLU redshift and distance from it in arcsec
    source is the source API data from get_source_api
    """

    if source['redshift']:
        z0 = source['redshift'] # True redshift

        if len(source['annotations'])>=1: # check that annotations exist
            # Now search through all crossmatches and find the closest match in redshift and fetch it's seperation in arcsec

            # Now check if CLU annotations exist
            for s in source['annotations']: # Looking through annotations
                if s['origin']=='clu:clu': # found CLU xmatch
                    z_list, z_sep_list = [], []
                    try:
                        for val in zip(s['data']['CLU_z'], s['data']['CLU_distance_arcsec']):
                            # append those redshift and seperations
                            z_list.append(val[0])
                            z_sep_list.append(val[1])

                        # Convert those redshifts and seperations into arrays
                        z_list, z_sep_list = np.array(z_list), np.array(z_sep_list)

                        # estimate the differences in redshift from that of z0
                        delta_z = abs(z0 - z_list)

                        # The smallest difference will likely be the galaxy and thus take that seperation
                        min_to_max_z = np.argsort(delta_z)

                        z_list_sorted, z_sep_sorted = z_list[min_to_max_z], z_sep_list[min_to_max_z]

                        # the first index is likely the cloest match based on the redshift
                        return (z0, z_sep_sorted[0], False)
                    except:
                        continue

                else:
                    continue

            else: # no CLU annotation exists
                #return (z0, False)
                return (z0, None, False)
        else:
            # No annotations in this source... just return redshift
            #return (z0, False)
            return (z0, None, False)

    # There is no redshift entry on Fritz -- now checking manually through annotations
    else:
        # Check if there are any annotations
        if source['annotations']:
            for s in source['annotations']:
                if s['origin']=='clu:clu': # found CLU xmatch
                        z_list, z_sep_list = [], []
                        try:
                            for val in zip(s['data']['CLU_z'], s['data']['CLU_distance_arcsec']):
                                # append those redshift and seperations
                                z_list.append(val[0])
                                z_sep_list.append(val[1])

                            # Convert those redshifts and seperations into arrays
                            z_list, z_sep_list = np.array(z_list), np.array(z_sep_list)
                            
                            # sort by the smallest seperation
                            min_to_max_z = np.argsort(z_sep_list)

                            z_list_sorted, z_sep_sorted = z_list[min_to_max_z], z_sep_list[min_to_max_z]

                            # If there's a negative redshift, likely in a very nearby galaxy (eg. M31, M33)
                            if np.any((z_list<0) & (z_list>-999)):
                                return (z_list_sorted[0], z_sep_sorted[0], True) #if nearby, don't want to exclude based on abs mag but further calculations may be inaccurate - raise flag instead
                            else:
                                # the first index is likely the cloest match based on the redshift
                                return (z_list_sorted[0], z_sep_sorted[0], False)
                        except:
                            #return (False, False)
                            return (None, None, False)

            else: # no annotations exist
                #return (False, False)
                return (None, None, False)

def peak_mag_CLU(photometry):
    """
    Return the peak apparent magnitude given the ztf name. This will return the first maximum detection.
    Input
    -----
    photometry (astropy.Table): extrapolated photometry of a given source

    Output
    ------
    peak_app_mag (float): peak apparent magnitude
    """
    if photometry: # if photometry table is not empty
        try:
            nan_index = photometry['mag'].data!=None # select index where no None's

            # Lightcurve parameters
            filt = photometry['filter'].data[nan_index]
            Mapp = photometry['mag'].data[nan_index]
            Mapp_err = photometry['mag_err'].data[nan_index]
            phase = photometry['mjd'].data[nan_index] - photometry['mjd'].data[nan_index][0]

            # Try to get the peak magnitude in r-band, if r-band does not exist, any other photometry is okay
            contain_red, contain_green, contain_i = False, False, False
            for f in filt:
                if f=="ztfr": # rband
                    contain_red = True
                elif f=='ztfg':
                    contain_green = True
                elif f=='ztfi':
                    contain_i = True

            if contain_red==True:
                Mapp_phot_band = Mapp[filt=='ztfr']
                filt_band = 'ztfr'
            elif contain_red==False and contain_green==True:
                Mapp_phot_band = Mapp[filt=='ztfg']
                filt_band = 'ztfg'
            elif contain_red==False and contain_green==False and contain_i==True:
                Mapp_phot_band = Mapp[filt=='ztfi']
                filt_band = 'ztfi'

            peak_app = np.argmin(Mapp_phot_band) # peak in apparent magnitude space

            return (Mapp_phot_band[peak_app], filt_band)
        except:
            return (None, None)
    else:
        return (None, None)

def CLU_distances(z, distance_cut=150.):
    """
    Returns a flag according to whether a source is within a distance limit and the distance in Mpc.
    Input:
    ---
    z (float): redshift of given ZTF source
    distance_cut (float): distance limit in Mpc
    
    Output:
    ---
    True if object distance < distance_cut
    False if object distance > distance_cut
    NaN if redshift is anomalous and could require further checking (eg. z=999)
    None if redshift (and thus, luminosity distance) is not available
    """
    if (z < 100.) and (z>-100.):
        dist = redshift_to_distance(z)

        if (dist <= distance_cut) or (z<=0.0331): # turns out z=0.0331 is like, 150.5 Mpc...
            return True, dist
        elif dist > distance_cut:
            return False, np.abs(dist)
        
    else:
        return None, None
    
def query_CLU_dict(ZTF_name, saved_at):
    """
    This function will append the parameters of a given ZTF souce to the
    global parameters (see query output).

    Query Input:
    -----------
    ZTF_name (str): obj_id of the given ZTF source
    saved_at (str): saved date of the given source

    Query Output:
    -------------
    obj_id (str): obj_id of a given ZTF source
    saved_at (str): date saved on fritz.scince
    ra, dec (float): RA, DEC coordinates of the source (mean position)
    z (float): redshift of CLU host galaxy
    clu_d (float): sky seperation of the transient with respect to the CLU host (arcsec)
    classification (str): the most recent classification on fritz.science
    classification_prob(float): probability of the most recent classification
    peak_app_mag (float): peak apparent magnitude of the lightcurve
    peak_abs_mag (float): peak absolute magnitude of the ligtcurve given the most recent redshift estimate
    luminous_event (str): if source is luminous or subluminous given a -17 luminosity cut
    """
    # Fetch source basic info
    source_api_data = get_source_api(ZTF_name, comments=False)

    # Fetch most recent classification
    clf = source_api_data['classifications']
    if clf:
        param_classification = clf[-1]['classification']
        param_prob = clf[-1]['probability']
    else:
        #param_classification = None
        param_classification = ''
        param_prob = None
    
    # fetch ra/dec
    param_ra = source_api_data['ra']
    param_dec = source_api_data['dec']
    
    # Fetch redshift and seperation from clu_galaxy
    redshift_and_sep_info = fetch_likely_CLU_z(source_api_data)
    param_z  = redshift_and_sep_info[0] # redshift
    param_clu_d = redshift_and_sep_info[1] # seperation from CLU host
    param_nearby = redshift_and_sep_info[2]
    
    # physical distance to each CLU transient: Delta_sep * (z/z_max)
    if param_clu_d and param_z:
        d0_kpc = param_clu_d  * (abs(param_z) / 0.05)
    else:
        #d0_kpc = False
        d0_kpc = None

    # photometry params: first, get photometry
    source_photometry = get_source_photometry(ZTF_name)
    source_photometry_ext = extrapolate_photometry(source_photometry)
    
    if param_z: # if param_z is not None
        luminosity_info = CLU_luminous(source_photometry_ext, param_z, param_nearby, luminosity_threshold=-17.0) # fetch luminosity boolean
        param_luminous_event = luminosity_info[0] # luminosity boolean
        param_peak_abs_mag = luminosity_info[1] # peak absolute magnitude
        param_within_vol, dist_Mpc = CLU_distances(param_z, distance_cut=150.)
    else: # if no redshift
        #param_luminous_event = False
        #param_peak_abs_mag = -99.0
        param_luminous_event = None
        param_peak_abs_mag = None
        param_within_vol = None
        dist_Mpc = None

    # one parameter to summarize if we should check a source for follow-up (in volume + subluminous)
    param_in_clu = (param_within_vol==True) and (param_luminous_event==False)
    
    # Fetch peak apparent magnitudes
    peak_apparent_magnitude, peak_filt = peak_mag_CLU(source_photometry_ext)
    if peak_apparent_magnitude!=None:
        param_peak_app_mag = peak_apparent_magnitude
    else:
        param_peak_app_mag = None

    N_photometry_r, N_photometry_g, N_photometry_i = get_source_phot_count(source_photometry) # number of rgi band det.
    last_mag, last_filt, last_det = get_source_lastphot(source_photometry) # last detection in days
    
    # get spectra info
    source_spec = get_source_spectra(ZTF_name)
    N_spec = count_fritz_spec(source_spec)
    N_SEDM, N_DBSP, N_LRIS = count_telescope(source_spec, inst='SEDM'), count_telescope(source_spec, inst='DBSP'), count_telescope(source_spec, inst='LRIS')

    # Extrapolate to main dict.
    db['source'].append({"ZTF_id": ZTF_name, "saved_date": saved_at,
                         "ra":param_ra, "dec":param_dec,
                         "z": param_z, "clu_d_host": param_clu_d, "nearby_flag": str(param_nearby),
                         "classification": param_classification, "classification_prob": param_prob,
                         "peak_abs_mag": param_peak_abs_mag, "in_clu":str(param_in_clu),
                         "luminous_event": str(param_luminous_event), "in_clu_vol": str(param_within_vol),
                         "peak_app_mag": param_peak_app_mag, 'peak_mag_filt':peak_filt, "clu_d_kpc":d0_kpc, "host_d_Mpc":dist_Mpc,
                         'N_photometry_r':N_photometry_r, 'N_photometry_g':N_photometry_g, 'N_photometry_i':N_photometry_i,
                         'last_mag': last_mag, 'last_filt':str(last_filt), 'last_det':last_det, 
                         'N_spec':N_spec, 'N_SEDM':N_SEDM, 'N_DBSP':N_DBSP, 'N_LRIS':N_LRIS})
    
def QUERYME(date_query='2021-01-01', existing_query_file=None, force_refresh=None):
    if date_query!='all':
        sources = get_all_sources(43, saved_after=date_query) # fetch all sources for ZTF CLU experiment
    else:
        sources = get_all_sources(43)

    # All souces (names & dates) of all saved ZTF CLU candidates
    names = np.array([s['obj_id'] for s in sources])
    dates = np.array([s['saved_at'] for s in sources])
    
    # start with existing file if specified
    if existing_query_file:
        existing_data = ascii.read(f"query_data/{existing_query_file}")
        
        # identify all candidates whose info needs updating and delete existing rows to be replaced
        if force_refresh:
            refresh_date = Time(force_refresh, format='isot')
            refresh = np.where(dates>=refresh_date)[0]
            dates = dates[refresh]
            names = names[refresh]
            
            existing_data.remove_rows(refresh[refresh<len(existing_data)])
        
        # just identify candidates not in the current catalog
        else:
            new_entries = np.empty(len(names), dtype=bool)
            for i, name in enumerate(names):
                if name in existing_data['ZTF_id']:
                    new_entries[i] = False
                else:
                    new_entries[i] = True

            names = names[new_entries]
            dates = dates[new_entries]
            
    list_threads = [] # list of threads for pooling
    buff = 1 # buffer time (NOTE: API sometimes crashes if request is to frequent) - start at 1 to account for query to load all sources in group
    for i in tqdm(range(len(names))):
        buff +=1
        if buff>=2: # feed to pool 2 sources a time - each source calls the API 3x and the rate limit is 5 calls per second
            time.sleep(2) # sleep for a second just to be safe
            buff = 0 # restart buffer

        t = threading.Thread(target=query_CLU_dict, args=[names[i], dates[i]])
        t.start()
        list_threads.append(t)

    for thread in list_threads:
        thread.join()
        
    # if started with existing file, append new/updated entries, sort by date, and save
    if existing_query_file:
        for i in range(len(db['source'])):
            source_dict = db['source'][i]
            row_mask = dict.fromkeys(source_dict)

            for key, val in source_dict.items():
                if val==None or val=='' or val=='None':
                    row_mask[key] = True
                else:
                    row_mask[key] = False
            
            existing_data.add_row(source_dict, mask=row_mask)
            
        existing_data.sort('saved_date')
        existing_data.write(f"query_data/CLU_Query_{today.year}_{today.month}_{today.day}.ascii", format='ascii')
        
    # if creating new query file, extrapolate dict to astropy Table, sort by date, and save
    else:
        param_names = ("ZTF_id", "saved_date", "ra", "dec", "z", 
                       "clu_d_host", "nearby_flag",
                       "classification", "classification_prob", 
                       "peak_app_mag", "peak_mag_filt", "peak_abs_mag", "in_clu",
                       "luminous_event", "in_clu_vol", "clu_d_kpc", "host_d_Mpc",
                       "N_photometry_r", "N_photometry_g", "N_photometry_i",
                       'last_mag', 'last_filt', 'last_det', 
                       'N_spec', "N_SEDM", "N_DBSP", "N_LRIS")
        
        cols = []
        for i, param_name in enumerate(param_names):
            param_vals = np.asarray([s[param_name] for s in db['source']])
            none_vals = np.isin(param_vals, [None, '', 'None']) # test if any of the None/no value flags exist, return True at index if it does
            
            # if there's a bad value, make a mask and save as astropy MaskedColumn
            if np.any(none_vals):
                col = MaskedColumn(param_vals, name=param_name, mask=none_vals)
            # if all the values are fine, save as regular astropy Column
            else:
                col = Column(param_vals, name=param_name)
                
            cols.append(col)
            
        data_table = Table(cols)
        data_table.sort('saved_date')
        data_table.write(f"query_data/CLU_Query_{today.year}_{today.month}_{today.day}.ascii", format='ascii')