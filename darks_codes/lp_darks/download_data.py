"""Downloads and organizes WFC3/UVIS raw dark or science data to be reduced.

Code to download darks or science data using astroquery instead of Quicklook.
To be used in the context of the STScI darks pipeline edited 
by Laura Prichard and implemented as outlined in darks_reduction_lp.ipynb
that this rountine should be with in lp_darks/.

For darks: 

    The code takes anneal cycle start and end dates in MJD for one cycle at a 
    time as inputs. For each anneal cycle, a directory is made that has folders 
    for the downloaded data, the raw dark copies to be processed and the CTE 
    correctd darks. Dark files are then selected using astroquery and there 
    is the option to download them (in case just the file organization is needed).
    Finally, downloaded data are copied over, from the nested directories done 
    as standard by astroquery, to a single directory ready for CTE correction. 

    Execution of this script will create the following file tree:
    ::
    <anneal_date:YYYYMMMDD>anneal_rawdarks_aquery<download_date:YYYYMMMDD>/
        ctecorr_darks/
        mastDownload/HST/  (created by the astroquery download command)    
            id**/          (creates ID specific folders with the files in)
        raw_darks/

For science data: 

    The code takes an HST/WFC3 program ID as an input and downloads 
    the raw science data for that program. For each proposal ID, a directory is 
    made and sub-directories are made for the downloaded data, the raw data 
    copies to be processed and the CTE correctd raw science data. Science data 
    files are then selected using astroquery and there is the option to download 
    them (in case just the file organization is needed).
    Finally, downloaded data are copied over, from the nested directories done 
    as standard by astroquery, to a single directory ready for CTE correction. 

    Execution of this script will create the following file tree:
    ::
    PID<proposal_id>_rawdata_aquery<download_date:YYYYMMMDD>/
        ctecorr_sci/       (for CTE corrected data)
        mastDownload/HST/  (created by the astroquery download command)    
            id**/          (creates ID specific folders with the files in)
        raw_sci/           (for raw data copied from the download directory into a single directory)
        calwf3_sci/        (for processing science data with new darks with calwf3)


Authors
-------
    Laura Prichard, 2019

Use
---

    The routine should be run for ONE anneal cycle at a time and can be 
    run from the command line using:
    ::

        python download_data.py [-t|--type] [-s|--data_start] [-e|--data_end] 
        [-p|--proposal_id] [-r|--raw_dir] [-d|--download] [-l|--dload_date]

    ``-t --type`` (Required) -- the type of data to be downloaded, either ``dark`` 
        or raw darks or ``science`` for raw science data. ``-s --anneal_start`` 
        and ``-e --anneal_end`` must be set if --type==``dark``. ``-p --proposal_id`` 
        should be set if --type==``science``.

    ``-s --data_start`` (Required if --type='dark') -- anneal start date of ONE 
        cycle in MJD listed to at least 6 d.p. as a float, available from 
        http://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/performance/monitoring/complete-anneal-history/_documents/uvis-anneal.txt. Optional if type=``science``, science data start date if a section of a proposal`s data is to be downloaded rather than the whole program.

    ``-e --data_end`` (Required if --type='dark') --  anneal end date of ONE 
        cycle in MJD listed to at least 6 d.p. as a float, available from 
        http://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/performance/monitoring/complete-anneal-history/_documents/uvis-anneal.txt. Optional if type=``science``, science data end date if a section of a proposal`s data is to be downloaded rather than the whole program. All files up to but NOT INCLUDING the end date are selected.

    ``-p --proposal_id`` (Required if --type='science') -- Proposal ID of raw 
        science data to be downloaded. Can be a string or float, must be set 
        if --type==``science``.

    ``-r --raw_dir`` (Required) -- Path to a raw file directory where the 
        directory for each anneal cycle should be made. Must be a string and have 
        no trailing "/".

    ``-d --download`` (Optional) -- Flag. The value is ``True`` (i.e.
        download data) if provided and ``False`` (i.e. don't download) 
        if not provided. The code does file organization in addition 
        to the download so the download can be switched off if needed.

    ``-l --dload_date`` (Required if ``-d --download`` not set) -- String 
        of the date the data was downloaded if just file organization is 
        required. This is to save re-downloadig the data which was timestamped.
        Should be in the format ``YYYYMonDD`` e.g. ``2019Aug17``.

"""


from astropy.time import Time
import argparse
import os
from pdb import set_trace as st
from astroquery.mast import Observations
import glob
import shutil


def copy_data(paths):
    """Copies the downloaded files out of the nested astroquery directories 
    and into one combined directory (raw downloaded data, RWD_DIR). Here, 
    the CTE correction code copies over only the darks with the right exposure 
    times to be CTE corrected. Or calwf3 is run on all the raw science data files.
    """

    # Move to the download directory
    os.chdir(paths['DLD_DIR'])

    # Copies data out of the astroquery sub-directories and into the RWD_DIR directory if empty
    if glob.glob(os.path.join(paths['RWD_DIR'], '*.fits')):
        print('Directory not empty!! Not copying files')
    else:
        print('Copying downloaded raws from astroquery subdirectories {}'.format(paths['DLD_DIR']))
        n=0
        for i, idob in enumerate(glob.glob('*')):
            f = os.path.join(paths['DLD_DIR'], idob, idob +'_raw.fits')
            shutil.copy(f, paths['RWD_DIR'])
            n+=1
        print('Copied {} files to combined raw data directory {}'.format(n, paths['RWD_DIR']))


def astroquery_darks(anneal_start, anneal_end, paths, download):
    """Retreives darks within an anneal cycle using astroquery

    To get the right files, astroquery requires the following inputs:

        - intentType='calibration' – for the darks
        - instrument_name="WFC3/UVIS" – or whichever instrument is needed  
        - t_min – Anneal dates start and end defined above
        - target_name=DARK* – all files identified as a type of dark

    NO exposure time should be given (too stringent a cut). The astroquery 
    command returns more than the required files than would be collected 
    from MAST. However, the darks CTE correction code then selects only the 
    files it needs from these, so this is not a problem.

    Astroquery takes the start and end dates of the anneal cycle. The end 
    anneal date in the astroquery command below is slightly reduced so that 
    the file at the start of the next anneal is not included. 
    e.g. [58170.7838 (start of Feb 21 2018 anneal), 58201.3608 (start of next Mar 24 2018 anneal] 
    --> [58170.7838, 58201.360]"""

    # Round the next anneal start date down so that the file from the new date isn't included
    anneal_end_rnd = anneal_end // 0.001 * 0.001

    # Collecting the raw dark files for one anneal cycle using astroquery
    darkobs = Observations.query_criteria(intentType='calibration',instrument_name="WFC3/UVIS", t_min=[anneal_start,anneal_end_rnd], target_name="DARK*")
    darkprod = Observations.get_product_list(darkobs)
    rawdark = Observations.filter_products(darkprod, productSubGroupDescription="RAW")

    # If download option is set, move to the anneal directory, check if empty, then download
    if download==True:
        os.chdir(paths['ANN_DIR'])
        if os.path.exists(paths['DLD_DIR']):
            print('Download directory exists!! Not downloading files')
        else: 
            Observations.download_products(rawdark, mrp_only=False)  
            print('Download to {} complete'.format(paths['DLD_DIR']))


def astroquery_sci(data_start, data_end, proposal_id, paths, download):
    """Retreives raw science data for a given proposal ID using astroquery

    To get the right files, astroquery requires the following inputs:

        - intentType='science' – for the darks
        - instrument_name="WFC3/UVIS" – or whichever instrument is needed  
        - proposal_id -- user-defined input string of program/proposal ID
    """

    # Collecting the raw science files for a program ID using astroquery
    if data_start:
        # Round the data end date down so that the file from the new date isn't included
        data_end_rnd = data_end // 0.001 * 0.001
        sciobs = Observations.query_criteria(intentType='science', t_min=[data_start,data_end_rnd], instrument_name="WFC3/UVIS", proposal_id=proposal_id)
    else:
        sciobs = Observations.query_criteria(intentType='science', instrument_name="WFC3/UVIS", proposal_id=proposal_id)

    sciprod = Observations.get_product_list(sciobs)
    rawsci = Observations.filter_products(sciprod, productSubGroupDescription="RAW")

    # If download option is set, move to the download directory, check if empty, then download
    if download==True:
        os.chdir(paths['PID_DIR'])
        if os.path.exists(paths['DLD_DIR']):
            print('Download directory exists!! Not downloading files')
        else: 
            Observations.download_products(rawsci, mrp_only=False)  
            print('Download to {} complete'.format(paths['DLD_DIR']))


def make_dirs(dtype, raw_dir, download, dload_date, anneal_name, pid_name):
    """Defining the path names for the downloaded raw, raw copies, 
    and CTE corrected data
    
    Returns
    -------
    paths : dict
        A dictonary whose keys are path identifiers and whose
        values are strings containing absolute paths.
    """

    paths = {}

    # If no download date is provided, assume the data is being downloaded and set date to today
    if dload_date==None: 
        # Setting today's/download date for directory namess
        now = Time.now()
        dload_date = now.strftime('%Y%b%d')
    print('Download date: ', dload_date)
    
    # Setting the names of the raw, anneal, raw copied data, CTE corrected, and download directories to paths
    paths['RAW_DIR'] = raw_dir
    if dtype=='dark':
        paths['ANN_DIR'] = os.path.join(paths['RAW_DIR'], anneal_name + 'anneal_rawdarks_aquery' + dload_date)
        paths['RWD_DIR']  = os.path.join(paths['ANN_DIR'], 'raw_darks')
        paths['CTE_CORR_DIR'] = os.path.join(paths['ANN_DIR'], 'ctecorr_darks')
        paths['DLD_DIR'] = os.path.join(paths['ANN_DIR'], 'mastDownload', 'HST')  #Naming convention from MAST, don't need to make DLD_DIR
    if dtype=='science':
        paths['PID_DIR'] = os.path.join(paths['RAW_DIR'], pid_name + '_rawdata_aquery' + dload_date)
        paths['RWD_DIR']  = os.path.join(paths['PID_DIR'], 'raw_sci')
        paths['CTE_CORR_DIR'] = os.path.join(paths['PID_DIR'], 'ctecorr_sci')
        paths['DLD_DIR'] = os.path.join(paths['PID_DIR'], 'mastDownload', 'HST')  #Naming convention from MAST, don't need to make DLD_DIR
        paths['CW3_DIR'] = os.path.join(paths['PID_DIR'], 'calwf3_sci')
    
    # Making RAW_DIR if it doesn't exist
    if not os.path.exists(paths['RAW_DIR']): 
            os.makedirs(paths['RAW_DIR'], 0o774)
            print('Created raw data directory: ', paths['RAW_DIR'])
    
    if dtype=='dark':        
        # Making the subdirectories if they don't exist
        if os.path.exists(paths['ANN_DIR']):   
            print('Anneal directory already exists!!! Not creating sub-directories.')
        else:
            os.mkdir(paths['ANN_DIR'], 0o774)
            print('Created anneal cycle directory: {}'.format(paths['ANN_DIR']))

            os.mkdir(paths['RWD_DIR'], 0o774)
            print('Created raw darks directory for copied raw files: {}'.format(paths['RWD_DIR']))

            os.mkdir(paths['CTE_CORR_DIR'], 0o774)
            print('Created CTE corrected darks directory: {}'.format(paths['CTE_CORR_DIR']))

    if dtype=='science':
        # Making the subdirectories if they don't exist
        if os.path.exists(paths['PID_DIR']):   
            print('Proposal ID directory already exists!!! Not creating sub-directories.')
        else:
            os.mkdir(paths['PID_DIR'], 0o774)
            print('Created proposal ID directory: {}'.format(paths['PID_DIR']))

            os.mkdir(paths['RWD_DIR'], 0o774)
            print('Created raw science data directory for copied raw files: {}'.format(paths['RWD_DIR']))

            os.mkdir(paths['CTE_CORR_DIR'], 0o774)
            print('Created CTE corrected science data directory: {}'.format(paths['CTE_CORR_DIR']))

            os.mkdir(paths['CW3_DIR'], 0o774)
            print('Created directory for processing science data with calwf3.e software: {}'.format(paths['CW3_DIR']))

    return paths


def get_anneal_name(anneal_start):

    t = Time(anneal_start, format='mjd')
    anneal_name = t.strftime('%Y%b%d')
    print('Anneal start date: ', anneal_name)

    return anneal_name


def parse_args():

    type_help='Required, the type of data that is being downloaded. String, options ``dark`` or ``science``.'
    data_start_help = 'Required if type=``dark``, anneal start date of ONE cycle. Optional if type=``science``, science data start date if a section of a proposal`s data is to be downloaded rather than the whole program. In MJD listed to at least 6 d.p. String or float.'
    data_end_help = 'Required if type=``dark``, anneal end date of ONE cycle. Optional if type=``science``, science data end date if a section of a proposal`s data is to be downloaded rather than the whole program. In MJD listed to at least 6 d.p. String or float.'
    pid_help = 'Required if type=``science``, proposal ID of raw science data to be downloaded. String.'
    raw_dir_help = 'Path to a raw file directory where the directory for each anneal cycle should be made. Must be a string and have no trailing "/".'
    download_help = 'Optional flag, the data will be downloaded by astroquery if set and will not be if not set.'
    dload_date_help = 'String of the date the data was downloaded if just file organization is required. Should be in the format ``YYYYMonDD`` e.g. ``2019Aug17``.'
    
    parser = argparse.ArgumentParser()
    
    # Argument to output CTE corrected raw directory
    parser.add_argument('-t --type',
        dest='type',
        action='store',
        type=str,
        required=True,   
        help=type_help) 

    # Argument to output CTE corrected raw directory
    parser.add_argument('-s --data_start',
        dest='data_start',
        action='store',
        type=float,
        required=False,   
        help=data_start_help) 

    # Argument to input downloaded raw data directory
    parser.add_argument('-e --data_end',
        dest='data_end',
        action='store',
        type=float,
        required=False,   
        help=data_end_help) 

    # Argument to input downloaded raw data directory
    parser.add_argument('-p --proposal_id',
        dest='proposal_id',
        action='store',
        type=str,
        required=False,   
        help=pid_help) 

    # Argument to input downloaded raw data directory
    parser.add_argument('-r --raw_dir',
        dest='raw_dir',
        action='store',
        type=str,
        required=True,   
        help=raw_dir_help) 

    # Argument to software directory
    parser.add_argument('-d --download',
        dest='download',
        action='store_true',    #Set's to True if called otherwise, False as set below
        help=download_help)

    # Argument to input downloaded raw data directory
    parser.add_argument('-l --dload_date',
        dest='dload_date',
        action='store',
        type=str,
        required=False,   
        help=dload_date_help) 

    parser.set_defaults(download=False)

    # Parse args
    args = parser.parse_args()

    return args


def test_args(args):
    """Ensures that the command line arguments are of proper format. If
    they are not, an assertion error is raised.

    Parameters
    ----------
    args : obj
        The ``argparse`` object containing the command line arguments.
    """
    assert args.type is str, 'Data type eith ``dark`` or ``science`` to be downloaded, must be a string.'
    assert args.data_start is float, 'Invalid date format. Format must be a float MJD value.'
    assert args.data_end is float, 'Invalid date format. Format must be a float MJD value.'
    assert args.proposal_id is str, 'Proposal ID of raw data set to be downloaded must be a string.'
    assert args.raw_dir is str, 'Path to the raw data directory (raw_dir) must be a string.'
    assert args.download in [True, False], 'Invalid option for download argument. Valid options are True and False.'
    assert args.dload_date is str, 'Invalid option for dload_date argument. Should be a string in the format ``YYYYMonDD`` e.g. ``2019Aug17``.'
    if (args.download is False) and (dload_date is None):
        parser = argparse.ArgumentParser()
        parser.error("If --download not set, provide --dload_date.")

if __name__ == '__main__':

    # Read in and test arguments
    args = parse_args()
    # test_args(args)
    anneal_name=''
    pid_name=''

    if args.type=='dark':
        # Define the anneal name from the MJD input
        anneal_name = get_anneal_name(args.data_start)

        # Make a dictionary of paths and make if they don't already exist
        paths = make_dirs(args.type, args.raw_dir, args.download, args.dload_date, anneal_name, pid_name)

        # Retreive and download darks (if download is set) using astroquery
        astroquery_darks(args.data_start, args.data_end, paths, args.download)
    
    elif args.type=='science':
        # Define the program ID name
        pid_name = 'PID' + str(args.proposal_id)

        # Make a dictionary of paths and make if they don't already exist
        paths = make_dirs(args.type, args.raw_dir, args.download, args.dload_date, anneal_name, pid_name)

        # Retreive and download science data (if download is set) using astroquery
        astroquery_sci(args.data_start, args.data_end, args.proposal_id, paths, args.download)

    # Copy darks out of the astroquery nested download 
    # directories into a single directory ready for processing
    copy_data(paths)
