#! /usr/bin/env python

"""Retrieves WFC3/UVIS darks and performs the CTE correction.

This script retrieves WFC3/UVIS science data from a given directory, copies 
them into a specified cte correction directory, and uses the WFC3 CTE
correction code to CTE correct any images in that directory that do not
already have an accompanying ``RAC`` file.  ``RAC`` files are stored in
the same directory.

Authors
-------

    Edited by Laura Prichard, September 2019
    Matthew Bourque, June 2014

Use
---

    This code is intended to be run from the command line using the following:
    ::

        python cal_uvis_make_ctecorr_darks.py [-c|--ctecorr_dir] 
            [-r|--rwd_dir] [-s|--software_dir]

    ``-c --ctecorr_dir`` (Required) -- path to output CTE corrected data directory. 
        String, no trailing slash "/".

    ``-r --rwd_dir`` (Required) -- path to raw data directory (all files in one 
        folder as done in previous step). String, no trailing slash "/".

    ``-s --software_dir`` (Required) -- path to the compiled CTE correction code 
        ./wfc3uv_ctereverse.e. String, no trailing slash "/".

"""

import glob
import logging
import multiprocessing
import os
import shutil
import subprocess
import socket    #LP added
import datetime  #LP added
from pdb import set_trace as st  #LP added
import argparse  #LP added

import astropy.io.fits as fits

# LP comment
# from pyql.logging.logging_functions import configure_logging   #LP commented, for Quicklook
# from pyql.logging.logging_functions import log_info
# from pyql.logging.logging_functions import log_fail


# LP commented, alternative definition below in copy_dark_lp
# def copy_dark(src, dst):
#     """Copies a file from source to destination if it meets the header
#     value criteria of ``TARGNAME = 'DARK'`` or ``'DARK-NM'``, and
#     ``EXPTIME = 900``.

#     Parameters
#     ----------
#     src : str
#         The absolute path of the source location.
#     dst : str
#         The absolute path of the destination location.
#     """

#     if 'raw' in src:

#         # Get header information
#         exptime = fits.getval(src, 'EXPTIME')
#         targname = fits.getval(src, 'TARGNAME')

#         # Only copy appropriate files
#         if targname == 'DARK' or targname == 'DARK-NM':
#             if exptime in [900.0, 1800.0]:
#                 if not os.path.exists(dst):
#                     logging.info('\tCopying {} to {}'.format(src, dst))
#                     shutil.copyfile(src, dst)


def copy_dark_lp(src, dst):
    """Copies a file from source to destination if it meets the header
    value criteria of ``TARGNAME = 'DARK'`` or ``'DARK-NM'``, and
    ``EXPTIME = 900`` or ``EXPTIME = 1800``.

    Parameters
    ----------
    src : str
        The absolute path of the source location.
    dst : str
        The absolute path of the destination location.

    """

    # Check for all files with raw in the name
    for filepath in glob.glob(os.path.join(src,'*')):

        # Getting just the filename
        file = os.path.basename(filepath)

        # If that file is raw, then check header
        if 'raw' in file:

            # Get header information
            exptime = fits.getval(filepath, 'EXPTIME')
            targname = fits.getval(filepath, 'TARGNAME')

            # Only copy appropriate files
            if targname == 'DARK' or targname == 'DARK-NM':
                if exptime in [900.0, 1800.0]:
                    dstfile = os.path.join(dst, file)
                    if not os.path.exists(dstfile):
                        # print('Copying {} to {}'.format(filepath, dst+file))
                        print('\tLP log: Copying {} to {}'.format(filepath, dstfile)) 
                        shutil.copyfile(filepath, dstfile)


def get_files_to_correct(CTE_CORR_DIR):   #LP added CTE_CORR_DIR
    """Generates a list of RAW files that need to be CTE corrected.

    Determines the RAW files that need to be CTE corrected by checking
    if there is an accompanying ``RAC`` file in the same directory.  If
    there is no accompanying ``RAC`` file, the basename of the ``RAW``
    file path is added to the list of files to be CTE corrected.

    Parameters - LP added
    ----------
    CTE_CORR_DIR : str
        The absolute path of the CTE corrected data.

    Returns
    -------
    files_to_correct : list
        The list of RAW files (basenames) that require CTE correction.
    """

    files_to_correct = []

    raw_files = glob.glob(os.path.join(CTE_CORR_DIR, '*raw.fits'))  #LP removed the /
    rac_files = glob.glob(os.path.join(CTE_CORR_DIR, '*rac.fits'))  #LP removed the /

    raw_basenames = [os.path.basename(raw).split('_')[0] for raw in raw_files]
    rac_basenames = [os.path.basename(rac).split('_')[0] for rac in rac_files]
 
    for raw_basename in raw_basenames:
        if raw_basename not in rac_basenames:
            files_to_correct.append('{}_raw.fits'.format(raw_basename))

    return files_to_correct


def perform_correction(filelist, CTE_CORR_DIR, SOFTWARE_DIR):   #LP added CTE_CORR_DIR, SOFTWARE_DIR
    """Executes the CTE correction code on each file in filelist

    Before the CTE correction code is executed, the current working
    directory is switched to the CTE correction directory as to avoid
    argument errors in the code.

    Parameters
    ----------
    filelist : list
        The list of RAW files (basenames) that require CTE correction.

    LP added:
    CTE_CORR_DIR : str
        The absolute path of the CTE corrected data.

    SOFTWARE_DIR : str
        The absolute path of the wfc3uv_ctereverse.e executable.
    """

    job_list = []

    print('')

    if len(filelist) > 0:
        for rawfile in filelist:
            # LP comment
            # job_list.append('cd {}; ./wfc3uv_ctereverse.e {}'.format(
            #     CTE_CORR_DIR, rawfile))
            # LP edit, included SOFTWARE_DIR as input
            job_list.append('cd {}; {}/./wfc3uv_ctereverse.e {}'.format(
                CTE_CORR_DIR, SOFTWARE_DIR, rawfile))
        # print('Processing {} jobs'.format(str(len(job_list))))
        print('\tLP log: Processing {} jobs'.format(str(len(job_list))))
        pool = multiprocessing.Pool(processes=4)
        pool.map(run_process, job_list)
        pool.close()
        pool.join()
    else:
        # print('No files to correct')
        print('\tLP log: No files to correct')


def run_process(cmd):
    """Calls ``subprocess`` with the command ``cmd``.

    Parameters
    ----------
    cmd : str
        The subprocess command to execute.

    Returns
    -------
    subprocess.call() : obj
        The ``subprocess`` call.
    """

    return subprocess.call(cmd, shell=True)


# @log_fail
# @log_info
# def cte_correct():
#     """The main function of the ``cal_uvis_make_ctecorr_darks`` module.
#     See module docstrings for further details.
#     """

#     # List of proposals with UVIS darks to use
#     root_dirs = [

#     # SMOV
#     #'/grp/hst/wfc3a/Quicklook/11431/',
#     #'/grp/hst/wfc3a/Quicklook/11446/',

#     # Cycle 17
#     #'/grp/hst/wfc3a/Quicklook/11905/',
#     #'/grp/hst/wfc3a/Quicklook/11909/',

#     # Cycle 18
#     #'/grp/hst/wfc3a/Quicklook/12342/',
#     #'/grp/hst/wfc3a/Quicklook/12343/',

#     # Cycle 19
#     #'/grp/hst/wfc3a/Quicklook/12689/',
#     #'/grp/hst/wfc3a/Quicklook/12687/',
#     #'/grp/hst/wfc3a/Quicklook/13103/',

#     # Cycle 20
#     #'/grp/hst/wfc3a/Quicklook/13071/',
#     #'/grp/hst/wfc3a/Quicklook/13073/',
#     #'/grp/hst/wfc3a/Quicklook/13074/',
#     #'/grp/hst/wfc3a/Quicklook/13075/',
#     #'/grp/hst/wfc3a/Quicklook/13076/',

#     # Cycle 21
#     #'/grp/hst/wfc3a/Quicklook/13554/',
#     #'/grp/hst/wfc3a/Quicklook/13556/',
#     #'/grp/hst/wfc3a/Quicklook/13557/',
#     #'/grp/hst/wfc3a/Quicklook/13558/',

#     # Cycle 22
#     # '/grp/hst/wfc3a/Quicklook/14000/',
#     # '/grp/hst/wfc3a/Quicklook/14002/',
#     # '/grp/hst/wfc3a/Quicklook/14003/',
#     # '/grp/hst/wfc3a/Quicklook/14004/',

#     # Cycle 23
#     #'/grp/hst/wfc3a/Quicklook/14366/',
#     #'/grp/hst/wfc3a/Quicklook/14368/',
#     #'/grp/hst/wfc3a/Quicklook/14369/',
#     #'/grp/hst/wfc3a/Quicklook/14370/',
#     #'/grp/hst/wfc3a/unlooked/14366/',
#     #'/grp/hst/wfc3a/unlooked/14368/',
#     #'/grp/hst/wfc3a/unlooked/14369/',
#     #'/grp/hst/wfc3a/unlooked/14370/'

#     # Cycle 24
#     '/grp/hst/wfc3a/Quicklook/14529/',
#     '/grp/hst/wfc3a/Quicklook/14531/',
#     '/grp/hst/wfc3a/Quicklook/14532/',
#     '/grp/hst/wfc3a/Quicklook/14533/',

#     #Cycle 25
#     '/grp/hst/wfc3a/Quicklook/14978/',
#     '/grp/hst/wfc3a/Quicklook/14980/',
#     '/grp/hst/wfc3a/Quicklook/14981/',
#     '/grp/hst/wfc3a/Quicklook/14982/',

#     #Cycle 26
#     '/grp/hst/wfc3a/Quicklook/15567/',
#     '/grp/hst/wfc3a/Quicklook/15569/',
#     '/grp/hst/wfc3a/Quicklook/15570/',
#     '/grp/hst/wfc3a/Quicklook/15571/']

#     # Search for and copy and new files
#     for root in root_dirs:
#         print('\tLP log: ')
#         print('\tLP log: Scanning {} for new dark files'.format(root))
#         print('\tLP log: ')
#         if os.path.exists(root):
#             for path, subdirs, files in os.walk(root):
#                 for name in files:
#                     src = os.path.join(path, name)
#                     dst = os.path.join(CTE_CORR_DIR, os.path.basename(src))

#                     copy_dark(src, dst)

#     # Determine which files to CTE correct and perform correction
#     files_to_correct = get_files_to_correct()
#     perform_correction(files_to_correct)


# @log_fail
# @log_info
def cte_correct_lp(RWD_DIR, CTE_CORR_DIR, SOFTWARE_DIR):
    """LP edit to not rely on Quicklook, instead data is downloaded with 
    astroquery and files sorted accordinaly in download_data.py.
    A copy of the files is made subject to checks with copy_dark_lp, a 
    list of files to correct (without accompanying rac files) is made, 
    then the CTE correction is applied. 

    This is the main function of the ``cal_uvis_make_ctecorr_darks`` module.
    See module docstrings for further details.

    Parameters - LP added
    ----------
    RWD_DIR : str
        The absolute path of the raw data directory (all in one folder).

    CTE_CORR_DIR : str
        The absolute path of the CTE corrected data.

    SOFTWARE_DIR : str
        The absolute path of the wfc3uv_ctereverse.e executable.
    """

    # Check files and copy over 
    src = RWD_DIR   #Source, see notebook for file arranging
    dst = CTE_CORR_DIR   #Desitnation for where darks will be processes
    #Checks file headers and copies darks with 900s, 1800s exposure 
    # and targname of DARK*
    copy_dark_lp(src, dst)  

    # For files in the filepath, make a list of files to correct
    files_to_correct = get_files_to_correct(CTE_CORR_DIR)      #LP added CTE_CORR_DIR
    #Checks for the files with no CTE correction in the CTE_CORR_DIR directory, puts the format into a list 

    # perform CTE correction
    perform_correction(files_to_correct, CTE_CORR_DIR, SOFTWARE_DIR)  #LP added CTE_CORR_DIR, SOFTWARE_DIR


# LP added function to take user inputs
def parse_args():
    """Parses command line arguments.

    Returns
    -------
    args : obj
        An ``argparse`` object containing all of the added arguments.
    """

    ctecorr_dir_help = 'LP arg: Required user-defined path to the location of the raw CTE-corrected data. No trailing slash "/".'
    rwd_dir_help = 'LP arg: Required user-defined path to the location of the raw downloaded darks, all in one directory (see notebook for file organization if using astroquery). No trailing slash "/".'
    software_dir_help = 'LP arg: Required user-defined path to the compiled CTE correction code ./wfc3uv_ctereverse.e available here: http://www.stsci.edu/~jayander/X/EXPORT_WFC3UV_CTE/wfc3uv_ctereverse.F . No trailing slash "/".'  #Only use in CTE correction code
    
    parser = argparse.ArgumentParser()
    
    # Argument to output CTE corrected raw directory
    parser.add_argument('-c --ctecorr_dir',
        dest='CTE_CORR_DIR',
        action='store',
        type=str,
        required=True,   
        help=software_dir_help) 

    # Argument to input downloaded raw data directory
    parser.add_argument('-r --rwd_dir',
        dest='RWD_DIR',
        action='store',
        type=str,
        required=True,   
        help=software_dir_help) 

    # Argument to software directory
    parser.add_argument('-s --software_dir',
        dest='SOFTWARE_DIR',
        action='store',
        type=str,
        required=True,   
        help=software_dir_help)

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
    # ------------------------------------------
    # LP added tests
    # Check that input directories are strings
    assert args.CTE_CORR_DIR is str, 'Path to the raw CTE corrected data (ctecorr_dir) must be a string.'
    assert args.RWD_DIR is str, 'Path to the calibration file directory (rwd_dir) must be a string.'
    assert args.SOFTWARE_DIR is str, 'Path to the compiled CTE correction code ./wfc3uv_ctereverse.e (software_dir) must be a string.'
    # ------------------------------------------


if __name__ == '__main__':

    # module = os.path.basename(__file__).strip('.py')
    # configure_logging(module)

    # LP added
    now = datetime.datetime.now()   #LP changed name to now
    today = datetime.datetime.strftime(now, '%Y%m%d')  #LP changed name to now
    print('*****************************************************************************')
    print('LP log: CTE correction started at ', now.strftime("%Y-%m-%d %H:%M"))
    print('*****************************************************************************')

    args = parse_args()   # LP added
    # test_args(args)  # LP added

    cte_correct_lp(args.RWD_DIR, args.CTE_CORR_DIR, args.SOFTWARE_DIR)

    # LP added
    now = datetime.datetime.now()   #LP changed name to now
    today = datetime.datetime.strftime(now, '%Y%m%d')  #LP changed name to now
    print('*****************************************************************************')
    print('LP log: CTE correction ended at ', now.strftime("%Y-%m-%d %H:%M"))
    print('*****************************************************************************')

