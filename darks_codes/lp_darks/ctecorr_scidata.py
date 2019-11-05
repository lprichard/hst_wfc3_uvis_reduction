#! /usr/bin/env python

"""Retrieves WFC3/UVIS science data and performs the CTE correction.

This script retrieves WFC3/UVIS science data from a given directory, copies 
them into a specified cte correction directory, and uses the WFC3 CTE
correction code to CTE correct any images in that directory that do not
already have an accompanying ``RAC`` file.  ``RAC`` files are stored in
the same directory.

Authors
-------

    Laura Prichard, adapted for science data, September 2019
    Matthew Bourque [originally the cal_uvis_make_ctecorr_darks.py code], June 2014

Use
---

    This code is intended to be run from the command line using the following:
    ::

        python ctecorr_scidata.py [-c|--ctecorr_dir] [-r|--rwd_dir] 
            [-s|--software_dir]

    ``-c --ctecorr_dir`` (Required) -- path to output CTE corrected data directory. 
        String, no trailing slash "/".

    ``-r --rwd_dir`` (Required) -- path to raw data directory (all files in one 
        folder as done in previous step). String, no trailing slash "/".

    ``-s --software_dir`` (Required) -- path to the compiled CTE correction code 
        ./wfc3uv_ctereverse.e. String, no trailing slash "/".

"""

import glob
# import logging
import multiprocessing
import os
import shutil
import subprocess
import socket    #LP added
import datetime  #LP added
from pdb import set_trace as st  #LP added
import argparse  #LP added

import astropy.io.fits as fits


def copy_sci_data(src, dst):
    """Copies a file from source to destination if it is raw and meets 
    the header value criteria of ``FILETYPE = 'SCI*'``.

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
        file = filepath.split('/')[-1]

        # If that file is raw, then check header
        if 'raw' in file:

            # Get header information
            filetype = fits.getval(filepath, 'FILETYPE')

            # Only copy appropriate files
            if filetype == 'SCI':
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
    """

    job_list = []

    print('')

    if len(filelist) > 0:
        for rawfile in filelist:
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


def cte_correct_lp(RWD_DIR, CTE_CORR_DIR, SOFTWARE_DIR):
    """The main function of the ``ctecorr_scidata`` module.
    See module docstrings for further details.
    """

    # Check files and copy over 
    src = RWD_DIR   #Source, see notebook for file arranging
    dst = CTE_CORR_DIR   #Destination for where science data will be processes
    #Checks file headers and copies science data
    copy_sci_data(src, dst)  

    # For files in the filepath, make a list of files to correct
    files_to_correct = get_files_to_correct(CTE_CORR_DIR)      #LP added CTE_CORR_DIR
    #Checks for the files with no CTE correction in the CTE_CORR_DIR directory, puts the format into a list 

    # perform CTE correction
    perform_correction(files_to_correct, CTE_CORR_DIR, SOFTWARE_DIR)  #LP added CTE_CORR_DIR, SOFTWARE_DIR


# LP added function to take user inputs
def parse_args():

    ctecorr_dir_help = 'Required user-defined path to the location of the raw CTE-corrected data. String, no trailing slash "/".'
    rwd_dir_help = 'Required user-defined path to the location of the raw downloaded darks, all in one directory (see notebook for file organization if using astroquery). String, no trailing slash "/".'
    software_dir_help = 'Required user-defined path to the compiled CTE correction code ./wfc3uv_ctereverse.e available here: http://www.stsci.edu/~jayander/X/EXPORT_WFC3UV_CTE/wfc3uv_ctereverse.F . String, no trailing slash "/".'  #Only use in CTE correction code
    
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

    # LP added
    now = datetime.datetime.now()   #LP changed name to now
    print('*****************************************************************************')
    print('LP log: CTE correction started at ', now.strftime("%Y-%m-%d %H:%M"))
    print('*****************************************************************************')

    args = parse_args()   # LP added
    # test_args(args)  # LP added

    cte_correct_lp(args.RWD_DIR, args.CTE_CORR_DIR, args.SOFTWARE_DIR)

    # LP added
    now = datetime.datetime.now()   #LP changed name to now
    print('*****************************************************************************')
    print('LP log: CTE correction ended at ', now.strftime("%Y-%m-%d %H:%M"))
    print('*****************************************************************************')

