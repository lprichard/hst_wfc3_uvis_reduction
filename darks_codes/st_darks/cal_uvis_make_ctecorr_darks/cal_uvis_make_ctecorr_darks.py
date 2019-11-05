#! /usr/bin/env python

"""Retrieves WFC3/UVIS darks and performs the CTE correction.

This script retrieves WFC3/UVIS darks from Quicklook, copies them into
the directory ``/grp/hst/wfc3n/ctecorr_darks/``, and uses the WFC3 CTE
correction code to CTE correct any images in that directory that do not
already have an accompanying ``RAC`` file.  ``RAC`` files are stored in
the same directory.  Only 900s darks from the UVIS anneal and UVIS CCD
Daily Monitor programs that are not already in the destination directory
are copied.

******************************************
*****Edit by Laura Prichard, May 2019*****
Without Quicklook access, download the required files for the anneal cycle
using astroquery. Select on (tmin=) anneal start and end date (stop just short of 
the new anneal cycle), target_name="DARK*", intentType='calibration', and 
e.g., instrument_name="WFC3/UVIS". Store these in e.g.:
``[filepath]/darks_reduction/lp_darks/raw_files/Feb21_2018anneal/HST/``
Code then copies the files into: 
``[filepath]/darks_reduction/lp_darks/raw_files/Feb21_2018anneal/ctecorr_darks/``
******************************************

Authors
-------

    Matthew Bourque, June 2014

Use
---

    This script is intended to run daily via a cron job, or via the
    command line as such:
    ::

        python cte_correct.py
"""

import glob
import logging
import multiprocessing
import os
import shutil
import subprocess

import astropy.io.fits as fits

from pyql.logging.logging_functions import configure_logging
from pyql.logging.logging_functions import log_info
from pyql.logging.logging_functions import log_fail

CTE_CORR_LOC = '/grp/hst/wfc3q/ctecorr_darks/'


def copy_dark(src, dst):
    """Copies a file from source to destination if it meets the header
    value criteria of ``TARGNAME = 'DARK'`` or ``'DARK-NM'``, and
    ``EXPTIME = 900``.

    Parameters
    ----------
    src : str
        The absolute path of the source location.
    dst : str
        The absolute path of the destination location.
    """

    if 'raw' in src:

        # Get header information
        exptime = fits.getval(src, 'EXPTIME')
        targname = fits.getval(src, 'TARGNAME')

        # Only copy appropriate files
        if targname == 'DARK' or targname == 'DARK-NM':
            if exptime in [900.0, 1800.0]:
                if not os.path.exists(dst):
                    logging.info('\tCopying {} to {}'.format(src, dst))
                    shutil.copyfile(src, dst)


def get_files_to_correct():
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

    raw_files = glob.glob(os.path.join(CTE_CORR_LOC, '*raw.fits'))
    rac_files = glob.glob(os.path.join(CTE_CORR_LOC, '*rac.fits'))

    raw_basenames = [os.path.basename(raw).split('_')[0] for raw in raw_files]
    rac_basenames = [os.path.basename(rac).split('_')[0] for rac in rac_files]

    for raw_basename in raw_basenames:
        if raw_basename not in rac_basenames:
            files_to_correct.append('{}_raw.fits'.format(raw_basename))

    return files_to_correct


def perform_correction(filelist):
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

    logging.info('')

    if len(filelist) > 0:
        for rawfile in filelist:
            job_list.append('cd {}; ./wfc3uv_ctereverse.e {}'.format(
                CTE_CORR_LOC, rawfile))
        logging.info('Processing {} jobs'.format(str(len(job_list))))
        pool = multiprocessing.Pool(processes=4)
        pool.map(run_process, job_list)
        pool.close()
        pool.join()
    else:
        logging.info('No files to correct')


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


@log_fail
@log_info
def cte_correct():
    """The main function of the ``cal_uvis_make_ctecorr_darks`` module.
    See module docstrings for further details.
    """

    # List of proposals with UVIS darks to use
    root_dirs = [

    # SMOV
    #'/grp/hst/wfc3a/Quicklook/11431/',
    #'/grp/hst/wfc3a/Quicklook/11446/',

    # Cycle 17
    #'/grp/hst/wfc3a/Quicklook/11905/',
    #'/grp/hst/wfc3a/Quicklook/11909/',

    # Cycle 18
    #'/grp/hst/wfc3a/Quicklook/12342/',
    #'/grp/hst/wfc3a/Quicklook/12343/',

    # Cycle 19
    #'/grp/hst/wfc3a/Quicklook/12689/',
    #'/grp/hst/wfc3a/Quicklook/12687/',
    #'/grp/hst/wfc3a/Quicklook/13103/',

    # Cycle 20
    #'/grp/hst/wfc3a/Quicklook/13071/',
    #'/grp/hst/wfc3a/Quicklook/13073/',
    #'/grp/hst/wfc3a/Quicklook/13074/',
    #'/grp/hst/wfc3a/Quicklook/13075/',
    #'/grp/hst/wfc3a/Quicklook/13076/',

    # Cycle 21
    #'/grp/hst/wfc3a/Quicklook/13554/',
    #'/grp/hst/wfc3a/Quicklook/13556/',
    #'/grp/hst/wfc3a/Quicklook/13557/',
    #'/grp/hst/wfc3a/Quicklook/13558/',

    # Cycle 22
    # '/grp/hst/wfc3a/Quicklook/14000/',
    # '/grp/hst/wfc3a/Quicklook/14002/',
    # '/grp/hst/wfc3a/Quicklook/14003/',
    # '/grp/hst/wfc3a/Quicklook/14004/',

    # Cycle 23
    #'/grp/hst/wfc3a/Quicklook/14366/',
    #'/grp/hst/wfc3a/Quicklook/14368/',
    #'/grp/hst/wfc3a/Quicklook/14369/',
    #'/grp/hst/wfc3a/Quicklook/14370/',
    #'/grp/hst/wfc3a/unlooked/14366/',
    #'/grp/hst/wfc3a/unlooked/14368/',
    #'/grp/hst/wfc3a/unlooked/14369/',
    #'/grp/hst/wfc3a/unlooked/14370/'

    # Cycle 24
    '/grp/hst/wfc3a/Quicklook/14529/',
    '/grp/hst/wfc3a/Quicklook/14531/',
    '/grp/hst/wfc3a/Quicklook/14532/',
    '/grp/hst/wfc3a/Quicklook/14533/',

    #Cycle 25
    '/grp/hst/wfc3a/Quicklook/14978/',
    '/grp/hst/wfc3a/Quicklook/14980/',
    '/grp/hst/wfc3a/Quicklook/14981/',
    '/grp/hst/wfc3a/Quicklook/14982/',

    #Cycle 26
    '/grp/hst/wfc3a/Quicklook/15567/',
    '/grp/hst/wfc3a/Quicklook/15569/',
    '/grp/hst/wfc3a/Quicklook/15570/',
    '/grp/hst/wfc3a/Quicklook/15571/']

    # Search for and copy and new files
    for root in root_dirs:
        logging.info('')
        logging.info('Scanning {} for new dark files'.format(root))
        logging.info('')
        if os.path.exists(root):
            for path, subdirs, files in os.walk(root):
                for name in files:
                    src = os.path.join(path, name)
                    dst = os.path.join(CTE_CORR_LOC, os.path.basename(src))

                    copy_dark(src, dst)

    # Determine which files to CTE correct and perform correction
    files_to_correct = get_files_to_correct()
    perform_correction(files_to_correct)


if __name__ == '__main__':

    module = os.path.basename(__file__).strip('.py')
    configure_logging(module)

    cte_correct()
