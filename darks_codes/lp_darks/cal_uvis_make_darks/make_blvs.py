"""Generate bias-corrected, cr-flagged, postflash-corrected
``blv_tmp.fits`` files.

This module contains several functions that set header keyword
values, builds association tables, performs basic calibration steps
(via ``CALWF3``), and unit conversions to create ``blv_tmp.fits``
files.  The ``blv_tmp.fits`` files are bias and postflash-corrected,
and have their cosmic rays flagged. As a consequence of the execution
of ``CALWF3``, dark combined images are also created.

Output products are CR-flagged, bias/postflashed (if necessary)
corrected ``blv_tmp.fits`` files, accompanying ``.tra`` files, combined
``*drk.fits`` files, ``*drk_asn.fits`` association tables, and
``*drk_crj*.fits`` cosmic ray rejected images.

Authors
-------

    Matthew Bourque, July 2014

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.

Notes
-----

    As a consequence of the execution of ``CALWF3``, dark combined
    ``*drk.fits`` images and accompanying cosmic-ray-rejected iamges
    are also created.  This step cannot be turned off in ``CALWF3``.
    The resulting combined images are useless for the remainder of the
    algorithm.
"""

import glob
import os
import logging
import warnings

from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyDeprecationWarning
from wfc3tools.calwf3 import calwf3
from pdb import set_trace as st

warnings.filterwarnings('ignore', category=AstropyDeprecationWarning)


def buildasntable(proc_dir):
    """Build an association table named by rootname containing all of
    the raw files in the working directory.

    The resulting association table (``dark_asn.fits``) is used by
    ``CALWF3`` to perform image combinations and generate
    ``blv_tmp.fits`` files.

    Parameters
    ----------
    proc_dir : str
        The absolute path of the directory in which the files
        to be processed reside.
    """

    print('\tBuilding the association table.')

    # Gather RAW files to process
    os.chdir(proc_dir)
    filelist = glob.glob('*raw.fits')
    filelist.append('dark')

    # Build the columns of the association table
    row = [i+1 for i in range(len(filelist))]
    memtype = ['EXP-CRJ' for f in filelist]
    memprsnt = [True for f in filelist]
    xoffset = [0. for f in filelist]
    yoffset = [0. for f in filelist]
    xdelta = [0. for f in filelist]
    ydelta = [0. for f in filelist]
    rotation = [0. for f in filelist]
    scale = [1. for f in filelist]

    # Change the last row to have a memtype of PROD-CRJ instead of EXP-CRJ
    memtype[-1] = 'PROD-CRJ'

    # Build the table
    columns = [row, filelist, memtype, memprsnt, xoffset, yoffset, xdelta, ydelta, rotation, scale]
    names = ['row', 'MEMNAME', 'MEMTYPE', 'MEMPRSNT', 'XOFFSET', 'YOFFSET', 'XDELTA', 'YDELTA', 'ROTATION', 'SCALE']
    asn_table = Table(columns, names=names)

    # Save the table to a FITS file
    asn_table.write('dark_asn.fits')

    # Update header to have required keywords
    hdulist = fits.open('dark_asn.fits', mode='update')
    hdulist[0].header['INSTRUME'] = 'WFC3'
    hdulist[0].header['DETECTOR'] = 'UVIS'
    hdulist.close()


def run_calwf3(iref_dir, calib_dir):  #LP added iref_dir, calib_dir
    """Perform cosmic-ray rejection and image combination on the
    associtation table (``dark_asn.fits``).

    Parameters - LP added
    ----------
    iref_dir : str
        The absolute path of the IREF files used in reduction.
    calib_dir : str
        The absolute path of the STScI calibration files used in the reduction.
    """

    print('\tRunning CALWFC3 to create blv_tmp.fits files.')

    # LP added, for custom iref files
    print('LP log: iref_dir environment:', iref_dir + '/')
    os.environ['iref'] = iref_dir + '/'

    # LP added, as long file paths are truncated
    print('LP log: calib_dir environment:', calib_dir + '/')
    os.environ['calib_dir'] = calib_dir + '/'  #Get final name

    # Run calwf3, which will make blv_tmp.fits files in the directory
    calwf3(input='dark_asn.fits', save_tmp=True)
    # LP added
    print('LP log: calwf3 has finished!!')


def make_blvs_main(proc_dir, iref_dir, calib_dir):    #LP added iref_dir, calib_dir
    """Executes ``CALWF3`` on input RAW darks to generate
    bias-corrected, cr-flagged, postflash-corrected ``blv_tmp.fits``
    files.


    Parameters
    ----------
    proc_dir : str
        The absolute path of the directory in which the files to be
        processed reside.
    LP added:
    iref_dir : str
        The absolute path of the IREF files used in reduction.
    calib_dir : str
        The absolute path of the STScI calibration files used in the reduction.
    """

    # LP comment
    # proc_dir_short = proc_dir.split('uvis_darks/')[1]
    print('')
    print('')
    print('')
    print('---------- Making blv_tmp.fits files for {} ----------'.format(proc_dir))   #LP remove proc_dir_short
    print('')

    # Gather list of files
    list_of_raws = glob.glob(os.path.join(proc_dir, '*raw.fits'))

    if len(list_of_raws) > 0:

        # Build the association table.
        buildasntable(proc_dir)

        # Make the blv frames with calwf3.
        run_calwf3(iref_dir, calib_dir)    #LP added calib dir

    else:
        print('\tNo data for {}'.format(proc_dir))
