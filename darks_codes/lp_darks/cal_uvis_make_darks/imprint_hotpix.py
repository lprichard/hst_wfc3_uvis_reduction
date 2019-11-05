"""Imprints hot pixel flags from superdarks back onto the corresponding
``nn*blv_tmp.fits`` files.

In order to create a high-quality anneal-cycle-averaged 'masterdark',
the input ``nn*blv_tmp.fits`` files must have an accurate measure of
the hot pixels present around the time of observations.  Thus, those
pixels that are flagged as hot in the nominal 4-day superdarks are
imprinted back onto the ``DQ`` arrays of the ``nn*blv_tmp.fits`` files
that went into the creation of the superdark itself.

Output products are ``tmp_ext<ext>_<filename>.fits`` files, which are
temporary files containing a common hot pixel mask between the
``nn*blv_tmp.fits`` file and the superdark for the particular
extension.  These files are simply one-extension files that contain the
``DQ`` array holding the common hot pixel flags.

Authors
-------

    Matthew Bourque, 2016

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import glob
import logging
import multiprocessing
import os

from astropy.io import fits
import numpy as np


def process_masks(superdark_dir):
    """Perform the imprinting of the hot pixel mask for the group of
    CR-grown ``blv_tmp.fits`` files in the given superdark directory.

    Parameters
    ----------
    superdark_dir : str
        The path to the superdark directory to process.
    """

    superdark = glob.glob(os.path.join(superdark_dir, '*drk.fits'))[0]
    nnblvs = glob.glob(os.path.join(superdark_dir, 'nn*blv_tmp.fits'))

    for nnblv in nnblvs:
        for ext in [3, 6]:

            print('\tImprinting hot pixel from {} onto {} EXT {}'.format(
                os.path.basename(superdark), os.path.basename(nnblv), ext))

            superdark_data = fits.getdata(superdark, ext)
            hdulist_nnblv = fits.open(nnblv, mode='update')
            nnblv_data = hdulist_nnblv[ext].data
            hdulist_nnblv[ext].data[np.where((nnblv_data > 1) | (superdark_data > 1))] = 2
            hdulist_nnblv.close()


def imprint_hotpix_main(paths):
    """The main function of the ``imprint_hotpix`` module.

    Gathers a list of ``nn*blv_tmp.fits`` files that correspond to a
    superdark and calls functions to create and place the common hot
    pixel mask back onto the ``nn*blv_tmp.fits`` files.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    """

    print('')
    print('')
    print('')
    print('---------- Imprinting hot pixels ----------')
    print('')

    superdark_dirs = glob.glob(paths['superdark_dirs'])

    # LP halved the number of processes from 30 to 15 due to processign capacity
    pool = multiprocessing.Pool(processes=15)
    pool.map(process_masks, superdark_dirs)
    pool.close()
    pool.join()
