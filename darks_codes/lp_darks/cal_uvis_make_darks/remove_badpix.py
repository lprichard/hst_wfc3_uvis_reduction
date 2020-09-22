"""Removes pixels in the ``DQ`` array of ``*blv_tmp.fits`` files that
have flags other than cosmic ray flags (``DQ`` value of ``8192``).

Operates on all ``blv_tmp.fits`` files in the given directory. Non-CR
pixels are removed from both ``DQ`` FITS extensions (``3`` and ``6``).
This is done before the CR-grow step as to avoid growing the affected
pixels (and thus eventually oversubtracting in the dark reference file
product).

Authors
-------
    Benjamin Kuhn, February 2020
    Matthew Bourque, April 2016

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import glob
import logging
import os

from astropy.io import fits
import numpy as np


def remove_badpix_main(proc_dir):
    """The main function of the ``remove_badpix`` module.  See module
    docstring for further details.

    Parameters
    ----------
    proc_dir
        The absolute path of the directory to process.
    """

    # LP edit, not the same directory now
    # proc_dir_short = proc_dir.split('uvis_darks/')[-1]

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Removing non-cosmic ray flags for {} ----------'.format(proc_dir))   #LP edit: removed proc_dir_short
    logging.info('')

    blvs = glob.glob(os.path.join(proc_dir, '*blv_tmp.fits'))

    if len(blvs) > 0:

        for blv in blvs:

            # Get the data
            hdulist = fits.open(blv, mode='update')
            data_ext3 = hdulist[3].data
            data_ext6 = hdulist[6].data

            # Find ALL saturated pixels and CR pixels
            # Returns an array of 0's, 256's, 8192's, and 8448's
            newdq3 = data_ext3 & (256 | 8192)
            newdq6 = data_ext6 & (256 | 8192)

            # Replace everything == 256 and greater than 8192 with 8192
            # We're treating saturated pixels as CRs
            newdq3[np.where((newdq3 == 256) | (newdq3 > 8192))] = 8192
            newdq6[np.where((newdq6 == 256) | (newdq6 > 8192))] = 8192

            # Populated the 3rd and 6th extensions with new DQ array
            hdulist[3].data = newdq3
            hdulist[6].data = newdq6

            hdulist.close()

    else:
        logging.info('\tNo data for {}'.format(proc_dir))
