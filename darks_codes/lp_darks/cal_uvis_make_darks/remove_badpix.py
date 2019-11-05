"""Removes pixels in the ``DQ`` array of ``*blv_tmp.fits`` files that
have flags other than cosmic ray flags (``DQ`` value of ``8192``).

Operates on all ``blv_tmp.fits`` files in the given directory. Non-CR
pixels are removed from both ``DQ`` FITS extensions (``3`` and ``6``).
This is done before the CR-grow step as to avoid growing the affected
pixels (and thus eventually oversubtracting in the dark reference file
product).

Authors
-------

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

    print('')
    print('')
    print('')
    print('---------- Removing non-cosmic ray flags for {} ----------'.format(proc_dir))   #removed proc_dir_short
    print('')

    blvs = glob.glob(os.path.join(proc_dir, '*blv_tmp.fits'))

    if len(blvs) > 0:

        for blv in blvs:

            # Get the data
            hdulist = fits.open(blv, mode='update')
            data_ext3 = hdulist[3].data
            data_ext6 = hdulist[6].data

            # Set 256 pixels (saturated pixels) to 8192 to treat them as
            # cosmic rays
            data_ext3[np.where(data_ext3 == 256.0)] = 8192.0
            data_ext6[np.where(data_ext6 == 256.0)] = 8192.0

            # Replace everything under 8192 (cosmmic ray) with 0
            data_ext3[np.where((0.0 < data_ext3) & (data_ext3 < 8192.0))] = 0.0
            data_ext6[np.where((0.0 < data_ext6) & (data_ext6 < 8192.0))] = 0.0

            # Reset everything above 8192 (cosmic ray + others) back to 8192
            data_ext3[np.where(data_ext3 > 8192.0)] = 8192.0
            data_ext6[np.where(data_ext6 > 8192.0)] = 8192.0

            hdulist.close()

    else:
        print('\tNo data for {}'.format(proc_dir))
