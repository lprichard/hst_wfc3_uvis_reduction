"""Grow the cosmic ray mask by a 2-pixel radius.

The cosmic ray masks (i.e. pixels with a ``DQ`` flag of ``8192``) are
"grown" (i.e. expanded) by a 2-pixel radius for all ``blv_tmp.fits``
files.  This is done as a conservative measure because ``CALWF3`` fails
to mask the wings of some cosmic rays properly.

Output products are ``nn*blv_tmp.fits`` files, which are cr-grown
``blv_tmp.fits`` files.

Authors
-------

    - Matthew Bourque, 2013
    - John Biretta, 2012

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import glob
import itertools
import logging
import multiprocessing
import os
import shutil

from astropy.io import fits
import skimage.morphology as morph
from skimage.morphology import disk


def process_mask(args):
    """Perform the growing routine.

    Parameters
    ----------
    args : tuple
        The arguments given by ``multiprocessing.Pool.map``. The 0th
        element is the file to process, and the 1st element is the
        directory to which to save the CR-grown file.
    """

    # Parse args
    blv_file = args[0]
    crgrow_dir = args[1]

    # Grow the mask
    logging.info('\tGrowing {}'.format(blv_file))
    dq3 = fits.getdata(blv_file, ext=3)
    dq6 = fits.getdata(blv_file, ext=6)
    dq3_grown = morph.dilation(dq3.byteswap().newbyteorder('='), disk(2))
    dq6_grown = morph.dilation(dq6.byteswap().newbyteorder('='), disk(2))

    # Write grown mask to new FITS file
    hdulist = fits.open(blv_file)
    hdulist[3].data = dq3_grown
    hdulist[6].data = dq6_grown
    hdulist[3].header['EXTNAME'] = 'DQ'
    hdulist[6].header['EXTNAME'] = 'DQ'
    hdulist[3].header['EXTVER'] = '1'
    hdulist[6].header['EXTVER'] = '2'
    hdulist.writeto(os.path.join(crgrow_dir, 'nn{}'.format(os.path.basename(blv_file))))


def grow_cr_mask_main(crrej_dir):
    """The main function of the ``grow_cr_mask`` module.

    A new ``crgrow/`` directory is made that houses the affected
    ``blv_tmp.fits`` files and the resulting ``nn*blv_tmp.fits`` files.

    Parameters
    ----------
    ccrej_dir : str
        The path to the crrej directory to create.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Growing cosmic ray mask for {} ----------'.format(crrej_dir.split('uvis_darks/')[1]))
    logging.info('')

    # Prepare crgrow directory
    crgrow_dir = os.path.join(crrej_dir, 'crgrow/')
    os.mkdir(crgrow_dir, 0o774)
    logging.info('\tCreated directory {}'.format(crgrow_dir))
    logging.info('')
    blv_list = glob.glob(os.path.join(crrej_dir, '*blv_tmp.fits'))

    if len(blv_list) > 0:
        for blv_file in blv_list:
            shutil.move(blv_file, crgrow_dir)
            logging.info('\tMoved {} to {}'.format(blv_file, crgrow_dir))
        logging.info('')

        # Grow CR mask
        blv_list = glob.glob(os.path.join(crgrow_dir, '*blv_tmp.fits'))
        mp_args = list(zip(blv_list, itertools.repeat(crgrow_dir)))
        pool = multiprocessing.Pool(processes=30)
        pool.map(process_mask, mp_args)
        pool.close()
        pool.join()

    else:
        logging.info('\tNo data for {}'.format(crrej_dir))
