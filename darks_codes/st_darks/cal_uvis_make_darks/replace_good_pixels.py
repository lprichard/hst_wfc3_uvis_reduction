"""Replace the pixels in each superdark that fall below the hot pixel
threshold with the values of the masterdark image.

The masterdark image from which good pixels are grabbed depends on the
``mode``.  If the ``mode`` is ``dev``, then the supplied anneal cycle
masterdark is used.  If the ``mode`` is ``prod``, then the previous
anneal cycle masterdark is used.  The previous anneal masterdarks are
used in production in order to keep UVIS dark reference file delivery
times reasonable, as opposed to waiting for the current anneal cycle to
end before creating dark reference files.

Authors
-------

    Matthew Bourque, 2014

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import glob
import logging
import os
import shutil

from astropy.io import fits
import numpy as np


def get_masterdark(anneal_date, paths, ctecorr, mode):
    """Return the full path of the previous anneal masterdark.

    Parameters
    ----------
    anneal_date : str
        The anneal date, in the format of ``YYYY-MM-DD``, to be used in
        the filename of the masterdark.
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    ctecorr : bool
        The CTE correction switch.
    mode : str
        The processing mode. Can either be ``dev``, in which
        masterdarks are made from the supplied anneal, or ``prod``
        (i.e. "production") in which the masterdarks are made from
        the previous anneal's darks.

    Returns
    -------
    masterdark : str
        The full path to the previous anneal cycle masterdark.
    """

    if mode == 'dev':
        masterdark = glob.glob(os.path.join(paths['masterdark_create_dir'], 'masterdark_*.fits'))[0]

    elif mode == 'prod':

        # Get list of masterdarks in the masterdarks directory
        if ctecorr:
            masterdarks = glob.glob(os.path.join(paths['masterdark_pool'], 'masterdark_*_ctecorr.fits'))
        else:
            masterdarks = glob.glob(os.path.join(paths['masterdark_pool'], 'masterdark_????-??-??.fits'))

        # Insert dummy current masterdark
        if ctecorr:
            current_masterdark = 'masterdark_{}_ctecorr'.format(anneal_date)
        else:
            current_masterdark = 'masterdark_{}'.format(anneal_date)
        current_masterdark = os.path.join(paths['masterdark_pool'], current_masterdark)
        masterdarks.append(current_masterdark)

        # Sort the masterdarks
        masterdarks = sorted(masterdarks)

        # Get index of current masterdark
        current_masterdark_index = masterdarks.index(current_masterdark)

        # Determine previous anneal masterdark
        masterdark = masterdarks[current_masterdark_index - 1]

    return masterdark


def process_anneal(superdark_list, masterdark):
    """Replace the pixels in each superdark that fall below the hot
    pixel threshold with the values of the masterdark image.

    Parameters
    ----------
    superdark_list : list
        A list of full paths to the superdarks.
    masterdark : str
        The full path to the masterdark.
    """

    logging.info('\tReplacing pixels in superdarks with values from masterdark {}'.format(masterdark))

    # Replace good pixels in superdarks
    for superdark in superdark_list:

        # Open the masterdark and the superdark
        masterdark_hdulist = fits.open(masterdark, mode='readonly')
        superdark_hdulist = fits.open(superdark, mode='update')

        # Find the non-good pixels in the superdark
        bad_pixels_ext3 = np.where(superdark_hdulist[3].data != 0)
        bad_pixels_ext6 = np.where(superdark_hdulist[6].data != 0)

        # Assume the new data takes the form of the masterdark
        new_data_ext1 = masterdark_hdulist[1].data.astype(np.float32)
        new_data_ext2 = masterdark_hdulist[2].data.astype(np.float32)
        new_data_ext4 = masterdark_hdulist[4].data.astype(np.float32)
        new_data_ext5 = masterdark_hdulist[5].data.astype(np.float32)

        # For non-good pixels, replace the values with that of the superdark
        new_data_ext1[bad_pixels_ext3] = superdark_hdulist[1].data[bad_pixels_ext3].astype(np.float32)
        new_data_ext2[bad_pixels_ext3] = superdark_hdulist[2].data[bad_pixels_ext3].astype(np.float32)
        new_data_ext4[bad_pixels_ext6] = superdark_hdulist[4].data[bad_pixels_ext6].astype(np.float32)
        new_data_ext5[bad_pixels_ext6] = superdark_hdulist[5].data[bad_pixels_ext6].astype(np.float32)
        superdark_hdulist[1].data = new_data_ext1
        superdark_hdulist[2].data = new_data_ext2
        superdark_hdulist[4].data = new_data_ext4
        superdark_hdulist[5].data = new_data_ext5

        # Save the changes
        superdark_hdulist.close()


def replace_good_pixels_main(paths, anneal_date, superdark_list, ctecorr, mode):
    """The main function of the ``replace_good_pixels`` module.  Please
    see module docstrings for further details.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    anneal_date : str
        The anneal date, in the format of ``YYYYMMDD``, to be used in
        the filename of the masterdark.
    superdark_list : list
        A list of full paths to the superdarks.
    ctecorr : bool
        The CTE correction switch.
    mode : str
        The processing mode. Can either be ``dev``, in which
        masterdarks are made from the supplied anneal, or ``prod``
        (i.e. "production") in which the masterdarks are made from
        the previous anneal's darks.
    """

    logging.info('\tReplacing pixels in superdarks with values from masterdark')

    # Determine which masterdark to use
    masterdark = get_masterdark(anneal_date, paths, ctecorr, mode)

    if mode == 'prod':
        # Copy masterdark to the local masterdark directory
        masterdark_dst = os.path.join(paths['masterdark_create_dir'], os.path.basename(masterdark))
        shutil.copyfile(masterdark, masterdark_dst)

    process_anneal(superdark_list, masterdark)
