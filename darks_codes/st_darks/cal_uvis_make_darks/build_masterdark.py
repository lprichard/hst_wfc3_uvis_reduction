"""Create "masterdark" images, which are anneal-cycle-averaged darks.

Masterdarks are created by combining the nnblv files from the entire
anneal cycle.  These nnblv files have had hot pixels (that were
identified in the superdarks) flagged, and thus the only pixels that
are combined are those that were not flagged in the individual nnblv
or 4-day superdark as hot.

Since a running 4-day window is used during superdark creation, some
nnblv files are used in different superdark-creation processes, and
thus have different sets of hot pixels.  From this pool of redundant
nnblv files, this module chooses the "best" nnblv files to use in the
masterdark creation.  The "best" files are determined by which nnblvs
were observed closest to the middle of their respective 4-day window;
these nnblv files theoretically should have the most accurate measure
of hot pixels.

The chosen nnblv files are then average-combined to create the
"masterdark" image. Once the masterdark image is created, its values
are copied upon the "good" pixels of the 4-day superdarks (i.e. those
pixels that fall below the hot pixel threshold).  What is left are
4-day superdarks that have entire-anneal-cycled-averaged values and
accurate hot pixel flags.  These superdarks, after having their
headers prepped, can be delivered to the user community as dark
reference files.

Outputs include:

    1. masterdark.fits - The masterdark image
    2. Updated superdarks with their "good" pixels replaced with
       values from the masterdark image.

Authors
-------

    Matthew Bourque, 2016

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import datetime
import glob
import logging
import os
import shutil

from astropy.io import fits

from automated_scripts.cal_uvis_make_darks.combine import combine
from automated_scripts.cal_uvis_make_darks.electrons_normalize import electrons_normalize
from automated_scripts.cal_uvis_make_darks.replace_good_pixels import replace_good_pixels_main


def get_nnblvs_to_combine(paths):
    """Determine which nnblv files to use in masterdark creation.

    Each superdark directory contains nnblv files from a 4-day window
    that were used in the superdark creation, as well as the superdark
    itself.  The masterdark requires the combination of all nnblv files
    from the anneal cycle.  However, since a sliding 4-day window is
    used in the creation of the superdarks, each superdark directory
    contains redundant nnblv files with other superdark directories.
    Thus, each superdark directory is searched for the "best" nnblv
    files to use for the masterdark creation.  "Best" in this case
    is defined as those nnblv files whose observation times are
    closest to that of the superdark ``USEAFTER`` + 2 days.  In other
    words, nnblv files whose observations are closest to the middle of
    the 4-day window are identified.

    The full paths of those nnblv files that are identified as "best"
    are placed into the ``nnblv_list``.  This list is later used to
    copy the appropriate files to the masterdark directory for the
    creation of the masterdark.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.

    Returns
    -------
    nnblv_list : list
        A list of full paths to the nnblv files in the superdark
        directories to be used in the creation of the masterdark.
    """

    superdark_dirs = glob.glob(paths['superdark_dirs'])
    nnblv_filenames, nnblv_dirs, nnblv_date_indices = [], [], []
    master_nnblv_dict = {}

    for superdark_dir in superdark_dirs:

        # Get list of files relative to the superdark_dir
        local_superdark = glob.glob(os.path.join(superdark_dir, '*drk.fits'))[0]
        local_nnblvs = glob.glob(os.path.join(superdark_dir, 'nn*blv_tmp.fits'))

        # Get the superark USEAFTER + 2 days
        superdark_hdulist = fits.open(local_superdark)
        superdark_useafter = superdark_hdulist[0].header['USEAFTER']
        superdark_useafter = datetime.datetime.strptime(superdark_useafter, '%b %d %Y %H:%M:%S')
        superdark_index = superdark_useafter + datetime.timedelta(days=2)

        for local_nnblv in local_nnblvs:

            # Determine the EXPSTART of the nnblv file
            nnblv_hdulist = fits.open(local_nnblv)
            nnblv_dateobs = nnblv_hdulist[0].header['DATE-OBS']
            nnblv_timeobs = nnblv_hdulist[0].header['TIME-OBS']
            nnblv_expstart = '{} {}'.format(nnblv_dateobs, nnblv_timeobs)
            nnblv_expstart = datetime.datetime.strptime(nnblv_expstart, '%Y-%m-%d %H:%M:%S')

            # Determine how the nnblv EXPSTART compares to the superdark USEAFTER + 2 days
            nnblv_time_diff = superdark_index - nnblv_expstart
            nnblv_time_diff = nnblv_time_diff.days + (float(nnblv_time_diff.seconds) / (60. * 60. * 24.))
            nnblv_time_diff = abs(nnblv_time_diff)

            # Add this information to lists
            nnblv_filenames.append(os.path.basename(local_nnblv))
            nnblv_dirs.append(os.path.dirname(local_nnblv))
            nnblv_date_indices.append(nnblv_time_diff)

    for nnblv_filename, nnblv_dir, nnblv_date_index in zip(nnblv_filenames, nnblv_dirs, nnblv_date_indices):

        # If the nnblv is not already in the dictionary, add it.  Otherwise, update the
        # dictionary entry if the relative time between observations is shorter.
        if nnblv_filename not in list(master_nnblv_dict.keys()):
            master_nnblv_dict[nnblv_filename] = [nnblv_dir, nnblv_date_index]
        else:
            date_index_to_compare = master_nnblv_dict[nnblv_filename][1]
            if nnblv_date_index < date_index_to_compare:
                master_nnblv_dict[nnblv_filename] = [nnblv_dir, nnblv_date_index]

    # Build list of full paths to the nnblv files to copy
    nnblv_list = [os.path.join(record[1][0], record[0]) for record in list(master_nnblv_dict.items())]

    return nnblv_list


def build_masterdark_main(paths, anneal_date, ctecorr, mode):
    """The main function of the ``build_masterdark`` module.

    A masterdark image is created by combining all nnblv files from the
    anneal cycle.  nnblv files are appropriately chosen from the
    superdark directories based on their observation date relative to
    the middle of the 4-day window.  The nnblv files are then combined.
    The units of the masterdark are then changed to e-/s.  Lastly, the
    "good" pixels in the superdarks are replaced by values in the
    masterdark.  What is left are 4-day window superdarks with their
    "good" pixels replaced by entire-anneal-cycle-averaged masterdark
    values.  "Good" pixels are defined as those that fall below the hot
    pixel threshold.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    anneal_date: str
        The anneal date, in the format of ``YYYYMMDD``, to be used in
        the filename of the masterdark.
    ctecorr : bool
        The CTE correction switch.
    mode : str
        The processing mode.  Can either be ``dev``, in which
        masterdarks are made from the supplied anneal, or ``prod``
        (i.e. "production") in which the masterdarks are made from
        the previous anneal's darks.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Building masterdark ----------')
    logging.info('')

    # Make masterdark directory
    os.mkdir(paths['masterdark_create_dir'], 0o774)
    logging.info('\tCreated directory {}'.format(paths['masterdark_create_dir']))
    logging.info('')

    # Fill masterdark directory with nnblvs and 4-day superdarks
    nnblvs_list = get_nnblvs_to_combine(paths)
    superdark_list = glob.glob(os.path.join(paths['superdark_dirs'], 'd*drk.fits'))
    for nnblv in nnblvs_list:
        shutil.copy(nnblv, paths['masterdark_create_dir'])
        logging.info('\tCopied {} to {}'.format(nnblv, paths['masterdark_create_dir']))
    logging.info('')
    for superdark in superdark_list:
        shutil.copy(superdark, paths['masterdark_create_dir'])
        logging.info('\tCopied {} to {}'.format(superdark, paths['masterdark_create_dir']))
    logging.info('')

    # Combine data
    nnblv_list = glob.glob(os.path.join(paths['masterdark_create_dir'], 'nn*blv_tmp.fits'))
    if ctecorr:
        masterdark = os.path.join(paths['masterdark_create_dir'], 'masterdark_{}_ctecorr.fits'.format(anneal_date))
    else:
        masterdark = os.path.join(paths['masterdark_create_dir'], 'masterdark_{}.fits'.format(anneal_date))
    combine(nnblv_list, masterdark)

    # Set the masterdark permissions
    os.chmod(masterdark, 0o774)
    os.chown(masterdark, -1, 340)

    # Convert units
    electrons_normalize(masterdark)

    # Replace "good pixels" in superdark with masterdark values
    new_superdark_list = glob.glob(os.path.join(paths['masterdark_create_dir'], 'd*drk.fits'))
    replace_good_pixels_main(paths, anneal_date, new_superdark_list, ctecorr, mode)
