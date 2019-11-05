"""Copies UVIS darks observed within a given time range to the
appropriate ``<SHUTRPOS>/`` directory.

UVIS darks are copied from the ``ctecorr_dir`` directory in the
``paths `` dict.  Only files that meet the criteria of
``TARGNAME = {'DARK' || 'DARK-NM'}``, ``EXPTIME = {'900' || '1800'}``,
and ``FLASHLVL = 12`` are copied.  Files are copied to
``/grp/hst/wfc3k/uvis_darks/<outdir>/<SHUTRPOS>``, where ``outdir``
is the main directory as determined by the ``-d --outdir`` argument.
Any files listed in ``/grp/hst/wfc3k/uvis_darks/exclude.list``, which
is a list of anomalous UVIS dark files, are ignored and are not copied.

Authors
-------

    - Matthew Bourque, 2013
    - John Biretta, 2012

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.


Notes
-----

    The absolute paths to any new UVIS dark programs in the Quicklook
    filesystem should be added to the ``dark_roots`` list as they come
    in. This occurs at the beginning of new calibration proposal
    cycles.
"""

import datetime
import glob
import logging
import multiprocessing
import os
import shutil

from astropy.io import fits


def copy_file(file_to_copy):
    """Copy the given file.

    Parameters
    ----------
    file_to_copy : tuple
        The 0th element of the tuple is the source, and the 1st
        element is the destination.
    """

    src = file_to_copy[0]
    dst = file_to_copy[1]

    shutil.copyfile(src, dst)


def pick_files(anneal_date, endtime, file_path, postflash, paths):
    """Copies a file from ``file_path`` to the appropriate
    ``<SHUTERPOS>/`` directory if the file meets the header value
    criteria of ``TARGNAME`` = ``DARK`` or ``DARK-NM``, ``EXPTIME``
    = ``900`` or ``1800``, and ``FLASHLVL`` = ``12`` (if processing
    postflashed data).

    Parameters
    ----------
    anneal_date : str
        The lower-limit of the ``DATE-OBS`` / ``TIME-OBS`` to be used
        as the starting point for data to copy, in the format
        ``YYYYMMDD-HH:MM:SS``
    endtime : str
        The upper-limit of the`` DATE-OBS`` / ``TIME-OBS`` to be used
        as the ending point for data to copy, in the format.
        ``YYYYMMDD-HH:MM:SS``
    file_path : str
        The absolute path of the potential UVIS dark file to copy.
    postflash : bool
        ``True`` if postflash processing is turned on, ``False`` if
        postflash processing is turned off. Used to determine which
        files to copy.
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.

    Notes
    -----
    If ``_rac.fits`` files are copied, they are renamed to
    ``_raw.fits`` for simplicity purposes.
    """

    # Open file and get header information
    hdulist = fits.open(file_path, mode='readonly')
    dateobs = hdulist[0].header['DATE-OBS']
    timeobs = hdulist[0].header['TIME-OBS']
    targname = hdulist[0].header['TARGNAME']
    exptime = hdulist[0].header['EXPTIME']
    try:
        flashlvl = hdulist[0].header['FLASHLVL']
    except:
        flashlvl = 0
    try:
        shutrpos = hdulist[0].header['SHUTRPOS']
    except:
        shutrpos = 'NA'
    hdulist.close()

    # Convert dateobs and timeobs to datetime objects
    dateobs = datetime.datetime.strptime(dateobs, '%Y-%m-%d')
    timeobs = datetime.datetime.strptime(timeobs, '%H:%M:%S')
    obs = datetime.datetime.combine(dateobs.date(), timeobs.timetz())

    # Determine which FLASHLVL is appropriate
    if postflash:
        flashlvl_to_use = 12
    else:
        flashlvl_to_use = 0

    # Find the appropriate files to copy
    if obs > anneal_date and obs < endtime and 'DARK' in targname and \
       exptime in [900, 1800] and flashlvl == flashlvl_to_use:

        logging.info('\t\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            file_path,
            str(obs),
            targname,
            str(exptime),
            str(flashlvl),
            shutrpos))

        rootname = os.path.basename(file_path).split('_')[0]
        basename = '{}_raw.fits'.format(rootname)

        if not postflash:
            dst = os.path.join(paths['preproc'], basename)
            return (file_path, dst)
        else:
            if shutrpos == 'A':
                dst = os.path.join(paths['shutterA'], basename)
                return (file_path, dst)
            elif shutrpos == 'B':
                dst = os.path.join(paths['shutterB'], basename)
                return (file_path, dst)

    return None


def copy_darks_main(anneal_date, endtime, ctecorr, postflash, paths):
    """The main function of the ``copy_darks`` module.

    Creates the ``<SHUTERPOS>/`` directories, if they do not already
    exist, iterates over the files in the root path, and copies
    appropriate files.

    Parameters
    ----------
    anneal_date : str
        The lower-limit of the ``DATE-OBS`` / ``TIME-OBS`` to be used
        as the starting point for data to copy, in the format
        ``YYYYMMDD-HH:MM:SS``.
    endtime : str
        The upper-limit of the ``DATE-OBS`` / ``TIME-OBS`` to be used
        as the ending point for data to copy, in the format
        ``YYYYMMDD-HH:MM:SS``.
    ctecorr : bool
        ``True`` (if CTE-correction is turned on) or ``False`` (if
        CTE-correction is turned off). Used to determine which root
        path to use when copying files.
    postflash : bool
        ``True`` if postflash processing is turned on, ``False`` if
        postflash processing is turned off. Used to determine
        which files to copy.
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Retrieving Data ----------')
    logging.info('')

    # Create output directories if they do not exist
    if postflash:
        outdir_list = [paths['main'], paths['shutterA'], paths['shutterB']]
    else:
        outdir_list = [paths['main'], paths['preproc']]
    for outdir in outdir_list:
        if os.path.exists(outdir) == False:
            os.mkdir(outdir, 0o774)
            logging.info('\tCreated Directory {}'.format(outdir))

    if ctecorr == True:
        filetype = 'rac.fits'
    else:
        filetype = 'raw.fits'

    # Get list of files to exclude
    with open(paths['exclude_list'], 'r') as f:
        exclude_list = f.readlines()
    exclude_list = [item.strip() for item in exclude_list]

    # Walk through dark root and pick appropriate files to copy
    logging.info('')
    logging.info('\tScanning {} for data:'.format(paths['ctecorr_dir']))
    logging.info('')

    files_to_copy = []
    files_to_check = glob.glob(os.path.join(paths['ctecorr_dir'], '*_{}'.format(filetype)))
    for file_path in files_to_check:
        if os.path.basename(file_path) not in exclude_list:
            file_to_copy = pick_files(anneal_date, endtime, file_path, postflash, paths)
            if file_to_copy is not None:
                files_to_copy.append(file_to_copy)

    # Copy the files
    logging.info('')
    logging.info('\tCopying the files')
    pool = multiprocessing.Pool(processes=30)
    pool.map(copy_file, files_to_copy)
    pool.close()
    pool.join()
