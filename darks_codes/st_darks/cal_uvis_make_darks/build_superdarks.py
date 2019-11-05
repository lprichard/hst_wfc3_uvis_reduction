"""Create "superdarks" by average-combinging ``nn*blv_tmp.fits`` files.

This module and the functions within creates "superdarks", which are
average-combined ``nn*blv_tmp.fits`` files over 4-day windows.  The
resulting superdark(s) are then normalized to convert units from DN
to e-/s, and their headers are modiefied to reflect the changes made
and prepare them for delivery.

A ``dark_create/`` directory is created and all ``nn*blv_tmp.fits``
files are copied to there.  The ``nn*blv_tmp.fits`` files are then
split and copied into separate ``ith_superdark/`` directories based on
their ``DATE-OBS``; each ``ith_superdark/`` directory contains a
maximum 4-days worth of ``nn*blv_tmp.fits files``.  Some
``ith_superdark/`` directories may contain less than 4-days worth of
data if the processing date happens to be shortly after the end of a
4-day window.

Within each ``ith_superdark/`` directory, the ``nn*blv_tmp.fits`` files
are average-combined using ``numpy`` masked arrays, resulting in a
"superdark", whose name is based on the ``USEAFTER`` of the earliest
input ``nn*blv_tmp.fits files``, and the total number of input
``nn*blv_tmp.fits files``.

Each "superdark" is then "normalized" (i.e. multiplied by the gain
and divided by the exposure time) to convert the units from DN to
e-/s.  A hot pixel threshold is then employed, in which pixels
in the ``SCI`` extensions that fall below 0.015 e-/s are considered
"good" and are set to the frame's median value, and pixels that fall
above the threshold are flagged with 16 in the corresponding ``DQ``
array. The headers of each "superdark" are then modified to reflect the
recent changes and to prepare them for delivery.

The output product are ``d<useafter><num_files>_drk.fits files`` --
average-combined "superdarks", where ``<useafter>`` is the ``DATE-OBS``
of the earliest input dark, and ``<num_files>`` is the number of files
that went into the combination, placed in the appropriate
``dark_create/ith_superdark/`` directory.

Authors
-------

    - Matthew Bourque, 2016
    - John Biretta, 2012
    - Megan Sosey, pre-2012

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import glob
import itertools
import logging
import math
import multiprocessing
import os
import shutil
import time

from astropy.io import fits
import numpy as np

from automated_scripts.cal_uvis_make_darks.combine import combine
from automated_scripts.cal_uvis_make_darks.electrons_normalize import electrons_normalize


def build_superdark_dirs(args):
    """Create superdark directories and copy appropriate
    ``nn*blv_tmp.fits`` files into them.

    Before ``nn*blv_tmp.fits`` files are copied, their median values
    are checked; if the frame's median value is > 9, the file is not
    copied and thus not used in the creation of the superdark.

    Parameters
    ----------
    superdark_filelist : list
        A list of paths to the input ``nn*blv_tmp.fits`` files.
    paths : dict
        A dictonary whose keys are path identifiers and whose
        values are strings containing absolute paths.
    """

    # Parse arguments
    superdark_files = args[0]
    i = args[1]
    paths = args[2]

    # Build superdark directory
    superdark_dir = os.path.join(paths['superdark_create_dir'], '{}th_superdark'.format(str(i + 1)))
    os.mkdir(superdark_dir, 0o774)
    logging.info('')
    logging.info('\tCreated directory {}'.format(superdark_dir))
    logging.info('')
    for frame in superdark_files:
        keep = get_dark_median_stats(frame)
        if keep:
            shutil.copyfile(frame, '{}/{}'.format(superdark_dir, os.path.basename(frame)))
            logging.info('\tCopied {} to {}'.format(frame, superdark_dir))


def get_dark_median_stats(frame):
    """Determine if mean value of dark frame falls below "good"
    threshold of 8 e-/s.

    Dark frames that have median values that exceed this threshold are
    considered bad and are ignored during the creation of superdarks.

    Parameters
    ----------
    frame : str
        The path to the dark image.

    Returns
    -------
    keep : bool
        ``True`` if dark median <= 8 e-/s, ``False`` if dark median >
        8 e-/s.
    """

    # Open image
    hdulist = fits.open(frame)

    # Get the exptime
    exptime = hdulist[0].header['EXPSTART']

    # Calculate dark statistics
    dq_ext1 = hdulist[3].data
    dq_ext4 = hdulist[6].data
    good_data_ext1 = hdulist[1].data[np.where(dq_ext1 == 0)]
    good_data_ext4 = hdulist[4].data[np.where(dq_ext4 == 0)]
    cold_data_ext1 = good_data_ext1[np.where(good_data_ext1 <= 9.)]
    cold_data_ext4 = good_data_ext4[np.where(good_data_ext4 <= 9.)]
    dark_median_ext1 = np.median(cold_data_ext1) * 1.5 * 3600. / exptime
    dark_median_ext4 = np.median(cold_data_ext4) * 1.5 * 3600. / exptime

    # Determine if to keep or toss
    if dark_median_ext1 > 8. or dark_median_ext4 > 8.:
        return False
    else:
        return True


def get_superdark_filelist(frame_list, expstart_list):
    """Return a list of nnblv datasets that occur in a running 4-day
    window.

    Each item in the list is a list of paths to nnblv files that were
    observed +/- 2 days around the increasing expstart index.  The
    expstart + 4 days of the earliest nnblv file is chosen as the
    first window, meaning that the first four days of nnblv files will
    make up the first superdark.  The expstart index is then increased
    by 3 days as to avoid redundant superdarks (i.e. increasing by only
    two days will result in the same 4-day window and thus the same
    filelist). The expstart index is increased by a day as datasets are
    gathered, until the latest expstart - 1 day is reached.  In
    other words, the 4-day window runs from (first observation + 4
    days) to (last observation - 1 day).

    Each item in the list of nnblv datasets will be used to average-
    combine the files within to create superdarks.

    Parameters
    ----------
    frame_list : list
        A list of paths to the nnblv files.
    expstart_list : list
        A list of ``EXPSTART`` times corresponding to the nnblv files.

    Returns
    -------
    superdark_filelist : list
        Each item in the list is a list of nnblv files occuring
        within the particular 4-day window, to be used to create
        superdarks.

    Notes
    -----
    Due to the scheduling of individual dark observations, some cases
    exist in which the ``EXPSTART`` time of the first observation for a
    given 4-day window can be the same as another. This causes
    superdarks with redundant USEAFTERs. This function contains some
    logic that eliminates these redundant datasets.
    """

    # Initializations
    first_expstart = min(expstart_list)
    superdark_filelist = []

    # Build list of files taken within first four days
    first_superdark = []
    for expstart, frame in zip(expstart_list, frame_list):
        if expstart < (first_expstart + 4):
            first_superdark.append(frame)
    superdark_filelist.append(first_superdark)

    # Determine the third day and the last day of observation
    expstart_index = first_expstart + 3
    last_day = max(expstart_list)

    # Build lists of files within running 4-day window
    while expstart_index <= last_day - 1:
        dataset = []
        for expstart, frame in zip(expstart_list, frame_list):
            if (expstart - 2) < expstart_index <= (expstart + 2):
                dataset.append(frame)
        if len(dataset) >= 5:
            superdark_filelist.append(dataset)
        expstart_index += 1

    # Remove redundant filelists
    min_expstart_list, nfiles, good_expstarts, good_nfiles, bad_indices = [], [], [], [], []
    index = 0
    for filelist in superdark_filelist:
        expstarts = [fits.open(item)[0].header['EXPSTART'] for item in filelist]
        min_expstart_list.append(min(expstarts))
        nfiles.append(len(filelist))
    for expstart, nfile in zip(min_expstart_list, nfiles):
        if expstart in good_expstarts:
            nfiles_to_compare = good_nfiles[-1]
            if nfile >= nfiles_to_compare:
                bad_indices.append(index-1)
        else:
            good_expstarts.append(expstart)
            good_nfiles.append(nfile)
        index += 1
    superdark_filelist = [item for i, item in enumerate(superdark_filelist) if i not in bad_indices]

    return superdark_filelist


def get_superdark_filename(nnblv_list, superdark_dir):
    """Return the filename of the superdark based on the observation
    dates and number of files in the nnblv dataset.

    The superdark filename convention is
    ``d<useafter><nfiles><superdark_number>_drk.fits``, where
    ``useafter`` is in the form of ``YYMMDD`` and ``superdark_number``
    is an integer indicator of which 4-day window superdark is being
    created.

    Parameters
    ----------
    nnblv_list : list
        A list of paths to the nnblv files.
    superdark_dir : str
        The path to the superdark directory that holds the nnblv
        files.

    Returns
    -------
    outputname : str
        The filename of the superdark.
    """

    nfiles = len(nnblv_list)
    dateobs_list = [fits.getval(frame, 'DATE-OBS', 0) for frame in nnblv_list]
    useafter = time.strftime('%y%m%d', time.strptime(min(dateobs_list), '%Y-%m-%d'))

    if superdark_dir.endswith('/'):
        directory_indicator = superdark_dir.split('/')[-2].split('th_')[0]
    else:
        directory_indicator = superdark_dir.split('/')[-1].split('th_')[0]

    outputname = 'd{}{}{}_drk.fits'.format(str(useafter), str(nfiles), str(directory_indicator))

    return outputname


def process_superdark(args):
    """The worker function for processing a superdark using
    ``multiprocessing``.  Creates the superdark.

    Parameters
    ----------
    args : list
        A list of the multiprocessing arguments.  The 0th item in the
        list is the ``superdark_dir`` (the directory in which the
        input files exist and the output superdark will be written).
        The 1st element is the ``anneal_date``.  The 2nd item is the
        ``ctecorr`` switch.
    """

    # Parse args
    superdark_dir = args[0]
    anneal_date = args[1]
    ctecorr = args[2]

    nnblv_list = glob.glob(os.path.join(superdark_dir, 'nn*blv_tmp.fits'))
    superdark = get_superdark_filename(nnblv_list, superdark_dir)
    superdark = os.path.join(superdark_dir, superdark)
    combine(nnblv_list, superdark)

    # Normalize by gain and exposure time
    electrons_normalize(superdark)

    # Set everything below hot pixel threshold to median in SCI and ERR externsions
    nfiles = len(nnblv_list)
    set_good_pixels(superdark, nfiles)

    # Update headers
    update_header_for_crds(superdark, nnblv_list, anneal_date, ctecorr)


def set_good_pixels(superdark, nfiles):
    """Set all pixels below the hot pixel threshold to the frame's
    median value.

    The hot pixel threshold is 0.015 e-/sec.

    Parameters
    ----------
    superdark : str
        The path to the superdark to process.
    nfiles : int
        The number of files that went into creating the superdark.
    """

    logging.info('\tSetting "good" pixels to median value of frame.')

    # Open the superdark
    hdulist = fits.open(superdark, mode='update')
    hot_thresh = 0.015

    # Get the data
    ext1 = hdulist[1].data
    ext2 = hdulist[2].data
    ext3 = hdulist[3].data
    ext4 = hdulist[4].data
    ext5 = hdulist[5].data
    ext6 = hdulist[6].data

    # Determine which pixels fall below hot pixel threshold
    index1 = ext1 <= hot_thresh
    index4 = ext4 <= hot_thresh

    # Calculate median and standard deviation
    med_ext1 = np.median(ext1[index1])
    std_ext1 = np.std(ext1[index1])
    med_ext4 = np.median(ext4[index4])
    std_ext4 = np.std(ext4[index4])

    # Set "good" pixels to median value
    ext1[index1] = med_ext1
    ext4[index4] = med_ext4
    ext2[index1] = std_ext1/math.sqrt(nfiles)
    ext5[index4] = std_ext4/math.sqrt(nfiles)

    # Set DQ flags of non-"good" pixels to 16 (hot pixel)
    hotpix_ext1 = np.where(ext1 != med_ext1)
    hotpix_ext4 = np.where(ext4 != med_ext4)
    ext3[hotpix_ext1] = 16
    ext6[hotpix_ext4] = 16

    # Save the image
    hdulist.flush()
    hdulist.close()


def sort_nnblvs(paths):
    """Return lists of ``nn*blv_tmp.fits`` files and expstart times,
    sorted by ``expstart``.

    Parameters
    ----------
    paths : dict
        A dictonary whose keys are path identifiers and whose
        values are strings containing absolute paths.

    Returns
    -------
    frame_list : list
        A list of paths to the nnblv files.
    expstart_list : list
        A list of ``EXPSTART`` times corresponding to the nnblv files.
    """

    nnblv_list = glob.glob(paths['blvs_for_dark_create'])
    expstart_list = [fits.getval(frame, 'EXPSTART', 0) for frame in nnblv_list]
    frame_list = [frame for frame in nnblv_list]
    expstart_list, frame_list = (list(x) for x in zip(*sorted(zip(expstart_list, frame_list))))

    return frame_list, expstart_list


def update_header_for_crds(superdark, nnblv_list, anneal_date, ctecorr):
    """Prepare header of superdark for CRDS delivery.

    CRDS keywords are specified in TIR CDBS 2008-01.

    Parameters
    ----------
    superdark : str
        The path to the superdark to process.
    nnblv_list : list
        The list of nnblv files that went into the creation of the
        superdark.
    anneal_date : str
        The anneal date and time, in the format of
        ``YYYYMMDD-HH:MM:SS``, to be used as the lower date/time
        limit for darks to process.
    ctecorr : bool
        The CTE correction switch.


    Notes
    -----
    Unlike the ``USEAFTERs`` used in other areas of this pipeline,
    the ``USEAFTER`` in the header for CRDS delivery is of the
    format ``yyyy-mm-dd``.
    """

    logging.info('\tUpdating header for CRDS delivery.')

    # Gather list of DATE-OBS from nnblv_list
    dateobs_list = []
    for frame in nnblv_list:
        dateobs = [fits.getval(frame, 'DATE-OBS', 0), fits.getval(frame, 'TIME-OBS', 0)]
        dateobs_list.append(dateobs)
    dateobs_list.sort()

    # Open superdark
    hdulist = fits.open(superdark, mode='update', checksum='remove')

    # Determine PEDIGREE
    mindate = time.strftime("%d/%m/%Y", time.strptime(dateobs_list[0][0], "%Y-%m-%d"))
    maxdate = time.strftime("%d/%m/%Y", time.strptime(dateobs_list[-1][0], "%Y-%m-%d"))
    pedigree = 'INFLIGHT {} {}'.format(mindate, maxdate)

    # Determine USEAFTER.  The USEAFTER is the anneal date for the first superdark,
    # otherwise it is the DATE-OBS of the earliest observation
    if '1th_superdark' == os.path.dirname(superdark).split('/')[-1]:
        useafter = '{} {}'.format(time.strftime("%b %d %Y", time.strptime(str(anneal_date.date()), "%Y-%m-%d")), str(anneal_date.time()))
    else:
        useafter = '{} {}'.format(time.strftime("%b %d %Y", time.strptime(dateobs_list[0][0], "%Y-%m-%d")), (dateobs_list[0][1]))

    # Determine CAL_VER and OPUS_VER
    cal_ver = hdulist[0].header['CAL_VER']
    opus_ver = hdulist[0].header['OPUS_VER']

    # Update the header with modified/delete keywords
    hdulist[0].header['DESCRIP'] = 'DARK created from in-flight WFC3/UVIS frames ----------------------'
    hdulist[0].header['PEDIGREE'] = pedigree
    hdulist[0].header['USEAFTER'] = useafter
    hdulist[0].header['SAMP_SEQ'] = 'NONE'
    hdulist[0].header['SUBTYPE'] = 'FullImag'
    hdulist[0].header['APERTURE'] = 'UVIS'
    hdulist[0].header['CCDAMP'] = 'ANY'
    hdulist[3].header['BUNIT'] = 'UNITLESS'
    hdulist[6].header['BUNIT'] = 'UNITLESS'
    del hdulist[1].header['BUNIT']
    del hdulist[2].header['BUNIT']
    del hdulist[4].header['BUNIT']
    del hdulist[5].header['BUNIT']

    # Update the FILETYPE depending on if CTE correction is applied or not
    if ctecorr:
        hdulist[0].header['FILETYPE'] = 'CTEDARK'
    else:
        hdulist[0].header['FILETYPE'] = 'DARK'

    # Remove unnecessary history added by CTE-correction software and other things, but capture the CTE SOFTWARE VERSION
    if ctecorr:
        cte_ver = hdulist[0].header['HISTORY'][4].split()[1]
        del hdulist[0].header['HISTORY']

    # Update header with HISTORY
    hdulist[0].header.add_history('--------------------------------------------------------------')
    hdulist[0].header.add_history('This reference file was generated by the WFC3/UVIS dark')
    hdulist[0].header.add_history('calibration pipeline.  The following software versions were')
    hdulist[0].header.add_history('used:')
    hdulist[0].header.add_history(' ')
    hdulist[0].header.add_history('    CALWF3: {}'.format(cal_ver))
    hdulist[0].header.add_history('    OPUS: {}'.format(opus_ver))

    if ctecorr:
        hdulist[0].header.add_history('    CTE software: {}'.format(cte_ver))

    hdulist[0].header.add_history(' ')
    hdulist[0].header.add_history('The following WFC3/UVIS dark frames were used in the creation')
    hdulist[0].header.add_history('of this reference file:')
    hdulist[0].header.add_history(' ')

    for nnblv in nnblv_list:
        rootname = '{}'.format(os.path.basename(nnblv).split('_')[0].split('nn')[1])
        hdulist[0].header.add_history(rootname)

    if ctecorr:
        hdulist[0].header.add_history(' ')
        hdulist[0].header.add_history('Individual dark frames were preprocessed by the WFC3/UVIS CTE')
        hdulist[0].header.add_history('correction software (wfc3uv_ctereverse.F) to remove CTE')
        hdulist[0].header.add_history('effects as well as by CALWF3 to perform bias correction,')
        hdulist[0].header.add_history('postflash correction, and cosmic ray flagging.  A custom')
        hdulist[0].header.add_history('CRREJTAB was used during this step using crsigmas=40, 39, 38')
        hdulist[0].header.add_history('and crradius=0.5.  The flagged cosmic rays were then "grown"')
        hdulist[0].header.add_history('with a radius of 2 pixels in order to completely mask the')
        hdulist[0].header.add_history('cosmic rays before identifying hot pixels.')
        hdulist[0].header.add_history(' ')

    else:
        hdulist[0].header.add_history(' ')
        hdulist[0].header.add_history('Individual dark frames were preprocessed by CALWF3 to perform')
        hdulist[0].header.add_history('bias correction, postflash correction, and cosmic ray')
        hdulist[0].header.add_history('flagging.  A custom CRREJTAB was used during this step using')
        hdulist[0].header.add_history('crsigmas=40, 39, 38 and crradius=0.5.  The flagged cosmic rays')
        hdulist[0].header.add_history('were then "grown" with a radius of 2 pixels in order to')
        hdulist[0].header.add_history('completely mask the cosmic rays before identifying hot pixels.')
        hdulist[0].header.add_history(' ')

    hdulist[0].header.add_history('The dark frames were then split into groups of 4-day windows')
    hdulist[0].header.add_history('and combined, resulting in "superdarks" spanning 4 days.')
    hdulist[0].header.add_history('The ERR arrays of the superdarks were generated via standard')
    hdulist[0].header.add_history('error propagation (i.e. summing the SCI values in quadrature).')
    hdulist[0].header.add_history('The superdarks were then normalized (i.e. multiplied by the')
    hdulist[0].header.add_history('gain and divided by the exposure time) to convert units from')
    hdulist[0].header.add_history('DN to e-/s.  A hot pixel threshold was then employed in which')
    hdulist[0].header.add_history('pixels with values greater than 0.015 e-/s are flagged as 16')
    hdulist[0].header.add_history('(i.e. hot pixel) in the DQ array.  Superdarks are restricted')
    hdulist[0].header.add_history('to 4-day windows due to the continusously changing population')
    hdulist[0].header.add_history('of hot pixels (See the WFC3 Instrument Handbook Section 5.4.8')
    hdulist[0].header.add_history('for further details.')
    hdulist[0].header.add_history(' ')

    hdulist[0].header.add_history('4-day superdarks were created surrounding each individual dark')
    hdulist[0].header.add_history('frame as to provide the most accurate measure of hot pixels')
    hdulist[0].header.add_history('surrounding a given observation as possible.  The hot pixel')
    hdulist[0].header.add_history('measurements from each superdark were then imprinted back onto')
    hdulist[0].header.add_history('the DQ arrays of the individual input dark frames.  A larger')
    hdulist[0].header.add_history('set of individual dark frames (of order 100) were then')
    hdulist[0].header.add_history('combined, resulting in a single "masterdark".  Thus, each')
    hdulist[0].header.add_history('"good" pixel (i.e. non-hot pixels) in the masterdark is a')
    hdulist[0].header.add_history('measure of the dark current over an approximate one month')
    hdulist[0].header.add_history('period from a large dataset taken close in time to the current')
    hdulist[0].header.add_history('epoch.')
    hdulist[0].header.add_history(' ')

    hdulist[0].header.add_history('After the creation of the masterdark image, its values were')
    hdulist[0].header.add_history('copied onto the "good" pixels of the 4-day superdarks.  The')
    hdulist[0].header.add_history('resulting image is a reference file containing the "good"')
    hdulist[0].header.add_history('pixels populated from a masterdark image and hot pixel values')
    hdulist[0].header.add_history('from a 4-day superdark.')

    # Remove unnecessary comments and set the reference file creator
    del hdulist[0].header['COMMENT']
    hdulist[0].header.add_comment("= 'Reference file created by C. Martlin.'")

    # EXTVER in DQ arrays should be an integer
    hdulist[3].header['EXTVER'] = int(hdulist[3].header['EXTVER'])
    hdulist[6].header['EXTVER'] = int(hdulist[6].header['EXTVER'])

    # Save file
    hdulist.close()


def build_superdarks_main(anneal_date, ctecorr, paths):
    """The main function of the ``build_superdarks`` module.  See
    module docstrings for further details.

    Parameters
    ----------
    anneal_date : str
        The anneal date and time, in the format of
        ``YYYYMMDD-HH:MM:SS``, to be used as the lower date/time
        limit for darks to process.
    ctecorr : bool
        The CTE correction switch.
    paths : dict
        A dictonary whose keys are path identifiers and whose
        values are strings containing absolute paths.

    Notes
    -----
    ``multiprocessing`` is used to perform superdark creation for
    different epochs concurrently.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Building superdarks ----------')
    logging.info('')

    # Gather and sort list of nnblv frames by EXPSTART
    frame_list, expstart_list = sort_nnblvs(paths)

    # Gather list of running 4-day window nnblv datasets
    superdark_filelist = get_superdark_filelist(frame_list, expstart_list)

    # Copy data into appropriate superdark directories
    i = list(range(len(superdark_filelist)))
    mp_args = list(zip(superdark_filelist, i, itertools.repeat(paths)))
    pool = multiprocessing.Pool(processes=30)
    pool.map(build_superdark_dirs, mp_args)
    pool.close()
    pool.join()

    # Make superdarks
    superdark_dirs = glob.glob(paths['superdark_dirs'])
    mp_args = list(zip(superdark_dirs, itertools.repeat(anneal_date), itertools.repeat(ctecorr)))
    pool = multiprocessing.Pool(processes=30)
    pool.map(process_superdark, mp_args)
    pool.close()
    pool.join()
