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
e-/s.  

*****************************************************************
LP edit: the ``fitpix`` flag is set to ``True`` if provided and False otherwise.
It is advised in this version of the code. Instead of using a constant
0.015 e-/s threshold, instead, the number of hot pixels are measured closest 
to the read out with this threshold. Then for each bundle of 50 rows, a 
threshold is found that matches this average number of hot pixels per row. 
What results is a more complete identification of hot pixels that remains 
roughly constant all the way through the chip, as is expected. All pixels that 
remain below the adaptive threshold found as a function of row number are 
considered good and are set to the frame's median value, and pixels 
that fall above the threshold are flagged with 16 in the corresponding 
``DQ`` array.

If not set, i.e. ``False``, the original St method is used:
A hot pixel threshold is then employed, in which pixels
in the ``SCI`` extensions that fall below 0.015 e-/s are considered
"good" and are set to the frame's median value, and pixels that fall
above the threshold are flagged with 16 in the corresponding ``DQ``
array. WARNING, more hot pixels are found close to the read out than 
further away from it, often by a factor of >60%.
*****************************************************************

The headers of each "superdark" are then modified to reflect the
recent changes and to prepare them for delivery.

The output product are ``d<useafter><num_files>_drk.fits files`` --
average-combined "superdarks", where ``<useafter>`` is the ``DATE-OBS``
of the earliest input dark, and ``<num_files>`` is the number of files
that went into the combination, placed in the appropriate
``dark_create/ith_superdark/`` directory.

Authors
-------

    - Edited by Laura Prichard, 2019
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
from pdb import set_trace as st

from astropy.io import fits
import numpy as np

# LP edited module locations: put the directory above cal_uvis_make_darks/ in PYTHONPATH, 
# e.g., export PYTHONPATH="/user/lprichard/darks_red_ext/lp_darks:$PYTHONPATH"
from cal_uvis_make_darks.combine import combine
from cal_uvis_make_darks.electrons_normalize import electrons_normalize

# LP added modules
from astropy.stats import gaussian_fwhm_to_sigma, SigmaClip, sigma_clipped_stats
from photutils import detect_sources, detect_threshold, Background2D, MedianBackground, make_source_mask
import datetime
import pandas as pd


def smooth(y, box_pts):
    '''LP added function
    Smooths a linear array.

    Parameters
    ----------
    y : arr
        Linear array to be smoothed.
    box_pts : int
        Degree of smoothing.

    Returns
    -------
    y_smooth : arr
        Linear array of smoothed data.
    '''
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


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
    print('')
    print('\tCreated directory {}'.format(superdark_dir))
    print('')
    for frame in superdark_files:
        keep = get_dark_median_stats(frame)
        if keep:
            shutil.copyfile(frame, '{}/{}'.format(superdark_dir, os.path.basename(frame)))
            print('\tCopied {} to {}'.format(frame, superdark_dir))

    # LP added
    print('LP log: finished copying nnblv temps to superdark_dirs')

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

        LP added the 4th item, the ``fitpix`` switch, whether a function 
        is fitted to find a roughly constant value of hot pixels accross 
        the chip by adapting the hot pixel threshold by row number.
    """
    # LP checks
    now = datetime.datetime.now()
    print('LP log: starting build_superdarks.process_superdark at ', now.strftime("%Y-%m-%d %H:%M"))

    # Parse args
    superdark_dir = args[0]
    anneal_date = args[1]
    ctecorr = args[2]
    fitpix = args[3]    #LP added fitpix 

    nnblv_list = glob.glob(os.path.join(superdark_dir, 'nn*blv_tmp.fits'))
    superdark = get_superdark_filename(nnblv_list, superdark_dir)
    superdark = os.path.join(superdark_dir, superdark)

    print('LP log: combining nnblvs in: ', superdark_dir)
    combine(nnblv_list, superdark)

    # Normalize by gain and exposure time
    print('LP log: normalising electrons in: ', superdark_dir)
    electrons_normalize(superdark)

    # LP added - save output superdark combined from nnblvs before it's overwritten
    scomb_file = os.path.join(superdark_dir, 'superdark_nnblvcombined.fits')
    print('Copying combined nnblv superdark to: ', scomb_file)
    shutil.copyfile(superdark, scomb_file)

    # Set everything below hot pixel threshold to median in SCI and ERR externsions
    nfiles = len(nnblv_list)
    set_good_pixels(superdark, nfiles, fitpix)     #LP added fitpix

    # Update headers
    update_header_for_crds(superdark, nnblv_list, anneal_date, ctecorr)


def calc_hot_thresh(indata, step_size=50, inc=0.0001, lim='avg', ro_tot=7000, ro_avg=140, ht_thresh=0.015):
    '''LP added function

    Parameters
    ----------
    indata : arr
        Data array to find hot pixels in and calculate thresholds to find those pixels.
    step_size : int
        Number of rows to find hot pixels in at a time.
    inc : float
        Increment to decrease the hot pixel threshold by in order to increase the number of hot pixels.
    lim : str
        The limit to use to find hot pixels, eithere the total number withing the 50 row steps ``tot`` 
        or average number per row ``avg``. The average is more stable to outliers
    ro_tot : int
        The total value of hot pixels close to the read out per step_size.
    ro_avg : int
        The average number of hot pixels per row close the the read out. 
    ht_thresh : float
        The input starting hot pixel threshold to work downwards from to increase the number of 
        hotpixels found.

    Returns
    -------
    df : pandas data frame
        A table of starting row number (``row``), hot pixel threshold for the set of 
        row numbers starting at row (``hot_thresh``), the total number of hotpixels 
        in that group of rows with values above the hot pixel threshold (``tot_hotpix``),
        the average number of hotpixels per row in that group of rows with values above 
        the hot pixel threshold (``avg_hotpix``).

    '''
    # There are 2051 rows, the code runs from row 0 to 2050
    steps = np.arange(0, indata.shape[0], step_size)

    # Storing values
    hot_threshs = []
    rows=[]
    totals=[]
    averages=[]

    # Loop over each 50 rows to determine number of hot pixels
    # Lowers the threshold to match the edge hot pixel values then 
    # returns the new threshold to meet this count
    for step in steps[:-1]:
        # Getting chunk of 50 rows, or whatever the step_size values is
        dat = indata[step:step+step_size]

        # Setting starting hot pixel threshold
        ht = ht_thresh

        # Determining number of hotpixels at the starting threshold
        shp, ahp = count_hp(dat, htpix_thresh=ht)
        print('***************************************')
        print('For rows {} to {}'.format(step, step+step_size))
        print('Total no. of hotpix:', shp)
        print('Average no. of hotpix:', ahp)

        # Reducing the threshold to increase the number of hotpixels 
        # found if less than the amount close to the readout
        if lim=='tot':
            while shp<ro_tot:
                # print('+++++++++++++++++++++++++++++++++++++++++++')
                # print('Iterating down threshold value')
                ht -= inc
                # print('Hot threshold:', ht)
                shp, ahp = count_hp(dat, htpix_thresh=ht)
                # print('Total no. of hotpix:', shp)
                # print('Average no. of hotpix:', ahp)
        elif lim=='avg':
            while ahp<ro_avg:
                # print('+++++++++++++++++++++++++++++++++++++++++++')
                # print('Iterating down threshold value')
                ht -= inc
                # print('Hot threshold:', ht)
                shp, ahp = count_hp(dat, htpix_thresh=ht)
                # print('Total no. of hotpix:', shp)
                # print('Average no. of hotpix:', ahp)

        # If the hot threshold has gone down, set it to the one before it 
        # went above the maximum no, hot pixels at the edge closest to the 
        # read out, then store the final values
        if ht<0.015:
            ht += inc
            shp, ahp = count_hp(dat, htpix_thresh=ht)


        print('FINAL hot pixel threshold:', ht) 
        print('FINAL Total no. of hotpix:', shp)
        print('FINAL Average no. of hotpix:', ahp)
        print('***************************************')

        # Store values
        rows.append(step)
        hot_threshs.append(ht)
        totals.append(shp)
        averages.append(ahp)

    # Place in a data frame depending on whether it is a total measure of avegrae measure
    df = pd.DataFrame({})
    df['row']=rows
    df['hot_thresh']=hot_threshs
    df['tot_hotpix']=totals
    df['avg_hotpix']=averages
    print(df)
    
    return df


def count_hp(dat, htpix_thresh=0.015):
    '''LP added function:
    Counts the number of hotpixels in a given data set and a provided threshold.
    Then calculates and returns the total number and average number per row number.

    Parameters
    ----------
    dat : arr
        Data array to find hot pixels in
    htpix_thresh : float
        Pixel vaue to be considered hot, default is 0.015 which is the ST standard.

    Returns
    -------
    shp : int
        The sum of all the hot pixels within the given data set.
    ahp : int
        The median per row number of hot pixels within the given data set.
    '''
    
    # Determine which pixels fall below hot pixel threshold
    hotpix = dat > htpix_thresh

    nhp=[]
    # Count the number of hotpix as a function of row number
    for i in range(dat.shape[0]):
        nhp.append(len(np.where(hotpix[i]==True)[0]))
    
    # Sum the no. hotpix 
    shp = np.nansum(nhp)
    # Average no. hotpix across the rows
    ahp = np.nanmedian(nhp)
    
    return shp, ahp


def set_good_pixels(superdark, nfiles, fitpix):      #LP added fitpix
    """Set all pixels below the hot pixel threshold to the frame's
    median value.

    The hot pixel threshold is 0.015 e-/sec.

    Parameters
    ----------
    superdark : str
        The path to the superdark to process.
    nfiles : int
        The number of files that went into creating the superdark.
    LP added:
    fitpix : bool
        Turns on (``True``) or off (``False``) the option to create a hot 
        pixel threshold function by row number to find the number of hot pixels 
        to match that close to the read out, i.e. a ~constant number of hot 
        pixels across the chip as expected. Otherwise a constant value is used.

    """

    print('\tSetting "good" pixels to median value of frame.')

    # Open the superdark
    hdulist = fits.open(superdark, mode='update')
    hot_thresh = 0.015   #The ST hot pixel threshold value that is the starting value for the fitpix method

    # Get the data
    ext1 = hdulist[1].data
    ext2 = hdulist[2].data
    ext3 = hdulist[3].data
    ext4 = hdulist[4].data
    ext5 = hdulist[5].data
    ext6 = hdulist[6].data

    # ------------------------------------------------------------------------
    # LP edit, new method of finding pixels that aren't hot (index1/4) 
    if fitpix==True:

        print("")
        print("")
        print("")
        print("---------------- Finding hot pixels ----------------")
        print("")
        # Calculating the max no. hot pixels at the read out (ro), 
        # The total (ro_tot) and average (ro_avg) hot pixels per row
        # For a chunk of 50 rows at either side (there are 2051 rows), 
        # using the tested ST threshold, and the average number
        # There are 2051 rows in ext1, sszie sets the number of rows to consider at a time
        ssize=50    #This can be changed if required
        steps = np.arange(0, ext1.shape[0], ssize)

        print("LP log: finding hot pixel threshold for {} row chunks to match the average no. hot pixels close to read out".format(ssize))
        print("LP log: creating hot pixel threshold as a function of row number")
        print("LP log: identifying hot pixels according to the threshold function, mapping onto superdark")
        print("LP log: for ", superdark)

        # Extension 1 maximum values, first ~50 rows
        dat1 = ext1[steps[0]:steps[1]]
        ro_tot1, ro_avg1 = count_hp(dat1, htpix_thresh=hot_thresh)
        print('LP log: Total no. of hotpix rows {}-{} in extension 1: {}'.format(steps[0], steps[1], ro_tot1))
        print('LP log: Average no. of hotpix rows {}-{} in extension 1: {}'.format(steps[0], steps[1], ro_avg1))

        # Extension 4 maximum values, last ~50 rows
        dat4 = ext4[steps[-2]:steps[-1]]
        ro_tot4, ro_avg4 = count_hp(dat4, htpix_thresh=hot_thresh)
        print('LP log: Total no. of hotpix rows {}-{} in extension 4: {}'.format(steps[-2], steps[-1], ro_tot4))
        print('LP log: Average no. of hotpix rows {}-{} in extension 4: {}'.format(steps[-2], steps[-1], ro_avg4))

        # Calculating the hot pixel threshold per 50 rows for each extension
        # Hot pixel threshold increment
        inc = 0.0001
        # Determines whether the limit is the average per row or total numver
        lim = 'avg'
        # Determines number of rows to measure at a time
        df1 = calc_hot_thresh(ext1, step_size=ssize, inc=inc, lim=lim, ro_tot=ro_tot1, ro_avg=ro_avg1, ht_thresh=hot_thresh)
        df4 = calc_hot_thresh(ext4, step_size=ssize, inc=inc, lim=lim, ro_tot=ro_tot4, ro_avg=ro_avg4, ht_thresh=hot_thresh)

        print('LP log: Median hot pixels per row for extension 1:', np.nanmedian(df1['avg_hotpix']))
        print('LP log: Median hot pixels per row for extension 4:', np.nanmedian(df4['avg_hotpix']))

        # Calculating third order polynomial fits to the hot pixel thresholds for every 50 rows
        # Extension 1
        fit1 = np.polyfit(df1['row']+(0.5*ssize), df1['hot_thresh'], 3)
        p1 = np.poly1d(fit1)

        # Extension 4
        fit4 = np.polyfit(df4['row']+(0.5*ssize), df4['hot_thresh'], 3)
        p4 = np.poly1d(fit4)

        # Creating boolean index arrays for good pixels
        index1 = np.empty_like(ext1, dtype=bool)
        index4 = np.empty_like(ext4, dtype=bool)

        # Looping over each 50 rows, finding the average threshold from the polynomial function
        # Identifying hot pixels at that threshold, storing in the good pixel array (index1/4)
        for step in steps[:-1]:
            # For the final step which has one more pixel
            if step==steps[-2]:
                up = ssize + (ext1.shape[0] % ssize)
            else:
                up = ssize

            # Define rows over which the threshold function will be averaged over
            rows=np.arange(step, step+up, 1)

            # Set the hot pixel threshold as the average over the 50 row chunk
            hot_thresh1 = np.nanmedian(p1(rows))
            hot_thresh4 = np.nanmedian(p4(rows))
            print('Rows {} to {}:'.format(rows[0], rows[-1]))
            print('LP log: Hot pixel threshold extension 1:', hot_thresh1)
            print('LP log: Hot pixel threshold extension 4:', hot_thresh4)

            # Identify the good pixels in the arrray as below the threshold for each 50 rows
            index1[step:step+up] = ext1[step:step+up] <= hot_thresh1
            index4[step:step+up] = ext4[step:step+up] <= hot_thresh4

    else:

        print("")
        print("")
        print("")
        print("---------------- Finding hot pixels ----------------")
        print("")
        print("LP log: using constant hotpix threshold ", hot_thresh)

        # Determine which pixels fall below hot pixel threshold
        index1 = ext1 <= hot_thresh
        index4 = ext4 <= hot_thresh

    # ------------------------------------------------------------------------

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
    
    print("LP log: completed build_superdarks.set_good_pixels for ", superdark)
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

    print('\tUpdating header for CRDS delivery.')

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


def build_superdarks_main(anneal_date, ctecorr, paths, fitpix):  #LP added fitpix to build_superdarks_main
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
    LP added:
    fitpix : bool
        Turns on (``True``) or off (``False``) the option to create a hot 
        pixel threshold function by row number to find the number of hot pixels 
        to match that close to the read out, i.e. a ~constant number of hot 
        pixels across the chip as expected. Otherwise a constant value is used.


    Notes
    -----
    ``multiprocessing`` is used to perform superdark creation for
    different epochs concurrently.
    """

    print('')
    print('')
    print('')
    print('---------- Building superdarks ----------')
    print('')

    # Gather and sort list of nnblv frames by EXPSTART
    frame_list, expstart_list = sort_nnblvs(paths)

    # Gather list of running 4-day window nnblv datasets
    superdark_filelist = get_superdark_filelist(frame_list, expstart_list)

    # Copy data into appropriate superdark directories
    i = list(range(len(superdark_filelist)))
    mp_args = list(zip(superdark_filelist, i, itertools.repeat(paths)))
    # LP halved the number of processes from 30 to 15 due to processign capacity
    pool = multiprocessing.Pool(processes=15)
    pool.map(build_superdark_dirs, mp_args)
    pool.close()
    pool.join()

    # Make superdarks
    superdark_dirs = glob.glob(paths['superdark_dirs'])
    mp_args = list(zip(superdark_dirs, itertools.repeat(anneal_date), itertools.repeat(ctecorr), itertools.repeat(fitpix)))
    
    # LP halved the number of processes from 30 to 15 due to processign capacity
    pool = multiprocessing.Pool(processes=15)

    # LP added, 
    print('LP log: running build_superdark.process_superdark')
    pool.map(process_superdark, mp_args)
    pool.close()
    pool.join()
