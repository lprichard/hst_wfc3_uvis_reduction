"""Computes dark current and hot pixel statistics on ``blv_tmp.fits`` and
superdark reference files.

The amount of hot pixels (in % of chip) are computed by counting the
number of pixels that have values > 54 electrons/hour in the FITS ``SCI``
extensions of the ``blv_tmp.fits`` and superdark reference files [1].  The dark
current is calculated by computing the median value of all pixels that have
values equal to or below the 54 electrons/hour threshold.

Outputs include four statistic files, two for the number of hot pixels
(one for each chip) and two for the median dark current (one for each
chip):

1. ``hotpix_stat_ext<ext>_<anneal_date>.dat``: hot pixel statistics
   used to generate the hot pixel trending plots, placed in the
   ``/grp/hst/wfc3k/uvis_darks/plots/blvs`` directory for the ``blv_tmp.fits``
   plot, or ``/grp/hst/wfc3k/uvis_darks/plots/superdarks`` for the
   reference file plot.
2. ``midpt_stat_ext<ext>_<anneal_date>.dat``: median dark current
   statistics used to generate the dark current trending plots,
   placed in the ``/grp/hst/wfc3k/uvis_darks/plots/blvs`` directory for the
   ``blv_tmp.fits`` plot, or ``/grp/hst/wfc3k/uvis_darks/plots/superdarks`` for
   the reference file plot.

Authors
-------
    - Jennifer V. Medina, 2018
    - Matthew Bourque, 2013
    - John Biretta, 2012

Use
---
    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.

Notes
-----
[1] This would be equivalent to 9 counts for the un-normalized ``blv_tmp.fits``
files, and 0.015 electrons/s for the normalized superdark reference files (When
calculating using the nominal gain of 1.5 electrons/count for the UVIS detector
and a 900 second exposure time for the individual dark exposures)
"""

import glob
import logging
import os

from astropy.io import fits
import numpy as np


def calc_stats(stat_type, paths, anneal_date, ctecorr):
    """Builds a file within the `blvs` and `superdarks` stats subdirectories
    containing either dark current or hot pixel statistics, based on
    ``stat_type``.

    Creates 2 output statistic files (one under ``blvs`` and one under
    ``superdarks``), removing any already-existing ones, containing either dark
    current statistics (if ``stat_type`` is ``midpt``) or hot pixel statistics
    (if ``stat_type`` is ``hotpix``).

    Parameters
    ----------
    stat_type : str
        Can be either ``midpt`` to calculate dark current statistics
        or ``hotpix`` to calculate hot pixel statistics.
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the pipeline.
    anneal_date : str
        The anneal date, in the format of ``YYYYMMDD``, to be used in
        the filename of the output file.
    ctecorr : bool
        ``True`` (if CTE-correction is turned on) or ``False`` (if
        CTE-correction is turned off). Used to determine the output
        filename.

    Notes
    -----
    Statistic files are replaced if they already exist so that the
    files are essentially updated, in case the pipeline is run multiple
    times throughout an anneal cycle.
    """

    # The process is performed twice:
    # 1) Using the blv_tmp.fits files to calculate the statistics
    # 2) Using the superdark files to calculate the statistics
    for data_type in ['blvs', 'superdarks']:

        # Build path to stat files we're about to create
        if ctecorr:
            path_ext1 = os.path.join(paths['plot_'+data_type], '{}/{}_stat_ext1_{}_ctecorr.dat'.format(stat_type, stat_type, anneal_date))
            path_ext4 = os.path.join(paths['plot_'+data_type], '{}/{}_stat_ext4_{}_ctecorr.dat'.format(stat_type, stat_type, anneal_date))
        else:
            path_ext1 = os.path.join(paths['plot_'+data_type], '{}/{}_stat_ext1_{}.dat'.format(stat_type, stat_type, anneal_date))
            path_ext4 = os.path.join(paths['plot_'+data_type], '{}/{}_stat_ext4_{}.dat'.format(stat_type, stat_type, anneal_date))

        # Remove stat files if they already exist
        if os.path.exists(path_ext1):
            os.remove(path_ext1)
        if os.path.exists(path_ext4):
            os.remove(path_ext4)

        # LP added, check that directory actually exists
        if not os.path.exists(os.path.dirname(path_ext1)): 
            os.makedirs(os.path.dirname(path_ext1), 0o774)
            print('\tLP log: Created driectory {}'.format(os.path.dirname(path_ext1)))
        if not os.path.exists(os.path.dirname(path_ext4)): 
            os.makedirs(os.path.dirname(path_ext4), 0o774)
            print('\tLP log: Created driectory {}'.format(os.path.dirname(path_ext4)))

        # Calculate stats
        if stat_type == 'midpt':
            print('\tCalculating dark median statistics.')
            print('')
            print('\t\tFilename\tEXPSTART\tDark median ext 1\tdark median ext4')
            files = glob.glob(paths[data_type+'_for_dark'])
            for frame in files:
                expstart, stats_ext1, stats_ext4 = get_dark_median(frame, data_type)
                filename = os.path.basename(frame)
                print('\t\t{}\t{}\t{}\t{}'.format(filename, str(expstart), str(stats_ext1), str(stats_ext4)))
                write_to_output(filename, expstart, stats_ext1, path_ext1)
                write_to_output(filename, expstart, stats_ext4, path_ext4)
            print('')
            print('\tDark median statistic files written to {}'.format(os.path.join(paths['plot_'+data_type], 'midpt/')))
            print('')

        elif stat_type == 'hotpix':
            print('\tCalculating hot pixel statistics.')
            print('')
            print('\t\tFilename\tEXPSTART\tPercent hotpix ext 1\tPercent hotpix ext4')
            files = glob.glob(paths[data_type+'_for_hotpix'])
            for frame in files:
                expstart, stats_ext1, stats_ext4 = get_hotpix_stats(frame, data_type)
                filename = os.path.basename(frame)
                print('\t\t{}\t{}\t{}\t{}'.format(filename, str(expstart), str(stats_ext1), str(stats_ext4)))
                write_to_output(filename, expstart, stats_ext1, path_ext1)
                write_to_output(filename, expstart, stats_ext4, path_ext4)
            print('')
            print('\tHot pixel statistic files written to {}'.format(os.path.join(paths['plot_'+data_type], 'hotpix/')))
            print('')


def get_dark_median(frame, data_type):
    """Computes the median dark current for each ``SCI`` extension of
    the frame.

    The median value of the dark current is only computed for pixels
    with values less than or equal to the 54 electrons/hour threshold
    (which is converted to 9 counts for the un-normalized blv files and
    0.015 electrons/second for the normalized superdark reference files [1]).

    Parameters
    ----------
    frame : str
        The filename of the image to gather the statistics from.

    Returns
    -------
    expstart : str
        The ``EXPSTART`` time of the image [2], in MJD format.
    dark_median_ext1 : float
        The median dark current for ``SCI`` extension 1.
    dark_median_ext4 : float
        The median dark current for ``SCI`` extension 4.

    Notes
    -----
    [1] The UVIS detector has a nominal gain of ~1.5 electrons/count in each
    amplifier, and the exposure time for each dark is 900 seconds.
    [2] This ``EXPSTART`` is actually the elapsed time since data started being
    taken (MJD 54993, or 11 July, 2009).
    """

    # Open image
    hdulist = fits.open(frame)

    # Determine header information
    exptime = hdulist[0].header['EXPTIME']
    expstart = hdulist[0].header['EXPSTART'] - 54993

    # Calculate dark statistics
    dq_ext1 = hdulist[3].data
    dq_ext4 = hdulist[6].data
    good_data_ext1 = hdulist[1].data[np.where(dq_ext1 == 0)]
    good_data_ext4 = hdulist[4].data[np.where(dq_ext4 == 0)]
    if data_type == 'blvs':
        cold_data_ext1 = good_data_ext1[np.where(good_data_ext1 <= 9.)]
        cold_data_ext4 = good_data_ext4[np.where(good_data_ext4 <= 9.)]
        dark_median_ext1 = np.median(cold_data_ext1) * 1.5 * 3600. / exptime
        dark_median_ext4 = np.median(cold_data_ext4) * 1.5 * 3600. / exptime
    elif data_type == 'superdarks':
        cold_data_ext1 = good_data_ext1[np.where(good_data_ext1 <= 0.015)]
        cold_data_ext4 = good_data_ext4[np.where(good_data_ext4 <= 0.015)]
        dark_median_ext1 = np.median(cold_data_ext1) * 3600.
        dark_median_ext4 = np.median(cold_data_ext4) * 3600.

    return expstart, dark_median_ext1, dark_median_ext4


def get_hotpix_stats(frame, data_type):
    """Computes the number of hot pixels for each ``SCI`` extension of
    the frame.

    The number of hot pixels is a count of the pixels that have values
    above the 9 count (~56 electrons/hr, in this case) threshold.
    The number of hot pixels are computed as percentage of chip
    (i.e. (number of hot pixels / number of pixels
    in the chip) * 100).

    Parameters
    ----------
    frame : str
        The filename of the image to gather the statistics from.

    Returns
    -------
    expstart : str
        The ``EXPSTART`` time of the image, in MJD format. (The MJDs since
        the first day data started being collected)
    percentage_hot_ext1 : float
        The percentage of hot pixels for ``SCI`` extension 1.
    dark_median_ext4 : float
        The percentage of hot pixels for ``SCI`` extension 4.
    """

    # Open image
    hdulist = fits.open(frame)

    # Determine header information
    expstart = hdulist[0].header['EXPSTART'] - 54993
    npix = float(np.size(hdulist[1].data))

    # Calculate hotpix statistics
    dq_ext1 = hdulist[3].data
    dq_ext4 = hdulist[6].data
    good_data_ext1 = hdulist[1].data[np.where(dq_ext1 == 0)]
    good_data_ext4 = hdulist[4].data[np.where(dq_ext4 == 0)]
    if data_type == 'blvs':
        hot_data_ext1 = good_data_ext1[np.where(good_data_ext1 > 9.)]
        hot_data_ext4 = good_data_ext4[np.where(good_data_ext4 > 9.)]
    elif data_type == 'superdarks':
        hot_data_ext1 = hdulist[1].data[np.where(dq_ext1 == 16)]
        hot_data_ext4 = hdulist[4].data[np.where(dq_ext4 == 16)]

    hot_npix_ext1 = float(np.size(hot_data_ext1))
    hot_npix_ext4 = float(np.size(hot_data_ext4))
    percentage_hot_ext1 = (hot_npix_ext1 / npix) * 100.
    percentage_hot_ext4 = (hot_npix_ext4 / npix) * 100.

    return expstart, percentage_hot_ext1, percentage_hot_ext4


def write_to_output(filename, expstart, statvalue, outfile):
    """Writes the output file that will be used to make the monitoring bokeh
    plots in monitoring_bokeh_plots.py.

    Parameters
    ----------
    filename : str
        The filename of the image that the statistics come from.
    expstart : str
        The MJD corresponding to the ``EXPSTART`` of the image that
        the statistics come from.
    statvalue : float
        The statistic value.
    outfile : str
        Filename of the output file to write to.
    """

    with open(outfile, 'a') as statfile:
        statfile.write("%s\t%f\t\t%6.4f\n"%(filename, expstart, statvalue))
    os.chmod(outfile, 0o775)


def dark_stats_main(paths, anneal_date, ctecorr):
    """The main function of the ``dark_stats`` module.

    Calls the ``calc_stats`` function for both types of statistics to
    compute (``midpt`` for median dark current and ``hotpix`` for the
    number of hot pixels).

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    anneal_date : str
        The anneal date, in the format of ``YYYYMMDD``, to be used in
        the filename of the output file.
    ctecorr : bool
        ``True`` (if CTE-correction is turned on) or ``False`` (if
        CTE-correction is turned off).
    """

    print('')
    print('')
    print('')
    print('---------- Calculating statistics ----------')
    print('')

    calc_stats('midpt', paths, anneal_date, ctecorr)
    calc_stats('hotpix', paths, anneal_date, ctecorr)
