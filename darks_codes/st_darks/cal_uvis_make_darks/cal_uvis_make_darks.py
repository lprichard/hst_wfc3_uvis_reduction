#! /usr/bin/env python

"""Generates WFC3 UVIS dark reference files.

This script serves as a pipeline to create WFC3 UVIS dark reference
files and is a wrapper around several modules that perform subtasks of
the dark reference file creation algorithm.  UVIS dark reference files
are created and delivered to the Calibration Database System (CDBS) on
a weekly basis in order to provide up-to-date WFC3/UVIS calibration
files.  When executed, this script will retrieve appropriate WFC3/UVIS
DARK images from the Quicklook filesystem, perform calibration tasks
such as header adjustments and image combinations, and create
deliverable dark reference files, accompanied by several plots,
statistic files, and quick-look JPEG files for analysis.  A general
outline of the algorithm as well as quick descriptions of the modules
within the wrapper are given in NOTES.  Please refer to the
documentation of the individual functions/modules for further details.

Execution of this script will create the following file tree:
::

    post-anneal-<anneal_date>_<timestamp>/
        a/
            crrej/
                crgrow/
        b/
            crrej/
                crgrow/
        dark_create/
            1th_superdark/
            2th_superdark/
            ...
            nth_superdark/
        masterdark_create/
            batch1/ *
            batch2/ *
        deliver/
        report/

    * if necessary

The ``a/`` and ``b/`` directories serve as the areas to independentally
process data using shutter blade A and shutter blade B, respectively
(since they use different postflash calibration files). The ``crrej/``
directories serve as the processing area of cosmic ray rejection and
initial run-through of CALWF3.  The ``crgrow/`` directories serve as
the processing area of the growth of the CR mask.  The ``dark_create/``
directory holds a copy of the CR-grown ``nn*blv_tmp.fits`` files and is
used as the processing area of image combination for the building of
superdarks.  Since the data are split up in 4-day intervals, there will
be 1-7 ``*th_superdark/`` directories nested within the
``dark_create/`` directory.  Within these ``*th_superdark/``
directories are the CR-grown images for the particular 4-day
combination.  The ``masterdark_create/`` directory serves as the
processing area for the creation of entire-anneal-cycle masterdarks.
The ``deliver/`` directory holds a copy of the dark reference files to
be delivered.  The ``report/`` directory holds various output products
used for visual inpection, analysis, and quality assurance.

The following output products are produced from the execution of
this script:

    1. JPEG images of all input RAW darks, placed in the ``report/``
       directory
    2. JPEG images of all output dark reference files, placed in the
       ``report/`` directory
    3. ``hotpix_plot_ext<ext>_<current_date>.pdf``: hot pixel trending
       plots, placed in the ``report/`` directory as well as the
       ``/grp/hst/wfc3k/uvis_darks/plots/`` directory
    4. ``midpt_plot_ext<ext>_<current_date>.pdf``: median dark current
       trending plots, placed in the ``report/`` directory as well as
       the ``/grp/hst/wfc3k/uvis_darks/plots/`` directory.
    5. ``hotpix_stat_ext<ext>_<anneal_date>.dat``: hot pixel statistics
       used to generate the hot pixel trending plots, placed in the
       ``/grp/hst/wfc3k/uvis_darks/plots/`` directory.
    6. ``midpt_stat_ext_<ext>_anneal_date>.dat``: median dark current
       statistics used to generate the dark current trending plots,
       placed in the ``/grp/hst/wfc3k/uvis_darks/plots/`` directory.
    7. ``imstat_ext<ext>.dat``: Files containing image statistics of
       the dark reference files, placed in the
       ``report/`` directory.
    8. Intermediate files, such as ``*_blv_tmp.fits files`` (CR-flagged
       bias-corrected darks), ``nn*_blv_tmp.fits`` (CR-grown
       blv_tmp.fits files), association tables, etc.

Authors
-------
    Jennifer V. Medina, 2018
    Matthew Bourque, 2013

Use
---

    This script is intended to be executed weekly via a cron job or via
    the command line as such:
    ::

        python cal_uvis_make_darks.py [-a|--anneal_date] [-e|--endtime]
        [-c|--ctecorr] [-p|--postflash] [-m|--mode]

    ``-a --anneal_date`` (Optional) - The anneal date and time that
    serves as the lower limit for ``DATE-OBS``/``TIME-OBS`` during file
    retrieval. Only images that were observed after ``anneal_date`` are
    retrieved for processing.  The value must be in the format
    ``YYYYMMDD-HH:MM:SS``.  The default value is the latest anneal date
    listed in ``/grp/hst/wfc3b/calibration/history.txt``.

    ``-e --endtime`` (Optional) -- The date and time that serves as the
    upper limit for the ``DATE-OBS``/``TIME-OBS`` during file retrieval.
    Only images that were observed before the endtime are retrieved
    for processing.  The value must be in the format
    ``YYYYMMDD-HH:MM:SS``.  The default value is the current date and
    time.

    ``-c --ctecorr`` (Optional) -- Turns on or off the use of
    CTE-corrected darks during processing.  The value is ``True`` (i.e.
    use CTE-corrected darks) if provided and ``False`` (i.e. don't use
    CTE-corrected darks) if not provided.

    ``-p --postflash`` (Optional) -- Turns on or off the use of
    postflashed darks during processing.  The value is ``True`` (i.e.
    use postflashed darks) if provided and ``False`` (i.e. don't use
    postflashed darks) if not provided.

    ``-m --mode`` (Optional) -- The processing mode. Can either be
    ``dev`` (i.e. development) in which the masterdark for the given
    anneal cycle will be used to replace "good" pixels in the
    superdark, or ``prod`` (i.e. production) to use the previous
    anneal masterdark.

Notes
-----

    This script calls the following modules (in order of execution):

    ``get_anneal_date.py`` - Determines the latest (default) anneal
    date.

    ``FileIO.py`` -- A handler for various file input/output tasks.

    ``copy_darks.py`` -- Retrieves RAW dark files from WFC3 Quicklook
    (if ``ctecorr`` is ``False``) or
    ``/grp/hst/wfc3k/uvis_darks/ctecorr/`` (if ``ctecorr`` is ``True``)
    used to generate UVIS dark reference files. In order to retrieved,
    the raw files must have the following characteristics:
    1. ``DATE-OBS`` and ``TIME-OBS`` between ``anneal_date`` and
    ``enddtime``
    2. ``TARGNAME`` = ``dark``
    3. ``EXPTIME`` = ``900``
    4. ``FLASHLVL`` = ``12``
    Files are then copied to directories based on their ``SHUTRPOS``.

    ``preprocess_data.py`` -- Processes the darks before superdark and
    masterdark creation by performing various calibrations.

    ``build_superdark.py`` -- Combines ``nn*blv_tmp.fits`` files to
    create superdarks.

    ``electrons_normalize.py`` -- Converts the units of the given dark
    reference file from DN to e-/s.

    ``imprint_hotpix.py`` -- Copies the DQ flags of the superdarks onto
    the individual CR-grown darks.

    ``build_masterdark.py`` -- Combines ``nn*blv_tmp.fits`` files to
    create masterdarks.

    ``dark_stats.py`` -- Computes dark current and hot pixel statistics
    on ``blv_tmp.fits`` files.

    ``monitoring_plots.py`` -- Generates hot pixel and dark current
    plots.

    ``make_report.py`` -- Constructs a directory containing jpegs of
    darks, image statistics, and dark median/hot pixel plots.

References
----------

    WFC3/UVIS dark calibration documentation:

    - WFC3 ISR 2016-08: WFC3/UVIS Dark Calibration: Monitoring Results
      and Improvements to Dark Reference Files (M. Bourque, S. Baggett)
    - WFC3 ISR 2014-04: WFC3 Cycle 19 & 20 Dark Calibration: Part I
      (J. Biretta, M. Bourque)
    - WFC3 TIR 2014-01: WFC3 Cycle 19 & 20 Dark Calibration: Part II
      (J. Biretta)
"""

import argparse
import datetime
import logging
import os
import re
import warnings

from pyql.logging.logging_functions import configure_logging
from pyql.logging.logging_functions import log_info
from pyql.logging.logging_functions import log_fail

from automated_scripts.cal_uvis_make_darks.build_masterdark import build_masterdark_main
from automated_scripts.cal_uvis_make_darks.build_superdarks import build_superdarks_main
from automated_scripts.cal_uvis_make_darks.copy_darks import copy_darks_main
from automated_scripts.cal_uvis_make_darks.dark_stats import dark_stats_main
from automated_scripts.cal_uvis_make_darks import fileIO
from automated_scripts.cal_uvis_make_darks.get_anneal_date import get_anneal_date
from automated_scripts.cal_uvis_make_darks.imprint_hotpix import imprint_hotpix_main
from automated_scripts.cal_uvis_make_darks.make_report import make_report_main
#from automated_scripts.cal_uvis_make_darks.monitoring_plots import monitoring_plots_main
from automated_scripts.cal_uvis_make_darks.monitoring_bokeh_plots import monitoring_plots_bokeh_main
from automated_scripts.cal_uvis_make_darks.preprocess_data import preprocess_data_main
from automated_scripts.cal_uvis_make_darks.remove_checksum import remove_checksum_main

warnings.filterwarnings('ignore', category=UserWarning)


@log_fail
@log_info
def process_darks_main(anneal_date, endtime, ctecorr, postflash, mode):
    """The main function of the 11cal_uvis_make_darks11 pipeline.

    This function serves as a wrapper around all other function and
    module calls in the pipeline.  Please see the documentation of the
    individual functions or modules for further details.

    Parameters
    ----------
    anneal_date : str
        The anneal date and time, in the format of
        ``YYYYMMDD-HH:MM:SS``, to be used as the lower date/time
        limit for darks to process.
    endtime : str
        The date and time, in the format of ``YYYYMMDD-HH:MM:SS``,
        to be used as the upper date/time limit for darks to
        process.
    ctecorr : bool
        Turns on (``True``) or off (``False``) the use of
        CTE-corrected darks during processing.
    postflash : bool
        Turns on (``True``) or off (``False``) the processing of
        postflashed data.
    mode : str
        The processing mode. Can either be ``dev``, in which
        masterdarks are made from the supplied anneal, or ``prod``
        (i.e. "production") in which the masterdarks are made from
        the previous anneal's darks.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('Using anneal date: {}'.format(str(anneal_date)))

    # Determine default outdir directory name
    today = datetime.datetime.now()
    today = datetime.datetime.strftime(today, '%Y%m%d')
    anneal_date_str = datetime.datetime.strftime(anneal_date, '%Y%m%d')
    if ctecorr:
        outdir = '/grp/hst/wfc3k/uvis_darks/new_algorithm/post-anneal-{}_{}_ctecorr'.format(anneal_date_str, today)
    else:
        outdir = '/grp/hst/wfc3k/uvis_darks/new_algorithm/post-anneal-{}_{}'.format(anneal_date_str, today)

    # Prepare directories and fill with data
    paths = fileIO.make_path_dict(outdir)
    copy_darks_main(anneal_date, endtime, ctecorr, postflash, paths)

    # Preprocess the data
    preprocess_data_main(paths, postflash)

    # Build superdark
    build_superdarks_main(anneal_date, ctecorr, paths)

    # Imprint hot pixels from superdark back onto nn*blvs
    imprint_hotpix_main(paths)

    # Build masterdark
    build_masterdark_main(paths, str(anneal_date.date()), ctecorr, mode)

    # Copy masterdark to external masterdark directory
    fileIO.copy_masterdark(paths)

    # Copy dark reference files to deliver directory
    fileIO.prep_delivery(paths, ctecorr)

    # Move blv files to main directory for statistics
    fileIO.move_for_stats(paths, postflash)

    # Run statistics
    dark_stats_main(paths, str(anneal_date.date()), ctecorr)

    # Make dark median and hot pixel Bokeh plots
    monitoring_plots_bokeh_main(ctecorr, paths)

    # Make report
    make_report_main(paths, ctecorr)

    # Copy report directory to automated_outputs/
    fileIO.send_to_automated_outputs(paths)

    # Copy daily output plots to daily_outputs/ folder
    fileIO.send_to_daily_outputs(paths)

    # Remove checksum
    remove_checksum_main(paths)

    # Set permissions for files and directories
    fileIO.set_permissions(paths)


def parse_args():
    """Parses command line arguments.

    Returns
    -------
    args : obj
        An ``argparse`` object containing all of the added arguments.
    """

    # Create help string
    anneal_date_help = 'Date of WFC3/UVIS anneal (YYYYMMDD-HH:MM:SS)'
    endtime_help = 'End time of time range (YYYYMMDD-HH:MM:SS).  The default is the current time.'
    ctecorr_help = 'Turn on/off use of CTE corrected darks.  The default is False (i.e. off) if not supplied.'
    postflash_help = 'Turn on/off processing of postflash data.  The default is False (i.e. off) if not supplied.'
    mode_help = ('The processing mode.  Can either be "dev" for masterdark development, or "prod"'
                 ' (i.e. "production") for nominal delivery output products.  The default is "prod".')

    # Add time arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-a --anneal_date',
        dest='anneal_date',
        action='store',
        type=str,
        required=False,
        help=anneal_date_help,
        default=get_anneal_date())
    parser.add_argument('-e --endtime',
        dest='endtime',
        action='store',
        type=str,
        required=False,
        help=endtime_help)

    # Add CTE correction switch argument
    parser.add_argument('-c --ctecorr',
        dest='ctecorr',
        action='store_true',
        help=ctecorr_help)

    # Add postflash argument
    parser.add_argument('-p --postflash',
        dest='postflash',
        action='store_true',
        help=postflash_help)

    # Add mode argument
    parser.add_argument('-m --mode',
        dest='mode',
        action='store',
        type=str,
        required=False,
        help=mode_help,
        default='prod')

    # Set defaults
    parser.set_defaults(ctecorr=False, postflash=False)

    # Parse args
    args = parser.parse_args()

    return args


def test_args(args):
    """Ensures that the command line arguments are of proper format. If
    they are not, an assertion error is raised.

    Parameters
    ----------
    args : obj
        The ``argparse`` object containing the command line arguments.
    """

    # Build regular expression that reflects anneal_date/endtime format
    reg = '((?:(?:[1]{1}\\d{1}\\d{1}\\d{1})|(?:[2]{1}\\d{3}))(?:[0]?[1-9]|'
    reg += '[1][012])(?:(?:[0-2]?\\d{1})|(?:[3][01]{1})))(?![\\d])'
    reg += '(-)((?:(?:[0-1][0-9])|(?:[2][0-3])|(?:[0-9])):(?:[0-5][0-9])'
    reg += '(?::[0-5][0-9])?(?:\\s?(?:am|AM|pm|PM))?)'
    re_test = re.compile(reg, re.IGNORECASE | re.DOTALL)

    # Assert that anneal_date and endtime are of proper format
    assert bool(re_test.search(args.anneal_date)) is True, 'Invalid date format. Format must be (YYYYMMDD-HH:MM:SS)'
    if args.endtime is not None:
        assert bool(re_test.search(args.endtime)) is True, 'Invalid date format. Format must be (YYYYMMDD-HH:MM:SS)'

    # Assert that the ctecorr switch is True or False
    assert args.ctecorr in [True, False], 'Invalid option for ctecorr argument. Valid options are True and False.'

    # Assert that the postflash switch is True or False
    assert args.postflash in [True, False], 'Invalid option for postflash argument. Valid options are True and False.'

    # Assert that the mode is either dev or prod
    assert args.mode in ['dev', 'prod'], 'Invalid option for mode argument. Valid options are "dev" and "prod".'


if __name__ == '__main__':

    module = os.path.basename(__file__).replace('.py', '')
    configure_logging(module)
    logging.getLogger(module)

    args = parse_args()

    # Convert times to datetime objects
    args.anneal_date = datetime.datetime.strptime(args.anneal_date, '%Y%m%d-%H:%M:%S')
    if args.endtime is None:
        args.endtime = datetime.datetime.now()
    else:
        args.endtime = datetime.datetime.strptime(args.endtime, '%Y%m%d-%H:%M:%S')

    process_darks_main(args.anneal_date, args.endtime, args.ctecorr, args.postflash, args.mode)
