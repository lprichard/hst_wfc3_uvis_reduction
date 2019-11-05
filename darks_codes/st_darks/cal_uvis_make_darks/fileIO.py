"""Handles various file copying and moving between various directories.

This module contains several functions that copy files, move files,
set file permissions, create directories, etc.

Authors
-------

    Matthew Bourque, 2013
    Catherine Martlin, 2018
    Jennifer V. Medina, 2018

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import glob
import logging
import os
import shutil


def make_path_dict(parent_dir):
    """Builds a dictionary containing absolute paths to be used
    throughout the pipeline.

    This function serves as a single place for hard-coded paths that
    are needed for file operations throughout the pipeline.  Some of
    paths are actually stings to be used in ``glob.glob()`` statements,
    (i.e. they contain wildcards).

    Parameters
    ----------
    parent_dir : str
        The absolute path of the main directory as determined by
        the ``-d --outdir`` argument of the ``cal_uvis_make_darks``
        module.

    Returns
    -------
    paths : dict
        A dictonary whose keys are path identifiers and whose
        values are strings containing absolute paths.
    """

    paths = {}

    # Hard-coded paths
    paths['exclude_list'] = '/grp/hst/wfc3k/uvis_darks/exclude.list'
    paths['dark_crrej'] = '/grp/hst/wfc3k/uvis_darks/crr_for_dark.fits'
    paths['hotpix_crrej'] = '/grp/hst/wfc3k/uvis_darks/crr_for_hotpix.fits'
    paths['masterdark_pool'] = '/grp/hst/wfc3k/uvis_darks/masterdarks/'
    paths['plot_blvs'] = '/grp/hst/wfc3k/uvis_darks/plots/new_algorithm/blvs/'
    paths['plot_superdarks'] = '/grp/hst/wfc3k/uvis_darks/plots/new_algorithm/superdarks/'
    paths['daily_outputs'] = '/grp/hst/wfc3a/automated_outputs/cal_uvis_make_darks/daily_outputs/'
    paths['ctecorr_dir'] = '/grp/hst/wfc3q/ctecorr_darks/'
    paths['tmp_dir'] = '/grp/hst/wfc3k/uvis_darks/tmp/'

    # Variable paths
    paths['main'] = parent_dir
    paths['preproc'] = os.path.join(paths['main'], 'preproc')
    paths['preproc_crrej'] = os.path.join(paths['preproc'], 'crrej')
    paths['shutterA'] = os.path.join(paths['main'], 'a')
    paths['shutterB'] = os.path.join(paths['main'], 'b')
    paths['crrej_dirA'] = os.path.join(paths['shutterA'], 'crrej')
    paths['crrej_dirB'] = os.path.join(paths['shutterB'], 'crrej')
    paths['superdark_create_dir'] = os.path.join(paths['main'], 'superdark_create')
    paths['masterdark_create_dir'] = os.path.join(paths['main'], 'masterdark_create')
    paths['deliver_dir'] = os.path.join(paths['main'], 'deliver')
    paths['report_dir'] = os.path.join(paths['main'], 'report')
    paths['automated_outputs'] = os.path.join('/grp/hst/wfc3a/automated_outputs/cal_uvis_make_darks/', paths['main'].split('/')[-1])

    # Paths for glob statements
    paths['blvs_for_dark_create'] = os.path.join(paths['superdark_create_dir'], 'nn*blv_tmp.fits')
    paths['superdark_dirs'] = os.path.join(paths['superdark_create_dir'], '*th_superdark/')
    paths['blvs_for_dark'] = os.path.join(paths['main'], '*blv_tmp.fits')
    paths['blvs_for_hotpix'] = os.path.join(paths['superdark_create_dir'], 'nn*blv_tmp.fits')
    paths['superdarks_for_dark'] = os.path.join(paths['deliver_dir'], '*d??.fits')
    paths['superdarks_for_hotpix'] = os.path.join(paths['deliver_dir'], '*d??.fits')
    paths['darks_to_deliver'] = os.path.join(paths['deliver_dir'], '*d??.fits')  # fix redundancy later
    paths['report_files'] = os.path.join(paths['report_dir'], '*')
    paths['daily_plots'] = os.path.join(paths['report_dir'], '*.png')
    paths['daily_plots_bokeh'] = os.path.join(paths['report_dir'], '*.html')

    return paths


def copy_masterdark(paths):
    """Copies the masterdark to the "masterdark pool", an external
    directory that holds masterdarks from other anneal cycles.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    """

    masterdark = glob.glob(os.path.join(paths['masterdark_create_dir'], 'masterdark_*.fits'))

    if len(masterdark) > 1:
        logging.info('\t There are multiple masterdarks in this directory.')
        for mdark in masterdark:
            mdark_name = mdark.split("/")[-1]
            masterdark_dst = os.path.join(paths['masterdark_pool'], mdark_name)
            if os.path.exists(masterdark_dst):
                os.remove(masterdark_dst)
                logging.info('\t {} was removed in order to replace it.'.format(masterdark_dst))
            shutil.copyfile(mdark, masterdark_dst)
            logging.info('\tCopied {} to {}'.format(mdark, masterdark_dst))
    else:
        mdark_name = masterdark[0].split("/")[-1]
        masterdark_dst = os.path.join(paths['masterdark_pool'], os.path.basename(mdark_name))
        shutil.copyfile(masterdark, masterdark_dst)
        logging.info('\tCopied {} to {}'.format(masterdark, masterdark_dst))


def move_for_stats(paths, postflash):
    """Moves ``blv_tmp.fits`` and ``th_superdark.fits`` files from their
    respective preprocessing directories back to the main directory.

    The ``blv_tmp.fits`` files in the preprocessing directories are
    moved back into the main directory and are reunited for future
    processing.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    postflash : bool
        ``True`` if postflash processing is turned on, ``False`` if
        postflash processing is turned off.
    """

    blv_files = []
    if postflash:
        for shutter in ['a', 'b']:
            blv_files.extend(glob.glob(os.path.join(paths['main'], '{}/*blv_tmp.fits'.format(shutter))))
    else:
        blv_files.extend(glob.glob(os.path.join(paths['preproc'], '*blv_tmp.fits')))

    for blv_file in blv_files:
        shutil.move(blv_file, paths['main'])


def move_crgrow_files(paths, postflash):
    """Moves newly-created CR-grown files (i.e. ``nn*blv_tmp.fits``
    files) from each shutter directory to a common ``dark_create/``
    directory.

    The ``nn*blv_tmp.fits`` files are moved to a common directory in
    order for future image combination.

    Parameters
    ----------
    paths : dict
        A dictonary whose keys are path identifiers and whose values
        are strings containing absolute paths.
    postflash : bool
        ``True`` if postflash processing is turned on, ``False`` if
        postflash processing is turned off.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Preparing superdark directory ----------')
    logging.info('')

    os.mkdir(paths['superdark_create_dir'], 0o774)
    logging.info('\tCreated directory {}'.format(paths['superdark_create_dir']))
    logging.info('')

    nnblvs = []
    if postflash:
        for shutter in ['a', 'b']:
            nnblvs.extend(glob.glob(os.path.join(paths['main'], '{}/crrej/crgrow/nn*.fits'.format(shutter))))
    else:
        nnblvs.extend(glob.glob(os.path.join(paths['preproc'], 'crrej/crgrow/nn*.fits')))

    if len(nnblvs) > 0:
        for nnblv in nnblvs:
            shutil.move(nnblv, paths['superdark_create_dir'])
            logging.info('\tMoved {} to {}'.format(nnblv, paths['superdark_create_dir']))


def prep_delivery(paths, ctecorr):
    """Copies masterdark reference files to a ``deliver/`` directory.

    The masterdark reference files are placed in a separate directory
    for easy delivery.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    ctecorr : bool
        ``True`` (if CTE-correction is turned on) or ``False`` (if
        CTE-correction is turned off).

    Notes
    -----
    If ``ctecorr`` switch is on, the ``_drk.fits`` files are renamed to
    ``_drc.fits`` to indicate that cte correction was performed.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Preparing deliver directory ----------')
    logging.info('')

    os.mkdir(paths['deliver_dir'], 0o774)
    logging.info('Created directory {}'.format(paths['deliver_dir']))
    logging.info('')

    masterdarks = glob.glob(os.path.join(paths['masterdark_create_dir'], 'd*drk.fits'))

    for masterdark in masterdarks:
        masterdark_filename = os.path.basename(masterdark)
        if ctecorr:
            masterdark_filename = masterdark_filename.replace('drk', 'dkc')
        masterdark_deliver_location = os.path.join(paths['deliver_dir'], masterdark_filename)
        shutil.copyfile(masterdark, masterdark_deliver_location)
        logging.info('\tCopied {} to {}'.format(masterdark, masterdark_deliver_location))


def send_to_automated_outputs(paths):
    """Copies various output products to the ``automated_outputs/``
    directory.

    The output products copied include the contents of the ``report/``
    directory and are copied to
    ``/grp/hst/wfc3a/automated_outputs/cal_uvis_make/darks/<outdir>/``.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Copying output products to automated outputs ----------')
    logging.info('')

    os.mkdir(paths['automated_outputs'], 0o774)
    logging.info('\tCreated directory {}'.format(paths['automated_outputs']))
    logging.info('')

    report_files = glob.glob(paths['report_files'])

    for report_file in report_files:
        shutil.copyfile(report_file, os.path.join(paths['automated_outputs'], os.path.basename(report_file)))
        logging.info('\tCopied {} to {}'.format(report_file, paths['automated_outputs']))

    logging.info('')


def send_to_daily_outputs(paths):
    """Copies specific output products to the ``daily_outputs/``
    directory.

    The output products copied include the hot pixel and median dark
    current trending plots, which are copied to
    ``/grp/hst/wfc3a/automated_outputs/cal_uvis_make_darks/daily_outputs/``
    These output products are meant to be viewed daily as part of the
    WFC3 Quicklook operations in order to monitor the hot pixel and
    dark current trends.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Copying output products to daily outputs ----------')
    logging.info('')

    # Remove what is already in the daily_outputs folder
    existing_files = glob.glob(os.path.join(paths['daily_outputs'], '*'))
    for existing_file in existing_files:
        os.remove(existing_file)

    # Copy updated plots to daily_outputs folder
    plots_to_copy = glob.glob(paths['daily_plots'])
    plots_to_copy_bokeh = glob.glob(paths['daily_plots_bokeh'])
    for plot_to_copy in plots_to_copy:
        dst = os.path.join(paths['daily_outputs'], os.path.basename(plot_to_copy))
        shutil.copy(plot_to_copy, dst)
        logging.info('\tCopied {} to {}'.format(plot_to_copy, dst))

    for plot_to_copy_bokeh in plots_to_copy_bokeh:
        dst_bokeh = os.path.join(paths['daily_outputs'], os.path.basename(plot_to_copy_bokeh))
        shutil.copy(plot_to_copy_bokeh, dst_bokeh)
        logging.info('\tCopied {} to {}'.format(plot_to_copy_bokeh, dst_bokeh))

    logging.info('')


def set_permissions(paths):
    """Sets the permissions of all files and directories under the main
    directory and the ``automated_outputs`` directory to ``rwxrwxr--``.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Setting permissions ----------')
    logging.info('')

    for directory in [paths['main'], paths['automated_outputs']]:
        os.chmod(directory, 0o774)
        for root, subdirs, files in os.walk(directory):
            for subdir in subdirs:
                os.chmod(os.path.join(root, subdir), 0o774)
            for name in files:
                os.chmod(os.path.join(root, name), 0o774)
