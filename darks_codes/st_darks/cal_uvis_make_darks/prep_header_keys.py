"""Set various header keywords in preparation for ``CALWF3``.

This module updates/sets several header keywords for the raw dark
files, which must be done prior to the ``CALWF3`` processing/creation
of ``blv_tmp.fits`` files.  Aside from normal calibration switches,
particular attention is directed towards the ``CRREJTAB`` and
``FLSHFILE``.  Custom ``CRREJTABs`` are used (which live in
``/grp/hst/wfc3k/uvis_darks/``) for the purpose of dark current and hot
pixel measurements.  The ``FLSHFILE`` must also be set for postflash
correction (if necessary).  However, as to avoid hardcoding the file in
the script, the ``crds.bestrefs`` routine is used to determine the best
postflash reference file to use.

Authors
-------

    Matthew Bourque, July 2014

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.

Notes
-----

    In order for the ``crds.bestrefs.BestrefsScript`` to run, the
    following environments must be configured:
    ::

        setenv CRDS_PATH /grp/crds/cache
        setenv CRDS_SERVER_URL https://hst-crds.stsci.edu
"""

import glob
import logging
import os
import shutil

from astropy.io import fits
from crds.bestrefs import BestrefsScript


def process_directory(directory, paths, postflash):
    """Process and individual preprocessing directory which contains
    raw files to prepare header keywords for.  This could be an
    individual shutter directory if processing postflashed data.

    Parameters
    ----------
    dictionary : str
        The path to the directory to process.
    paths : dict
        A dictonary whose keys are path identifiers and whose values
        are strings containing absolute paths.
    postflash : bool
        ``True`` if postflash processing is turned on, ``False`` if
        postflash processing is turned off. Used to determine which
        files to copy.
    """

    # Update header keywords for each raw file in shutter directory
    filelist = glob.glob(os.path.join(directory, '*raw.fits'))

    if len(filelist) > 0:
        set_keywords(filelist, paths, postflash)
    else:
        logging.info('')
        logging.info('\tNo data for {}'.format(directory))

    # Create a crrej directory if it doesn't exist
    crrej_path = os.path.join(directory, 'crrej')
    os.mkdir(crrej_path, 0o774)
    logging.info('\tCreated directory {}'.format(crrej_path))

    if len(filelist) > 0:

        # Populate the ccrej directory with raw files
        logging.info('\tCopying files to crrej directory:')
        logging.info('')
        for rawfile in filelist:
            shutil.copy(rawfile, crrej_path)
            logging.info('\t\tCopied {} to {}'.format(rawfile, os.path.join(crrej_path, os.path.basename(rawfile))))
        logging.info('')

        # Update header keywords for raw files in the crrej directory
        crrej_filelist = glob.glob(os.path.join(crrej_path, '*raw.fits'))
        set_keywords(crrej_filelist, paths, postflash)

    else:
        logging.info('\tNo data for {}'.format(directory))


def set_flshfile(image_list):
    """Set the ``FLSHFILE`` header keyword with the best postflash
    reference file.

    Parameters
    ----------
    image_list : list
        The list of absolute paths to the images to update.
    """

    images = ' '.join(image_list)
    bestrefs_arg = "crds.bestrefs --files {} --types FLSHFILE --update-bestrefs --verbosity 0".format(images)
    script = BestrefsScript(argv=bestrefs_arg)
    script.run()


def set_keywords(image_list, paths, postflash):
    """Set various header keywords in preparation for running
    ``CALWF3``.

    Parameters
    ----------
    image_list : list
        The list of absolute paths to the images to update.
    paths : dict
        A dictonary whose keys are path identifiers and whose values
        are strings containing absolute paths.
    postflash : bool
        ``True`` if postflash processing is turned on, ``False`` if
        postflash processing is turned off. Used to determine which
        files to copy.

    Notes
    -----
    Two custom ``CCREJTABs`` are used, one for dark current
    measurements, and one for hot pixel measurements.
    """

    logging.info('\tSetting header keywords for {}'.format(os.path.dirname(image_list[0])))

    num_files = len(image_list)
    for image in image_list:

        # Open image and get the header
        hdulist = fits.open(image, 'update')
        header = hdulist[0].header

        # Set the header keywords
        header['CRSPLIT'] = num_files
        header['DQICORR'] = 'PERFORM'
        header['ATODCORR'] = 'OMIT'
        header['BLEVCORR'] = 'PERFORM'
        header['BIASCORR'] = 'PERFORM'
        header['CRCORR'] = 'PERFORM'
        header['EXPSCORR'] = 'OMIT'
        header['SHADCORR'] = 'OMIT'
        header['DARKCORR'] = 'OMIT'
        header['FLATCORR'] = 'OMIT'
        header['PHOTCORR'] = 'OMIT'
        header['DRIZCORR'] = 'OMIT'

        # Set the PCTECORR if necessary:
        try:
            header['PCTECORR'] = 'OMIT'
        except:
            pass

        # Set the FLSHCORR if necessary:
        if postflash:
            header['FLSHCORR'] = 'PERFORM'
        else:
            header['FLSHCORR'] = 'OMIT'

        # Set CRREJTAB based on if the image is in the crrej directory or not
        if 'crrej' in image:
            header['CRREJTAB'] = paths['hotpix_crrej']
        else:
            header['CRREJTAB'] = paths['dark_crrej']

        # Save changes
        hdulist.flush()
        hdulist.close()

    # Use bestref process to update FLSHFILE with best reference file
    if postflash:
        set_flshfile(image_list)


def prep_header_keys_main(paths, postflash):
    """The main function of the ``prep_header_keys`` module.

    Various header keyword values are set/updated for raw darks located
    in the data preprocessing directory.  ``crrej/`` subdirectories are
    created for each preprocessing, the corresponding raw darks are
    copied into them, and different ``CRREJTAB`` is applied for the
    purposes of measuring hot pixels.

    Parameters
    ----------
    paths : dict
        A dictonary whose keys are path identifiers and whose values
        are strings containing absolute paths.
    postflash : bool
        ``True`` if postflash processing is turned on, ``False`` if
        postflash processing is turned off. Used to determine which
        files to copy.

    Notes
    -----
    If there are no data for a particular preprocessing directory, then
    the process is skiped. Also, if the data are postflashed,
    preprocessing is split into two "shutter directories" (see NOTES in
    ``preprocess_data`` module).
    """

    logging.info('')
    logging.info('')
    logging.info('')
    logging.info('---------- Preparing Header Keywords ----------')
    logging.info('')

    # For postflashed data
    if postflash:
        for shutter_path in [paths['shutterA'], paths['shutterB']]:
            process_directory(shutter_path, paths, postflash)

    # For non-postflashed data
    else:
        process_directory(paths['preproc'], paths, postflash)
