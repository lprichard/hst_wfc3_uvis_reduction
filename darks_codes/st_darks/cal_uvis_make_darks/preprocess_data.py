"""Preprocess the darks before superdark and masterdark creation by
performing various calibration steps.

Before the darks can be combined, they must be bias-corrected,
postflash-corrected (if necessary), CR-flagged, and have the CR mask
'grown'. These steps are performed by the various modules that this
module imports from, namely:

- ``prep_header_keys.py`` -- Copies files to a ``crrej`` directory and
  sets the various header keywords parameters before processing.

- ``make_blvs.py`` -- Creates bias-corrected, postflash-corrected (if
  necessary), CR-flagged ``blv_tmp.fits`` files by running the raw
  darks through ``CALWF3``.

- ``remove_badpix.py`` -- Removes bad pixels (i.e. ``DQ`` flags of
  ``4`` and ``512``) from the ``DQ`` arrays of ``blv_tmp.fits`` files.

- ``grow_cr_mask.py`` -- Masks the outer edges of cosmic rays by
  'growing' the CR mask.

See the docstrings of the individual modules for further details.

Authors
-------

    Matthew Bourque, April 2015

Use
---
    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file calibration pipeline.

Notes
-----

    Since each shutter blade requires its own postflash reference file
    (i.e. ``FLSHFILE``), data from each shutter blade must be processed
    independently if postflash correction is needed.  To speed this
    process up, multiprocessing is used to calibrate each shutter
    independently and concurrently.
"""

from automated_scripts.cal_uvis_make_darks import fileIO
from automated_scripts.cal_uvis_make_darks.make_blvs import make_blvs_main
from automated_scripts.cal_uvis_make_darks.remove_badpix import remove_badpix_main
from automated_scripts.cal_uvis_make_darks.grow_cr_mask import grow_cr_mask_main
from automated_scripts.cal_uvis_make_darks.prep_header_keys import prep_header_keys_main


def preprocess_data_main(paths, postflash):
    """The main function of the ``preprocess_data`` module.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    postflash : bool
        ``True`` if postflash processing is turned on, ``False`` if
        postflash processing is turned off.
    """

    # Prep header keywords
    prep_header_keys_main(paths, postflash)

    if postflash:

        # Separate data by shutter to process independentaly
        shutter_dirs = [paths['shutterA'], paths['shutterB']]
        crrej_dirs = [paths['crrej_dirA'], paths['crrej_dirB']]
        for directory, crrej_dir in zip(shutter_dirs, crrej_dirs):

            # Make blv files
            make_blvs_main(directory)
            make_blvs_main(crrej_dir)

            # Remove bad pixels
            remove_badpix_main(directory)
            remove_badpix_main(crrej_dir)

            # Grow CR Mask
            grow_cr_mask_main(crrej_dir)

    else:
        directory = paths['preproc']
        crrej_dir = paths['preproc_crrej']

        # Make blv files
        make_blvs_main(directory)
        make_blvs_main(crrej_dir)

        # Remove bad pixels
        remove_badpix_main(directory)
        remove_badpix_main(crrej_dir)

        # Grow CR Mask
        grow_cr_mask_main(crrej_dir)

    # Move the crgrow files from the two shutters into a common directory.
    fileIO.move_crgrow_files(paths, postflash)
