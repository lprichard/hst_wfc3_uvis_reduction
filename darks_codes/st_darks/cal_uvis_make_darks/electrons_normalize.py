"""Convert an image from units of ``DN`` to ``electrons/s``.

The given image is multiplied by the gain and normalized by exposure
time to convert its units from ``DN`` to ``e-/s``.  The ``SCI`` and
``ERR`` exentions (extensions 1,2,4, and 5) are affected, while the
``DQ`` extensions are untouched.

Authors
-------

    - Matthew Bourque, 2014
    - John Biretta, 2012

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import logging

from astropy.io import fits


def apply_norm(frame, operand, norm_value, region):
    """Perform the calculation to apply image normalization.

    The given region in the given extension (``ext``) of the given
    frame is either multiplied or divided by (as determined by
    ``operand``) the ``norm_value``.

    Parameters
    ----------
    frame : array
        The data from the particular FITS extension.
    operand : str
        Either ``*`` for multiplication or ``/`` for division.
    norm_value : float
        The value to multiply or divide by.
    region : str
        Either ``regionAorC`` to apply normalization to amps A or C,
        ``regionBorD`` to apply normalization to amps B or D, or
        ``None`` to apply normalization to the entire chip.
    """

    # Determine region
    if region == 'regionAorC':
        x1 = 0
        x2 = 2048
    elif region == 'regionBorD':
        x1 = 2048
        x2 = 4096
    elif region == 'None':
        x1 = 0
        x2 = 4096

    # Apply gain to specific region
    if operand == '*':
        frame[0:2051, x1:x2] = frame[0:2051, x1:x2] * norm_value
    elif operand == '/':
        frame[0:2051, x1:x2] = frame[0:2051, x1:x2] / norm_value


def electrons_normalize(superdark):
    """The main function of the ``electrons_normalize`` module.

    The image given by outputfinal is converted from units of ``DN``
    to ``e-/s``.  This is done by multiplying the image by the gain
    and dividing the image by the exposure time.  The gain and exposure
    time are grabbed from the header.

    Parameters
    ----------
    superdark : str
        The path to the superdark to process.
    """

    logging.info('\tConverting {} to electrons.'.format(superdark))

    # Open the image and get the data
    hdulist = fits.open(superdark, 'update')
    sci1 = hdulist[1].data
    err2 = hdulist[2].data
    sci4 = hdulist[4].data
    err5 = hdulist[5].data

    # Find gains and exposure time
    gain = {}
    gain['A'] = hdulist[0].header['ATODGNA']
    gain['B'] = hdulist[0].header['ATODGNB']
    gain['C'] = hdulist[0].header['ATODGNC']
    gain['D'] = hdulist[0].header['ATODGND']
    exptime = hdulist[0].header['EXPTIME']

    # Multiply each "half" of the extensions by the appropriate gain.
    logging.info('\tMultiplying each quadrant by its gain.')
    apply_norm(sci1, '*', gain['C'], 'regionAorC')
    apply_norm(err2, '*', gain['C'], 'regionAorC')
    apply_norm(sci1, '*', gain['D'], 'regionBorD')
    apply_norm(err2, '*', gain['D'], 'regionBorD')
    apply_norm(sci4, '*', gain['A'], 'regionAorC')
    apply_norm(err5, '*', gain['A'], 'regionAorC')
    apply_norm(sci4, '*', gain['B'], 'regionBorD')
    apply_norm(err5, '*', gain['B'], 'regionBorD')

    # Normalizing the gain to 1 is not necessary since calwf3
    # doesn't look at this keyword. It already assumes the units
    # of the darks are e-/sec and will use the gains in CCDTAB to
    # reconvert the darks to DNs. But we do it for consistency.
    logging.info('\tNormalizing the SCI and ERR extensions (1, 2, 4, 5) ' + \
        'by the integration time.')
    apply_norm(sci1, '/', exptime, 'None')
    apply_norm(err2, '/', exptime, 'None')
    apply_norm(sci4, '/', exptime, 'None')
    apply_norm(err5, '/', exptime, 'None')

    # Update necessary keywords
    for ext in range(7):
        hdulist[ext].header['CCDGAIN'] = 1.0
    hdulist[0].header['EXPTIME'] = 1.0
    hdulist[0].header['TEXPTIME'] = 1.0
    hdulist.close()
