"""Average-combine an input list of darks to create "superdarks" or
"masterdarks".

An input list of ``nn*blv_tmp.fits`` files are combined using ``numpy``
masked arrays, where the masks are determined by the ``DQ`` arrays of
the input files.  While the ``SCI`` arrays are a simple
``np.ma.mean()`` operation, the ``ERR`` arrays are calculated via
nominal error propagation (i.e. summing the errors in quadrature)[1].
The ``DQ`` arrays are zeroed out.

Authors
-------

    Matthew Bourque, 2016

Use
---

    This module is intended to be called by ``build_superdarks.py`` and
    ``build_masterdark.py`` as part of the UVIS dark reference file
    pipeline.

References
----------

    1. Bevington, Phillip R. and Robinson, D. Keith, 'Data Reduction
    and error analysis for the physical sciences', 2003.
"""

import logging

from astropy.io import fits
import numpy as np


def combine(nnblv_list, outputname):
    """Perform image combination of individual dark frames.

    Parameters
    ----------
    nnblv_list : list
        The list of paths to input nnblv files to process.
    outputname : str
        The name of the output file that will be written.
    """

    print('')
    print('\tPerforming image combination for {}.'.format(outputname))

    # Get header of first nnblv file
    hdulist_nnblv = fits.open(nnblv_list[0], mode='readonly')

    # Create masked arrays
    masked_nnblvs_ext1, masked_nnblvs_ext2, masked_nnblvs_ext4, masked_nnblvs_ext5 = [], [], [], []
    for nnblv in nnblv_list:
        with fits.open(nnblv, mode='readonly') as hdulist:

            # For chip 2
            mask_ext3 = np.zeros(hdulist[3].data.shape)
            mask_ext3[np.where(hdulist[3].data != 0)] = 1
            masked_nnblv_ext1 = np.ma.masked_array(hdulist[1].data, mask=mask_ext3)
            masked_nnblvs_ext1.append(masked_nnblv_ext1)
            masked_nnblv_ext2 = np.ma.masked_array(hdulist[2].data, mask=mask_ext3)
            masked_nnblvs_ext2.append(masked_nnblv_ext2)

            # For chip 1
            mask_ext6 = np.zeros(hdulist[6].data.shape)
            mask_ext6[np.where(hdulist[6].data != 0)] = 1
            masked_nnblv_ext4 = np.ma.masked_array(hdulist[4].data, mask=mask_ext6)
            masked_nnblvs_ext4.append(masked_nnblv_ext4)
            masked_nnblv_ext5 = np.ma.masked_array(hdulist[5].data, mask=mask_ext6)
            masked_nnblvs_ext5.append(masked_nnblv_ext5)

    # Average combine SCI arrays
    comb_ext1 = np.ma.mean(masked_nnblvs_ext1, axis=0).data.astype(np.float32)
    comb_ext4 = np.ma.mean(masked_nnblvs_ext4, axis=0).data.astype(np.float32)

    # Create empty DQ arrays
    comb_ext3 = np.zeros((2051, 4096)).astype(np.int16)
    comb_ext6 = np.zeros((2051, 4096)).astype(np.int16)

    # Propoagate uncertainties for ERR arrays
    weight_image_ext1 = np.zeros((2051, 4096))
    weight_image_ext4 = np.zeros((2051, 4096))
    for nnblv in masked_nnblvs_ext1:
        mask = nnblv.mask
        weight_image_ext1[np.where(mask == False)] += 1.0
    for nnblv in masked_nnblvs_ext4:
        mask = nnblv.mask
        weight_image_ext4[np.where(mask == False)] += 1.0

    masked_nnblvs_ext2_squared = [(item * (1/weight_image_ext1))**2 for item in masked_nnblvs_ext2]
    masked_nnblvs_ext5_squared = [(item * (1/weight_image_ext4))**2 for item in masked_nnblvs_ext5]
    comb_ext2 = np.sqrt(np.ma.sum(masked_nnblvs_ext2_squared, axis=0)).data.astype(np.float32)
    comb_ext5 = np.sqrt(np.ma.sum(masked_nnblvs_ext5_squared, axis=0)).data.astype(np.float32)

    # Build FITS file
    hdu0 = fits.PrimaryHDU(header=hdulist_nnblv[0].header)
    hdu1 = fits.ImageHDU(comb_ext1, header=hdulist_nnblv[1].header)
    hdu2 = fits.ImageHDU(comb_ext2, header=hdulist_nnblv[2].header)
    hdu3 = fits.ImageHDU(comb_ext3, header=hdulist_nnblv[3].header)
    hdu4 = fits.ImageHDU(comb_ext4, header=hdulist_nnblv[4].header)
    hdu5 = fits.ImageHDU(comb_ext5, header=hdulist_nnblv[5].header)
    hdu6 = fits.ImageHDU(comb_ext6, header=hdulist_nnblv[6].header)
    hdulist = fits.HDUList([hdu0, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6])

    # Save FITS file
    hdulist.writeto(outputname)
