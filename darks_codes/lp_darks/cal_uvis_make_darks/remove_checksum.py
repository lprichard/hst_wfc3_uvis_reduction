"""This module is to ensure all the headers of our superdarks to be
delivered are devoid of the ``CHECKSUM`` keyword.

Authors
-------

    Catherine Martlin, 2016.

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline. It should
    be the last step to ensure the file isn't opened and used again
    which would, possibly, cause the ``CHECKSUM`` keyword to be put in
    again.
"""

import glob
import logging

from astropy.io import fits


def remove_checksum_main(paths):
    """Final preperation of the header of superdark for CRDS
    delivery.

    Parameters
    ----------
    paths : dict
        A dictonary whose keys are path identifiers and whose
        values are strings containing absolute paths.

    Notes
    -----
    This should be the last thing run on the files before delivery.

    References
    ----------
    CRDS keywords are specified in ``TIR CDBS 2008-01``.
    """

    print('---------- Removing CHECKSUM ----------')
    print('')

    superdark_delivery_list = glob.glob(paths['darks_to_deliver'])
    for superdark in superdark_delivery_list:
        hdulist = fits.open(superdark, mode='update', checksum='remove')
        hdulist.close()
