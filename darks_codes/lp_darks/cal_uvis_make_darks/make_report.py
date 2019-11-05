"""Construct a directory containing output products for quality
control.

This module creates a ``report/`` directory under the main processing
directory that contains JPEGs of input darks, JPEGs of output dark
reference files, data files with image statistics, and dark current &
hot pixel plots.  These output products should be analyized and
reviewed before delivery of the reference files.

Output products include:

1. A ``report/`` directory
2. JPEGs of input raw darks
3. JPEGs of dark reference files
4. ``imstat_ext<ext>.dat`` - image statistics of the dark reference
   files for each ``SCI`` extentions (1 and 4)
5. a copy of ``midpt_plot_ext<ext>_<proc_date>.png``, the dark current
   plots
6. a copy of ``hotpix_plot_ext<ext>_<proc_date>.png``, the hot pixel
   plots.
7. a copy of ``midpt_plot_ext<ext>_<proc_date>.html``, the dark current
   plots
8. a copy of ``hotpix_plot_ext<ext>_<proc_date>.html``, the hot pixel
   plots.

Authors
-------

    Matthew Bourque, 2014
    Catherine Martlin, 2018
    Jennifer V. Medina, 2018

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import copy
import datetime
import glob
import logging
import os
import shutil

from astropy.io import ascii
from astropy.io import fits
import numpy as np
from PIL import Image


def copy_plots_to_report(paths, ctecorr):
    """Copy the dark current and hot pixel plots to the ``report/``
    directory.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    ctecorr : bool
        The CTE correction switch. If ``True`` (on), the plots
        containing the CTE corrected data are copied.
    """

    today = datetime.datetime.now()
    today = datetime.datetime.strftime(today, '%Y-%m-%d')

    for data_type in ['blvs', 'superdarks']:
        for plot_type in ['midpt', 'hotpix']:
            for ext in ['1', '4']:

                if ctecorr:
                    filename_bokeh = '{}_{}_plot_ext{}_{}_ctecorr.html'.format(data_type, plot_type, ext, today)
                else:
                    filename_bokeh = '{}_{}_plot_ext{}_{}.html'.format(data_type, plot_type, ext, today)

                # Copy bokeh plot
                src = os.path.join(paths['plot_'+data_type], plot_type, filename_bokeh)
                dst = os.path.join(paths['report_dir'], filename_bokeh)
                shutil.copyfile(src, dst)
                print('\tCopied {} to {}'.format(src, dst))


def create_imstat_files(paths):
    """Create text files containing image statistics of dark reference
    files.

    Statistics are computed are dark reference files located in the
    ``deliver/`` directory and are stored in text files placed in the
    ``report/`` directory.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    """

    images = glob.glob(os.path.join(paths['deliver_dir'], '*d??.fits'))

    for ext in [1, 4]:
        print('')
        print('\tGathering extension {} statistics.'.format(ext))

        # Initialize table
        table_data = {'image': [], 'npix': [], 'midpt': [], 'stddev': [], 'min': [], 'max': []}

        # Get stats
        for image in images:
            data = fits.getdata(image, ext)
            table_data['image'].append(os.path.basename(image))
            table_data['npix'].append(data.size)
            table_data['midpt'].append(np.median(data))
            table_data['stddev'].append(np.std(data))
            table_data['min'].append(data.min())
            table_data['max'].append(data.max())

        # Write out table
        table_file = os.path.join(paths['report_dir'], 'imstat_ext{}.dat'.format(ext))
        names = ['image', 'npix', 'midpt', 'stddev', 'min', 'max']
        ascii.write(table_data, table_file, names=names)

        print('\t\tstatistics written to {}'.format(os.path.join(paths['report_dir'], 'imstat_ext{}.dat'.format(ext))))


def make_jpeg(frame, paths):
    """Create a JPEG image from a FITS file.

    Parameters
    ----------
    frame : str
        The path to the FITS file.
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    """

    for ext in [1, 4]:

        # Get the image data
        hdulist = fits.open(frame)
        data = hdulist[ext].data

        # Clip the top and bottom 1% of pixels
        sorted_data = copy.copy(data)
        sorted_data = sorted_data.ravel()
        sorted_data.sort()
        top = sorted_data[int(len(sorted_data) * 0.99)]
        bottom = sorted_data[int(len(sorted_data) * 0.01)]
        top_index = np.where(data > top)
        data[top_index] = top
        bottom_index = np.where(data < bottom)
        data[bottom_index] = bottom

        # Scale the data
        data = data - data.min()
        data = (data / data.max()) * 255.
        data = np.flipud(data)
        data = np.uint8(data)

        # Write the images to a JPEG
        image = Image.fromarray(data)
        jpeg_output = os.path.join(paths['report_dir'], '{}_ext{}.jpg'.format(os.path.basename(frame).split('.')[0], str(ext)))
        image.save(jpeg_output)
        print('\t\tSaved JPEG to {}'.format(jpeg_output))


def make_report_main(paths, ctecorr):
    """The main function of the ``make_report`` module.

    A new ``report/`` directory is created containing products to
    review for quality assurance before dark reference file delivery.

    Parameters
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    ctecorr : bool
        The CTE correction switch. If ``True`` (on), the plots
        containing the CTE corrected data are copied.
    """

    print('')
    print('')
    print('')
    print('---------- Making report directory ----------')
    print('')

    # LP added, check that directory actually exists, make the parent
    if not os.path.exists(os.path.dirname(paths['report_dir'])): 
        os.makedirs(os.path.dirname(paths['report_dir']), 0o774)
        print('\tLP log: Created driectory {}'.format(os.path.dirname(paths['report_dir'])))

    # Create report directory
    os.mkdir(paths['report_dir'], 0o774)
    print('\tCreated directory {}'.format(paths['report_dir']))
    print('')

    # Copy dark median and hot pixel plots to report directory
    copy_plots_to_report(paths, ctecorr)

    # Create dark reference file statistics file using imstat
    create_imstat_files(paths)

    # Make JPEGs of dark reference files
    # For blv_tmp.fits files
    print('')
    print('\tCreating JPEGs.')
    print('')
    blv_files = glob.glob(paths['blvs_for_dark'])
    for blv_file in blv_files:
        make_jpeg(blv_file, paths)

    # For masterdarks
    masterdark_files = glob.glob(paths['darks_to_deliver'])
    for masterdark_file in masterdark_files:
        make_jpeg(masterdark_file, paths)
