"""Generates hot pixel and dark current ``bokeh`` plots.

This script produces several ``bokeh`` plots that describe the WFC3/
UVIS hot pixel and dark current populations. It is based on the
original ``matplotlib`` plotting in ``monitoring_plots.py``.

Authors
-------
    - Catherine A. Martlin, November 2015; March 2019.
    - Jennifer V. Medina, August 2018.
    - Matthew Bourque

Use
---

    This module is imported and used by ``cal_uvis_make_darks.py`` as
    such:
    ::

        from cal_uvis_make_darks.monitoring_bokeh_plots import monitoring_plots_bokeh_main
        monitoring_plots_bokeh_main()
"""

import datetime
import glob
import logging
import subprocess
import numpy as np

from astropy.io import ascii
from astropy.time import Time

from bokeh.plotting import figure, output_file, show, ColumnDataSource, save
from bokeh.models import HoverTool

# LP commented, not called anywhere
# MAIN_DIR = '/grp/hst/wfc3k/uvis_darks/plots/new_algorithm'

def convert_dates_for_axis(time_ext_x):
    """Converts the x axis from days since data began being taken
    (June 11th 2009) to ``datetime`` objects for plotting.

    Parameters
    ----------
    time_ext_x : list
        A list of dates to convert.

    Returns
    -------
    converted_date : list
        A list of dates converted to ``datetime`` objects.
    """
    d = datetime.date(year=2009, month=6, day=11)
    start_date = d.toordinal()

    dates_1 = []
    for val in range(len(time_ext_x)):
        time_1 = np.asarray(time_ext_x)
        date1 = float(time_ext_x[val])
        dates_1.append(float(start_date)+date1)

    # Intrepret date string via strptime and convert to datetime
    converted_date = []
    for i in range(len(dates_1)):
        dt = dates_1[i]
        int_dt = int(dt)
        temp_1_greg = datetime.date.fromordinal(int_dt)
        converted_date.append(temp_1_greg)

    return converted_date


def get_data(ext, data_type, plot_type, ctecorr, paths):
    """Retrieves data from each ``d???_cy??_ext?.dat`` file.

    Parameters
    ----------
    ext : int
        The extension (i.e. ``1`` or ``4``).
    data_type : str
        The data type used to calculate the statistics (i.e. ``blvs`` or
        ``superdarks``).
    plot_type : str
        The plot type (i.e. ``hotpix`` or ``midpt``).

    Returns
    -------
    time : array
        An array of time values (in MJD)
    ydata : array
        An array of values for the given ``data_type`` and ``plot_type``
    """
    # Determine which data files to read in
    if ctecorr:
        data_files = glob.glob('{}/{}/{}_stat_ext{}_????-??-??_ctecorr.dat'.format(paths['plot_'+data_type], plot_type, plot_type, ext))
    else:
        data_files = glob.glob('{}/{}/{}_stat_ext{}_????-??-??.dat'.format(paths['plot_'+data_type], plot_type, plot_type, ext))

    filename, x, y = [], [], []

    for data_file in data_files:
        data = ascii.read(data_file, names=['filename', 'time', 'ydata'])
        for name, xdata, ydata in zip(data['filename'], data['time'], data['ydata']):
            filename.append(str(name))
            x.append(float(xdata))
            y.append(float(ydata))

    x = np.array(x)
    y = np.array(y)

    return filename, x, y


def get_history_list(hist_type, paths):
    """Gathers a list of anneal dates or sic&dh failures (depending on
    the ``hist_type``) by grepping the ``history.txt`` file in
    ``/grp/hst/wfc3b/calibration/``.
    
    Parameters
    ----------
    hist_type : str
        The type of history to grep. Can either be ``anneal`` or
        ``SIC&DH``.
    LP added:
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.

    Returns
    -------
    hist_list : list
        A list of times for the given ``hist_type``.
    """
    # Build grep command
    grep_command = 'grep -E "54.*' + hist_type + '|55.*' + hist_type + \
                   '|56.*' + hist_type + '|57.*' + hist_type + '|58.*' + hist_type + \
                   '|59.*' + hist_type + '" ' + \
                   paths['history']    #LP added paths

    # Execute grep command and capture results
    grep_results = subprocess.check_output(grep_command, shell=True)
    grep_results = grep_results.split(b'\n')

    # Build anneal list
    hist_list = [result[0:9].decode() for result in grep_results]
    hist_list = [item for item in hist_list if item]
    hist_list = [(float(hist)) for hist in hist_list]

    return hist_list


def plot_sicdh_failures(p, sicdh_list_in_seconds):
    """Plots the SIC&DH Failures.

    Parameters
    ----------
    p : obj
        The ``bokeh`` plot object for the plot.
    sicdh_list_in_seconds : list
        A list of times for sic&dh failures (in seconds).
    """

    p.ray(x = sicdh_list_in_seconds[0], y=0.0, length=0, angle=1.57079633, color='red', legend='SIC&DH Lockup')

    for sicdh in sicdh_list_in_seconds:
        p.ray(x=sicdh, y=0.0, length=0, angle=1.57079633, color='red')


def shade_anneals(anneal_list, anneal_in_seconds, xmax, ymin, ymax, p):
    """Adds a shaded region between every other set of anneals.

    Parameters
    ----------
    anneal_list : list
        A list of anneal dates.
    anneal_in_seconds : list
        A list of anneal dates (in seconds).
    xmax : float
        The maximum x axis value.
    ymin : float
        The minimum y axis value.
    ymax : float
        The maximum y axis value.
    p : obj
        The ``bokeh`` plot object for the plot.

    Returns
    -------
    shaded_regions : obj
        ``bokeh`` plot objects to shade the anneal regions.
    """    
    anneal_begin = list(range(0, len(anneal_list), 2))
    anneal_end = list(range(1, len(anneal_list), 2))

    source_anneal = ColumnDataSource(data=dict(start=anneal_begin, end=anneal_end))
    TOOLTIPS_anneal = [
                ('anneal start', '@start'),
                ('anneal end', '@end'),
                ]

    for i, j in zip(anneal_begin, anneal_end):
        leftside = anneal_in_seconds[i]
        rightside = anneal_in_seconds[j]
        shaded_regions = p.quad(top=ymax, bottom=ymin, left=leftside, right=rightside,
                                source=source_anneal, color="#C0C0C0")
    if len(anneal_list) & 1 == 1:
        leftside = float(anneal_in_seconds[-1])
        rightside = leftside + (float(anneal_in_seconds[1]) - float(anneal_in_seconds[0]))
        shaded_regions = p.quad(top=ymax, bottom=ymin, left=leftside, right=rightside,
                                source=source_anneal, color="#C0C0C0")

    return shaded_regions, TOOLTIPS_anneal


def unix_time_millis(dt):
    """Must transform the ``datetime`` of the sicdh failure lines and
    the anneal shading into seconds since epoch to give the glyphs for
    the units they want for the x-axis ``datetime`` plotting of the
    main ratios.

    Parameters
    ----------
    dt : list
        A list of ``datetime`` objects to transform

    Returns
    second_dt : list
        A list of `datetime`` objects transformed into seconds.`
    """

    epoch = datetime.datetime.utcfromtimestamp(0)
    second_dt = []

    for val in range(len(dt)):
        y = dt[val]
        second_dt.append((y - epoch).total_seconds() * 1000.0)

    return second_dt


def monitoring_plots_bokeh_main(ctecorr, paths):
    """The main function of the ``bokeh_monitoring_dark_plots`` module.
    See module docstrings for further details.

    Parameters
    ----------
    ctecorr : bool
        A bool object. If True, plots and files will be labeled with
        ``CTE-Corrected`` and ``ctecorr`` (respectively).

    paths : dict
        A dictonary whose keys are path identifiers and whose
        values are strings containing absolute paths.
    """

    # Read data fron anneal_list.dat and sich_list.dat
    anneal_list = get_history_list('anneal', paths)  #LP added paths
    sicdh_list = get_history_list('SIC&DH', paths)   #LP added paths

    for data_type in ['blvs', 'superdarks']:
        for ext in [1, 4]:
            for plot_type in ['hotpix', 'midpt']:

                # Read in data
                filenames, times, ydata = get_data(ext, data_type, plot_type, ctecorr, paths)

                # Set title:
                if plot_type == 'hotpix':
                    if ext == 1:
                        if not ctecorr:
                            plot_title = ('WFC3/UVIS Hot Pixel Evolution (Chip 2)')
                            label_title = 'Percent of Hot Pixels'
                        elif ctecorr:
                            plot_title = ('WFC3/UVIS Hot Pixel Evolution (Chip 2 CTE-Corrected)')
                            label_title = 'Percent of Hot Pixels'
                    elif ext == 4:
                        if not ctecorr:
                            plot_title = ('WFC3/UVIS Hot Pixel Evolution (Chip 1)')
                            label_title = 'Percent of Hot Pixels'
                        elif ctecorr:
                            plot_title = ('WFC3/UVIS Hot Pixel Evolution (Chip 1 CTE-Corrected)')
                            label_title = 'Percent of Hot Pixels'
                elif plot_type == 'midpt':
                    if ext == 1:
                        if not ctecorr:
                            plot_title = ('WFC3/UVIS Dark Current Evolution (Chip 2)')
                            label_title = 'Median Dark Current'
                        elif ctecorr:
                            plot_title = ('WFC3/UVIS Dark Current Evolution (Chip 2 CTE-Corrected)')
                            label_title = 'Median Dark Current'
                    elif ext == 4:
                        if not ctecorr:
                            plot_title = ('WFC3/UVIS Dark Current Evolution (Chip 1)')
                            label_title = 'Median Dark Current'
                        elif ctecorr:
                            plot_title = ('WFC3/UVIS Dark Current Evolution (Chip 1 CTE-Corrected)')
                            label_title = 'Median Dark Current'

                # Set axis limits
                xmax = max(times) + 20.
                ymin = 0.0
                if plot_type == 'hotpix':
                    ymax = 7.4
                elif plot_type == 'midpt':
                    ymax = 12.5

                # Get shaded anneals:
                anneals_list = Time(anneal_list, format='mjd')
                anneal_list_dates = anneals_list.datetime
                anneal_in_seconds = unix_time_millis(anneal_list_dates)

                # Get the SIC&DH Failures:
                sicdh_list = Time(sicdh_list, format='mjd')
                sicdh_list_dates = sicdh_list.datetime
                sicdh_in_seconds = unix_time_millis(sicdh_list_dates)

                # Get value of post-flash start:
                postflash_start = datetime.datetime(2012, 11, 8, 0, 0)
                postflash_date = Time(postflash_start, format='datetime')
                postflash_date = postflash_date.datetime
                epoch = datetime.datetime.utcfromtimestamp(0)
                postflash_seconds = (postflash_date - epoch).total_seconds() * 1000.0

                # Convert x axis to datetime dates:
                converted_dates = convert_dates_for_axis(times)


                # Create plots

                # Create source for hover tools:
                time_array = []
                for converted_date in converted_dates:
                    t = str(converted_date)
                    time_array.append(t)

                source = ColumnDataSource(data=dict(name=filenames, time=converted_dates, greg_time=time_array, ydata=ydata))

                TOOLTIPS = [
                ('filename', '@name'),
                ('USEAFTER', '@greg_time'),
                ('y', '@ydata'),
                ]

                # Set general plotting parameters
                tools = "pan,wheel_zoom,box_zoom,reset,save"

                # Plots
                p = figure(x_axis_type="datetime", y_range=(ymin, ymax), title=plot_title+' using '+data_type, tools=tools, plot_width=850, plot_height=650)#, tooltips=TOOLTIPS)

                # Plot shaded anneals:
                shaded_regions, TOOLTIPS_anneal = shade_anneals(anneal_list, anneal_in_seconds, xmax, ymin, ymax, p)

                # Plot the SIC&DH Failures:
                plot_sicdh_failures(p, sicdh_in_seconds)

                # Set y label:
                if plot_type == 'hotpix':
                    p.yaxis.axis_label = ('Number of Hot Pixels (% of Chip)')
                elif plot_type == 'midpt':
                    p.yaxis.axis_label = ('Median Dark Current (e-/hr)')

                # Set x label:
                p.xaxis.axis_label = ('Date')

                # Creating bokeh plot
                if data_type == 'blvs':
                    r1 = p.circle('time', 'ydata', color='blue', alpha=0.5, source=source)
                elif data_type == 'superdarks':
                    r2 = p.square('time', 'ydata', color='blue', alpha=0.5, source=source)

                r3 = p.ray(x = postflash_seconds, y=0.0, length=0, angle=1.57079633, line_width=3, color='green', legend = 'Post-Flash Start')
                p.legend.location = 'top_left'

                # Adding hovertools to the superdark bokeh plots
                if data_type == 'superdarks':
                    hover = HoverTool(tooltips=TOOLTIPS, renderers=[r2])
                    p.add_tools(hover)

                elif data_type == 'blvs':
                    pass

                # Save figure
                today = datetime.datetime.now()
                today = datetime.datetime.strftime(today, '%Y-%m-%d')

                if ctecorr == True:
                    savefig_html = '{}/{}/{}_{}_plot_ext{}_{}_ctecorr.html'.format(paths['plot_'+data_type], plot_type, data_type, plot_type, ext, today)
                else:
                    savefig_html = '{}/{}/{}_{}_plot_ext{}_{}.html'.format(paths['plot_'+data_type], plot_type, data_type, plot_type, ext, today)

                output_file(savefig_html)
                save(p)
                print('Saved figure to ' + savefig_html)
