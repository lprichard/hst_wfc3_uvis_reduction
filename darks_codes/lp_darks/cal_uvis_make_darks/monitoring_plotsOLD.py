"""Generate hot pixel and median dark current plots.

As part of the WFC3/UVIS monitoring, hot pixel and median dark current
plots are generated showing the evolution over time since the
installation of WFC3.  Anneal cycles are plotted as shaded gray
regions. SIC&DH failures are plotted as red vertical lines.  A green
vertical line is plotted to represent the start of post-flashing. Plots
are saved to the plotting directory
(``/grp/hst/wfc3k/uvis_darks/plots/``).

Output products include:

1. ``hotpix_plot_ext<ext>_<today>_<recent?>_<ctecorr?>.png`` - hot
   pixel plots for each ``SCI`` extension.
2. ``midpt_plot_ext<ext>_<today>_<recent?>_<ctecorr?>.png`` - median
   dark current plots for each ``SCI`` extension.

Plots are written to ``/grp/hst/wfc3k/uvis_darks/plots/``.  Note that
``<recent>`` and ``<ctecorr>`` are optional.

Authors
-------

    - Matthew Bourque, 2013
    - John Biretta, 2012

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import datetime
import glob
import logging
import subprocess

from astropy.io import ascii
import numpy as np

# Turn off DISPLAY invoke in matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class MonitoringPlot(paths):
    """Class for WFC3/UVIS monitoring plots.
    
    This class allows the user to plot UVIS hot pixel and dark current
    data over time since the start of WFC3.  The plots contain lines
    representing UVIS anneals, post-flash start, and SIC&DH failures.

    Parameters - LP added
    ----------
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.
    """

    def __init__(self, paths):
        """Initialize monitoring plot attributes.

        # LP edit, added paths to class to call on history.txt file
        This method will set anneal and sic&dh dates as class
        attributes, as well as set basic plot and axis parameters.
        """
        self.anneal_list = get_history_list('anneal', paths)
        self.sicdh_list = get_history_list('SIC&DH', paths)
        self.set_plotting_parameters()
        self.set_axis()

    def set_plotting_parameters(self):
        """Set basic ``matplotlib`` parameters."""

        plt.rcParams['legend.fontsize'] = 10
        plt.minorticks_on()
        plt.grid(True)

    def set_axis(self):
        """Set basic axis parameters."""

        self.ax = plt.subplot(111)
        for tick in self.ax.get_xticklabels():
            tick.set_fontsize(12)
        for tick in self.ax.get_yticklabels():
            tick.set_fontsize(12)

    def plot_rec(self, x, y, output):
        """Create a monitoring plot.

        Parameters:
            x : numpy array
                An array containing the times since first observation.
            y : numpy array
                An array containing either the hotpix or midpt stats,
                depending on the plot_type.
            output: string
                The path to the output save location.

        Returns:
            nothing

        Outputs:
            A saved png of the monitoring plot, written to what is
            specified in output.
        """
        # Plot sicdh failures
        plt.axvline(x=self.sicdh_list[0], linewidth=0.9, color='red', label='SIC&DH Lockup')
        for sicdh in self.sicdh_list:
            plt.axvline(x=sicdh, linewidth=0.9, color='red')

        # Shade anneals
        anneal_begin = list(range(0, len(self.anneal_list), 2))
        anneal_end = list(range(1, len(self.anneal_list), 2))
        for begin, end in zip(anneal_begin, anneal_end):
            plt.axvspan(self.anneal_list[begin], self.anneal_list[end], facecolor='0.5', alpha=0.5)
        if len(self.anneal_list) & 1 == 1: # For last anneal to plot end
            plt.axvspan(self.anneal_list[-1], max(x), facecolor='0.5', alpha=0.5)

        # Plot postflash start
        plt.axvline(x=1246, linewidth=3.0, color='green', label='Post-flash Start')

        # Get recent data
        lengthx = len(x)
        start_recx = len(x) - 750
        x_recent = x[start_recx:lengthx]
        lengthy = len(y)
        start_recy = len(y) - 750
        y_recent = y[start_recy:lengthy]

        # Plot data
        self.ax.plot(x_recent, y_recent, 'kx')

        # Plot legend
        self.ax.legend(loc='upper left')

        # Convert time labels to year
        ticks, tick_labels = get_year_labels(x)
        self.ax.set_xticks(ticks)
        self.ax.set_xticklabels(tick_labels)

        # Save plot
        plt.savefig(output)
        print('\tSaved monitoring plot to {}'.format(output))

        # Clear the figure
        plt.clf()

    def plot(self, x, y, output):
        """Create a monitoring plot.

        Parameters
        ----------
        x : array
            An array containing the times since first observation.
        y : array
            An array containing either the ``hotpix`` or ``midpt``
            stats, depending on the ``plot_type``.
        output : str
            The path to the output save location.
        """

        # Plot sicdh failures
        plt.axvline(x=self.sicdh_list[0], linewidth=0.9, color='red', label='SIC&DH Lockup')
        for sicdh in self.sicdh_list:
            plt.axvline(x=sicdh, linewidth=0.9, color='red')

        # Shade anneals
        anneal_begin = list(range(0, len(self.anneal_list), 2))
        anneal_end = list(range(1, len(self.anneal_list), 2))
        for begin, end in zip(anneal_begin, anneal_end):
            plt.axvspan(self.anneal_list[begin], self.anneal_list[end], facecolor='0.5', alpha=0.5)
        if len(self.anneal_list) & 1 == 1:  # For last anneal to plot end
            plt.axvspan(self.anneal_list[-1], max(x), facecolor='0.5', alpha=0.5)

        # Plot postflash start
        plt.axvline(x=1246, linewidth=3.0, color='green', label='Post-flash Start')

        # Plot data
        self.ax.plot(x, y, 'kx')

        # Plot legend
        self.ax.legend(loc='upper left')

        # Convert time labels to year
        ticks, tick_labels = get_year_labels(x)

        self.ax.set_xticks(ticks)
        self.ax.set_xticklabels(tick_labels)

        # Save plot
        plt.savefig(output)
        print('\tSaved monitoring plot to {}'.format(output))

        # Clear the figure
        plt.clf()


def get_data(plot_type, ext, ctecorr, paths):
    """Read in data from each ``d???_cy??_ext?.dat`` file.

    Parameters
    ----------
    plot_type : str
        The plot type (can be ``hotpix`` or ``midpt``).
    ext : int
        The FITS SCI extension (can be ``1`` or ``4``).
    ctecorr : bool
        ``True`` (if CTE-correction is turned on) or ``False`` (if
        CTE-correction is turned off).
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.

    Returns
    -------
    x : array
        An array containing the times since first observation.
    y : array
        An array containing either the ``hotpix`` or ``midpt`` stats,
        depending on the ``plot_type``.

    Notes
    -----
    If ``ctecorr`` is on, then all data files are read in. If
    ``ctecorr`` is off, then only non-ctecorr data files are read in.
    """

    # Determine which data files to read in
    if ctecorr:
        data_files = glob.glob('{}/{}/{}_stat_ext{}_????-??-??_ctecorr.dat'.format(paths['plot_dir'], plot_type, plot_type, ext))
    else:
        data_files = glob.glob('{}/{}/{}_stat_ext{}_????-??-??.dat'.format(paths['plot_dir'], plot_type, plot_type, ext))

    x, y = [], []

    for data_file in data_files:
        data = ascii.read(data_file, names=['nnblv', 'time', 'ydata'])
        for xdata, ydata in zip(data['time'], data['ydata']):
            x.append(float(xdata))
            y.append(float(ydata))

    x = np.array(x)
    y = np.array(y)

    return x, y


def get_history_list(hist_type, paths):
    """Return a list of anneal dates or SIC&DH failures (depending on
    the ``hist_type``) by grepping the ``history.txt`` file.

    Parameters
    ----------
    hist_type : str
        The type of history to grep for (can be ``anneal`` or
        ``SIC&DH``).
    LP added:
    paths : dict
        The dictionary containing the absolute paths of directories
        used throughout the the pipeline.

    Returns
    -------
    hist_list : list
        A list of anneal dates or SIC&DH failure dates, depending
        on the ``hist_type``.
    """

    # Build grep command
    history_file = paths['history']
    grep_command = 'grep -E "54.*{}|55.*{}|56.*{}|57.*{}|58.*{}" {}'.format(
        hist_type, hist_type, hist_type, hist_type, hist_type, history_file)

    # Execute grep command and capture results
    grep_results = subprocess.check_output(grep_command, shell=True)
    grep_results = grep_results.split(b'\n')

    # Build anneal list, subtracting out zero point
    hist_list = [result[0:9].decode() for result in grep_results]
    hist_list = [item for item in hist_list if item]
    hist_list = [(float(hist) - 54993) for hist in hist_list]

    return hist_list


def get_year_labels(x):
    """Returns two arrays - one of the placement of year
        tick marks and one of the labels for the years.

    Parameters:
        x : integer
            the number of days since observations began

    Returns:
        ticks : numpy array
            An array of where the tick labels should go

        tick_marks: numpy array
            An array of year labels

    Outputs:
        nothing

    Notes:

    """

    # Convert time labels to year
    days = max(x) - 205 # started on 205th day of 2009
    num_years = int(days / 365.25)
    ticks = np.asarray([205 + 365 * n for n in range(num_years + 1)])
    ticks[3:] += 1 # account for leap year in 2012
    ticks[7:] += 1 # account for leap year in 2016
    tick_labels = [str(2010 + n) for n in range(len(ticks))]

    return ticks, tick_labels


def monitoring_plots_main(ctecorr, paths):
    """The main function of the ``monitoring_plots`` module.

    Hot pixel and median dark current plots are generated showing the
    evolution since the installation of WFC3.  See module documentation
    for further details.

    Parameters
    ----------
    ctecorr : bool
        ``True`` (if CTE-correction is turned on) or ``False`` (if
        CTE-correction is turned off).
    paths : dict
        The dictionary containing the absolute paths of directories
    """

    # Determine today for output name
    today = datetime.datetime.now()
    today = datetime.datetime.strftime(today, '%Y-%m-%d')

    # FITS extensions to iterate over
    exts = [1, 4]
    chips = [2, 1]

    # Plot the nominal hot pixel evolution
    for ext, chip in zip(exts, chips):
        x, y = get_data('hotpix', ext, ctecorr, paths)
        mp = MonitoringPlot(paths)  #LP edit to call on history.txt later
        mp.ax.axis([0, max(x) + 20, 0, 7])
        mp.ax.set_title('WFC3/UVIS Hot Pixel Evolution (Chip {})'.format(chip))
        mp.ax.set_ylabel('Number of Hot Pixels (% of chip)', size=12)
        output = '{}/hotpix/hotpix_plot_ext{}_{}.png'.format(paths['plot_dir'], ext, today)

        # Convert time labels to year
        ticks, tick_labels = get_year_labels(x)
        mp.ax.set_xticks(ticks)
        mp.ax.set_xticklabels(tick_labels)

        if ctecorr == True:
            output = output.replace('.png', '_ctecorr.png')
        mp.plot(x, y, output)

    # Plot the recent hot pixel evolution
    for ext, chip in zip(exts, chips):
        x, y = get_data('hotpix', ext, ctecorr, paths)
        # Get recent data
        lengthx = len(x)
        start_recx = len(x) - 500
        x_recent = x[start_recx:lengthx]

        mp = MonitoringPlot(paths)  #LP edit to call on history.txt later
        mp.ax.axis([min(x_recent) - 150, max(x_recent) + 150, 3.2, 6])
        mp.ax.set_title('Recent WFC3/UVIS Hot Pixel Evolution (Chip {})'.format(chip))
        mp.ax.set_ylabel('Number of Hot Pixels (% of chip)', size=12)
        output = '{}/hotpix/hotpix_plot_ext{}_{}_recent.png'.format(paths['plot_dir'], ext, today)

        # Convert time labels to year
        ticks, tick_labels = get_year_labels(x)

        mp.ax.set_xticks(ticks)
        mp.ax.set_xticklabels(tick_labels)

        if ctecorr:
            output = output.replace('_recent.png', '_ctecorr_recent.png')
        mp.plot_rec(x, y, output)

    # Plot the nominal dark current evolution
    for ext, chip in zip(exts, chips):
        x, y = get_data('midpt', ext, ctecorr, paths)
        mp = MonitoringPlot(paths)  #LP edit to call on history.txt later
        mp.ax.axis([0, max(x) + 20, 0, 12])
        mp.ax.set_title('WFC3/UVIS Dark Current Evolution (Chip {})'.format(chip))
        mp.ax.set_ylabel('Median Dark Current (e-/hr)', size=12)
        output = '{}/midpt/midpt_plot_ext{}_{}.png'.format(paths['plot_dir'], ext, today)

        # Convert time labels to year
        ticks, tick_labels = get_year_labels(x)

        mp.ax.set_xticks(ticks)
        mp.ax.set_xticklabels(tick_labels)

        if ctecorr:
            output = output.replace('.png', '_ctecorr.png')
        mp.plot(x, y, output)

    # Plot the recent dark current evolution
    for ext, chip in zip(exts, chips):
        x, y = get_data('midpt', ext, ctecorr, paths)
        # Get recent data
        lengthx = len(x)
        start_recx = len(x) - 500
        x_recent = x[start_recx:lengthx]

        mp = MonitoringPlot(paths)  #LP edit to call on history.txt later
        mp.ax.axis([min(x_recent) - 75, max(x_recent) + 75, 4, 11])
        mp.ax.set_title('Recent WFC3/UVIS Dark Current Evolution (Chip {})'.format(chip))
        mp.ax.set_ylabel('Median Dark Current (e-/hr)', size=12)
        output = '{}/midpt/midpt_plot_ext{}_{}_recent.png'.format(paths['plot_dir'], ext, today)

        # Convert time labels to year
        ticks, tick_labels = get_year_labels(x)
        mp.ax.set_xticks(ticks)
        mp.ax.set_xticklabels(tick_labels)

        if ctecorr == True:
            output = output.replace('_recent.png', '_ctecorr_recent.png')
        mp.plot_rec(x, y, output)
