"""Determine the current latest anneal date from the
``/grp/hst/wfc3a/automated_outputs/ql_wfc3_make_anneal_list/anneal_dates.txt`` file.

Authors
-------

    Matthew Bourque, 2013
    Catherine Martlin, 2018

Use
---

    This module is intended to be called by ``cal_uvis_make_darks.py``
    as part of the UVIS dark reference file creation pipeline.
"""

import datetime
import subprocess


def get_anneal_date():
    """Determines the most recent anneal date and time.

    This function determines the date and time from the latest anneal,
    based on information from 
    ``/grp/hst/wfc3a/automated_outputs/ql_wfc3_make_anneal_list/anneal_dates.txt``.
    If today's date is within 5 days after the latest anneal, then the
    second-to-last anneal date is returned.  Conversely, if today's
    date is beyond 5 days after the lastest anneal, then the latest
    anneal date is returned.  The 5 day buffer is used to allow for
    enough post-anneal data to be available before processing.

    Returns
    -------
    last_anneal_date OR second_to_last_anneal_date : str
        The most recent anneal date, in the format
        ``YYYYMMDD-HH:MM:SS``.
    """

    #Open file and read lines into a list:
    anneal_info = []
    anneal_file = "/grp/hst/wfc3a/automated_outputs/ql_wfc3_make_anneal_list/anneal_dates.txt"
    for line in open(anneal_file):
        anneal_info.append(line)


    # Parse list
    anneal_times = [anneal.split()[2] for anneal in anneal_info]
    anneal_months = [anneal.split()[3] for anneal in anneal_info]
    anneal_days = [anneal.split()[4] for anneal in anneal_info]
    anneal_years = [anneal.split()[1].split('.')[0] for anneal in anneal_info]

    # Determine lastest anneal dates:
    last_anneal_date = '{} {} {} {}'.format(
        anneal_months[-1],
        anneal_days[-1],
        anneal_years[-1],
        anneal_times[-1])
    second_to_last_anneal_date = '{} {} {} {}'.format(
        anneal_months[-2],
        anneal_days[-2],
        anneal_years[-2],
        anneal_times[-2])

    # Convert latest anneal dates to datetime objects:
    last_anneal_date = datetime.datetime.strptime(
        last_anneal_date, '%b %d %Y %H:%M:%S')
    second_to_last_anneal_date = datetime.datetime.strptime(
        second_to_last_anneal_date, '%b %d %Y %H:%M:%S')

    # Determine correct anneal date to used based on current date
    now = datetime.datetime.now()
    five_days = datetime.timedelta(days=5)
    if last_anneal_date + five_days < now:
        last_anneal_date = datetime.datetime.strftime(
            last_anneal_date, '%Y%m%d-%H:%M:%S')
        return last_anneal_date
    elif second_to_last_anneal_date < now < last_anneal_date:
        second_to_last_anneal_date = datetime.datetime.strftime(
            second_to_last_anneal_date, '%Y%m%d-%H:%M:%S')
        return second_to_last_anneal_date