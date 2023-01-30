#!/usr/bin/env python
# Copyright (c) 2019-2023, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

import datetime as dt
import glob
import os
import re


def get_file_date_regex(filename):
    """Get regex from filename"""
    # regex for ""%Y-%m-%d--%H:%M:%S"
    reg0 = r"\d{4}.\d{2}.\d{2}..\d{2}.\d{2}.\d{2}"
    # regex for "%Y%m%d%H%M%S"
    reg1 = r"\d{14}"
    match = re.search(reg0, os.path.basename(filename))
    return reg1 if match is None else reg0


def get_datetime_from_filename(filename, regex):
    """Get datetime from filename"""
    fmt = "%Y%m%d%H%M%S"
    match = re.search(regex, os.path.basename(filename))
    match = "".join(re.findall(r"[0-9]+", match.group()))
    return dt.datetime.strptime(match, fmt)


def create_xpol_filelist(path_glob, starttime, endtime):
    """Create filelist from path_glob and filename dates"""
    file_names = sorted(glob.glob(path_glob))
    regex = get_file_date_regex(file_names[0])
    for fname in file_names:
        time = get_datetime_from_filename(fname, regex)
        if time >= starttime and time < endtime:
            yield fname


def create_cpol_filelist(inpath, starttime, endtime, moments):
    """Create filelist from path_glob and filename dates"""
    for mom in moments:
        pglob = inpath + f'{mom}_*'
        file_names = sorted(glob.glob(pglob))
        regex = get_file_date_regex(file_names[0])
        for fname in file_names:
            time = get_datetime_from_filename(fname, regex)
            if time >= starttime and time < endtime:
                yield fname
