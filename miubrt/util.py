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


def polyarea(p):
    """Calculate area of polygon using shoelace formula.

    https://en.wikipedia.org/wiki/Shoelace_formula

    .. math::

        A = \frac{1}{2} \\sum_{i=1}^{n} x_i(y_{i+1} - y_{i-1})
    """
    area = 0.0
    x = p[:, 0]
    y = p[:, 1]
    for i in range(len(x)):
        area += x[i - 1] * (y[i] - y[i - 2])
    return abs(area) / 2.0


