#!/usr/bin/env python
# Copyright (c) 2023, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

import os

import contourpy
import numpy as np

from .util import polyarea


class FilteredContourGenerator(contourpy.Mpl2014ContourGenerator):
    def create_contour(self, level):
        return self.lines(level)

    def create_filled_contour(self, lower_level, upper_level):
        return self.filled(lower_level, upper_level)

    def filled(self, lower_level, upper_level):
        min_area = float(os.environ.get("MPL_MIN_CONTOUR_AREA", "0"))
        # get vertices/kinds
        vertices, kinds = super().filled(lower_level, upper_level)
        # return early if nothing to do
        if min_area == 0.:
            return vertices, kinds
        ncseg = []
        nckind = []
        # iterate over vertices/kinds
        for k, (seg, kind) in enumerate(zip(vertices, kinds)):
            nseg = []
            nkind = []
            # iterate over single contours
            start = np.where(kind == 1)[0]
            stop = np.where(kind == 79)[0]
            for s0, s1 in zip(start, stop):
                # calculate area and extend
                if polyarea(seg[s0: s1]) >= min_area:
                    nseg.extend(seg[s0: s1])
                    nkind.extend(kind[s0: s1])
            # combine again
            if nseg:
                ncseg.append(np.vstack(nseg))
                nckind.append(nkind)
        return ncseg, nckind

    def lines(self, level):
        min_area = float(os.environ["MPL_MIN_CONTOUR_AREA"])
        # get vertices/kinds
        vertices, kinds = super().lines(level)
        # return early if nothing to do
        if min_area == 0.:
            return vertices, kinds

        ncseg = []
        nckind = []
        # iterate over vertices/kinds
        for k, (seg, kind) in enumerate(zip(vertices, kinds)):
            # calculate area and append
            if polyarea(seg) >= min_area:
                ncseg.append(seg)
                nckind.append(kind)
        return ncseg, nckind


# overwrite contour algorithm
contourpy._class_lookup["mpl2014"] = FilteredContourGenerator
