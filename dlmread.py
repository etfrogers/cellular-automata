#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 13:29:48 2017

@author: user
"""

import csv
import numpy as np

def readnfile(fname, dl):
    arr = [];
    with open(fname) as f:
        reader = csv.reader(f, delimiter=dl)
        for row in reader:
            if row != []:
                arr.append([float(num) for num in row])
    return np.array(arr)