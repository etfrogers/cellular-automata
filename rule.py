#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 19:09:12 2017

@author: Edward Rogers
"""

class Rule:
    stateFrom = []
    stateTo = []
    Nneighbours = 0;
    prob = -1;
    
    
    def __init__(self, sFrom, sTo, N, prob=1):
        self.stateFrom = sFrom
        self.stateTo = sTo
        self.Nneighbours = N
        self.prob = prob