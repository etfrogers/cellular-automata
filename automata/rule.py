#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 19:09:12 2017

@author: Edward Rogers
"""
from scipy import signal
import numpy as np
from functools import reduce

class Rule:
    
    def __init__(self, sFrom, sTo, sNb, N, prob=1):
        self.stateFrom = sFrom
        self.stateTo = sTo
        try:
            _ = [a for a in sNb]
            self.stateNbours = sNb
        except TypeError:
            self.stateNbours = [sNb]            
        #assert(len(N) == 1 or len(N) == 2)
        self.Nneighbours = N
        
        self.prob = prob
        
    def apply(self, newState, refState, nhood):
        # newState is the current state of the layout and should be modified 
        # by this function
        # refState is the state of the layout at the start of the timestep and 
        # should be used to calculate which changes to make
        
        fromVals = reduce(np.logical_or, [refState == st for st in self.stateNbours]) 
        nbours = signal.convolve(fromVals.astype(int), nhood, 'same')
        #print(refState)
        #print('neighouurs')
        #print(nbours)
        validInds = np.logical_and(nbours>=min(self.Nneighbours), nbours<=max(self.Nneighbours))
        validInds = np.logical_and(validInds, refState==self.stateFrom) 
        if self.prob == 1:
            newState[validInds] = self.stateTo
        else:
            raise NotImplementedError('Probability not yet implemented')
        #print(newState)