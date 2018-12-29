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

    def __init__(self, from_, to, neighbours, n, prob=1):
        self.state_from = from_
        self.state_to = to
        try:
            _ = [a for a in neighbours]
            self.state_neighbours = neighbours
        except TypeError:
            self.state_neighbours = [neighbours]
            # assert(len(n) == 1 or len(n) == 2)
        self.n_neighbours = n

        self.prob = prob

    def apply(self, new_state, ref_state, nhood):
        # new_state is the current state of the layout and should be modified
        # by this function
        # ref_state is the state of the layout at the start of the timestep and
        # should be used to calculate which changes to make

        # may need to add wrapping

        from_vals = reduce(np.logical_or, [ref_state == st for st in self.state_neighbours])
        nbours = signal.convolve(from_vals.astype(int), nhood, 'same')
        # print(ref_state)
        # print('neighouurs')
        # print(nbours)
        valid_inds = np.logical_and(nbours >= min(self.n_neighbours), nbours <= max(self.n_neighbours))
        valid_inds = np.logical_and(valid_inds, ref_state == self.state_from)
        if self.prob < 1:
            valid_inds = np.logical_and(valid_inds, np.random.random(ref_state.shape) < self.prob)
        new_state[valid_inds] = self.state_to
        # print(new_state)
