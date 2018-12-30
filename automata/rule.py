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
        self.state_neighbours = neighbours
        self.n_neighbours = n

        self.prob = prob

    def apply(self, new_state, ref_state, nhood):
        # new_state is the current state of the layout and should be modified
        # by this function
        # ref_state is the state of the layout at the start of the timestep and
        # should be used to calculate which changes to make

        # may need to add wrapping

        valid_inds = self.get_valid_inds(nhood, ref_state, self.state_from, self.state_neighbours, self.n_neighbours)
        valid_inds = self.apply_prob(ref_state, valid_inds)
        new_state[valid_inds] = self.state_to

    def apply_prob(self, ref_state, valid_inds):
        if self.prob < 1:
            valid_inds = np.logical_and(valid_inds, np.random.random(ref_state.shape) < self.prob)
        return valid_inds

    @staticmethod
    def get_valid_inds(nhood, ref_state, state_from, state_neighbours, n_neighbours):
        try:
            _ = [a for a in state_neighbours]
            state_neighbours = state_neighbours
        except TypeError:
            state_neighbours = [state_neighbours]
        from_vals = reduce(np.logical_or, [ref_state == st for st in state_neighbours])
        nbours = signal.convolve(from_vals.astype(int), nhood, 'same')
        valid_inds = np.logical_and(nbours >= min(n_neighbours), nbours <= max(n_neighbours))
        valid_inds = np.logical_and(valid_inds, ref_state == state_from)
        return valid_inds


class MultiRule(Rule):
    def __init__(self, from_, to, neighbours, n, prob=1, combine_function=None):
        super().__init__(from_, to, neighbours, n, prob)
        self.combine_function = combine_function

    def apply(self, new_state, ref_state, nhood):
        valid_inds = [self.get_valid_inds(nhood, ref_state, state_from, state_neighbours, n_neighbours) for
                      state_from, state_neighbours, n_neighbours in
                      zip(self.state_from, self.state_neighbours, self.n_neighbours)]
        valid_inds = reduce(self.combine_function, valid_inds)
        valid_inds = self.apply_prob(ref_state, valid_inds)
        new_state[valid_inds] = self.state_to
