#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:38:45 2017

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt

plt.ion()


class Automata:
    # variables defined here are static
    nhood = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])

    def __init__(self, size, states):
        # variables defined here are instance variables
        assert (len(size) == 2)  # only implemented 2d automata
        self.layout = np.zeros(size) + states[0]
        self.states = states
        self.rules = []
        self._fig_handle = None
        self._im_handle = None
        self.do_pause = True

    def add_rule(self, rule):
        self.rules.append(rule)

    def evolve(self, n=1):
        for _ in range(n):
            self._evolve()

    def _evolve(self):
        old_state = self.layout
        new_state = np.copy(old_state)

        for rule in self.rules:
            rule.apply(new_state, old_state, self.nhood)
        # print('out')
        # print(new_state)
        self.layout = new_state

    def show(self, label=''):
        clims = (min(self.states), max(self.states))

        if self._fig_handle is None:
            self._fig_handle = plt.figure(1)
            self._im_handle = plt.imshow(self.layout, cmap="gray", clim=clims, animated=True, interpolation='none')
        else:
            plt.figure(self._fig_handle.number)
            self._im_handle.set_array(self.layout)
        plt.title(label)
        plt.draw()
        if self.do_pause:
            plt.pause(0.001)

    def randomise(self):
        self.layout = np.random.choice(self.states, self.layout.shape)

    @property
    def max_nbours(self):
        return np.sum(np.sum(self.nhood))

    @property
    def fig_handle(self):
        return self._fig_handle
