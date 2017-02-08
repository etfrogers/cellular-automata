#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:38:45 2017

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt


class Automata:
    rules = []
    layout = []
    states = []
    
    def __init__(self, size, states):
        self.layout = np.zeros(size)
        self.states = states
        
    def addRule(self, rule):
        self.rules.append(rule)
        
    def evolve(self):
        oldState = self.layout
        newState = oldState
        
        self.layout = newState
    
    def show(self):
        clims = (min(self.states), max(self.states))
        plt.imshow(self.layout, cmap="gray", clim=clims)
        
    def randomise(self):
        self.layout = np.random.choice(self.states, self.layout.shape)