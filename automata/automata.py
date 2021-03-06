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
    
    #variables defined here are static
    nhood = np.array([[1,1,1],[1,0,1],[1,1,1]])
    
    def __init__(self, size, states):
        #variables defined here are instance variables
        assert(len(size) == 2) #only implmented 2d automata
        self.layout = np.zeros(size) + states[0]
        self.states = states
        self.rules = []
        self._figHandle = None
        self._imHandle = None
        self.doPause = True
        
    def addRule(self, rule):
        self.rules.append(rule)
        
    def evolve(self):
        oldState = self.layout
        newState = np.copy(oldState)
        
        for rule in self.rules:
            rule.apply(newState, oldState, self.nhood)
        #print('out')
        #print(newState)
        self.layout = newState
    
    def show(self, label =''):
        clims = (min(self.states), max(self.states))
        
        if self._figHandle == None:
            self._figHandle = plt.figure(1)
            self._imHandle = plt.imshow(self.layout, cmap="gray", clim=clims, animated=True, interpolation='none')
        else:
            plt.figure(self._figHandle.number)
            self._imHandle.set_array(self.layout)            
        plt.title(label)
        plt.draw()
        if self.doPause:
            plt.pause(0.001)
        
    def randomise(self):
        self.layout = np.random.choice(self.states, self.layout.shape)
        
    @property 
    def maxNbours(self):
        return np.sum(np.sum(self.nhood))
    @property 
    def figHandle(self):
        return self.figHandle
    

    
    