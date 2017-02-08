#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:40:08 2017

@author: user
"""

from rule import Rule
from automata import Automata


gol = Automata([50,50], [0,1])

gol.addRule(Rule(0,1,1,[3]))
gol.addRule(Rule(1,0,1,[0,1]))
gol.addRule(Rule(1,0,1,[4,8]))

gol.layout[1,2] =1
gol.layout[2,2] =1
gol.layout[3,2] =1
gol.randomise()
gol.show()


for ii in range(0,10):
    
    gol.evolve()
    gol.show()
    