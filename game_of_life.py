#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:40:08 2017

@author: user
"""

from rule import Rule
from automata import Automata


gol = Automata([50,50], [0,1])

gol.addRule(Rule(0,1,3))
gol.addRule(Rule(1,0,[0,2]))
gol.addRule(Rule(1,0,[4,9]))

gol.randomise()
gol.show()