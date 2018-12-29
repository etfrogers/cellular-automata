#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:40:08 2017

@author: user
"""

from .automata.rule import Rule
from .automata.automata import Automata

gol = Automata([300, 300], [0, 1])

gol.add_rule(Rule(0, 1, 1, [3]))
gol.add_rule(Rule(1, 0, 1, [0, 1]))
gol.add_rule(Rule(1, 0, 1, [4, 8]))

gol.layout[1, 2] = 1
gol.layout[2, 2] = 1
gol.layout[3, 2] = 1
gol.randomise()
gol.show()

for ii in range(0, 50):
    gol.evolve()
    gol.show(str(ii))
