#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:28:22 2026

@author: syedtariqishtiaq
"""

# Testing out things that come up in the code

x = 3
lst1 = [1] + [x]*4
print(lst1)


# Partial functions
from functools import partial

def f(a,b,c,x):
    return 100*a + 10*b + 2*c + x
g = partial(f, 2, 3, 4)

# So g is a partial function of f, calling f using a=2, b=3
# and c=4. 
print(g(3))

def g1(v):
    return 100*2 + 10*3 + 2*4 + v

print(g1(69))
print(g(69) == g1(69))




import os
os.chdir('/Users/syedtariqishtiaq/Desktop/UCL Coding Project 2026/DOFitting')
print(os.getcwd())















