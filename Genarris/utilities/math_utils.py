
"""
Created on Sun Apr  1 16:29:48 2018

@author: timcrose
"""

import math, random

def mean(lst):
    return sum(lst) / float(len(lst))

def randrange_float(start, stop, step, num_decimal_places=4):
    return round(random.randint(0, int((stop - start) / step)) * step + start, num_decimal_places)
