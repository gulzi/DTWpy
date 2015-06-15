#!/usr/bin/env python
# -*- coding: utf-8 -*-


import patterns
import numpy as np
from memprof import *
def __derivation(ts):
        drts = []

        for i in range(1, len(ts) - 1):
            derived = (ts[i] - ts[i - 1]) + ((ts[i + 1] - ts[i - 1]) / 2.0) / 2.0
            drts.append(derived)

        return drts

def ddtw(x, y, dist=lambda a, b: abs(a - b), windowtype=None, windowsize=None, pattern="symmetric1", normalized=False, dist_only=False):
    x = __derivation(x)
    y = __derivation(y)
    len_x, len_y = len(x), len(y)
    window = ''
    if windowtype == "scband":
        window = __getSCwindow(windowsize, len_x, len_y)
    elif windowtype == "itakura":
        window = __getItakurawindow(len_x, len_y)
    elif windowtype == "paliwal":
        window = __getpaliwalwindow(windowsize, len_x, len_y)
    else:
        window = [(i, j) for i in range(len_x) for j in range(len_y)]   
    D = __cost_matrix(x, y, window, dist, pattern, normalized)
    if dist_only:        
        return D[len_x - 1][ len_y - 1]
    else:
        path = __backtrack(D, len_x - 1, len_y - 1)
    
    return (window,D[len_x - 1][ len_y - 1], path)

def dtw(x, y, dist=lambda a, b: abs(a - b), windowtype=None, windowsize=None, pattern="symmetric1", normalized=False,dist_only=False,cost=False):
    len_x, len_y = x.size, y.size
    np.empty([len_x, len_y], dtype=float)
    if windowtype == "scband":
        window = __getSCwindow(windowsize, len_x, len_y)
    elif windowtype == "itakura":
        window = __getItakurawindow(len_x, len_y)
    elif windowtype == "paliwal":
        window = __getpaliwalwindow(windowsize, len_x, len_y)
    else:
        window = [(i, j) for i in np.arange(len_x) for j in np.arange(len_y)] 
    D = __cost_matrix(x, y, window, dist, pattern, normalized)
    if dist_only:        
        return D[len_x - 1][ len_y - 1]
    else:
        path = __backtrack(D, len_x - 1, len_y - 1)
    if cost == False:
        return (D[len_x - 1][ len_y - 1], path)
    return (window,D[len_x - 1][ len_y - 1], path)
    

def fastdtw(x, y, radius=1, dist=lambda x, y:abs(x - y), pattern="symmetric1", normalized=False):
    min_time_size = radius + 2

    if len(x) < min_time_size or len(y) < min_time_size:
        return constrained_dtw(x, y, dist, pattern, normalized)

    x_shrinked = __reduce_by_half(x)
    y_shrinked = __reduce_by_half(y)
    _, path = fastdtw(x_shrinked, y_shrinked, radius, dist, pattern, normalized)
    window = __expand_window(path, len(x), len(y), radius)
    return constrained_dtw(x, y, dist, pattern, normalized, window)    

def constrained_dtw(x, y, dist, pattern, normalized, window=None):
    len_x, len_y = len(x), len(y)
    if window is None:
        window = [(i, j) for i in range(len_x) for j in range(len_y)]
        
    D = __cost_matrix(x, y, window, dist, pattern, normalized)
    path = __backtrack(D, len_x - 1, len_y - 1)
    
    return (D[len_x - 1][ len_y - 1], path)

def __cost_matrix(ts_x, ts_y, window, dist, pattern, normalized):
    if pattern == "symmetric1":
        return patterns.symmetric1(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "symmetric2" or pattern == "symmetricP0":
        return patterns.symmetric2(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "asymmetric":
        return patterns.asymmetric(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "symmetricP1":
        return patterns.symmetricP1(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "asymmetricP0":
        return patterns.asymmetricP0(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "asymmetricP1":
        return patterns.asymmetricP1(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "symmetricP2":
        return patterns.symmetricP2(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "asymmetricP2":
        return patterns.asymmetricP2(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "symmetricP05":
        return patterns.symmetricP05(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "asymmetricP05":
        return patterns.asymmetricP05(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "asymmetricItakura":
        return patterns.asymmetricItakura(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIa":
        return patterns.typeIa(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIb":
        return patterns.typeIb(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIc":
        return patterns.typeIc(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeId":
        return patterns.typeId(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIas":
        return patterns.typeIas(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIbs":
        return patterns.typeIbs(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIcs":
        return patterns.typeIcs(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIds":
        return patterns.typeIds(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIIa":
        return patterns.typeIIa(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIIb":
        return patterns.typeIIb(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIIc":
        return patterns.typeIIc(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIId":
        return patterns.typeIId(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIIIc":
        return patterns.typeIIIc(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "typeIVc":
        return patterns.typeIVc(ts_x, ts_y, window, dist, pattern, normalized)
    elif pattern == "mori2006":
        return patterns.mori2006(ts_x, ts_y, window, dist, pattern, normalized)
    
    
    
def __backtrack(D, max_x, max_y):
    path = []
    
    i, j = max_x, max_y
    path.append((i, j))
    while i > 0 or j > 0:
        diag_cost = float('inf')
        left_cost = float('inf')
        down_cost = float('inf')
        
        if (i > 0) and (j > 0): 
            diag_cost = D[i - 1][ j - 1]
        if i > 0: 
            left_cost = D[i - 1][ j ]
        if j > 0:             
            down_cost = D[i][ j - 1]
        
        if (diag_cost <= left_cost and diag_cost <= down_cost):
            i, j = i - 1, j - 1            
        elif (left_cost < diag_cost and left_cost < down_cost):
            i = i - 1
        elif (down_cost < diag_cost and down_cost < left_cost):
            j = j - 1
        elif i <= j:
            j = j - 1
        else:
            i = i - 1
        
        path.append((i, j))

    path.reverse()
    return path



def __reduce_by_half(x):
    return [(x[i // 2] + x[1 + i // 2]) / 2 for i in range(0, len(x), 2)]


def __expand_window(path, len_x, len_y, radius):
    path_ = set(path)
    for i, j in path:

        for a, b in ((i + a, j + b) for a in range(-radius, radius + 1) for b in range(-radius, radius + 1)):
            path_.add((a, b))

    window_ = set()
    for i, j in path_:
        for a, b in ((i * 2, j * 2), (i * 2, j * 2 + 1), (i * 2 + 1, j * 2), (i * 2 + 1, j * 2 + 1)):
            window_.add((a, b))

    window = []
    start_j = 0
    for i in range(0, len_x):
        new_start_j = None
        for j in range(start_j, len_y):
            if (i, j) in window_:
                window.append((i, j))
                if new_start_j is None:
                    new_start_j = j
            elif new_start_j is not None:
                break
        start_j = new_start_j

    return window

def __getSCwindow(windowsize, len_x, len_y):
    window = []
    for i in range(len_x):
        for j in range(len_y):
            if abs(i - j) < windowsize:
                window.append((i, j))
    return window

def __getItakurawindow(len_x, len_y):
    window = []
    for i in range(len_x):
        for j in range(len_y):
            if  ((j < 2 * i) and (i <= 2 * j) and (i >= len_x - 1 - 2 * (len_y - j)) and (j > len_y - 1 - 2 * (len_x - i))):                
                window.append((i, j))
    return window

def __getpaliwalwindow(windowsize, len_x, len_y):
    window = []
    for i in range(len_x):
        for j in range(len_y):
            if abs(i*len_y/len_x - j) < windowsize:
                window.append((i, j))
    return window

