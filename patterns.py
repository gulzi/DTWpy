#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def symmetric1(ts_x, ts_y, window, dist, pattern, normalized):    
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:        
        dt = dist(ts_x[i], ts_y[j])
        if i == j == 0:
            cost_matrix[i][j] = dt        
        else:
            cost_matrix[i][ j] = min(cost_matrix[i][ j - 1] + dt, cost_matrix[i - 1][ j] + dt, cost_matrix[i - 1][ j - 1] + dt)
    if normalized:return cost_matrix
    return cost_matrix
    

def symmetric2(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.empty([ts_x.size, ts_y.size])
    cost_matrix[:] = float('inf')
    for i, j in np.nditer(window):                
        dt = dist(ts_x[i], ts_y[j])
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][ j] = np.minimum(cost_matrix[i][ j - 1] + dt, cost_matrix[i - 1][ j] + dt, cost_matrix[i - 1][ j - 1] + 2 * dt)
            
    if normalized:cost_matrix[i][ j] = (cost_matrix[i][ j] / (len(ts_x) + len(ts_y)))
    return cost_matrix

def asymmetric(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:        
        dt = dist(ts_x[i], ts_y[j])
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:                
            cost_matrix[i][ j] = min(cost_matrix[i - 1][ j - 2] + dt, cost_matrix[i - 1][ j] + dt, cost_matrix[i - 1][j - 1] + dt)
            
    if normalized:cost_matrix[i][j] = (cost_matrix[i][j] / (len(ts_x)))
    return cost_matrix
    
def asymmetricP0(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:        
        dt = dist(ts_x[i], ts_y[j])
        if i == j == 0:
            cost_matrix[i][j] = dt 
        elif i == 0:
            cost_matrix[i][j] = cost_matrix[i][j - 1] + dt
        elif j == 0:
            cost_matrix[i][j] = cost_matrix[i - 1][ j] + dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i][j - 1], cost_matrix[i - 1][ j] + dt, cost_matrix[i - 1][ j - 1] + dt)
    if normalized:cost_matrix[i][ j] = (cost_matrix[i][j] / (len(ts_x)))
    return cost_matrix

def symmetricVelichkoZagoruyko(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:        
        dt = dist(ts_x[i], ts_y[j])
        if i == j == 0:
            cost_matrix[i][ j] = dt  # 0#2 * dt #float('inf')  #
        else:
            cost_matrix[i][ j] = min(cost_matrix[i ][ j - 1] + dt, cost_matrix[i - 1][ j - 1] + dt * .001, cost_matrix[i - 1][ j])
    
    return cost_matrix

def symmetricP1(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:        
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        if i - 1 < 0:
            dt2 = float('inf')
        if j - 1 < 0:
            dt1 = float('inf')
        if i == j == 0:
            cost_matrix[i][ j] = dt  # 0#2 * dt #float('inf')  #
        else:
            cost_matrix[i][ j] = min(cost_matrix[i - 1][ j - 2] + 2 * dt1 + dt, cost_matrix[i - 1][ j - 1] + 2 * dt, cost_matrix[i - 2][ j - 1] + 2 * dt2 + dt)
    if normalized:cost_matrix[i][ j] = (cost_matrix[i][ j] / (len(ts_x) + len(ts_y)))
    return cost_matrix

def asymmetricP1(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:        
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        if i - 1 < 0:
            dt2 = float('inf')
        if j - 1 < 0:
            dt1 = dist(ts_x[i], ts_y[j - 1])
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 1][ j - 2] + dt1 / 2 + dt / 2, cost_matrix[i - 1][ j - 1] + dt, cost_matrix[i - 2][ j - 1] + dt2 + dt)
    if normalized:cost_matrix[i][ j] = cost_matrix[i][ j] / (len(ts_x))
    return cost_matrix

def symmetricP2(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        dt3 = dist(ts_x[i - 1], ts_y[j - 2])
        dt4 = dist(ts_x[i - 2], ts_y[j - 1])
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        if  i - 2 < 0:
            dt4 = float('inf')
        if j - 2 < 0:
            dt3 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 2][ j - 3] + dt1 * 2 + dt + dt3 * 2, cost_matrix[i - 1][ j - 1] + 2 * dt, cost_matrix[i - 3][j - 2] + 2 * dt4 + 2 * dt2 + dt)
    if normalized:cost_matrix[i][j] = cost_matrix[i, j] / (len(ts_x) + len(ts_y))
    return cost_matrix
    
def asymmetricP2(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        dt3 = dist(ts_x[i - 1], ts_y[j - 2])
        dt4 = dist(ts_x[i - 2], ts_y[j - 1])
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        if  i - 2 < 0:
            dt4 = float('inf')
        if j - 2 < 0:
            dt3 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][ j] = min(cost_matrix[i - 2][ j - 3] + dt1 * (2 / 3) + dt * (2 / 3) + dt3 * (2 / 3), cost_matrix[i - 1][ j - 1] + dt, cost_matrix[i - 3][ j - 2] + dt4 + dt2 + dt)
    if normalized:cost_matrix[i][j] = cost_matrix[i][j] / (len(ts_x))
    return cost_matrix
    
def symmetricP05(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        dt3 = dist(ts_x[i], ts_y[j - 2])
        dt4 = dist(ts_x[i - 2], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        if  i - 2 < 0:
            dt4 = float('inf')
        if j - 2 < 0:
            dt3 = float('inf')
                    
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 1][ j - 3] + dt3 * 2 + dt1 + dt,
                                   cost_matrix[i - 1][ j - 2] + dt1 * 2 + dt,
                                   cost_matrix[i - 1][ j - 1] + 2 * dt,
                                   cost_matrix[i - 2][ j - 1] + 2 * dt2 + dt,
                                   cost_matrix[i - 3][ j - 1] + 2 * dt4 + dt2 + dt)
    if normalized:cost_matrix[i][ j] = cost_matrix[i, j] / (len(ts_x) + len(ts_y))
    return cost_matrix

def asymmetricP05(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        dt3 = dist(ts_x[i], ts_y[j - 2])
        dt4 = dist(ts_x[i - 2], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        if  i - 2 < 0:
            dt4 = float('inf')
        if j - 2 < 0:
            dt3 = float('inf')                    
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 1][ j - 3] + dt3 * (1 / 3) + dt1 * (1 / 3) + dt * (1 / 3) ,
                                   cost_matrix[i - 1][ j - 2] + dt1 * (1 / 2) + dt * (1 / 2),
                                   cost_matrix[i - 1][ j - 1] + dt,
                                   cost_matrix[i - 2][ j - 1] + dt2 + dt,
                                   cost_matrix[i - 3][ j - 1] + dt4 + dt2 + dt)
    if normalized:cost_matrix[i, j] = cost_matrix[i][ j][0] / (len(ts_x))
    return cost_matrix

def asymmetricItakura(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])       
        dt2 = dist(ts_x[i - 1], ts_y[j])                
        if j - 1 < 0 :
            dt2 = float('inf')                            
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][ j] = min(cost_matrix[i - 1][ j - 2] + dt ,
                                   cost_matrix[i - 1][ j - 1] + dt,
                                   cost_matrix[i - 2][ j - 1] + dt2 + dt,
                                   cost_matrix[i - 2][ j - 2] + dt2 + dt)    
    return cost_matrix
    
def typeIa(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 2][ j - 1] + dt2,
                                   cost_matrix[i - 1][j - 1] + dt,
                                   cost_matrix[i - 1][ j - 2] + dt1)    
    return cost_matrix

def typeIb(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 2][ j - 1] + dt2 + dt ,
                                   cost_matrix[i - 1][j - 1] + dt,
                                   cost_matrix[i - 1][ j - 2] + dt1 + dt)
    return cost_matrix

def typeIc(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 2][ j - 1] + dt2 + dt ,
                                   cost_matrix[i - 1][ j - 1] + dt,
                                   cost_matrix[i - 1][j - 2] + dt1)
    if normalized:cost_matrix[i][j] = cost_matrix[i, j] / (len(ts_x))
    return cost_matrix

def typeId(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 2][ j - 1] + dt2 * 2 + dt  ,
                                   cost_matrix[i - 1][j - 1] + dt * 2,
                                   cost_matrix[i - 1][ j - 2] + dt1 * 2 + dt)
    if normalized:cost_matrix[i][j] = cost_matrix[i][ j] / (len(ts_x) + len(ts_y))
    return cost_matrix

def typeIas(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][ j] = min(cost_matrix[i - 2][ j - 1] + dt2 / 2 + dt / 2 ,
                                   cost_matrix[i - 1][j - 1] + dt,
                                   cost_matrix[i - 1][ j - 2] + dt1 / 2 + dt / 2)
    if normalized:cost_matrix[i][j] = cost_matrix[i][ j] / (len(ts_x))
    return cost_matrix

def typeIbs(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][ j] = min(cost_matrix[i - 2][ j - 1] + dt2 + dt ,
                                   cost_matrix[i - 1][ j - 1] + dt,
                                   cost_matrix[i - 1][ j - 2] + dt1 + dt)
    if normalized:cost_matrix[i][j] = cost_matrix[i][ j] / (len(ts_x))
    return cost_matrix

def typeIcs(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 2][ j - 1] + dt2 + dt,
                                   cost_matrix[i - 1][j - 1] + dt,
                                   cost_matrix[i - 1][j - 2] + dt1 / 2 + dt / 2)
    if normalized:cost_matrix[i][j] = cost_matrix[i][j] / (len(ts_x))
    return cost_matrix

def typeIds(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])
        dt1 = dist(ts_x[i], ts_y[j - 1])
        dt2 = dist(ts_x[i - 1], ts_y[j])
        
        if i - 1 < 0 :
            dt2 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 2][ j - 1] + dt2 * 1.5 + dt * 1.5 ,
                                   cost_matrix[i - 1][ j - 1] + dt * 2,
                                   cost_matrix[i - 1][j - 2] + dt1 * 1.5 + dt * 1.5)
    if normalized:cost_matrix[i][j] = cost_matrix[i][j] / (len(ts_x) + len(ts_y))
    return cost_matrix

def typeIIa(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])        
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][ j] = min(cost_matrix[i - 1][ j - 1] + dt,
                                    cost_matrix[i - 2][j - 1] + dt,
                                   cost_matrix[i - 1][ j - 2] + dt)
    if normalized:cost_matrix[i][j] = cost_matrix[i][ j] / (len(ts_x))
    return cost_matrix
def typeIIb(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])        
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 1][ j - 1] + dt,
                                    cost_matrix[i - 2][ j - 1] + dt * 2 ,
                                   cost_matrix[i - 1][ j - 2] + dt * 2)
    if normalized:cost_matrix[i][j] = cost_matrix[i][j] / (len(ts_x))
    return cost_matrix
def typeIIc(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])        
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 1][ j - 1] + dt,
                                    cost_matrix[i - 1][ j - 2] + dt,
                                   cost_matrix[i - 2][ j - 1] + dt * 2)
    if normalized:cost_matrix[i][j] = cost_matrix[i][j] / (len(ts_x))
    return cost_matrix

def typeIId(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])        
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 1][ j - 1] + dt * 2,
                                    cost_matrix[i - 1][ j - 2] + dt * 3 ,
                                   cost_matrix[i - 2][ j - 1] + dt * 3)
    if normalized:cost_matrix[i][ j] = cost_matrix[i][ j] / (len(ts_x) + len(ts_y))
    return cost_matrix
def typeIIIc(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])        
        dt1 = dist(ts_x[i - 1], ts_y[j])
        
        if i - 1 < 0 :
            dt1 = float('inf')
        if i == j == 0:
            cost_matrix[i][ j] = dt  # 0#2 * dt #float('inf')  #
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 1][ j - 2] + dt,
                                    cost_matrix[i - 1][ j - 1] + dt,
                                    cost_matrix[i - 2][j - 1] + dt1 + dt,
                                   cost_matrix[i - 2][j - 2] + dt1 + dt)
    if normalized:cost_matrix[i][j] = cost_matrix[i][j] / (len(ts_x))
    return cost_matrix

def typeIVc(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])        
        dt1 = dist(ts_x[i - 1], ts_y[j])
        dt2 = dist(ts_x[i - 2], ts_y[j])
        
        if i - 1 < 0 :
            dt1 = float('inf')
        if i - 2 < 0 :
            dt2 = float('inf')
        if i == j == 0:
            cost_matrix[i][ j] = dt
        else:
            cost_matrix[i][ j] = min(cost_matrix[i - 1][ j - 1] + dt ,
                                    cost_matrix[i - 1][j - 2] + dt,
                                    cost_matrix[i - 1][ j - 3] + dt ,
                                    cost_matrix[i - 2][ j - 1] + dt1 + dt,
                                    cost_matrix[i - 2][ j - 2] + dt1 + dt,
                                    cost_matrix[i - 2][j - 3] + dt1 + dt,
                                    cost_matrix[i - 3][ j - 1] + dt2 + dt1 + dt,
                                    cost_matrix[i - 3][j - 2]+ dt2 + dt1 + dt,
                                   cost_matrix[i - 3][ j - 3] + dt2 + dt1 + dt )
    if normalized:cost_matrix[i][j] = cost_matrix[i][j] / (len(ts_x))
    return cost_matrix

def mori2006(ts_x, ts_y, window, dist, pattern, normalized):
    cost_matrix = np.zeros([len(ts_x), len(ts_y)])
    cost_matrix[:] = float('inf')
    for i, j in window:
        dt = dist(ts_x[i], ts_y[j])        
        dt1 = dist(ts_x[i - 1], ts_y[j])
        dt2 = dist(ts_x[i], ts_y[j - 1])
        if i - 1 < 0 :
            dt1 = float('inf')
        if j - 1 < 0 :
            dt2 = float('inf')
        if i == j == 0:
            cost_matrix[i][j] = dt
        else:
            cost_matrix[i][j] = min(cost_matrix[i - 2][ j - 1] + dt1 * 2 + dt,
                                    cost_matrix[i - 1][ j - 1] + dt * 3,
                                   cost_matrix[i - 1][ j - 2] + dt2 * 3 + dt * 3 )
    if normalized:cost_matrix[i][j] = cost_matrix[i][ j] / (len(ts_y))
    return cost_matrix
