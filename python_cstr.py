# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-20 15:01:08
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-04-27 15:23:47

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate 
import csv

Qout = []
with open('outflow.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=', ')
    for row in reader:
        print(', '.join(row))
        #Qout.append(', '.join(row))

def tank(t, y):
    # Read from user dictionary
    Cin = 10
    Qin = 5
    Qout = sim._model.getNodeResult('Tank', 1)
    V = sim._model.getNodeResult('Tank', 3)
    k = 0.01
    # Setup variables for ODE solver
    C = y[0]
    n = len(y)
    dCdt = np.zeros((n,1))
    dCdt[0] = (Qin*Cin - Qout*C)/V - k*C
    return dCdt