#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 18:33:23 2018

@author: luiggi
"""

import numpy as np

class Line():
    
    def __init__(self, N, L):
        self.__N = N
        self.__NB = 2
        self.__NI = N - self.__NB
        self.__L = L
        self.__x = np.zeros(N)
        
    def N(self):
        return self.__N
    
    def NB(self):
        return self.__NB
    
    def NI(self):
        return self.__NI

    def L(self):
        return self.__L
    
    def x(self):
        return self.__x
    
    def homogeneous(self):
        h = self.__L / (self.__N - 1.0)
        for i in range(self.__NI):
            self.__x[i] = (i+1) * h
        self.__x[self.__NI] = 0
        self.__x[-1] = self.__L
        
if __name__ == '__main__':

    N = 11
    L = 1.0
    linea = Line(N, L)
    linea.homogeneous()    
    
    print(linea.N(), linea.NB(), linea.NI(), linea.L())
    print(linea.x())
    
    import matplotlib.pyplot as plt
    plt.plot(linea.x(),'o-')
    plt.grid()
    plt.show()    
        
        