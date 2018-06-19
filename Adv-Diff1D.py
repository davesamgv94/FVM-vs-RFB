#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 20:03:13 2018

@author: 
"""

import numpy as np
import Knots as kn
import RadialBasisFunctions as RBF
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla

def fillGrammMatrix(G, line, rbf, vel,D):
    N = line.N()
    NI = line.NI()
    x = line.x()
    
# ----- Wl matrix
    for j in range(N):
        for i in range(NI):
            G[i,j] = ( vel * rbf.d1x(x[i], x[j]) - \
                           D * rbf.d2x(x[i], x[j]) )
# ----- Wb matrix
    for j in range(N):
        for i in range(NI,N):
            G[i,j] = rbf.mq(x[i], x[j])

def analyticSol(x):
    return (np.exp(rho * u * x / D) - 1) / (np.exp(rho * u * L / D) - 1) * (phiB - phiA) + phiA

            
# Datos del Problema (Variables)
D =  0.1
rho = 1.0 # kg/m^3
L = 1.0 # m
u = 0.1 # m/s
phiA = 1 #
phiB = 0 #
N = 50 # Número de nodos

max_iter = 500
vel = 0.1
linea = kn.Line(N, L)
linea.homogeneous()
    
phi_a = analyticSol(linea.x())

G = np.eye(N)

u = np.zeros(N)
lam = np.zeros(N)

# Boundary conditions
u[linea.NI()] = phiA
u[-1] = phiB

c = 1.0/np.sqrt(N)
rbf = RBF.Multiquadric(c*5)
    
fillGrammMatrix(G, linea, rbf, vel,D)

print(G)

lam = spla.bicgstab(G,u)
#lam = np.linalg.solve(G,u)
u = RBF.evaluate(linea, lam[0], rbf)

zzz = linea.x()

#plt.plot(linea.x(),u,'s')
plt.title('Solución de $\partial(p u \phi)/\partial x= \partial (\Gamma \partial\phi/\partial x)/\partial x$ con RBF')
plt.plot(zzz[0:N-2],phi_a[0:N-2], '-s', label = 'Sol. analítica') 
plt.plot(linea.x(),u,'s', label = 'Sol. RBF')
plt.legend()
plt.grid()
plt.show()
