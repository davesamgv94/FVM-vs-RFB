#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 20:03:13 2018

@author: 
"""
#Librerias de Python
import numpy as np
import Knots as kn
import RadialBasisFunctions as RBF
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla
import time
#Ingresar el kernel a resolver
kernel = input("Escribe el Kernel que deseas que resuelva como Multicuadratico(MQ) o Multricuadratico Inverso (IMQ) = ")
#Inicio del tiempo a resolver
t1 =  time.time()
#Funcion que realiza la matriz de Gramm
def fillGrammMatrix(G, line, rbf, vel,D):
    """
    Esta funcion realiza el llenado de la matriz de Grimm con los siguientes 
    parametros:
    D: Matriz 2-D con unos en la diagonal y ceros en otra parte,
    line: division de la linea (1D) de la cantidad de nodos
    rbf: kernel a escoger en RBF
    vel: variable de velocidad
    D:   equivalente a Gamma en FVM 
    """
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
            if kernel == "MQ":
                G[i,j] = rbf.mq(x[i], x[j])
            elif kernel == "IMQ":
                G[i,j] = rbf.inmq(x[i], x[j])
                
#Solucion Analitica
def analyticSol(x):
    return (np.exp(rho * u * x / D) - 1) / (np.exp(rho * u * L / D) - 1) * (phiB - phiA) + phiA

            
# Datos del Problema (Variables)
D =  0.1 #Gamma
rho = 1.0 # kg/m^3
L = 1.0 # m
u = 0.1 # m/s
phiA = 1 #
phiB = 0 #
N = 50 # Número de nodos
vel = 0.1 # m/s
# declaracion de la cantidad de nodos en una dimension
linea = kn.Line(N, L)
linea.homogeneous()
# Para la solucion Analitica      
u_a = analyticSol(linea.x())
#Matriz identidad
G = np.eye(N)
#Creacion de matrices ceros
u = np.zeros(N)
lam = np.zeros(N)

# Condiciones de Frontera
u[linea.NI()] = phiA
u[-1] = phiB

#Condicionador para los kernels
c = 1.0/np.sqrt(N)
#Seleccion del Kernel a utilizar
if kernel == "MQ":
    rbf = RBF.Multiquadric(c*5)            
elif kernel == "IMQ":
    rbf = RBF.InMultiquadric(c*5)

#Llenado de la matriz de Gramm
fillGrammMatrix(G, linea, rbf, vel,D)
#Impresion de la matriz de GRIMM
print(G)
#Solucion con el metodo de GMRES
lam = spla.gmres(G,u)
#Solucion final ya con RBF
u = RBF.evaluate(linea, lam[0], rbf,kernel)

#Parametro en x
zzz = linea.x()
#Error analitica vs RBF
error = np.absolute(u_a - u)
#Termino del tiempo que tardo en resolver e impresion
t2 = time.time()
print("Tiempo Realizado con GMRES de la libreria de Python en 50 Nodos")
print('GMRES: ' + str((t2 - t1)) + "\n")
print('||Error|| = ', np.linalg.norm(error))


#Impresion de la solucion junto con la graficas

plt.title('Solución de $\partial(p u \phi)/\partial x= \partial (\Gamma \partial\phi/\partial x)/\partial x$ con RBF')
plt.plot(zzz,u_a, 'bs', label = 'Sol. analítica') 
plt.plot(linea.x(),u,'r^', label = 'Sol. RBF')
plt.xlabel('$x$ [m]')
plt.ylabel('$U$ [...]')
plt.legend()
plt.grid()
plt.show()
