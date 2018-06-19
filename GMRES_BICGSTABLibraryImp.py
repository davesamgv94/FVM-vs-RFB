# -*- coding: utf-8 -*-
"""
Created on Thu May 31 00:30:23 2018

@author: saiby
"""

#Librerias de Python
import numpy as np

#Algoritmo de GMRES desarrollado
def GMRes(A, b, x0, e, nmax_iter, restart=None):
    """
    Esta funcion realiza el calculo del algoritmo de GMRES con los siguientes 
    parametros A: la matriz real o compleja NxN del sistema lineal,
    b: Lado derecho del sistema lineal. Tiene forma (N,) o (N, 1)
    e: Tolerancia para lograr. El algoritmo termina cuando el residual relativo o absoluto 
    está por debajo de e.
    nmax_iter: Número máximo de iteraciones (ciclos de reinicio). La iteración se detendrá 
    después de los pasos del maxiter aunque no se haya alcanzado la tolerancia especificada
    restart: Número de iteraciones entre reinicios. Los valores más grandes aumentan 
    el costo de la iteración, pero pueden ser necesarios para la convergencia.
    """  
    r = b - np.asarray(np.dot(A, x0)).reshape(-1)
    x = []
    q = [0] * (nmax_iter)
    x.append(r)
    q[0] = r / np.linalg.norm(r)
    h = np.zeros((nmax_iter + 1, nmax_iter))

    for k in range(min(nmax_iter, A.shape[0])):
        y = np.asarray(np.dot(A, q[k])).reshape(-1)

        for j in range(k+1):
            h[j, k] = np.dot(q[j], y)
            y = y - h[j, k] * q[j]
        h[k + 1, k] = np.linalg.norm(y)
        if (h[k + 1, k] != 0 and k != nmax_iter - 1):
            q[k + 1] = y / h[k + 1, k]

        b = np.zeros(nmax_iter + 1)
        b[0] = np.linalg.norm(r)

        result = np.linalg.lstsq(h, b)[0]

        x.append(np.dot(np.asarray(q).transpose(), result) + x0)

    return x

#Algoritmo de BICGSTAB desarrollado
def bicgstab(A,b,x0,e,nmax_iter,restart=True):
    """
    Esta funcion realiza el calculo del algoritmo de BICGSTAB con los siguientes 
    parametros A: la matriz real o compleja NxN del sistema lineal,
    b: Lado derecho del sistema lineal. Tiene forma (N,) o (N, 1)
    e: Tolerancia para lograr. El algoritmo termina cuando el residual relativo o absoluto 
    está por debajo de e.
    nmax_iter: Número máximo de iteraciones (ciclos de reinicio). La iteración se detendrá 
    después de los pasos del maxiter aunque no se haya alcanzado la tolerancia especificada
    restart: Número de iteraciones entre reinicios. Los valores más grandes aumentan 
    el costo de la iteración, pero pueden ser necesarios para la convergencia..
    """  
    
    r = b - np.dot(A,x0)
    rc = r
    p = rc-r
    v = rc-r
    rho = 1
    alpha = 1
    omega = 1
    k = 0
    snormv = np.zeros(nmax_iter)
    while (k < nmax_iter) :
        rhop1 = float(np.dot(rc.T,r))
        beta = (rhop1/rho)*(alpha/omega)
        p = r + beta*(p - omega*v)
        v = np.dot(A,p)
        alpha = rhop1/float(np.dot(rc.T,v))
        s = r - alpha*v
        t = np.dot(A,s)
        omega = float(np.dot(t.T,s))/float(np.dot(t.T,t))
        x0 += alpha*p + omega*s
        snorm = np.linalg.norm(s)
        snormv[k] = snorm
        k += 1
        
        if snorm < e:
            break
        r = s-omega*t
        rho = rhop1
    
    if restart:
        return x0, snorm, k, snormv
    else:
        return x0
