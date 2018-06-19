"""
Created on Tue Mar  6 15:11:05 2018

@author: Dr. Luis Miguel de la Cruz
adaptado por: Samuel Garcia

"""
### Librerias de Python ###
import numpy as np
import time
from pandas import DataFrame

# Importar cada uno los Scripts donde se encuentra los diferentes metodos
# ya se de los coeficientes, parte difusiva o de adveción, creación de la matriz, etc.
from Mesh import Mesh
from Coefficients import Coefficients
from Diffusion import Diffusion1D
from Advection import Advection1D
from Matrix import Matrix

# Declaración de la función Cronometro para la toma del tiempo de Ejecución
def crono(f):
 	"""
 	Regresa el tiempo que toma en ejecutarse la funcion.
 	"""
 	def eTime(A,b):
 		t1 = time.time()
 		f(A,b)
 		t2 = time.time()
 		return 'Elapsed time: ' + str((t2 - t1)) + "\n"
 	return eTime

### Impresion de Portada en Consola####
def decorate(f):
    def nicePrint(**kargs):
        line = '-' * 70
        print('.'+ line + '.')
        print('|{:^70}|'.format('Tarea 1'))
        print('.'+ line + '.')
        print('|{:^70}|'.format(' Autor: Luis M. De la Cruz, Adaptaciones: Samuel Garcia'))
        print('.'+ line + '.')
        f(**kargs)
        print('.'+ line + '.')
    return nicePrint
 
#### Impresion de Datos Ingresados con Formato #####   
@decorate
def printData(**kargs):
	for (key,value) in kargs.items():
		print('|{:^70}|'.format('{0:>15s} = {1:10.5e}'.format(key, value)))

# Calculo el error porcentual y agrego al DataFrame
# una columna con esos datos llamada 'Error %'
def printFrame(d):
    d['Error %'] = d['Error'] / d['Analytic'] 
    print(DataFrame(d))
    print('.'+ '-'*70 + '.')

def calcError(phiA, phiN):
    return np.absolute(phiA - phiN)
        
# Parte Principal del Script donde se invica cada uno de los metodos en los demas scripts
if __name__ == '__main__':
 
    Coefficients.alloc(5)
    m = Mesh(nodes = 5)
    d = Diffusion1D(m.volumes())
    ma = Matrix(m.volumes())
    a = Advection1D(m.volumes())

    print(m.delta(), d.aP(), a.aP(), ma.mat(), sep='\n')

    printData(nvx =5, nx = 6, longitud = 1.3)

