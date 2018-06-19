"""
Created on Tue Mar  6 15:11:05 2018

@author: Dr. Luis Miguel de la Cruz
adaptado por: Samuel Garcia

"""
# Libreria de Python
import numpy as np

# Declaración de la Clase Coeficientes
class Coefficients():
    """
    Esta clase define los arreglos principales para los coeficientes del
    metodo de Volumen Finito. Los arreglos son definidos como variables de
    clase para que sean compartidos por todos los objetos de esta clase.
    """
    # Declaración de Cada uno de las Coeficientes    
    __aP = None
    __aE = None
    __aW = None
    # Declaración de los Coeficientes para la parte de Upwind II y Quick 
    __aEE = None
    __aWW = None
    __Su = None
    # Declaración de la varieable de Numero de Volumenes
    __nvx = None
    # Declaración de la varieable de Tamaño de los Volumenes 
    __delta = None

    #### Inicialización de cada uno de las variables##########
    # Se realiza de la siguiente manera con el fin de realizar
    # Con Programación Orientada a Objetos
    def __init__(self, nvx = None, delta = None):
        Coefficients.__nvx = nvx
        Coefficients.__delta = delta


    @staticmethod
    def alloc(n):
        if Coefficients.__nvx:
            nvx = Coefficients.__nvx
        else:
            nvx = n
        Coefficients.__aP = np.zeros(nvx)
        Coefficients.__aE = np.zeros(nvx)
        Coefficients.__aW = np.zeros(nvx)
        Coefficients.__aEE = np.zeros(nvx)
        Coefficients.__aWW = np.zeros(nvx)
        Coefficients.__Su = np.zeros(nvx)
    
    def setVolumes(self, nvx):
        Coefficients.__nvx = nvx
        
    def setDelta(self, delta):
        Coefficients.__delta = delta
        
    def aP(self):
        return Coefficients.__aP

    def aE(self):
        return Coefficients.__aE
    
    def aW(self):
        return Coefficients.__aW
    
    def aEE(self):
        return Coefficients.__aEE
    
    def aWW(self):
        return Coefficients.__aWW
    
    def Su(self):
        return Coefficients.__Su
    
    ####################################################
    
    ##### Metodos de Ajuste en Fronteras###############
    ##### Frontera Tipo Dirichlet #####################
    @staticmethod
    def bcDirichlet(wall, phi):
        """Metodo que ajusta las fronteras"""
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        Su = Coefficients.__Su

        if wall == 'LEFT_WALL':
            aP[1] += aW[1]
            Su[1] += 2 * aW[1] * phi
        elif wall == 'RIGHT_WALL':
            aP[-2] += aE[-2]
            Su[-2] += 2 * aE[-2] * phi
            
    @staticmethod
    def bcDirichlet_LUD(phiA, phiB):
        """Metodo que ajusta las fronetras para tipo Upwind II"""
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        aWW = Coefficients.__aWW
        aEE = Coefficients.__aEE
        Su = Coefficients.__Su

        #if wall == 'LEFT_WALL':
        aP[1] += aW[1] + 3 * aWW[1]
        Su[1] += (2 * aW[1] + 4 * aWW[1]) * phiA
        aW[2] -= aWW[2]  # condición del segundo nodo (requerida en métodos de orden mayor a 1)
        Su[2] += (2 * aWW[2]) * phiA
        #elif wall == 'RIGHT_WALL':
        aP[-2] += aE[-2] + 3 * aEE[1]
        Su[-2] += (2 * aE[-2] + 4 * aEE[1]) * phiB
        aE[2] -= aEE[2]  # condición del penúltimo nodo (requerida en métodos de orden mayor a 1)
        Su[-3] += (2 * aEE[1]) * phiB
        
    
    @staticmethod        
    def bcDirichlet_QuicK(phiA, phiB, delta, gamma, DE, DW, rho, u):
        """Metodo que ajsuta las fronetras para tipo Quick"""
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        aWW = Coefficients.__aWW
        Su = Coefficients.__Su
        DA = gamma / delta  # Malalasekera pg. 160 half cell approximation      
        DB = gamma /delta

        """ Nodo primero """
        aE[1] = ((1/3) * DA) - ((3/8) * rho * u) + DE[0]
        Sp = -(((8/3) * DA)  +((2/8) * rho * u) + rho * u)
        aP[1] = aE[1] - Sp
        Su[1] = (-Sp) * phiA
        
        """ Segundo Nodo """
        aW[2] = ((7 / 8) * rho * u) + ((1 / 8) * rho * u) + DW[0]
        Sp2 = (1/4) * rho * u
        aP[2] =  aW[2] + aE[2] - Sp2
        Su[2] = -(Sp2) * phiA

        """ Ultimo nodo """
        aW[-2] = ((1/3)* DB) + ((6/8)* rho * u) + DW[-1]
        Spl = -(((8 / 3) * DB) - rho * u)
        aP[-2] = aW[-2] + aWW[-2] - Spl
        Su[-2] = (-Spl) * phiB
    ###################################################
    
    ##### Frontera Tipo Neumman #######################
    @staticmethod
    def bcNeumman(wall, flux):
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        Su = Coefficients.__Su
        dx = Coefficients.__delta

        if wall == 'LEFT_WALL':
            aP[1] -= aW[1]
            Su[1] -= aW[1] * flux * dx
        elif wall == 'RIGHT_WALL':
            aP[-2] -= aE[-2]
            Su[-2] += aE[-2] * flux * dx  
    #####################################################        
    
    #### Calculo de los Coeficientes Su y Sp ###########
    def setSu(self, q):
        Su = Coefficients.__Su
        dx = Coefficients.__delta
        Su += q * dx
        
    def setSp(self, Sp):
        aP = Coefficients.__aP
        dx = Coefficients.__delta
        aP -= Sp * dx
    
    # Impresión de los Coeficientes
    def printCoefficients(self):
        print('aP = {}'.format(self.__aP), 
              'aE = {}'.format(self.__aE), 
              'aP = {}'.format(self.__aW),
              'Su = {}'.format(self.__Su), sep='\n')

    # Liberación de la Memoria (Coeficientes)
    def cleanCoefficients(self):
        Coefficients.__aP[:] = 0.0
        Coefficients.__aE[:] = 0.0
        Coefficients.__aW[:] = 0.0
        Coefficients.__Su[:] = 0.0
        Coefficients.__aEE[:] = 0.0
        Coefficients.__aWW[:] = 0.0

        
# Parte Principal del Script donde se invica la clase e Impresión de Pruebas
if __name__ == '__main__':
        
    coef1 = Coefficients(6, 0.25)
    coef1.alloc(6)
    coef1.setSu(100)
    coef1.setSp(-2)
    
    print('-' * 20)  
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.Su(), sep = '\n')
    print('-' * 20)  

    ap = coef1.aP()
    ap[2] = 25
    print(ap, coef1.aP(),sep='\n')
    print('-' * 20)  

    ae = coef1.aE()
    aw = coef1.aW()
    su = coef1.Su()
    ae.fill(5)
    aw.fill(5)
    ap.fill(10)
    coef1.setSp(-2)
    coef1.bcDirichlet('LEFT_WALL', 2)
    coef1.bcNeumman('RIGHT_WALL', 1)
    coef1.printCoefficients()
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.Su(), sep = '\n')
    print('-' * 20)  

