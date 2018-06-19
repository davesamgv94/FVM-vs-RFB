"""
Created on Tue Mar  6 15:11:05 2018

@author: Dr. Luis Miguel de la Cruz
adaptado por: Samuel Garcia

"""
# Libreria de Python
import numpy as np
# Importar Script donde se encuentra la parte
# de los Coeficienetes
from Coefficients import Coefficients

# Declaración de la Clase Advección
class Advection1D(Coefficients):
    
    """
    Esta clase realiza el calculo en los arreglos para la parte Advectiva 
    con los arreglos principales de los coeficientes del metodo de Volumen Finito. .
    """  
    # Declaración de variables
    def __init__(self, nvx = None, rho = None, dx = None):
        super().__init__(nvx)
        self.__nvx = nvx
        self.__rho = rho
        self.__dx = dx
        self.__u = np.zeros(nvx-1)

    def __del__(self):
        del(self.__nvx)
        del(self.__rho)
        del(self.__dx)
        del(self.__u)

    def setU(self, u):
        if type(u) == float:
            self.__u.fill(u)
        else:
            self.__u = u

    def u(self):
        return self.__u
    
    def calcCoef(self, typeAp = '', phiA = 0, phiB = 0):
        aP = self.aP()
        aE = self.aE()
        aW = self.aW()
        aEE = self.aEE()
        aWW = self.aWW()                
        u = self.__u
        rho = self.__rho
        
        # Calculos de Coeficientes con Diferentes Metodos
        
        ##########################################################################
		 # Diferencias Centrales
        if typeAp == 'Centradas':
            for i in range(1,self.__nvx-1):
                #
                CE = - rho * u[i] * 0.5
                CW =   rho * u[i-1] * 0.5             
                aE[i] += CE 
                aW[i] += CW
                aP[i] += CE + CW + rho * (u[i] - u[i-1])
		 ##########################################################################
		 # Upwind
        elif typeAp == 'Upwind':
            for i in range(1,self.__nvx-1):
                CE = max((-u[i],0)) 
                CW = max((u[i-1],0))
                aE[i] += CE 
                aW[i] += CW
                aP[i] += CE + CW + rho * (u[i] - u[i-1])
        ###########################################################################
        #UpwindII
        elif typeAp == 'UpwindII' :
            for i in range(1,self.__nvx-1):

                CE = max((rho*u[i]*0.5,0))
                CW = max((rho*u[i-1]*0.5,0))
                CEe = -max((-rho*u[i]*0.5,0))
                CWw = -max((-rho*u[i-1]*0.5,0))
                aE[i] += -3*CEe-CWw
                aW[i] += CE+3*CW
                aEE[i] += CEe
                aWW[i] += -CW
                aP[i] += 3*CE-3*CWw + rho * (u[i] - u[i-1])
                
        ###########################################################################
        # Quick
        elif typeAp == 'Quick':   ### pagina 177 pdf
            for i in range(2,self.__nvx-2):
                """
                CE = Ce
                CW = cw
                CEe = Ce*
                CWw = Cw*
                """
                CE = max((rho*u[i]*0.125,0))     # 1/8 = 0.125
                CW = max((rho*u[i-1]*0.125,0))
                CEe = -max((-rho*u[i]*0.125,0))
                CWw = -max((-rho*u[i-1]*0.125,0))
                aE[i] += -3*CE - 6*CEe - CWw
                aW[i] += 6*CW + CE + 3*CWw
                aEE[i] += CEe
                aWW[i] += -CW
                aP[i] += 6*CE - 3*CW + 3*CEe - 6*CWw + rho * (u[i] - u[i-1])
                               
                
# Parte Principal del Script donde se invica la clase e Impresión de Pruebas          
if __name__ == '__main__':
    
    nx = 5
    u = np.sin(np.linspace(0,1,nx))
#    u = np.ones(nx)
    print('-' * 20)  
    print(u)
    print('-' * 20)  

    af1 = Advection1D(6, 1, 1)
    af1.alloc(6)
    af1.setU(u)
    print(af1.u())
    print('-' * 20)  

    af1.calcCoef('Quick´')
    print(af1.aP(), af1.aE(), af1.aW(), af1.Su(), af1.aWW(), af1.aEE(), sep = '\n')
    print('-' * 20)  

#    af1.bcDirichlet('LEFT_WALL', 2)
#    af1.bcDirichlet('RIGHT_WALL', 1)
#    print(af1.aP(), af1.aE(), af1.aW(), af1.Su(), af1.aWW(), af1.aEE(), sep = '\n')
#    print('-' * 20)  



