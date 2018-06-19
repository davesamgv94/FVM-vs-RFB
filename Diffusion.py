"""
Created on Tue Mar  6 15:11:05 2018

@author: Dr. Luis Miguel de la Cruz
adaptado por: Samuel Garcia

"""
# Importar Script donde se encuentra la parte
# de los Coeficienetes
from Coefficients import Coefficients

# Declaración de la Clase Difusion

class Diffusion1D(Coefficients):
    
    """
    Esta clase realiza el calculo en los arreglos para la parte Difusiva con una Gamma dada
    con los arreglos principales de los coeficientes del metodo de Volumen Finito.
    """ 
    ##### Declaaración de las Varibles de Numero y Tamaño de Volumenes
    # en conjunto con la variable Gamma para la parte Difusiva
    ### Inicialización####
    def __init__(self, nvx = None, Gamma = None, dx = None):
        super().__init__(nvx, dx)
        self.__nvx = nvx
        self.__Gamma = Gamma
        self.__dx = dx

    def __del__(self):
        del(self.__Gamma)
        del(self.__dx)
    
    # Calculo de los Coeficientes con la Gamma####
    def calcCoef(self):
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
                
        aE += self.__Gamma / self.__dx
        aW += self.__Gamma / self.__dx
        aP += aE + aW
        
        return aE, aW, aP
 
#        for i in range(self.__nvx):
#            aE[i] += self.__Gamma / self.__dx
#            aW[i] += self.__Gamma / self.__dx
#            aP[i] += aE[i] + aW[i]


# Parte Principal del Script donde se invica la clase e Impresión de Pruebas
if __name__ == '__main__':
    
    df1 = Diffusion1D(5, 5, 1)
    df1.alloc(5)
    df1.calcCoef()
    df1.setSu(100)

    print('-' * 20)  
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  
