__author__ = 'Edoardo Mancini', 'Diego Chiocci'

from matplotlib.pyplot import vlines
import numpy as np
from typing import List

class PlyDef : 
    '''
    Definition of the ply properties of a UD compostite ply starting from the fiber and material properties and the respective percentages. Different methods were developed to calculate the composite properties.
# References: 
    - Johnes, Mechanics of Composite Materials, 2nd edition
    - Barbero, Introduction to Composite Materials Design, 3rd edition 

'''
    def __init__(self, fiber_properties: list[int] or np.ndarray[int], 
                fiber_name: str, matrix_properties: list[int] or np.ndarray[int], 
                matrix_name: str, fiber_frac: float=0, 
                matrix_frac: float=0, mass_or_vol_frac: str='vol', 
                draping: str='UD', grams_per_square_meter: int=0, compute_cured_thickness: bool=False,  
                ) :
        '''
# DESCRIPTION:
Initialization of the class "ply" creating a composite layer from its constituents 

# INPUTS:

    Required
    -   fiber_properties            : [GPa, adimensional, kg/m3] - list or numpy array of fiber mechanical properties in this order: E, ni, rho.
    -   fiber_name                  : name of the fiber used
    -   matrix_properties           : [GPa, adimensional, kg/m3] - list or numpy array of matrix mechanical properties in this order: E, ni, rho.
    -   matrix_name                 : name of the fiber used

    Optional 
    
    -   fiber_frac                  : fiber volume or mass fraction. 
    -   matrix_frac                 : matrix volume or mass fraction. It is required to provide either this quantity or the former. Also, both can be provided but in this case the values must agree
    -   mass_or_vol_frac            : whether the fraction is of mass or volume. Default "vol"
    -   draping                     : type of draping. Default "UD"
    -   grams_per_square_meter      : grams per square meter of the fiber. Default 0
    -   compute_cured_thickness     : whether to compute or not the cured ply thickness. Default False.

# OUTPUTS: 
    -   Class object with all the mechanical properties (E_i, n_ij, G_ij with i,j = 1, 2, 3 with i != j), and the other properties. In such a way that "output_name".properties returns the specific ply property

# Example:
    import PyComp_local as comp
    fiber_props = [513, .3, 1910]
    fiber_name = 'M55J/Toray'
    fiber_mass_frac = 1 - .31

    matrix_props = [513, .3, 1910]
    matrix_name = 'M55J/Toray'
    gpsqm = 160

    ExamplePly1 = comp.PlyDef(fiber_props, fiber_name, matrix_props, matrix_name, fiber_frac = fiber_mass_frac, mass_or_vol_frac='wgt')

    ''' 
        # Checks on the inputs types
        if isinstance(fiber_properties, list) is False and isinstance(fiber_properties, np.ndarray) is False:
            raise Exception('"fiber_properties" must either be a list or a np.array')
        if isinstance(fiber_name, str) is False:
            raise Exception('"fiber_name" must be a string')
        if isinstance(matrix_properties, list) is False and isinstance(matrix_properties, np.ndarray) is False:
            raise Exception('"matrix_properties" must either be a list or a np.array')
        if isinstance(matrix_name, str) is False:
            raise Exception('"matrix_name" must be a string')
        if isinstance(fiber_frac, float) is False and fiber_frac != 0 and fiber_frac != 1:
            raise Exception('"fiber_frac"  must be either float, 0 or 1.')
        if isinstance(matrix_frac, float) is False and matrix_frac != 0 and matrix_frac != 1:
            raise Exception('"matrix_frac" must be either float, 0 or 1.')
        if isinstance(mass_or_vol_frac, str) is False:
            raise Exception('"mass_or_vol_frac" must be a string')
        if isinstance(draping, str) is False:
            raise Exception('"draping" must either be a string')
        if isinstance(grams_per_square_meter, int) is False and isinstance(grams_per_square_meter, float) is False:
            raise Exception('"grams_per_square_meter" must be a float or an integer')
        if grams_per_square_meter < 0: 
            raise Exception('grams_per_square_meter must be positive')
        if isinstance(compute_cured_thickness, bool) is False:
            raise Exception('"compute_cured_thickness" must be a boolean')

        # Additional checks
        if fiber_frac == 0 and matrix_frac == 0 :
            raise Exception('Either fiber_vol_frac or matrix_vol_frac must be different from 0')
        if fiber_frac < 0:
            raise Exception('Error, fiber_frac must be positive.')
        if matrix_frac <0:
            raise Exception('Error, matrix_frac must be positive.')
        if fiber_frac > 1:
            raise Exception('Error, fiber_frac must be smaller or equal than 1.')
        if matrix_frac > 1:
            raise Exception('Error, matrix_frac must be smaller or equal than 1.')
        if fiber_frac != 0 and matrix_frac != 0 : 
            if 1 - fiber_frac != matrix_frac : 
                raise Exception('You provided both the fiber and the matrix volume fractions but their sum is different from 1. Please provide only one of the two or input coherent values')
        if mass_or_vol_frac != 'vol' and mass_or_vol_frac != 'wgt' :
            raise Exception('The value of the variable "mass_or_vol_frac" is not acceptable. Accepted values are: "vol" or "wgt".')
        if draping != 'UD' : 
            raise Exception('The value of the variable "draping" is not acceptable. Accepted values are "UD".')

        self.name = f'{fiber_name} - {matrix_name}'
        self.fiber_E, self.fiber_ni, self.fiber_rho = fiber_properties
        self.fiber_G = self.fiber_E / (2 * (1 + self.fiber_ni))
        
        self.matrix_E, self.matrix_ni, self.matrix_rho = matrix_properties
        self.matrix_G = self.matrix_E / (2 * (1 + self.matrix_ni))

        if mass_or_vol_frac == 'vol' :
            if fiber_frac == 0 :
                self.matrix_vol_frac = matrix_frac
                self.fiber_vol_frac = 1 - self.matrix_vol_frac
            else: 
                self.fiber_vol_frac = fiber_frac
                self.matrix_vol_frac = 1 - self.fiber_vol_frac
            
            self.rho = self.fiber_rho * self.fiber_vol_frac + self.matrix_rho * self.matrix_vol_frac

            self.fiber_wgt_frac = self.fiber_vol_frac * self.fiber_rho / self.rho 
            self.matrix_wgt_frac = self.matrix_vol_frac * self.matrix_rho / self.rho

        else : 
            if fiber_frac == 0 : 
                self.matrix_wgt_frac = matrix_frac
                self.fiber_wgt_frac = 1 - self.matrix_wgt_frac
            else: 
                self.fiber_wgt_frac = fiber_frac
                self.matrix_wgt_frac = 1 - self.fiber_wgt_frac    

            self.rho = self.fiber_rho * self.matrix_rho / (self.fiber_wgt_frac * self.matrix_rho + self.matrix_wgt_frac * self.fiber_rho)
            
            self.fiber_vol_frac = self.fiber_wgt_frac * self.rho / self.fiber_rho 
            self.matrix_vol_frac = self.matrix_wgt_frac * self.rho / self.matrix_rho        
        
        if grams_per_square_meter == 0:
            self.grm2 = 'NA'
            self.cured_thickness = 'NA'
        if grams_per_square_meter != 0 : 
           self.grm2 = grams_per_square_meter
           if compute_cured_thickness is True :
               self.cured_thickness = grams_per_square_meter / (1000 * (self.fiber_rho / 1000) * self.fiber_vol_frac)
           else :
               self.cured_thickness = 'NA'

    def ROM(self, print_cntrl:bool=False) :
        '''
# DESCRIPTION:
The class "ply" is used to generate a ply starting from its properties 

# INPUTS:

    Required
    - None
    Optional   
    - print_cntrl: wheather to print or not the computed ply properties
# OUTPUTS: 
    - None
# Example:
    import PyComp_local as comp
    fiber_props = [513, .3, 1910]
    fiber_name = 'M55J/Toray'
    fiber_mass_frac = 1 - .31

    matrix_props = [513, .3, 1910]
    matrix_name = 'M55J/Toray'
    gpsqm = 160

    ExamplePly1 = comp.PlyDef(fiber_props, fiber_name, matrix_props, matrix_name, fiber_frac = fiber_mass_frac, mass_or_vol_frac='wgt')

    ''' 
        if isinstance(print_cntrl, bool) is False:
            raise Exception('Error. The variable "print_cntrl" must be boolean.')
        
        self.E1 = self.fiber_E * self.fiber_vol_frac + self.matrix_E * self.matrix_vol_frac
        self.E2 = self.fiber_E * self.matrix_E / (self.matrix_vol_frac * self.fiber_E + self.fiber_vol_frac * self.matrix_E)
        self.E3 = self.E2

        self.G12 = self.fiber_G * self.matrix_G / (self.matrix_vol_frac * self.fiber_G + self.fiber_vol_frac * self.matrix_G)
        self.G13 = self.G12
        self.G23 = self.matrix_G

        self.ni12 = self.fiber_ni * self.fiber_vol_frac + self.matrix_ni * self.matrix_vol_frac            
        self.ni13 = self.ni12
        self.ni23 = self.E2 / (self.G23 * 2) - 1
        self.ni21 = self.ni12 * self.E2 / self.E1
        self.ni31 = self.ni21
        self.ni32 = self.ni23
        print('\033[33m','WARNING:The value "ni23" and "G23" were set to the matrix value. I a more precise value is needed use another method.')
        print('\033[37m',' ')

        self.mech_props = [self.name, self.E1, self.E2, self.E3, self.ni12, self.ni13, self.ni23, self.G12, self.G13, self.G23, self.rho, self.cured_thickness]

        if print_cntrl is True: 
            self.print_properties()
    def SSP_calc(self, csi_E: float=2) :

            # Halphin_Tsai method used, SSP too complicated 
            self.E1 = self.fiber_E * self.fiber_vol_frac + self.matrix_E * self.matrix_vol_frac
            self.ni12 = self.fiber_ni * self.fiber_vol_frac + self.matrix_ni * self.matrix_vol_frac            

            eta_E = (self.fiber_E / self.matrix_E - 1) / (self.fiber_E / self.matrix_E + csi_E)
            self.E2 = self.matrix_E * (1 + csi_E * eta_E * self.fiber_vol_frac) / (1 - eta_E * self.fiber_vol_frac)
            self.E3 = self.E2

            eta_ssp_23 = (3 - self.matrix_ni + self.matrix_G / self.fiber_G) / (4 * (1 - self.matrix_ni))
            self.G23 = self.matrix_G * (self.fiber_vol_frac + eta_ssp_23 * (1 - self.fiber_vol_frac))   
            eta_s = .5 * (1 + self.matrix_G / self.matrix_E)
            self.G12 = (self.fiber_vol_frac + eta_s * (1 - self.fiber_vol_frac)) / (eta_s * (1 - self.fiber_vol_frac) + self.fiber_vol_frac * self.matrix_G / self.fiber_G)
            self.G13 = self.G12

            self.ni13 = self.ni12
            self.ni23 = self.E2 / (self.G23 * 2) - 1
            self.ni21 = self.ni12 * self.E2 / self.E1
            self.ni31 = self.ni21
            self.ni32 = self.ni23

            # E1 = self.E1
            # E2 = self.E2
            # E3 = self.E3
            # ni12 = self.ni12
            # ni13 = self.ni13
            # ni21 = self.ni21
            # ni23 = self.ni23
            # ni32 = self.ni32
            # ni31 = self.ni31
            # G12 = self.G12
            # G23 = self.G23
            # G13 = self.G13

            # return(E1, E2, E3, ni12, ni13, ni21, ni23, ni31, ni32, G12, G23, G13)
            self.mech_props = [self.name, self.E1, self.E2, self.E3, self.ni12, self.ni13, self.ni23, self.G12, self.G13, self.G23, self.rho, self.cured_thickness]

    def Halphin_Tsai_calc(self, csi_E: float=2, csi_G: str='def') :
           
            if isinstance(csi_E, float) is False and isinstance(csi_E, int) is False:
                raise ValueError('"csi_E" must either be a float or an integer')
            if csi_E < 0 : 
                raise ValueError('"csi_E" must be > 0.')
            if isinstance(csi_G, float) is False and isinstance(csi_G, int) is False and csi_G != 'def':
                raise ValueError('"csi_G" must either be a float or an integer')
            if csi_G == 'def' :
                pass 
            elif csi_G < 0 : 
                raise ValueError('"csi_G" must be > 0.')
            # if isinstance(csi_ni23, float) is False and isinstance(csi_ni23, int) is False:
            #     raise ValueError('"csi_ni23" must either be a float or an integer')
            # if csi_ni23 < 0 : 
            #     raise ValueError('"csi_ni23" must be > 0.')

            # Computed with ROM 
            # if csi_E == 'def' :
            #     csi_E = 1 + 40 * self.fiber_vol_frac ** 40
            # elif csi_E == 1 :
            #     pass
            self.E1 = self.fiber_E * self.fiber_vol_frac + self.matrix_E * self.matrix_vol_frac
            self.ni12 = self.fiber_ni * self.fiber_vol_frac + self.matrix_ni * self.matrix_vol_frac
            

            eta_E = (self.fiber_E / self.matrix_E - 1) / (self.fiber_E / self.matrix_E + csi_E)
            self.E2 = self.matrix_E * (1 + csi_E * eta_E * self.fiber_vol_frac) / (1 - eta_E * self.fiber_vol_frac)

            if csi_G == 'def' :
                csi_G = 1 + 40 * self.fiber_vol_frac ** 10
            elif csi_G == 1 :
                pass
            
            eta_G = (self.fiber_G / self.matrix_G - 1) / (self.fiber_G / self.matrix_G + csi_G)
            self.G12 = self.matrix_G * (1 + csi_G * eta_G * self.fiber_vol_frac) / (1 - eta_E * self.fiber_vol_frac)
            
            #SSP for G23
            eta_23 = (3 - 4 * self.matrix_ni + self.matrix_G / self.fiber_G) / (4 * (1 - self.matrix_ni))
            self.G23 = self.matrix_G * (self.fiber_vol_frac + eta_23*(1 - self.fiber_vol_frac)) / (eta_23 * (1 - self.fiber_vol_frac) + (self.fiber_vol_frac * self.matrix_G / self.fiber_G)) 
            
            #calcolo della ni23
            # eta_ni23=(self.fiber_ni/self.matrix_ni - 1)/(self.fiber_ni/self.matrix_ni +csi_ni23)
            # self.ni23=self.matrix_ni * (1 + csi_ni23 * eta_ni23 * self.fiber_vol_frac)/(1- eta_ni23 * self.fiber_vol_frac)


            self.G13 = self.G12
            self.E3 = self.E2
            self.ni13 = self.ni12
            # self.ni23 = self.E2 / (self.G23 * 2) - 1 
            self.ni23= self.ni12
            self.ni21 = self.ni12 * self.E2 / self.E1
            self.ni31 = self.ni21
            self.ni32 = self.ni23
            
            self.mech_props = [self.name, self.E1, self.E2, self.E3, self.ni12, self.ni13, self.ni23, self.G12, self.G13, self.G23, self.rho, self.cured_thickness,csi_G]
            
    def PMM(self) : 
        
            mu_0 = self.matrix_G
            mu_1 = self.fiber_G
            ni_0 = self.matrix_ni
            ni_1 = self.fiber_ni
            f = self.fiber_vol_frac

            # lambda_0 = (self.matrix_ni * self.matrix_E) / (1 + self.matrix_ni) * (1 - self.matrix_ni * 2)
            lambda_0 = (self.matrix_E * self.matrix_ni) / ((1 + self.matrix_ni) * (1 - self.matrix_ni * 2))
            # lambda_1 = self.fiber_ni * self.fiber_E / (1 + self.fiber_ni) * (1 - self.fiber_ni * 2)
            
            a = mu_1 - mu_0 - 2 * mu_1 * ni_0 + 2 * mu_0 * ni_1
            b = - mu_0 * ni_0 + mu_1 * ni_1 + 2 * mu_0 * ni_0 * ni_1 - 2 * mu_1 * ni_0 * ni_1 
            c = (mu_0 - mu_1) * (mu_1 - mu_0 + mu_1 * ni_1 - mu_0 * ni_0 + 2 * mu_0 * ni_1 - 2 * mu_1 * ni_0 + 2 * mu_0 * ni_0 * ni_1 - 2 * mu_1 * ni_0 * ni_1)
            g = (2 - 2 * ni_0)
            
            S3 = 0.49247 - (0.47603 * f) - 0.02748 * (f ** 2)
            S6 = 0.36844 - (0.14944 * f) - 0.27152 * (f ** 2)
            S7 = 0.12346 - (0.32035 * f) + 0.23517 * (f ** 2)

            D = (a * (S3 ** 2)) / (2 * c * mu_0**2) - (a * S6 * S3 )/ (g * c * (mu_0**2)) + (a  * ((S6 ** 2) - (S7 ** 2)))/ (2 * (mu_0 ** 2) * (g ** 2) * c) + (S3 * ((b ** 2) - (a ** 2))) / (2 * mu_0 * (c ** 2) ) + (S6 *( (a ** 2) - (b ** 2)) + S7 * ( a * b + (b**2))) / (2* mu_0 * g * (c ** 2)) + ((a ** 3) - 2 * (b ** 3) - 3 * a * (b ** 2)) / (8 * (c ** 3)) 
            C11 = lambda_0 + 2 * mu_0 - f * ((S3 ** 2) / (mu_0 ** 2) - (2 * S6 * S3) / ((mu_0 ** 2) * g) - (a * S3) / (mu_0 * c) + ((S6 ** 2) - (S7 ** 2)) / ((mu_0 ** 2) * (g ** 2)) + ( a* S6 + b * S7) / (mu_0 * g * c) + ((a ** 2) - (b ** 2)) / (4 * (c **2))) / D
            C12 = lambda_0 + f * b * (S3 / (2 * c * mu_0) - (S6 - S7) / (2 * c * mu_0 * g) - (a + b) / (4 * (c ** 2))) / D
            C23 = lambda_0 + f * ( (a * S7) / (2 * mu_0 * g * c) - (b * a + (b ** 2)) / (4 * (c ** 2))) / D
            C22 = lambda_0 + 2 * mu_0 - f * (- (a * S3) / (2 * mu_0 * c) + (a * S6) / (2 * mu_0 * g * c) + ((a ** 2) - (b ** 2)) / (4 * (c ** 2))) / D
            C44 = mu_0 - f * (- (2 * S3) / mu_0 + ((mu_0 - mu_1) ** - 1) + 4 * S7 / (mu_0 * (2 - 2 * ni_0)))**-1
            C66 = mu_0 - f * (- S3 / mu_0 + ((mu_0 - mu_1) ** -1)) ** -1

            C55 = C66
            C13 = C12 
            C31 = C13 
            C21 = C12
            C32 = C23
            C33 = C22 
            C = np.array([[C11, C12, C13, 0, 0, 0], [C21, C22, C23, 0, 0, 0], [C31, C32, C33, 0, 0, 0], [0, 0, 0, C44, 0, 0], [0, 0, 0, 0, C55, 0 ], [0, 0, 0, 0, 0, C66]]) 
            self.C = C
            # S = np.reverse(C)
            # S=np.linalg.inv(C)

            # self.E1 = 1 / S[0, 0]
            # self.E2 = 1 / S[1, 1]
            # self.E3 = 1 / S[2, 2]
            # self.ni12 = - S[1, 0] * self.E1
            # self.ni23 = - S[2, 1] * self.E2 
            # self.G12 = 1 / S[3, 3]
            # self.G13 = 1 / S[4, 4]
            # self.G23 = 1 / S[5, 5]


            # self.E1 = C[0,0]- (2 * C[0,1]**2)/(C[1,1]+C[1,2])
            # self.E2 = ((2 * C[0,0] * C[1,1] + 2 * C[0,0] * C[1,2] - 4 * C[0,1]**2)*(C[1,1] - C[1,2] + 2 * C[3,3]))/(3 * C[0,0] * C[1,1] + C[0,0] * C[1,2] + 2 * C[0,0] * C[3,3] - 4 * C[0,1]**2)
            # self.E3 = self.E2
            # self.ni12 = C[0,1]/(C[1,1] + C[1,2])
            # self.ni13 = self.ni12
            # self.ni23 = (C[0,0] * C[1,1] + 3 * C[0,0] * C[1,2] - 2 * C[0,0] * C[3,3] - 4 * C[0,1]**2)/(3 * C[0,0] * C[1,1] + C[0,0] * C[1,2] + 2 * C[0,0] * C[3,3] - 4 * C[0,1]**2)
            # self.G12 = C[5,5]
            # self.G13 = self.G12
            # self.G23 = C[1,1]/4 - C[1,2]/4 + C[3,3]/2

            self.E1 = C11 - ((2 * (C12**2))/(C22 + C23))
            self.E2 = ((2 * C11 * C22 + 2 * C11 * C23 - 4 * (C12**2)) * (C22 - C23 + 2 * C44))/(3 * C11 * C22 + C11 * C23 + 2 * C11 * C44 - 4 * (C12**2))
            self.E3 = self.E2
            self.ni12 = (C12)/(C22 + C23)
            self.ni13 = self.ni12
            self.ni23 = (C11 * C22 + 3 * C11 * C23 - 2 * C11 * C44 - 4 * (C12**2))/(3 * C11 * C22 + C11 * C23 + 2 * C11 * C44 - 4 * (C12**2)) 
            self.G12 = C66
            self.G13 = self.G12
            self.G23 = (C22/4) - (C23/4) + (C44/2)

            # self.ni13 = self.ni12
            # # self.ni23 = self.E2 / (self.G23 * 2) - 1
            # self.ni21 = self.ni12 * self.E2 / self.E1
            # self.ni31 = self.ni21
            # self.ni32 = self.ni23
            self.mech_props = [self.name, self.E1, self.E2, self.E3, self.ni12, self.ni13, self.ni23, self.G12, self.G13, self.G23, self.rho, self.cured_thickness]

    def print_properties(self) :
        print('E1 = ', str(self.E1), ' GPa')
        print('E2 = ', str(self.E2), ' GPa')
        print('E3 = ', str(self.E3), ' GPa')
        print('ni12 = ', str(self.ni12))
        print('ni13 = ', str(self.ni13))
        print('ni23 = ', str(self.ni23))
        print('G12 = ', str(self.G12), ' GPa')
        print('G23 = ', str(self.G23), ' GPa')
        print('G13 = ', str(self.G13), ' GPa')
        print('rho = ', str(self.rho), ' kg/m^3')
        print('Ply thickness = ', str(self.cured_thickness), ' mm')

    def error_percent(self, data: list[float] or np.ndarray[float]):
        self.errorE1 = np.abs((data[0] - self.E1*(10**3)))/data[0]
        self.errorE2 = np.abs((data[1] - self.E2*(10**3)))/data[1]
        self.errorE3 = np.abs((data[2] - self.E3*(10**3)))/data[2]
        self.errorG12 = np.abs((data[3] - self.G12*(10**3)))/data[3]
        self.errorG23 = np.abs((data[4] - self.G23*(10**3)))/data[4]
        self.errorG13 = np.abs((data[5] - self.G13*(10**3)))/data[5]
        self.errorni12 = np.abs((data[6] - self.ni12))/data[6]
        self.errorni13 = np.abs((data[7] - self.ni13))/data[7]
        self.errorni23 = np.abs((data[8] - self.ni23))/data[8]
        self.errorrho = np.abs((data[9] - self.rho*((1/1000)*(1/10**9))))/data[9]

        self.ERRORI = [self.errorE1, self.errorE2, self.errorE3, self.errorG12, self.errorG23, self.errorG13, self.errorni12, self.errorni13, self.errorni23, self.errorrho]

        