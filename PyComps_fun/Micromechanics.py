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
    import PyComps as comp
    fiber_props = [513, .3, 1910]
    fiber_name = 'M55J/Toray'
    fiber_mass_frac = 1 - .31

    matrix_props = [5, .3, 1156]
    matrix_name = 'Ex-1515'
    gpsqm = 160

    example_ply_1 = comp.PlyDef(fiber_props, fiber_name, matrix_props, matrix_name, fiber_frac = fiber_mass_frac, mass_or_vol_frac='wgt')

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
ROM method computes ply equivalent properties following the Rule-Of-Mixtures approach

# INPUTS:

    Required
    - None
    Optional   
    - print_cntrl: wheather to print or not the computed ply properties
# OUTPUTS: 
    - None
# Example:
    import PyComps as comp
    fiber_props = [513, .3, 1910]
    fiber_name = 'M55J/Toray'
    fiber_mass_frac = 1 - .31

    matrix_props = [5, .3, 1156]
    matrix_name = 'Ex-1515'
    gpsqm = 160

    example_ply_2 = comp.PlyDef(fiber_props, fiber_name, matrix_props, matrix_name, fiber_frac=fiber_mass_frac, grams_per_square_meter=gpsqm, compute_cured_thickness=True, mass_or_vol_frac='wgt')
    example_ply_2.ROM(print_cntrl=True)
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
        self.ni23 = 'NA'
        self.ni21 = self.ni12 * self.E2 / self.E1
        self.ni31 = self.ni21
        self.ni32 = self.ni23
        print('\033[35m','Note: "ni23" and "G23" are not computed but set to the matrix value. \n For precise values please use another method.')
        print('\033[37m',' ')

        self.mech_props = [self.name, self.E1, self.E2, self.E3, self.ni12, self.ni23, \
                           self.ni13, self.G12, self.G23, self.G13, self.rho, self.cured_thickness]

        if print_cntrl is True: 
            self.print_properties()
 
    def Halphin_Tsai(self, csi_E: float=2, csi_G: str='def', print_cntrl:bool=False) :
        '''
# DESCRIPTION:
Halphin_Tsai method computes ply equivalent properties following the Halphin-Tsai approach

# INPUTS:
    Required
    - None
    Optional   
    - csi_E: value of csi used for the computation of the elastic modulus. Default is 2.
    - csi_G: value of csi used for the computation of the shear modulus. Default is 'def' which means: (1 + 40 * fiber_vol_frac ** 10).
    - print_cntrl: wheather to print or not the computed ply properties

# OUTPUTS: 
    - None

# Example:
    import PyComps as comp
    fiber_props = [513, .3, 1910]
    fiber_name = 'M55J/Toray'
    fiber_mass_frac = 1 - .31

    matrix_props = [5, .3, 1156]
    matrix_name = 'Ex-1515'
    gpsqm = 160

    example_ply_1 = comp.PlyDef(fiber_props, fiber_name, matrix_props, matrix_name, fiber_frac = fiber_mass_frac, mass_or_vol_frac='wgt', grams_per_square_meter=gpsqm, compute_cured_thickness=True)
    example_ply_1.Halphin_Tsai(print_cntrl=True)
    ''' 
        if isinstance(print_cntrl, bool) is False:
            raise Exception('Error. The variable "print_cntrl" must be boolean.')
                   
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
        self.G23 = self.matrix_G * (self.fiber_vol_frac + eta_23*(1 - self.fiber_vol_frac)) \
              / (eta_23 * (1 - self.fiber_vol_frac) + (self.fiber_vol_frac * self.matrix_G / self.fiber_G)) 
        
        self.G13 = self.G12
        self.E3 = self.E2
        self.ni13 = self.ni12
        self.ni23= self.ni12
        self.ni21 = self.ni12 * self.E2 / self.E1
        self.ni31 = self.ni21
        self.ni32 = self.ni23
        
        self.mech_props = [self.name, self.E1, self.E2, self.E3, self.ni12, self.ni23, \
                           self.ni13, self.G12, self.G23, self.G13, self.rho, self.cured_thickness]
        
        print('\033[35m','Note: "ni23" is not computed but set equal to ni12. \n For a more precise value please use another method.')
        print('\033[37m',' ')

        if print_cntrl is True: 
            self.print_properties()

    def PMM(self, print_cntrl:bool=False) : 
        '''
# DESCRIPTION:
Following the Periodic Microstructure Model method, the method computes ply equivalent properties.

# INPUTS:
    Required
    - None
    Optional   
    - print_cntrl: wheather to print or not the computed ply properties

# OUTPUTS: 
    - None

# Example:
    import PyComps as comp
    fiber_props = [513, .3, 1910]
    fiber_name = 'M55J/Toray'
    fiber_mass_frac = 1 - .31

    matrix_props = [5, .3, 1156]
    matrix_name = 'Ex-1515'
    gpsqm = 160

    example_ply_1 = comp.PlyDef(fiber_props, fiber_name, matrix_props, matrix_name, fiber_frac = fiber_mass_frac, mass_or_vol_frac='wgt', grams_per_square_meter=gpsqm, compute_cured_thickness=True)
    example_ply_1.PMM(print_cntrl=True)
    '''

        if isinstance(print_cntrl, bool) is False:
            raise Exception('Error. The variable "print_cntrl" must be boolean.')     
           
        mu_0 = self.matrix_G
        mu_1 = self.fiber_G
        ni_0 = self.matrix_ni
        ni_1 = self.fiber_ni
        f = self.fiber_vol_frac

        lambda_0 = (self.matrix_E * self.matrix_ni) / ((1 + self.matrix_ni) * (1 - self.matrix_ni * 2))
        
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

        self.E1 = C11 - ((2 * (C12**2))/(C22 + C23))
        self.E2 = ((2 * C11 * C22 + 2 * C11 * C23 - 4 * (C12**2)) * (C22 - C23 + 2 * C44))/(3 * C11 * C22 + C11 * C23 + 2 * C11 * C44 - 4 * (C12**2))
        self.E3 = self.E2
        self.ni12 = (C12)/(C22 + C23)
        self.ni13 = self.ni12
        self.ni23 = (C11 * C22 + 3 * C11 * C23 - 2 * C11 * C44 - 4 * (C12**2))/(3 * C11 * C22 + C11 * C23 + 2 * C11 * C44 - 4 * (C12**2)) 
        self.G12 = C66
        self.G13 = self.G12
        self.G23 = (C22/4) - (C23/4) + (C44/2)


        self.mech_props = [self.name, self.E1, self.E2, self.E3, self.ni12, self.ni23, \
                           self.ni13, self.G12, self.G23, self.G13, self.rho, self.cured_thickness]
        
        if print_cntrl is True: 
            self.print_properties()

    def print_properties(self):
        '''
# DESCRIPTION:
Print ply properties 

# INPUTS:
    Required
    - None
    Optional   
    
# OUTPUTS: 
    - None

# Example:
    import PyComps as comp
    fiber_props = [513, .3, 1910]
    fiber_name = 'M55J/Toray'
    fiber_mass_frac = 1 - .31

    matrix_props = [5, .3, 1156]
    matrix_name = 'Ex-1515'
    gpsqm = 160

    example_ply_1 = comp.PlyDef(fiber_props, fiber_name, matrix_props, matrix_name, fiber_frac = fiber_mass_frac, mass_or_vol_frac='wgt', grams_per_square_meter=gpsqm, compute_cured_thickness=True)
    example_ply_1.PMM()
    example_ply_1.print_properties()

    '''
        check_ply_properties = True
        try:
            print('E1 = ', str(np.round(self.E1, 4)), ' GPa')
            print('E2 = ', str(np.round(self.E2, 4)), ' GPa')
            print('E3 = ', str(np.round(self.E3, 4)), ' GPa')
            print('ni12 = ', str(np.round(self.ni12, 4)))
            print('ni13 = ', str(np.round(self.ni13, 4)))
            if self.ni23 != 'NA':
                print('ni23 = ', str(np.round(self.ni23, 4)))
            else: 
                print('ni23 = NA')
            print('G12 = ', str(np.round(self.G12, 4)), ' GPa')
            print('G23 = ', str(np.round(self.G23, 4)), ' GPa')
            print('G13 = ', str(np.round(self.G13, 4)), ' GPa')
        except AttributeError:
            check_ply_properties = False
        if check_ply_properties is False:    
            raise Exception('Error. E1 is not defined. Use one of the available method to compute the ply properties before printing.')
        print('rho = ', str(np.round(self.rho, 4)), ' kg/m^3')
        print('Ply thickness = ', str(np.round(self.cured_thickness, 4)), ' mm')

    def error_percent(self, data: list[float or int] or np.ndarray[float or int], print_cntrl:bool= False):
        ''''
# DESCRIPTION:
Compute the error percent between the computed ply properties and those from an external model
# INPUTS:
    Required
    - data_in: vector of material data for comparison. Expected format is: [E1, E2, E3, ni12, ni23, ni13, G12, G23, G13, rho]
    Optional   
    - print_cntrl, if "True" prints the % error. Default is false.
# OUTPUTS: 
    - None

# Example:
    import PyComps as comp
    fiber_props = [513, .3, 1910]
    fiber_name = 'M55J/Toray'
    fiber_mass_frac = 1 - .31

    matrix_props = [5, .3, 1156]
    matrix_name = 'Ex-1515'
    gpsqm = 160

    example_ply_1 = comp.PlyDef(fiber_props, fiber_name, matrix_props, matrix_name, fiber_frac = fiber_mass_frac, mass_or_vol_frac='wgt', grams_per_square_meter=gpsqm, compute_cured_thickness=True)
    example_ply_1.PMM()
    data_in = [295, 17, 17, .25, .25, .27, 6.8, 6.2, 6.8, 1575]
    example_ply_1.error_percent(data_in, print_cntrl=True)

    '''

        print('\033[35m', 'Note: Elastic properties units are in GPa.')
        print('Density units are in kg/m^3.', '\033[37m')

        if isinstance(data, list) is False and isinstance(data, np.ndarray) is False:
            raise Exception('Input data must be either a list or a numpy array.')
        if len(data) != 10: 
            raise Exception('Input data lenght must be 10')
        for i in data:
            if isinstance(i, float) is False and isinstance(i, int) is False:
                raise Exception('Input data values must be either float or integer.')
        if isinstance(print_cntrl, bool) is False:
            raise Exception('Error. The variable "print_cntrl" must be boolean.')     
        
        check_ply_properties = True
        try:
            self.E1 
        except AttributeError:
            check_ply_properties = False
        if check_ply_properties is False:    
            raise Exception('Error. E1 is not defined. Use one of the available method to compute the ply properties before computing the error.')

        self.error_percent_E1 = np.abs((data[0] - self.E1)) / data[0] * 100
        self.error_percent_E2 = np.abs((data[1] - self.E2)) / data[1] * 100
        self.error_percent_E3 = np.abs((data[2] - self.E3)) / data[2] * 100
        self.error_percent_ni12 = np.abs((data[3] - self.ni12)) / data[3] * 100
        if self.ni23 == 'NA':
            self.error_percent_ni23 = 'NA'
        else:
            self.error_percent_ni23 = np.abs((data[4] - self.ni23)) / data[4] * 100
        self.error_percent_ni13 = np.abs((data[5] - self.ni13)) / data[5] * 100
        self.error_percent_G12 = np.abs((data[6] - self.G12)) / data[6] * 100
        self.error_percent_G23 = np.abs((data[7] - self.G23)) / data[7] * 100
        self.error_percent_G13 = np.abs((data[8] - self.G13)) / data[8] * 100
        self.error_percent_rho = np.abs((data[9] - self.rho)) / data[9] * 100

        self.error_E1 = self.error_percent_E1 / 100 
        self.error_E2 = self.error_percent_E2 / 100 
        self.error_E3 = self.error_percent_E3 / 100 
        self.error_ni12 = self.error_percent_ni12  / 100
        self.error_ni13 = self.error_percent_ni13  / 100
        if self.ni23 == 'NA':
            self.error_ni23 = 'NA'
        else: 
            self.error_ni23 = self.error_percent_ni23 / 100
        self.error_G12 = self.error_percent_G12 / 100
        self.error_G23 = self.error_percent_G23 / 100
        self.error_G13 = self.error_percent_G13 / 100
        self.error_rho = self.error_percent_rho / 100

        self.errors = [self.error_E1, self.error_E2, self.error_E3, self.error_G12, \
                       self.error_G23, self.error_G13, self.error_ni12, self.error_ni23, self.error_ni13, self.error_rho]
        self.errors_percent = [self.error_percent_E1, self.error_percent_E2, self.error_percent_E3, \
                               self.error_percent_G12, self.error_percent_G23, self.error_percent_G13, \
                                self.error_percent_ni12, self.error_percent_ni23, self.error_percent_ni13]
        if print_cntrl is True:
            print('error E1 % = ', str(np.round(self.error_percent_E1, 4)))
            print('error E2 % = ', str(np.round(self.error_percent_E2, 4)))
            print('error E3 % = ', str(np.round(self.error_percent_E3, 4)))
            print('error ni12 % = ', str(np.round(self.error_percent_ni12, 4)))    
            if self.ni23 != 'NA':
                print('error ni23 % = ', str(np.round(self.error_percent_ni23, 4)))
            else: 
                print('error ni23 %  = NA')
            print('error ni13 % = ', str(np.round(self.error_percent_ni13, 4)))
            print('error G12 % = ', str(np.round(self.error_percent_G12, 4)))
            print('error G23 % = ', str(np.round(self.error_percent_G23, 4)))
            print('error G13 % = ', str(np.round(self.error_percent_G13, 4)))
            print('error rho % = ', str(np.round(self.error_percent_rho, 4)))
