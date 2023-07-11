__author__ = 'Edoardo Mancini'
import numpy as np
import math
from typing import List
from math import sin as sin
from math import cos as cos
from numpy.linalg import inv
import matplotlib.pyplot as plt

class Laminate: 
    '''
# References: 
    - Johnes, Mechanics of Composite Materials, 2nd edition
    - Barbero, Introduction to Composite Materials Design, 3rd edition 
    - ANSYS, theory reference for mechanical applications, Release 12.0
# Description:
    The class creates a laminate object, inizitalized according to the inputed stack-up. Additional methods allow the user to compute 
    composite properties such as equivalent stiffness, stresses and deformations

'''
# Initialization method 
    def __init__(self, *varargs:list[str, list or np.ndarray, list or np.ndarray], mech_prop_units:str='GPa', hide_text:bool=False):
        '''
# DESCRIPTION:
    In the init method the single plies are defined and their properties exctracted 

    # INPUTS:

    ## Required
        -   *varargs: variable number of three-elements lists. Each list represents a block of plies to be stack-up on to of the previous. 
        Again, each of these block is identified with a three-elements list containing: the material name (str), the properties (11 numerical elements list) 
        and the stackup (numerical list of ply orientations - angles in deg).
        mech_props vect/list:
        - E1: GPa / MPa
        - E2: GPa / MPa
        - E3: GPa / MPa
        - ni13: adim
        - ni12: adim
        - ni23: adim
        - G12: GPa / MPa
        - G23: GPa / MPa
        - G13: GPa / MPa
        - rho: kg/m^3
        - ply thickness: mm  
        # note: If E3, ni13 and ni23 are not used their value can be set to 0.
        
        ####     Attention: 
            the plyblocks must be input from bottom to top
            
        ####     Attention: 
            this function assumes that each ply is orthotropic

    ## Optional
        - mech_prop_units: choose between providing the mechanical properties in MPa or GPa. Default is GPa.

    # OUTPUTS:
        None
    # ATTRIBUTES:
        self.Q          # list of CLT plies' stiffness matrices in the ply reference system
        self.Qstar          # list of FSDT plies' stiffness matrices in the ply reference system
        self.Qs         # Q + Qstar   
        self.Qbars          # list of CLT plies' stiffness matrices in the laminate reference system
        self.Qbars_star         # list of FSDT plies' stiffness matrices in the laminate reference system
        self.density_avg        # average laminate density
        self.thickness_vector           # list of plies' thickness bot to top
        self.thickness          # overall laminate thickness
        self.plies          # list of different plies types (one for each input triplette)
        self.stackup            # list of vocabuliers with one element per ply containing the ply name, number, angle and, material 
        self.stackup_long           # extended version of stackup including also the plies' mechanical properties
        self.stackup_legend         # legend for the stackup list 
        self.stackup_long_legend         # legend for the stackup long list
        self.no_of_plies            # number of plies in the laminate
        self.A          # A matrix 
        self.B          # B matrix 
        self.D          # D matrix
        self.H          # H matrix
        self.stiff          # assembled stiffness matrix CLT (A, B, D )
        self.stiff_full          # assembled stiffness matrix FSDT (A, B, D, H)
        self.inv_stiff          # inverse of the stiffness matrix CLT (A, B, D )
        self.inv_stiff_full          # inverse of the stiffness matrix FSDT (A, B, D, H)

    # Example:

        import PyComp as comp
        ply_name = 'Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
  

    '''         
        # USER WARNINGS        
        if hide_text is False:
            print('\033[35m', 'Assumed input density\'s units are kg/m3', '\033[37m',' ')
            print('\033[35m', 'Assumed input ply thickness\' units are mm', '\033[37m',' ')
            print('\033[35m', 'Assumed input angle\'s units are degs', '\033[37m',' ')

        if mech_prop_units != 'MPa' and mech_prop_units != 'GPa' :
            raise Exception('The variable "mech_prop_units" must be either equal to "MPa" or to "GPa".')

        ## INITIALIZATION
        # A block it is defined as a set of n layers of the same material with any orientation
        blocks = varargs
        stackup = []
        stackup_list = []
        stackup_list_long = []
        plies = []
        # input checks, data massaging and storage performed on each of the input blocks
        for i, ply_block in enumerate(blocks):

            ## CHECKS

            # The function check that the data were provided in the correct format [name, mech_props, orientations]
            if len(ply_block) != 3:
                raise Exception('Error in ply_block no. '+ str(i + 1) \
                    + '. Wrong ply_block format. Ply tuple must have three elements.')
            if isinstance(ply_block[0], str) is False:
                raise Exception('Error in ply_block no. '+ str(i + 1) \
                    + '. Wrong ply_block format. The firse element (ply_block name) of ply_block np. '+ str(i + 1) + ' must be a string.')
            if isinstance(ply_block[1], list) is False and isinstance(ply_block[1], np.ndarray) is False:
                raise Exception('Error in ply_block no. '+ str(i + 1) \
                    + '. Wrong ply_block format. The second element (mechanical properties) of the ply_block no. ' \
                        + str(i + 1) + ' must be either a list or a numpy array.')
            if len(ply_block[1]) != 11:
                raise Exception('Error in ply_block no. '+ str(i + 1) \
                    + '. Wrong ply_block format. The lenght of the second element (mechanical properties) of the ply_block no. ' + str(i + 1) + ' is not 12.')
            for k in ply_block[1]:
                if isinstance(k, float) is False and isinstance(k, int) is False:
                    raise Exception('Error in ply_block no. ' + str(i + 1) \
                        + ' one or more elements of the mech_props input is neither a float nor an integer')
            #### Include check on the vector components nature (int or float)
            if isinstance(ply_block[2], list) is False and isinstance(ply_block[2], np.ndarray) is False:
                raise Exception('Error in ply_block no. '+ str(i + 1) \
                    + '. Wrong ply_block format. The third element (stackup) of the ply_block no. ' + str(i + 1) + ' must be either a list or a numpy array.')
            if len(ply_block[2]) < 1:
                raise Exception('Error in ply_block no. '+ str(i + 1) \
                    + '. Wrong ply_block format. The third element (stackup) of the ply_block no. ' + str(i + 1) + ' must have at least one element.')
            for kk in ply_block[2]:
                if isinstance(kk, float) is False and isinstance(kk, int) is False:
                    raise Exception('Error in ply_block no. ' + str(i + 1) \
                        + ' one or more elements of the angle input is neither a float nor an integer')
            
            ## DATA STORAGE AND MASSAGING
            # round the value of the ply mechanical properties
            for idx, kkk in enumerate(ply_block[1]):
                ply_block[1][idx] = self.__round(kkk, digits=4)

            # creation of a dictionary containing only the plies properties
            ply_list_elements = {'name': ply_block[0], 'ni12': ply_block[1][3], 'ni13': ply_block[1][4], 'ni23': ply_block[1][5], \
                        'density': ply_block[1][9], 'thickness': ply_block[1][10]}
            if mech_prop_units == 'GPa':
                if hide_text is False:
                    print('\033[35m', 'Assumed input elastic modulus\' units are GPa', '\033[37m',' ')
                ply_list_elements.update({ 'E1': ply_block[1][0] * 1000, 'E2': ply_block[1][1] * 1000, 'E3': ply_block[1][2] * 1000,\
                        'G12': ply_block[1][6] * 1000, 'G13': ply_block[1][7] * 1000, 'G23': ply_block[1][8] * 1000})
            else:
                ply_list_elements.update({'E1': ply_block[1][0], 'E2': ply_block[1][1], 'E3': ply_block[1][2],\
                        'G12': ply_block[1][6], 'G13': ply_block[1][7], 'G23': ply_block[1][8]})
                if hide_text is False:
                    print('\033[35m', 'Assumed input elastic modulus\' units are MPa', '\033[37m',' ')
            plies.append(ply_list_elements)

            # the code is already looping over the input blocks (each block is a set of same-material plies with diffrent orientations)  
            # this second loop sift through the plies in each block to store their information in dictionaries and lists
            for j in range(len(ply_block[2])):
                # stackup = [[ply_1_mech_props, ply_1_orientation],....., [ply_n_mech_props, ply_n_orientation]] to be used later
                stackup.append([ply_block[1], ply_block[2][j]]) 
                # in the latter two dictionaries are created "ply_long" and "ply" 
                # the first includes all the properties of each ply, the second only name, thickness and angle
                ply_long = {'name': ply_block[0], 'ply no. ': len(stackup_list) + 1, 'ni12': ply_block[1][3], 'ni13': ply_block[1][4], 'ni23': ply_block[1][5], \
                        'density': ply_block[1][9], 'thickness': ply_block[1][10], 'angle': ply_block[2][j]}
                # material properties are always stored in MPa therefore if the input is in GPa the properties have to be rescaled                
                if mech_prop_units == 'GPa':
                    ply_long.update({'E1': ply_block[1][0] * 1000, 'E2': ply_block[1][1] * 1000, 'E3': ply_block[1][2] * 1000,\
                        'G12': ply_block[1][6] * 1000, 'G13': ply_block[1][7] * 1000, 'G23': ply_block[1][8] * 1000})
                # in case the input is in MPa there is no need to rescale 
                else:
                    ply_long.update({'E1': ply_block[1][0], 'E2': ply_block[1][1], 'E3': ply_block[1][2],\
                        'G12': ply_block[1][6], 'G13': ply_block[1][7], 'G23': ply_block[1][8]})
                # the second dictionary is easier and does not require any "if" loop
                ply = {'name': ply_block[0], 'ply no. ': len(stackup_list) + 1, 'thickness': ply_block[1][10], 'angle': ply_block[2][j]}
                
                # here two list are created. Their lenght is equal to the number of plies while their elements are dictionaries containing the plies properties
                # short version
                stackup_list.append(ply)
                # long version
                stackup_list_long.append(ply_long)
        self.plies = plies
        # additional class propeties containing the description of the content of the former dictionaries
        self.stackup_legend = ['ply name', 'thickness [mm]', 'ply angle [deg]']
        self.stackup_long_legend = ['ply name', 'E1 [MPa]', 'E2 [MPa]', 'E3 [MPa]', 'ni12', \
            'ni13', 'ni23', 'G12 [MPa]', 'G13 [MPa]', 'G23 [MPa]','thickness [mm]', 'ply angle [deg]']
        
        # Definition of the number of plies as a class property
        self.no_of_plies = len(stackup)

        # Definition of the stackup lists as class properties
        self.stackup = stackup_list
        self.stackup_long = stackup_list_long

        ## MACROMECHANICS - Now the code proceeds with the calculatin of the laminate mechanical properties

        # Vectors initialization
        Qs = []
        Qbars = []
        Qbars_star = []
        dens = []
        thick = []
        z_coord = [0]
        
        # Looping on the plies to get the matrix of Q matrices, a matrix of Q start matrices, the oredered list of densities and the ordered list of thicknesses
        for i in range(self.no_of_plies): 
            # Qs matrices generation
            temp_Qs = self.__createQs(stackup[i][0], units = mech_prop_units)
            # Qs matrices rotation
            temp_Q_bar, temp_Q_star_bar = self.__rotation(temp_Qs[0], temp_Qs[1], stackup[i][1])
            # Definition of the z coordinate of each ply starting from the 0 point at the bottom. ATTENTION: the z_coordinate vector starts from 0 so it will be one point longer than the thickness vector
            temp_z_coordinate = z_coord[-1] + stackup[i][0][10] 
        # Vectors generation
            # Qs vector of matrices creation
            Qs.append(temp_Qs)
            # Qbars vector of matrices creation
            Qbars.append(temp_Q_bar)
            # Qbars star vector of matrices creation
            Qbars_star.append(temp_Q_star_bar)
            # 0-at-the-bottom z coordinate vector creation
            z_coord.append(temp_z_coordinate)
            # Density vector creation
            dens.append(stackup[i][0][9])
            # Thickness vector creation
            thick.append(stackup[i][0][10])

        # Calculation of the overall thickness and creation of the respective class property
        thickness = np.sum(thick)

        # Definition of the mid plane and coordinates rescaling
        z_mid = thickness / 2 
        z_coord = np.array(z_coord)
        mid_z_coord = z_coord - z_mid
        
        Q = []
        Qstar = []
        for i in Qs:
            Q.append(i[0])
            Qstar.append(i[1])
        self.Q = Q
        self.Qstar = Qstar
        self.Qs = Qs
        self.Qbars = Qbars
        self.Qbars_star = Qbars_star
        self.__density_vector = dens
        self.density_avg = np.average(self.__density_vector)
        self.thickness_vector = thick
        self.__mid_z_coord = mid_z_coord
        self.thickness = thickness

        self.__compute_laminate_stiffeness_matrices()
# Service class method (internal start with "__")
    # performs a matrix rotation
    def __rotation_matrices(self, angle:int) -> np.ndarray:
        '''
# DESCRIPTION:
    Internal method to caltulate the rotation matrix T.

    # INPUTS:

    ## Required
        -   angle: rotation angle in degs
    # OUTPUTS: 
        -   T: numpy array representing the rotation matrix T which brings a value from the reference system of the lamniate to the lamina's

    ### ATTENTION: the angle is expected in degs. 
# Example:
    NA - internal method
    
    '''     
        
        # Change from degrees to radiants
        ang = angle * math.pi * 2 / 360

        # Definition of the inverse rotation matrix
        T = np.array([[cos(ang) ** 2, sin(ang) ** 2, 2 * sin(ang) * cos(ang)], 
                    [sin(ang) ** 2, cos(ang) ** 2, - 2 * sin(ang) * cos (ang)], 
                    [- sin(ang) * cos(ang), sin(ang) * cos(ang), cos(ang) ** 2 - sin(ang) ** 2]]
                    )
        T_star = np.array([[cos(ang), sin(ang)], [-sin(ang), cos(ang)]]) ## new

        return (T, T_star)
    # rotates Q and Q_star
    def __rotation(self, Q:np.ndarray, Q_star:np.ndarray, angle:int) -> tuple[np.ndarray, np.ndarray]:
        '''
# DESCRIPTION:
    Internal method to rotate the Q and Qstar matrices

    # INPUTS:

    ## Required
        -   Q: Q matrix 
        -   Q_star: Qstar matrix 
        -   angle: rotation angle
    # OUTPUTS: 
        -   [Q, Qstar]: tuple containing two numpy array. The first is Q_bar(rotated), the second is Qstar_bar(rotated).

    ### ATTENTION: no checks are performed, inputs are check in the external methods. 
# Example:
    NA - internal method
    
    '''       
        
        # Definition of the inverse rotation matrix
        T, T_star = self.__rotation_matrices(angle)
        T_inv = np.linalg.inv(T)
        # Traspose of the rotation matrix
        T_inv_trasp = np.transpose(T_inv)

        # Same operations for the "star" matrices (the ones associated with the out-of-plane shear)
        T_star_inv = np.linalg.inv(T_star)
        T_star_inv_trasp = np.transpose(T_star_inv)

        # Rotation of Q
        Q_bar = T_inv.dot(Q.dot(T_inv_trasp))
        # Rotation of Q start
        Q_star_bar = T_star_inv.dot(Q_star.dot(T_star_inv_trasp))

        # Removal of numerical zeros (values below the threshold)
        thrshld = 1e-8
        idx = np.where(abs(Q_bar) < thrshld)
        Q_bar[idx] = 0
        idx_star = np.where(abs(Q_star_bar) < thrshld)
        Q_star_bar[idx_star] = 0
        
        return(Q_bar, Q_star_bar)
    # creates a list of the different kinds of plies in the laminate 
    # e.g. if the laminate has 10 plies of a specific material and 5 of another the list will contain 2 elements (kind1 and kin2)
    def __diffent_plies_list_generation(self):
        different_plies_list = []
        plies_name_list = []
        # two lists are created in the for loop. The first, plies_name_lists, which contains all the plies names in order. 
        # The second, different_plies_lists, contains with the name of the diffrent plies in the laminate repeated only once independently 
        # on their disposisition and their repetitions
        for i in range (len(self.stackup)):
            temp = (self.stackup[i].get('name'))
            plies_name_list.append(temp)
            if temp not in different_plies_list:
                different_plies_list.append(temp)
            del temp ###
        self.__different_plies_list = different_plies_list
    # rounds a number
    def __round(self, input:float or int, digits:int = 3, lower_limit:float = 1e-22) -> float or int: 
        '''
# DESCRIPTION:
    this method is used to round the output values. Independetly on the number of zeros the final digits are always the same. e.g. 10.3454 -> 10.34, 1.3454 -> 1.34, 0.0013454 -> 0.00134

    # INPUTS:

    ## Required
        -   input: number to round.
    ## Optinal
        - digits: number of non zero digits after the ".".
        - lower_limit: minimum value below which the number is not rounded but returned as provided. Needed to avoid an infinite loop when one of the input is 0.
    # OUTPUTS: 
        -   output: rounded number. 

    # Example:
    NA - internal 
'''
        #inputs check
        if isinstance(input, int) is False and isinstance(input, float) is False:
            raise Exception('The input is neither a float nor an integer.')
        if isinstance(digits, int) is False:
            raise Exception('The vairble "digits" is not an integer.')
        
        # take the absolute value
        number = abs(input)
        # check if the input is zero. In case it is the rounding is skipped otherwise it will cause an infinite loop. 
        if number < lower_limit: 
            output = input
            return output
        else:
            check = number
            # define an increment variable
            k = 0
            # check how many zeros follow the "." by multipling the number by 10 k times, until it is equal or higher than 1. 
            while check <= 1: 
                check = check * 10 
                k = k + 1
            # if k != 0 the number is smaller than 0. Perform the round operation with an increased number of digits
            if k > 0: 
                output = round(input, k + digits - 1)
            # if k == 0 the number was already larger than 0. Perform the normal round operation 
            else:
                output = round(input, digits)

            return output
    # rounds matrices
    def __round_matrices(self, input:np.ndarray or list, digits:int = 3, lower_limit: float = 1e-22) -> np.ndarray or list: 
        '''
# DESCRIPTION:
    implements the __round function on matrices or vectors

    # INPUTS:

    ## Required
        -   input: array/list to round.
    ## Optinal
        - digits: number of non zero digits after the ".".
        - lower_limit: minimum value below which the number is not rounded but returned as provided. Needed to avoid an infinite loop when one of the input is 0. 

    # OUTPUTS: 
        -   output: rounded array/list. 

    # Example:
    NA - internal 
'''
        # inputs check
        if isinstance(digits, int) is False:
            raise Exception('The vairble "digits" is not an integer.')
        
        if isinstance(input, np.ndarray) is False and isinstance(input, list) is False:
            raise Exception('Error. The input is neither a numpy array nor a list.')
        # array input
        elif isinstance(input, np.ndarray) is True:
            array = input
            # matrix input
            try :
                # initialize the output array
                output = np.zeros((input.shape[0], input.shape[1]))
                # loop over the array
                for idx1, vector in enumerate (array):
                    for idx2, element in enumerate (vector):
                        # discern between positive and negative values
                        if element >= 0:
                            output[idx1, idx2] = self.__round(element, digits=digits, lower_limit=lower_limit)
                        else: 
                            output[idx1, idx2] = - self.__round(abs(element), digits=digits, lower_limit=lower_limit)
            # vector input with only 1D
            except IndexError:
                # initialize the output array
                output = np.zeros((input.shape[0]))
                # loop over the array (only 1D)
                for idx1, vector in enumerate (array):
                    # discern between positive and negative
                    if vector >= 0:
                        output[idx1] = self.__round(vector, digits=digits, lower_limit=lower_limit)
                    else:
                        output[idx1] = - self.__round(abs(vector), digits=digits, lower_limit=lower_limit)
            
        # list inputs
        else: 
            array = input 
            output = []
            # more than 1D
            try:
                # loop
                for idx1, vector in enumerate (array):
                    output_temp = []
                    for idx2, element in enumerate (vector):
                        # discern positive and negative
                        if element >= 0:
                            output_temp.append(self.__round(element, digits=digits, lower_limit=lower_limit))
                        else: 
                            output_temp.append(- self.__round(abs(element), digits=digits, lower_limit=lower_limit))
                    # stack the rounded value
                    output.append(output_temp)
            # 1D
            except TypeError:
                # loop
                for idx1, vector in enumerate (array):
                    output_temp = 0
                    # discer positive and negative
                    if vector >= 0:
                        output_temp = (self.__round(vector, digits=digits, lower_limit=lower_limit))
                    else: 
                        output_temp = (- self.__round(abs(vector), digits=digits, lower_limit=lower_limit))
                    # stack the value
                    output.append(output_temp)
        # return output
        return output
    # create the laminae's Q matrices from the mechanical properties 
    def __createQs(self, mech_props:list, units:str) -> tuple[np.ndarray, np.ndarray]:
        '''
# DESCRIPTION:
    Internal method to compute the Q matrix for a ply starting from the its mechanical properties 

    # INPUTS:

    ## Required
        -   mech_props: mechanical properties of the ply 
        -   units: either GPa or MPa
    # OUTPUTS: 
        -   [Q, Qstar]: tuple containing two numpy array. The first is the Q matrix, the second the Qstar matrix (out-of-plane shear).

    ### ATTENTION: no checks are performed, inputs are check in the external methods. 
    # Example:
        NA - internal method
    
    '''       
        # Exctraction of the single properties from the list
        if units == 'MPa' :
            E1 = mech_props[0]
            E2 = mech_props[1]
            ni12 = mech_props[3]
            G12 = mech_props[6]
            G13 = mech_props[7]
            G23 = mech_props[8]
            delta = 1 - ni12 ** 2 * E2 / E1

        elif units == 'GPa' :
            E1 = mech_props[0] * 1000
            E2 = mech_props[1] * 1000
            ni12 = mech_props[3] 
            G12 = mech_props[6] * 1000
            G13 = mech_props[7] * 1000
            G23 = mech_props[8] * 1000
            delta = 1 - ni12 ** 2 * E2 / E1

        # Definition of the matrix terms
        Q11 = E1 / delta
        Q12 = ni12 * E2 / delta 
        Q21 = Q12 
        Q22 = E2 / delta
        Q16 = 0
        Q61 = Q16
        Q26 = 0
        Q62 = Q26
        Q66 = G12
        Q44 = G23
        Q55 = G13
        Q45 = 0
        Q54 = Q45

        # Building of matrix Q and Q*
        return (np.array([[Q11, Q12, Q16], [Q21, Q22, Q26], [Q61, Q62, Q66]]), np.array([[Q44, Q45], [Q54, Q55]]))
    # compute the laminate stiffness matrices (single -A, B, D, H - and global - 6x6 or 8x8 matrices including all the single -)
    def __compute_laminate_stiffeness_matrices(self):
        '''
# DESCRIPTION:
    Internal method to compute stiif, stiff_full, inv_stiff, inv_stiff_full, A, B, D, H, a, b, d, h

    # INPUTS:

    ## Required
        -   mech_props: mechanical properties of the ply 
        -   units: either GPa or MPa
    # OUTPUTS: 
        -   [Q, Qstar]: tuple containing two numpy array. The first is the Q matrix, the second the Qstar matrix (out-of-plane shear).

    ### ATTENTION: no checks are performed, inputs are check in the external methods. 
    # Example:
        NA - internal method
    
        '''

                # calculation of the laminate mechanical properties
        A11 = 0; A12 = 0; A16 = 0; A21 = 0; A22 = 0; A26 = 0; A61 = 0; A62 = 0; A66 = 0
        B11 = 0; B12 = 0; B16 = 0; B21 = 0; B22 = 0; B26 = 0; B61 = 0; B62 = 0; B66 = 0        
        D11 = 0; D12 = 0; D16 = 0; D21 = 0; D22 = 0; D26 = 0; D61 = 0; D62 = 0; D66 = 0
        H11 = 0; H12 = 0; H21 = 0; H22 = 0

        z_cd = self.__mid_z_coord
        for i in range(self.no_of_plies): 
            temp_A = self.Qbars[i] * (z_cd[i + 1] - z_cd[i])
            A11 = A11 + temp_A[0, 0]
            A12 = A12 + temp_A[0, 1]
            A16 = A16 + temp_A[0, 2]
            A21 = A21 + temp_A[1, 0]
            A22 = A22 + temp_A[1, 1]
            A26 = A26 + temp_A[1, 2]
            A61 = A61 + temp_A[2, 0]
            A62 = A62 + temp_A[2, 1]
            A66 = A66 + temp_A[2, 2]
            
            temp_B = .5 * self.Qbars[i] * (z_cd[i + 1] ** 2 - z_cd[i] ** 2)
            B11 = B11 + temp_B[0, 0]
            B12 = B12 + temp_B[0, 1]
            B16 = B16 + temp_B[0, 2]
            B21 = B21 + temp_B[1, 0]
            B22 = B22 + temp_B[1, 1]
            B26 = B26 + temp_B[1, 2]
            B61 = B61 + temp_B[2, 0]
            B62 = B62 + temp_B[2, 1]
            B66 = B66 + temp_B[2, 2]
            
            temp_D = (1 / 3) * self.Qbars[i] * (z_cd[i + 1] ** 3 - z_cd[i] ** 3)
            D11 = D11 + temp_D[0, 0]
            D12 = D12 + temp_D[0, 1]
            D16 = D16 + temp_D[0, 2]
            D21 = D21 + temp_D[1, 0]
            D22 = D22 + temp_D[1, 1]
            D26 = D26 + temp_D[1, 2]
            D61 = D61 + temp_D[2, 0]
            D62 = D62 + temp_D[2, 1]
            D66 = D66 + temp_D[2, 2]
            
            temp_H = (5 / 4) * self.Qbars_star[i] * (self.thickness_vector[i] - 4 / self.thickness ** 2 * \
                (self.thickness_vector[i] * ((z_cd[i + 1] + z_cd[i]) / 2) ** 2 + self.thickness_vector[i] ** 3 / 12))
            H11 = H11 + temp_H[0, 0]
            H12 = H12 + temp_H[0, 1]
            H21 = H21 + temp_H[1, 0]
            H22 = H22 + temp_H[1, 1]

        A = np.array([[A11, A12, A16], [A21, A22, A26], [A61, A62, A66]]) 
        B = np.array([[B11, B12, B16], [B21, B22, B26], [B61, B62, B66]]) 
        D = np.array([[D11, D12, D16], [D21, D22, D26], [D61, D62, D66]]) 
        H = np.array([[H11, H12], [H21, H22]]) 
        
        thrshld = 1e-8
        idx_A = np.where(abs(A)< thrshld)
        A[idx_A] = 0
        idx_B = np.where(abs(B) < thrshld)
        B[idx_B] = 0
        idx_D = np.where(abs(D) < thrshld)
        D[idx_D] = 0
        idx_H = np.where(abs(H) < thrshld)
        H[idx_H] = 0

        # define the single stiffness matrices
        self.A = self.__round_matrices(A) 
        self.B = self.__round_matrices(B) 
        self.D = self.__round_matrices(D) 
        self.H = self.__round_matrices(H)

        # generate the 6x6 stiffness matrix (A, B, D) - (N, M) = (self.stiff).dot(epsilon0, k) 
        self.stiff = np.vstack((np.hstack((self.A, self.B)), np.hstack((self.B, self.D))))
        # invert the stiffness matrix
        self.inv_stiff = self.__round_matrices(inv(self.stiff))

        # generate the 8x8 stiffness matrix (A, B, D, H) - (N, M, V) = (self.stiff_full).dot(epsilon0, k, gamma_out_of_plane) 
        # increase the matrix size from 6x6 to 8x8
        zeros = np.zeros((6))
        zeros2 = np.zeros((8))
        # add two 1x6 zero coulums at the end of the matrix - stiff_full_temp is a 6x8 matrix
        stiff_full_temp = np.c_[np.c_[self.stiff, zeros], zeros]
        # add two 1x8 zero rows at the end of the matrix - stiff_full_temp2 is an 8x8 matrix
        stiff_full_temp2 = np.vstack((np.vstack((stiff_full_temp, zeros2)), zeros2))
        # loop over H and change the zeros in postions [6, 6], [6, 7], [7, 6] and [7, 7] with the matrix H 
        # (rememeber that the first position is [0, 0] so [7, 7] is the last element of an 8x8 matrix) 
        for i in range(2):
            for j in range (2):
                stiff_full_temp2[i + 6, j + 6] = self.H[i, j]
        # declare the class property stiff_full and inv_stiff_full (perfom the inversion at the same time)
        self.stiff_full = stiff_full_temp2
        self.inv_stiff_full = self.__round_matrices(inv(self.stiff_full))

        # comunte the inverse of the signle stifness matrices, store and round them
        self.a = self.inv_stiff[:3, :3] 
        self.b = self.inv_stiff[:3, 3:6]
        self.d = self.inv_stiff[3:6, 3:6]
        self.h = self.inv_stiff_full[6:8, 6:8]
   # compute the hygrothermal stresses in the laminate
    def __ply_hygrothermal_stresses(self, cntrl_case:int, T_in:int or float = 0, m_in: int or float = 0):

        #### WORK IN PROGRESS #####
        dummy_var = cntrl_case
    # calculate the ply stress for a given loading condition    
    def __ply_stresses(self, N_in:list or np.array, M_in:list or np.array, V_in:list or np.array, \
        cntrl_external_actions:int = 0, T_in: int or float = 0, m_in: int or float = 0):
        '''
# DESCRIPTION:
    Internal method to compute the ply stresses in the laminate reference system

    # INPUTS:

    ## Required
        -   N_in: input forces per unit lenght. Three elements vector or list [Nx, Ny, Nxy] - [N/m]. 
        -   M_in: input moments per unit lenght. Three elements vector or list [Mx, My, Mxy] - [N].
    ## Optional
        -   V_in: input out-of-plane forces per unit lenght. Three elements vector or list [Vx, Vy] - [N/m].
        -   cntrl_external_actions: defines which actions are applied on the laminate. Case 0: only mechanical actions. /
        Case 1: mechanical and thermal. Case 2: mechanical and mositure. Case 3: mechanical, thermal and moisture. 
        -   T_in: laminate temperature [K].
        -   m_in: laminate moisture absorption [%of_mass].
    # OUTPUTS: 
        -   None: the stresses are stored as class properties

    # Example:
        NA - internal method

    ### ATTENTION: checks are not performed in internal methods
    
    '''     
        # Note: there is no check because the check were already preformed in the calling function. N_in and M_in can only be lists or arrays
        
        ## OPERATIONS
        # change the format of the inputs if given as lists or rename the input variable (in case they are already arrays)
        if isinstance(N_in, list) is True:
            N = np.array(N_in)
        else: 
            N = N_in
        if isinstance(M_in, list) is True:
            M = np.array(M_in)   
        else:
            M = M_in
    ## NEW
        if isinstance(V_in, str) is False:
            if isinstance(V_in, list) is True:
                V = np.array(V_in)   
            else:
                V = V_in
            skipV = 0
        else: 
            skipV = 1
    
            # input check
        cntrl_e_a = cntrl_external_actions
        if cntrl_e_a != 0 and cntrl_e_a != 1 and cntrl_e_a != 2 and cntrl_e_a != 3:
            raise Exception('the variable "cntrl_external_actions" must be either equal to 0, 1, 2 or 3.')
        # computation of Mt and Nt, to be added to N and M 
        # no temperature nor moisture effects
        elif cntrl_e_a == 0:
            Nt = np.array([0, 0, 0])
            Mt = np.array([0, 0, 0])
        # mechanical and thermal
        elif cntrl_e_a == 1:
            raise Exception('this portion of the code is still under development and will be available in the next releases')
            # Nt, Mt = self.__ply_hygrothermal_stresses(cntrl_external_actions, Tsf=T_in)
        # mechanical and hygro
        elif cntrl_e_a == 2:
            raise Exception('this portion of the code is still under development and will be available in the next releases')
            # Nt, Mt = self.__ply_hygrothermal_stresses(cntrl_external_actions, msf=m_in)
        # mechanical and hygrothermal
        else:
            raise Exception('this portion of the code is still under development and will be available in the next releases')
            # Nt, Mt = self.__ply_hygrothermal_stresses(cntrl_external_actions, Tsf=T_in, msf=m_in)
    ##
       
        # comupute the deformation with the inverted version of (N, M) = [[A, B], [B, D]] (eps0, k)
        # - theory reminder - Classical Lamination Theory assumes that the deformation of the object is due to two contributions: 
        # the deformation of the mid line (neutral axis) epsilon0 and the stretching or compression due to the curvature (k*z)
        N = N + Nt
        M = M + Mt
        deformation = self.inv_stiff.dot(np.hstack((N, M)))

        # exctract epsilon0 and k from the deformation vector
        self.epsilon0 = deformation[0:3]
        self.k0 = deformation[3:]

        Qbars = np.array(self.Qbars)
        Qbar_star = np.array(self.Qbars_star) ## new
        sigma_top = np.zeros((self.no_of_plies, 3))
        sigma_bot = np.zeros((self.no_of_plies, 3))
        def_top = np.zeros((self.no_of_plies, 3))
        def_bot = np.zeros((self.no_of_plies, 3))

        if skipV == 0:
            shear_stiff = self.H ## new
            self.inv_shear_stiff = inv(shear_stiff) ## new
            shear_deformation_avrg = self.inv_shear_stiff.dot(V) ## new
            self.gamma_oop_avrg = shear_deformation_avrg  ## new
            tau_oop_top = np.zeros((self.no_of_plies, 2)) ## new
            tau_oop_bot = np.zeros((self.no_of_plies, 2)) ## new
            tau_oop_mid = np.zeros((self.no_of_plies, 2)) ## new

        for i in range(self.no_of_plies): 
            
            # The stress is linear in the ply therefore I only need the stress at the top and bottom of each ply
            # The stress is computed through the deformation. Hence, I need the deformation at the top and the bottom
            # The deformation is epsiolon0 + k * z so I need the z coordinate of the top and bottom since e0 and k are already known
            # The deformation is computed from the center line of the laminate and so is the z coordinate. Therefore, negative values of z are allowed (below the mid line)
            # Hence, k0 is multiplied by z_top and z_bottom. z_bottom is the thickness of the laminate up to the i-th ply minus half of the laminate thickness (move the zero to the center of the laminate)
            # z_top is z_bottom plus the ply thickness
            z_top = (sum(self.thickness_vector[:i]) - self.thickness / 2 + self.thickness_vector[i])
            def_ply_top = (self.epsilon0 + self.k0 * z_top)
            # retrievement of the stresses 
            sigma_top_temp_x = (Qbars[i].dot(def_ply_top))[0]          
            sigma_top_temp_y = (Qbars[i].dot(def_ply_top))[1]          
            sigma_top_temp_xy = (Qbars[i].dot(def_ply_top))[2]

            sigma_top_temp = np.hstack((sigma_top_temp_x, sigma_top_temp_y, sigma_top_temp_xy))          
            sigma_top[i, :] = sigma_top_temp 
            def_top_temp = np.hstack((def_ply_top[0], def_ply_top[1], def_ply_top[2])) 
            def_top[i, :] = def_top_temp

            z_bot = (sum(self.thickness_vector[:i]) - self.thickness / 2)
            def_ply_bot = (self.epsilon0 + self.k0 * z_bot)
            sigma_bot_temp_x = (Qbars[i].dot(def_ply_bot))[0]
            sigma_bot_temp_y = (Qbars[i].dot(def_ply_bot))[1]
            sigma_bot_temp_xy = (Qbars[i].dot(def_ply_bot))[2]

            # stresses matrix creation
            sigma_bot_temp = np.hstack((sigma_bot_temp_x, sigma_bot_temp_y, sigma_bot_temp_xy))          
            sigma_bot[i, :] = sigma_bot_temp
            def_bot_temp = np.hstack((def_ply_bot[0], def_ply_bot[1], def_ply_bot[2])) 
            def_bot[i, :] = def_bot_temp

            if skipV == 0:
                # gamma_oop_top = self.gamma_oop_avrg * 5/4 * (1 - (z_top / (self.thickness / 2)) ** 2)
                gamma_oop_top = self.gamma_oop_avrg * 5/6 * (1 - (z_top / (self.thickness / 2)) ** 2)
                tau_oop_top_temp_yz = (Qbar_star[i].dot(gamma_oop_top))[0] ## new
                tau_oop_top_temp_xz = (Qbar_star[i].dot(gamma_oop_top))[1] ## new
                tau_top_temp = np.hstack((tau_oop_top_temp_yz, tau_oop_top_temp_xz)) ## new
                tau_oop_top[i, :] = tau_top_temp ## new
                
                # gamma_oop_bot = self.gamma_oop_avrg * 5/4 * (1 - (z_bot / (self.thickness / 2)) ** 2)
                gamma_oop_bot = self.gamma_oop_avrg * 5/6 * (1 - (z_bot / (self.thickness / 2)) ** 2)
                tau_oop_bot_temp_yz = (Qbar_star[i].dot(gamma_oop_bot))[0] ## new
                tau_oop_bot_temp_xz = (Qbar_star[i].dot(gamma_oop_bot))[1] ## new
                tau_bot_temp = np.hstack((tau_oop_bot_temp_yz, tau_oop_bot_temp_xz)) ## new
                tau_oop_bot[i, :] = tau_bot_temp ## new
                
        if (self.no_of_plies % 2 ) == 0 and skipV == 0: 
            # gamma_oop_mid = self.gamma_oop_avrg * 5/4 
            gamma_oop_mid = self.gamma_oop_avrg * 5/6 
            mid_index = int(self.no_of_plies / 2) - 1
            tau_oop_mid_temp_yz = (Qbar_star[mid_index].dot(gamma_oop_mid))[0] ## new
            tau_oop_mid_temp_xz = (Qbar_star[mid_index].dot(gamma_oop_mid))[1] ## new
            tau_mid_temp = np.hstack((tau_oop_mid_temp_yz, tau_oop_mid_temp_xz)) ## new
            tau_oop_mid = tau_mid_temp ## new
            # operation that will be performed in the latter for the top and bottom of each ply
            # are performed here because the mid value is a signle value and not a vector 
            self.tau_oop_mid = tau_oop_mid  ## new
            self.tau_oop_mid = tau_oop_mid  ## new
            T, T_star = self.__rotation_matrices(self.stackup[mid_index].get('angle')) ## new
            tau_lcl_mid = T_star.dot(self.tau_oop_mid)
            self.tau_oop_mid_ply_ref = tau_lcl_mid
        else:
            self.tau_oop_mid = 'NA. V_in not provided or the number of plies is odd. The maximum is at the top or bottom of one of the plies.'  ## new
            self.tau_oop_top_ply_ref = 'NA. V_in not provided or the number of plies is odd. The maximum is at the top or bottom of one of the plies.'

        # Interlaminar Shear Stress
        #ISS = sigma_top - sigma_bot

        # store of the stresses as class properties
        self.sigma_top = sigma_top
        self.sigma_bot = sigma_bot
        self.def_top = def_top
        self.def_bot = def_bot
        self.sigma_avrg = (sigma_top + sigma_bot) / 2
        self.def_avrg = (def_top + def_bot) / 2 
        #self.ply_ISS = ISS

        # So far the function provides the stresses in the global reference system (laminate reference system)
        # Nevertheless to verify the plies the stresses must be provided in the local reference systems and so the vectors shall be rotated       
        sigma_lcl_top = np.zeros((self.sigma_top.shape[0], self.sigma_top.shape[1]))
        sigma_lcl_bot = np.zeros((self.sigma_top.shape[0], self.sigma_top.shape[1]))
        sigma_lcl_avrg = np.zeros((self.sigma_top.shape[0], self.sigma_top.shape[1]))
        def_lcl_top = np.zeros((self.def_top.shape[0], self.def_top.shape[1]))
        def_lcl_bot = np.zeros((self.def_top.shape[0], self.def_top.shape[1]))
        def_lcl_avrg = np.zeros((self.def_top.shape[0], self.def_top.shape[1]))

        if skipV == 0:
            self.tau_oop_top = tau_oop_top  ## new
            self.tau_oop_bot = tau_oop_bot  ## new
            self.tau_oop_avrg = (tau_oop_top + tau_oop_bot) / 2  ## new
            tau_lcl_top = np.zeros((self.tau_oop_top.shape[0], self.tau_oop_top.shape[1])) ## new
            tau_lcl_bot = np.zeros((self.tau_oop_bot.shape[0], self.tau_oop_bot.shape[1])) ## new
            tau_lcl_avrg = np.zeros((self.tau_oop_bot.shape[0], self.tau_oop_bot.shape[1])) ## new
        
        # rotation of the stress vectors in the reference system of the ply 
        for i, element in enumerate(self.stackup):
            T, T_star = self.__rotation_matrices(element.get('angle')) ## new
            sigma_lcl_top[i, :] = T.dot(self.sigma_top[i, :])
            sigma_lcl_bot[i, :] = T.dot(self.sigma_bot[i, :])
            sigma_lcl_avrg[i, :] = T.dot(self.sigma_avrg[i, :])
            def_lcl_top[i, :] = T.dot(self.def_top[i, :])
            def_lcl_bot[i, :] = T.dot(self.def_bot[i, :])
            def_lcl_avrg[i, :] = T.dot(self.def_avrg[i, :])
            if skipV == 0:
                tau_lcl_top[i, :] = T_star.dot(self.tau_oop_top[i, :]) ## new
                tau_lcl_bot[i, :] = T_star.dot(self.tau_oop_bot[i, :]) ## new
                tau_lcl_avrg[i, :] = T_star.dot(self.tau_oop_avrg[i, :]) ## new

        # stresses in the local reference system
        self.sigma_bot_ply_ref = sigma_lcl_bot
        self.sigma_top_ply_ref = sigma_lcl_top
        self.sigma_avrg_ply_ref = sigma_lcl_avrg
        self.def_bot_ply_ref = def_lcl_bot
        self.def_top_ply_ref = def_lcl_top
        self.def_avrg_ply_ref = def_lcl_avrg
        
        if skipV == 0: 
            self.tau_oop_top_ply_ref = tau_lcl_top ## new
            self.tau_oop_bot_ply_ref = tau_lcl_bot ## new
            self.tau_oop_avrg_ply_ref = tau_lcl_avrg ## new
        else: 
            self.tau_oop_bot = 'NA. V_in not provided. Out-of-plane stresses not computed'
            self.tau_oop_top = 'NA. V_in not provided. Out-of-plane stresses not computed'
            self.tau_oop_avrg = 'NA. V_in not provided. Out-of-plane stresses not computed'
            self.tau_oop_bot_ply_ref = 'NA. V_in not provided. Out-of-plane stresses not computed'
            self.tau_oop_top_ply_ref = 'NA. V_in not provided. Out-of-plane stresses not computed'
            self.tau_oop_avrg_ply_ref = 'NA. V_in not provided. Out-of-plane stresses not computed'
    # add the strength properties of the laminae to the laminate properties
    def __add_max_strenght_properties(self, *varargs:list[list or np.ndarray], cntrl_shear:bool=False):
        '''
# DESCRIPTION:
    Internal method to exctract and store the strenght properties

    # INPUTS:

    ## Required
        -   *varargs: n block each containing the strenght values of the laminate. The number of block n is the number of different types of plies 
        in the laminate: how many different types of plies are present indepentenly on how many times such plies are repeated. 
        e.g. for a laminate with 10 plies of a given material and 2 plies of another one the code requires 2 blocks. 

    # OUTPUTS:
        -   None: the stresses are stored as class properties

    # Example:
        NA - internal method

### ATTENTION: this is the only internal method which performs checks because the check operation is not straight forward. For code readibility this part was removed from the main method.
    
    '''            
        # varargs = [[Xt, Xc, Yt, Yc, S12]]
        strenghts_vect = varargs[0]

        self.__diffent_plies_list_generation()

        ## CHECKS
        # check for the input to have the same lenght as the number of different plies: the code requires the strenght of each of the different plies  
        # in the laminate. In seek of clarity, if the laminate is composed by 3 blocks of different materials
        # with 10 plies each only 3 set of strenght characteristics are required 
        # because only three different constituents are present
        if len(strenghts_vect) != len(self.__different_plies_list) :
            raise Exception("The code is expecting " + str(len(self.__different_plies_list)) + " lists with 5 constants each instead " \
                + str(len(strenghts_vect)) + " were provided.")

        # input check on the strenght vectors
        for i3, strength_i in enumerate(strenghts_vect):
            # The function check that the data were provided in the correct format [name, mech_props, orientations]
            if cntrl_shear is False: 
                if len(strength_i) != 5:
                    raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. Wrong format. Vector must have five elements.')
            else: 
                if len(strength_i) != 7:
                    raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. Wrong format. Vector must have seven elements.')
            if isinstance(strength_i, list) is True: 
                for k in strength_i:
                    if isinstance(k, int) is False and isinstance(k, float) is False:
                        raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. Some of the input values are neither float nor integers.')
            elif isinstance(strength_i, np.ndarray) is True:
                for k in strength_i:
                    if k.dtype != 'int32' and k.dtype != 'float64':
                        raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. Strenght is a numpy array but its elements are neither "int32" nor "float64".')
            else:
                raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. Wrong format. Input must be either a list or a vector.') 

        for idx1, i2 in enumerate(self.plies):
            i2.update({'Xt': strenghts_vect[idx1][0]})
            i2.update({'Xc': strenghts_vect[idx1][1]})
            i2.update({'Yt': strenghts_vect[idx1][2]})
            i2.update({'Yc': strenghts_vect[idx1][3]})
            i2.update({'S12': strenghts_vect[idx1][4]})
            if cntrl_shear is True: 
                i2.update({'S13': strenghts_vect[idx1][5]})
                i2.update({'S23': strenghts_vect[idx1][6]})

        # addition of the strenght properties to the dictionaries in the "stackup_long" class property
        for i4 in self.stackup_long:
            for j, ply_name in enumerate(self.__different_plies_list):
                if i4.get('name') == ply_name:
                    i4.update({'Xt': strenghts_vect[j][0]})
                    i4.update({'Xc': strenghts_vect[j][1]})
                    i4.update({'Yt': strenghts_vect[j][2]})
                    i4.update({'Yc': strenghts_vect[j][3]})
                    i4.update({'S12': strenghts_vect[j][4]})
                    if cntrl_shear is True:
                        i4.update({'S13': strenghts_vect[j][5]})
                        i4.update({'S23': strenghts_vect[j][6]})
        
        # legend update 
        if cntrl_shear is True:
            self.stackup_long_legend.append(['Xt [MPa]', 'Xc [MPa]', 'Yt [MPa]', 'Yc [MPa]', 'S12 [MPa]', 'S13 [MPa]', 'S23 [MPa]' ])
        else:
            self.stackup_long_legend.append(['Xt [MPa]', 'Xc [MPa]', 'Yt [MPa]', 'Yc [MPa]', 'S12 [MPa]'])
    # add the strain properties of the laminae to the laminate properties
    def __add_max_strain_properties(self, *varargs:list[list or np.ndarray], cntrl_shear:bool=False):
        '''
# DESCRIPTION:
    Internal method to exctract and store the strenght properties

    # INPUTS:

    ## Required
        -   *varargs: n block each containing the strenght values of the laminate. The number of block n is the number of different types of plies 
        in the laminate: how many different types of plies are present indepentenly on how many times such plies are repeated. 
        e.g. for a laminate with 10 plies of a given material and 2 plies of another one the code requires 2 blocks. 

    # OUTPUTS:
        -   None: the stresses are stored as class properties

    # Example:
        NA - internal method

    ### ATTENTION: this is the only internal method which performs checks because the check operation is not straight forward. For code readibility this part was removed from the main method.
        
    '''            
        # varargs = [epsx_t, epsx_c, epsy_t, epsy_c, gamma_12]
        strains_vect = varargs[0]

        self.__diffent_plies_list_generation()

        ## CHECKgamma_
        # check for the input to have the same lenght as the number of different plies: the code requires the strenght of each of the different plies  
        # in the laminate. In seek of clarity, if the laminate is composed by 3 blocks of different materials
        # with 10 plies each only 3 set of strenght characteristics are required 
        # because only three different constituents are present
        if len(strains_vect) != len(self.__different_plies_list) :
            raise Exception("The code is expecting " + str(len(self.__different_plies_list)) + " lists with 5 constants each instead " \
                + str(len(strains_vect)) + " were provided.")

        # input check on the strenght vectors
        for i3, strains_i in enumerate(strains_vect):
            # The function check that the data were provided in the correct format [name, mech_props, orientations]
            if cntrl_shear is False:
                if len(strains_i) != 5:
                    raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. Wrong format. Vector must have five elements.')
            else:
                if len(strains_i) != 7:
                    raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. Wrong format. Vector must have five elements.')
            if isinstance(strains_i, list) is True: 
                for k in strains_i:
                    if isinstance(k, int) is False and isinstance(k, float) is False:
                        raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. gamma_ome of the input values are neither float nor integers.')
            elif isinstance(strains_i, np.ndarray) is True:
                for k in strains_i:
                    if k.dtype != 'int32' and k.dtype != 'float64':
                        raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. gamma_trenght is a numpy array but its elements are neither "int32" nor "float64".')
            else:
                raise Exception('Error in the strenght properties of ply no. '+ str(i3 + 1) + '. Wrong format. Input must be either a list or a vector.') 

        for idx1, i2 in enumerate(self.plies):
            i2.update({'epsx_t': strains_vect[idx1][0]})
            i2.update({'epsx_c': strains_vect[idx1][1]})
            i2.update({'epsy_t': strains_vect[idx1][2]})
            i2.update({'epsy_c': strains_vect[idx1][3]})
            i2.update({'gamma_12': strains_vect[idx1][4]})
            if cntrl_shear is True:
                i2.update({'gamma_13': strains_vect[idx1][5]})
                i2.update({'gamma_23': strains_vect[idx1][6]})

        # addition of the strenght properties to the dictionaries in the "stackup_long" class property
        for i4 in self.stackup_long:
            for j, ply_name in enumerate(self.__different_plies_list):
                if i4.get('name') == ply_name:
                    i4.update({'epsx_t': strains_vect[j][0]})
                    i4.update({'epsx_c': strains_vect[j][1]})
                    i4.update({'epsy_t': strains_vect[j][2]})
                    i4.update({'epsy_c': strains_vect[j][3]})
                    i4.update({'gamma_12': strains_vect[j][4]})
                    if cntrl_shear is True: 
                        i4.update({'gamma_13': strains_vect[j][5]})
                        i4.update({'gamma_23': strains_vect[j][6]})
        
        # legend update 
        if cntrl_shear is True:
            self.stackup_long_legend.append(['epsx_t [mm/mm]', 'epsx_c [mm/mm]', 'epsy_t [mm/mm]', 'epsy_c [mm/mm]', \
                                             'gamma_12 [mm/mm]', 'gamma_13 [mm/mm]', 'gamma_23 [mm/mm]'])
        else:
            self.stackup_long_legend.append(['epsx_t [mm/mm]', 'epsx_c [mm/mm]', 'epsy_t [mm/mm]', 'epsy_c [mm/mm]', 'gamma_12 [mm/mm]'])
# Callable class method (user accessible) 
    # compute the laminate equivalent properties        
    def calc_equivalent_properties(self, print_cntrl:bool=False, method='Barbero', disp_waring=True):
        '''
# DESCRIPTION:
    this method is used to compute the equivalent properties of the laminate if regarded as a homogeneous object. 
    Results are stored as class properties

    # INPUTS:

    ## Required
        -   None
    ## Optinal
        - print_cntrl: if True the equivalent properties are also printed. Default is true.
        - method: Barbero or ANSYS. Default is Barbero.

    # OUTPUTS: 
        -   None

    # Example:
        import PyComp as comp

        ply_name = 'Toray T300 - Epoxy 8552'
        ply_mech_props = [133.15, 16.931, 16.931, .264, .4361, .264, 5.8944, 5.7868, 5.8944, 1556.9, .275]
        ply_stkup = [0, 90, 45, -45, 45, -45, 45, 90, 0, 45]
        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)

        laminate.calc_equivalent_properties(print_cntrl=True, method='Barbero')
        '''  
        # internal variables are created from the class properties defined in __init__
        
        A = self.A
        B = self.B
        D = self.D
        H = self.H
        t = self.thickness
        a = self.a 
        b = self.b 
        h = self.h 
        d = self.d

        self.rn = self.__round((A[0, 2] / A[0, 0] ** 2) + (A[1, 2] / A[1, 1]) ** 2)
        self.rm = self.__round((D[0, 2] / D[0, 0] ** 2) + (D[1, 2] / D[1, 1]) ** 2)
        sum = 0
        for vector in B:
            for element in vector: 
                sum = sum + element ** 2
        self.rb = self.__round((3 * (A[0, 0] + A[1, 1] + A[2, 2]) * self.thickness) ** - 1 * (sum) ** .5)

        if disp_waring is True:
            if A[0, 2] > 1e-5 or A[1, 2] > 1e-5 or A[2, 0] > 1e-5 or A[2, 1] > 1e-5 : 
                print('\033[33m', 'One or more terms between A16, A26 are different from 0.',  \
    'The laminate is not orthotropic,', '\033[41m', 'rn = ', str(self.rn), '\033[49m', '. Check the resuls carefully: \
    equivalent properties are computed but the tensile-shear coupling is not considered.', '\033[37m',' ')
                print('')
            if D[0, 2] > 1e-5 or D[1, 2] > 1e-5 or D[2, 0] > 1e-5 or D[2, 1] > 1e-5 : 
                print('\033[33m','One or more terms between D16, D26 are different from 0.', '\033[33m', \
    'The laminate is not orthotropic,', '\033[41m', 'rm = ', str(self.rm), '\033[49m', '. Check the resuls carefully: \
    equivalent properties are computed but the general flexural-torsional coupling is not considered.', '\033[37m',' ')
                print('')
            if np.where(B > 1e-5)[0].shape[0] > 0: 
                print('\033[33m','One or more terms in the B are different from 0.', '\033[33m', \
    'The laminate is not orthotropic,', '\033[41m', 'rb = ', str(self.rb) , '\033[49m', '. Check the resuls carefully: \
    equivalent properties are computed but the general tensile-flexural coupling is not considered.', '\033[37m',' ')
        ## CHECKS
        if isinstance(print_cntrl, bool) is False:
            raise Exception('Error. The variable "print_cntrl" should be a boolean.')
        if isinstance(method, str) is False:
            raise Exception('Error. The variable "method" should be a string')
        if method != 'Barbero' and method != 'ANSYS':
            raise Exception('Error. Tha variable "method" should be either equal to "Barbero" or "ANSYS".')
        
        ## CALCULATIONS
        if method == 'Barbero':
            # check ref.2 sec 6.4 laminate moduli 
            self.Ex_eq = self.__round((A[0, 0] * A[1, 1] - A[0, 1] ** 2) / (t * A[1, 1]))
            self.Ey_eq = self.__round((A[0, 0] * A[1, 1] - A[0, 1] ** 2) / (t * A[0, 0]))
            self.G_eq = self.__round(A[2, 2] / t)
            self.ni_eq = self.__round(A[0, 1] / A[1, 1])

            self.G23_eq = self.__round(H[0, 0] / t) 
            self.G13_eq = self.__round(H[1, 1] / t) 
            
            self.Ex_flex_eq = self.__round(12 * (D[0, 0] * D[1, 1] - D[0, 1] ** 2) / (t ** 3 * D[1, 1]))
            self.Ey_flex_eq = self.__round(12 * (D[0, 0] * D[1, 1] - D[0, 1] ** 2) / (t ** 3 * D[0, 0]))
            self.G_flex_eq = self.__round(12 * D[2, 2] / t ** 3)
            self.ni_flex_eq = self.__round(D[0, 1] / D[1, 1])

            self.eq_props_vector = np.array([self.G_flex_eq, self.Ex_flex_eq, self.Ey_flex_eq, self.G_eq, self.Ex_eq, self.Ey_eq,self.ni_eq, self.G23_eq, \
                                            self.G13_eq, self.ni_flex_eq])
            self.eq_props_vector_legend = ['G_flex_eq', 'Ex_flex_eq', 'Ey_flex_eq', 'G_eq', 'Ex_eq', 'Ey_eq', 'ni_eq', 'G23_eq', 'G13_eq', 'ni_flex_eq']     
        elif method == 'ANSYS':
            a_star = a * t
            b_star = b * t ** 2 /2
            d_star = d * t ** 3 / 12
            h_star = h * t
            
            self.Ex_eq = 1 / a_star[0, 0]
            self.Ey_eq = 1 / a_star[1, 1]
            self.G_eq = 1 / a_star[2, 2]
            # self.ni_eq = self.__round(A[0, 1] / A[1, 1])

            self.G23_eq = h_star[0, 0] 
            self.G13_eq = h_star[1, 1] 
            
            self.Ex_flex_eq = 1 / d_star[0, 0]
            self.Ey_flex_eq = 1 / d_star[1, 1]
            self.G_flex_eq = 1 / d_star[2, 2]

            self.eq_props_vector = np.array([self.G_flex_eq, self.Ex_flex_eq, self.Ey_flex_eq, self.G_eq, self.Ex_eq, self.Ey_eq, self.G23_eq, \
                                            self.G13_eq])
            self.eq_props_vector_legend = ['G_flex_eq', 'Ex_flex_eq', 'Ey_flex_eq', 'G_eq', 'Ex_eq', 'Ey_eq', 'G23_eq', 'G13_eq']     

        if print_cntrl is True:
            print('G_eq_flex: ' + str(self.G_flex_eq) + ' MPa')
            print(' ')
            print('Ex_eq_flex: ' + str(self.Ex_flex_eq) + ' MPa')
            print(' ')
            print('Ey_eq_flex: ' + str(self.Ey_flex_eq) + ' MPa')
            print(' ')
            if method == 'Barbero':
                print('ni_eq_flex: ' + str(self.ni_flex_eq))
                print(' ')
            print('G_eq: ' + str(self.G_eq) + ' MPa')
            print(' ')
            print('Ex_eq: ' + str(self.Ex_eq) + ' MPa')
            print(' ')
            print('Ey_eq: ' + str(self.Ey_eq) + ' MPa')
            print(' ')
            print('G23_eq: ' + str(self.G23_eq))
            print(' ')
            print('G13_eq: ' + str(self.G13_eq))
            print(' ')
            if method == 'Barbero':
                print('ni_eq: ' + str(self.ni_eq))
                print(' ')
    # print the laminate stiffness matrix A, B, D, H
    def print_stifness_matrices(self):
        '''
# DESCRIPTION:
this method is used to print the stiffeness matrices A, B, D and, H

    # INPUTS:

    ## Required
        -   None
    # OUTPUTS: 
        -   None: values are printed

    # Example:
        import PyComp as comp
        ply_name = 'ANSYS Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
        laminate.print_stifness_matrices()
    '''         
        # no check needed because A, B, D and H are defined in __init__
        print('A :' + str(self.A))
        print(' ')
        print('B :' + str(self.B))
        print(' ')
        print('D :' + str(self.D))
        print(' ')
        print('H :' + str(self.H))
    # print the lamninate full stiffness matrix (still A, B, D, H but in a 6x6 or 8x8 format)
    def print_lam_stifness_matrix(self, FSDT = False):
        '''
# DESCRIPTION:
this method is used to print the assembled stiffeness matrix 

    # INPUTS:

    ## Required
        -   None

    ## optional
        -   FSDT: if the stiffness matrix has to be print in the 6x6 version (A, B, D) or in the ANSYS version (A, B, D, H)
    # OUTPUTS: 
        -   None: values are printed

    # Example:
        import PyComp as comp
        ply_name = 'ANSYS Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
        laminate.print_lam_stifness_matrix()
    '''         
        if isinstance(FSDT, bool) is False: 
            raise Exception('The variable "with_H" must be boolean (either True or False).')
        # wheather to print the stiffness matrix (A, B, D) or the stiffness matrix in the ANSYS form (A, B, D, H)
        if FSDT is False:
            print(self.stiff)
        else: 
            print(self.stiff_full)
    # compute the hygrothermal properties of the laminate
    def compute_hygrothermal_properties(self, cntrl_case:int, *varargs:list[float or int]): 
            #### WORK IN PROGRESS ####
        dummy_var = cntrl_case
    # calculate the stres state for a given loading condition
    def calculate_stress_state(self, N_in:list or np.array, M_in:list or np.array, V_in:str or list or np.array = 'None',\
        print:bool = True, print_shear:bool=False, length:int=3, \
        cntrl_external_actions:int = 0, T_in:float or int = 0, m_in:float or int = 0, \
            font_title:int=18, figwidth:int = 20, figheight:int=12, font_ax:int=15, title:str = 'None', ):
        '''
# DESCRIPTION:
Computation of the stress state in a point of the laminate under given external conditions

    # INPUTS:

    ## Required
        -   N_in: input forces per unit lenght. Three elements vector or list [Nx, Ny, Nxy].
        -   M_in: input moments per unit lenght. Three elements vector or list [Mx, My, Mxy]. 

    ## Optional
        -   V_in: input out-of-plane shera forces per unit lenght. Two elements vector or list [Vx, Vy]. 
        To be provided for the computation of out-of-plane shear stresses by default this operation is skipped. Default is "None".
        -   print: if enabled(=True) prints the laminate's stress and strain plots. Default is True.
        -   print_shear: if enabled(=True) prints the laminate's shear stress and strain plots. Default is False.
        -   length: number of points per ply in which stresses and strains are computed. Default is 20. 
        -   cntrl_external_action: Discriminate beteween four different loading cases. Depending on the case different external actions are considered.
        Case "0" only mechanical inputs. Case "1" mechanical and thermal inputs. Case "2" mechanical and hygrometric inputs. Case "3" all three. Default is 0.
        -   T_in: when the former is 1 or 3 the point temperature of the laminate is needed. Default is 0.
        -   m_in: when cntrl_external_actions is 2 or 3 the moisture absorption is needed. Default is 0.
        -   font_title: font sinze for the plots' title. Defaul is 18.
        -   figwidth: width of the output picture. Defaul is 20.
        -   figheight: height of the output picture. Defaul is 12.
        -   font_ax: font sinze for the plots' axis and labels. Defaul is 15.
    # OUTPUTS:
        -   None
    # ATTRIBUTES:
        self.sigma_top          # stress at the top of each ply
        self.sigma_bot          # stress at the bottom of each ply
        self.sigma_avrg          # average stress in the middle of each ply
        self.sigma_top_ply_ref          # stress at the top of each ply with respect to the ply reference system
        self.sigma_bot_ply_ref          # stress at the bottom of each ply with respect to the ply reference system
        self.sigma_avrg_ply_ref          # stress in the middle of each ply with respect to the ply reference system
        self.sigma_out          # through-the-thickness stress vectors (x, y, xy)
        self.tau_oop_top          # out-of-plane shear  stress at the top of each ply
        self.tau_oop_bot          # out-of-plane shear  stress at the bottom of each ply
        self.tau_oop_avrg            # average out-of-plane shear  stress in the middle of each ply
        self.tau_oop_mid            # out-of-plane shear  stress in the middle of the laminate
        self.tau_oop_top_ply_ref          # out-of-plane shear  stress at the top of each ply with respect to the ply reference system
        self.tau_oop_bot_ply_ref          # out-of-plane shear  stress at the bottom of each ply with respect to the ply reference system
        self.tau_oop_avrg_ply_ref          # average shear stress in the middle of each ply with respect to the ply reference system
        self.tau_oop_out          # through-the-thickness out-of-plane shear stress vectors (yz, xz)
        self.def_top          # deformation at the top of each ply
        self.def_bot          # deformation at the bottom of each ply
        self.def_average          # deformation in the middle of each ply
        self.def_top_ply_ref          # deformation at the top of each ply with respect to the ply reference system
        self.def_bot_ply_ref          # deformation at the bottom of each ply with respect to the ply reference system
        self.def_average_ply_ref          # averagedeformation in the middle of each ply with respect to the ply reference system
        self.def_out          # through-the-thickness deformation vectors (x, y, xy)
        self.def_out_percent          # percentual through-the-thickness deformation vectors (x, y, xy)
        
        self.epsilon0           # laminate mid plane deformation vector
        self.k0         # laminate mid plane curvature vector
        self.gamma_oop_avrg          # out-of-plane average shear deformation 
        self.gamma_out          # through-the-thickness out-of plane shear deformation vectors (yz, xz)
        self.gamma_out_percent          # percentual through-the-thickness out-of plane shear deformation vectors (yz, xz) 


    # Example:
        import PyComp as comp
        ply_name = 'Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
        N = [230, .10, -2.5]
        M = [-160, .012, -0.3]
        V = [.0005, 3.5]
        laminate.calculate_stress_state(N, M, V, print=True, print_shear=True)
    '''           
        if isinstance(N_in, list) is True:
            for i in N_in:
                if isinstance(i, float) is False and isinstance(i, int) is False or isinstance(i, bool) is True:
                    raise Exception('Error. Some elements in N_in are neither floats nor integers.')
        elif isinstance(N_in, np.ndarray) is True:
            for i in N_in:
                if i.dtype != 'int32' and i.dtype != 'float64':
                    raise Exception('Error. N_in is a numpy array but its elements are neither "int32" nor "float64".')
        else:
            raise Exception('Error. The input N_in must be either a list or a numpy array.')
        if isinstance(M_in, list) is True:
            for i in M_in:
                if isinstance(i, float) is False and isinstance(i, int) is False or isinstance(i, bool) is True:
                    raise Exception('Error. Some elements in M_in are neither floats nor integers.')
        elif isinstance(M_in, np.ndarray) is True:
            for i in M_in:
                if i.dtype != 'int32' and i.dtype != 'float64':
                    raise Exception('Error. M_in is a numpy array but its elements are neither "int32" nor "float64".')
        else:
            raise Exception('Error. The input M_in must be either a list or a numpy array.')
        if len(N_in) != 3:
                raise Exception('Error. N must have 3 elements')
        if len(M_in) != 3:
                raise Exception('Error. M must have 3 elements')      
        if isinstance(V_in, str) is False:
            calc_shear = True
            if isinstance(V_in, list) is True:
                for i in V_in:
                    if isinstance(i, float) is False and isinstance(i, int) is False or isinstance(i, bool) is True:
                        raise Exception('Error. Some elements in V_in are neither floats nor integers.')
            elif isinstance(V_in, np.ndarray) is True:
                for i in V_in:
                    if i.dtype != 'int32' and i.dtype != 'float64':
                        raise Exception('Error. V_in is a numpy array but its elements are neither "int32" nor "float64".')
            else:
                raise Exception('Error. The input V_in must be either a list or a numpy array.')       
            
            if len(V_in) != 2:
                    raise Exception('Error. V must have 2 elements')
        else:
            calc_shear = False
            if V_in != 'None':
                raise Exception('Error. The input V_in must be either a list or a numpy array.') 
        if isinstance(cntrl_external_actions, int) is False or isinstance(cntrl_external_actions, bool) is True:
            raise Exception('Error. The variable cntrl_external_actions must be an integer.')
        if cntrl_external_actions != 0 and cntrl_external_actions != 1 and cntrl_external_actions != 2 and cntrl_external_actions != 3:
            raise Exception('Error. The variable "cntrl_external_actions" must be either equal to 0, 1, 2 or 3.')
        if isinstance(T_in, float) is False and isinstance(T_in, int) is False or isinstance(T_in, bool) is True:
            raise Exception('Error. The variable T is neither a float nor an integer.')
        if isinstance(m_in, float) is False and isinstance(m_in, int) is False or isinstance(m_in, bool) is True:
            raise Exception('Error. The variable m is neither a float nor an integer.')
        # exctraction of the stresses (saved as class properties)
        if cntrl_external_actions == 0: 
            self.__ply_stresses(N_in, M_in, V_in) 
        elif cntrl_external_actions == 1: 
            raise Exception('this portion of the code is still under development and will be available in the next releases')
            # self.__ply_stresses(N_in, M_in, V_in, T_in=T_in, cntrl_external_actions=1) 
        elif cntrl_external_actions == 2: 
            raise Exception('this portion of the code is still under development and will be available in the next releases')
            # self.__ply_stresses(N_in, M_in, V_in, m_in=m_in, cntrl_external_actions=2) 
        else: 
            raise Exception('this portion of the code is still under development and will be available in the next releases')
            # self.__ply_stresses(N_in, M_in, V_in, T_in=T_in, m_in=m_in, cntrl_external_actions=3) 

        # internal variables
        sigma_top = self.sigma_top
        sigma_bot = self.sigma_bot
        # font size in the printed graphs
        font1 = font_title
        font2 = font_ax

        ## VARIABLES INITIALIZATION
        # number of points in the ply: so far only the top and bottom values of stress in the ply where computed
        # to have a better plot, the stress at some internal points is calculated
        # the number of such points is length 
        length = length
        # sigma_temp is an internal variable used in the second loop
        sigma_temp = np.zeros((length, len(self.thickness_vector)))
        # sigma is the varible to plot. It has nx3 dimensions where n is the number of plies time the variable lenght(points in each ply)
        sigma = np.zeros((length * len(self.thickness_vector), 3))
        def_temp = np.zeros((length * len(self.thickness_vector), 3))
        zed = np.zeros((length * len(self.thickness_vector)))
        ply_top_coord = np.zeros((len(self.thickness_vector)))

        # loop over the plies  
        for i in range (3):
            for j in range(len(self.stackup)) :
                
                # overall z vector (abscissa of the graph) - created only the first time beacause the geometry does not change
                # generated only once
                if i == 0:
                    zed[length * j:length * j + length] = np.linspace(0, self.thickness_vector[j], length) + np.sum(self.thickness_vector[:j])   
                
                # z vector. Used to compute the stressess at the internal points
                z = np.linspace(0, self.thickness_vector[j], length)
                # thickness of the j-th ply
                z1 = self.thickness_vector[j]
                # stress in the j-th ply (top and bottom) along the i-th direction
                sigma0 = sigma_bot[j, i]
                sigma1 = sigma_top[j, i]

                # through the thickness stress profile (linear in the ply)
                sigma_temp = sigma0 + (sigma1 - sigma0) / z1 * z

                # descibed in the initialization. Stress in n (=length) point of the ply
                sigma[length * j:length * j + length, i] = sigma_temp
                # top abscissa of each ply
                ply_top_coord[j] = sum(self.thickness_vector[:j + 1]) 

            
            # laminate point deformation
            deformation = self.epsilon0[i] + self.k0[i] * (zed - self.thickness / 2)
            def_temp[:, i] = deformation
            # storage of sigma values into class properties
            if i == 0:
                self.sigmax = sigma[:, i] 
                self.zed = zed
            elif i == 1:
                self.sigmay = sigma[:, i]
            elif i == 2:
                self.sigmaxy = sigma[:, i]
            else: 
                raise Exception('i assumed a forbidden value')

        # accessible only if print is "True"
        self.sigma_out = sigma
        self.zed_out = zed
        self.def_out = def_temp
        self.def_out_percent = def_temp * 100

        if print is True: ## new

            # loop over the three directions (sigma1, sigma2, and shear)
            for i in range(3) :
                
                deformation_precent = def_temp[:, i] * 100

            # graph definition (2 graphs: one for the deformation and another for the stresses) 
                fig, ax = plt.subplots(1, 2)
                if title != 'None':
                    fig.suptitle(title, fontsize = font_title)
                ax[0].axhline(0, color = 'red')
                ax[1].axhline(0, color = 'red')

                for kk in range(len(self.stackup)):
                    # horizontal lines to indicate the plies  
                    ax[0].axhline(ply_top_coord[kk], color = 'red')
                    ax[1].axhline(ply_top_coord[kk], color = 'red')
                    #x_coord0 = min(deformation) - max(deformation) / 100
                    #x_coord1 = min(sigma[:, i]) - max(sigma[:, i]) / 100

                # plots 
                ax[1].plot(sigma[:, i], zed)
                ax[1].set_title('sigma' + str(i + 1), fontsize = font1) 
                ax[0].plot(deformation_precent, zed)
                ax[0].set_title('epsilon' + str(i + 1), fontsize = font1)
                ax[1].set_xlabel('Stress [MPa]', fontsize = font2)
                #ax[1].set_ylabel('z [mm]')
                ax[0].set_xlabel('Deformation [%]', fontsize = font2)
                ax[0].set_ylabel('z [mm]', fontsize = font2)
                ax[0].tick_params(axis = 'x', labelsize = font2)
                ax[0].tick_params(axis = 'y', which = 'both', labelsize = font2)
                ax[1].tick_params(axis = 'x', labelsize = font2)
                ax[1].tick_params(axis = 'y', which = 'both', labelsize = font2)
                fig.set_figwidth(figwidth)
                fig.set_figheight(figheight)

                delta_def = max(deformation_precent) - min(deformation_precent)
                mean_def = (max(deformation_precent) + min(deformation_precent)) / 2
                delta_sigma =  max(sigma[:, i]) - min(sigma[:, i])
                mean_sigma = (max(sigma[:, i]) + min(sigma[:, i])) / 2
                ax[0].axvline(mean_def, color = 'black', ls = '--', lw = .5)
                ax[1].axvline(mean_sigma, color = 'black', ls = '--', lw = .5)

                # Definition of boundaries for plots
                # Check if there is an appereciable difference between max and min
                if  delta_def <= .0004:
                    # if not set the boundaries as +/- .04 about the max value (which is equal or almost equal to the min value)
                    boundary0_min = -.04 + max(deformation_precent) 
                    boundary0_max = .04 + max(deformation_precent)
                else: 
                    # if yes check the sign of the maximum
                    if max(deformation_precent) <= 0: 
                        # if it is negative take a 10% margin by multipling it by .9 
                        boundary0_max = max(deformation_precent) * .9
                    else:
                        # else multiply it by 1.1 to take the same margin
                        boundary0_max = max(deformation_precent) * 1.1
                    # to have a sym graph take the interval from before 
                    increment = abs(boundary0_max - max(deformation_precent))
                    boundary0_min = min(deformation_precent) - increment

                if delta_sigma <= .1:
                    boundary1_min = - 5 + max(sigma[:, i]) 
                    boundary1_max = 5 + max(sigma[:, i]) 
                else: 
                    if max(sigma[:, i]) <= 0:
                        boundary1_max = max(sigma[:, i]) * .9
                    else:
                        boundary1_max = max(sigma[:, i]) * 1.1
                    increment2 = abs(boundary1_max - max(sigma[:, i]))
                    boundary1_min = min(sigma[:, i]) - increment2

                max_x_0 = boundary0_max
                min_x_0 = boundary0_min
                max_x_1 = boundary1_max
                min_x_1 = boundary1_min

                ax[0].set_xlim(min_x_0, max_x_0)
                ax[1].set_xlim(min_x_1, max_x_1)
                x_coord0 = max_x_0 + (max_x_0 - min_x_0) * .01
                #print('x_coord0 ' + str(x_coord0))
                x_coord1 = max_x_1 + (max_x_1 - min_x_1) * .01 
                #print('x_coord1 ' + str(x_coord1))
                for j in range(len(self.stackup)) :
                    #ply_top_coord = sum(self.thickness_vector[:j + 1]) 
                    y_coord = ply_top_coord[j] - self.thickness_vector[j] / 2
                    ax[0].text(x_coord0, y_coord, 'ply ' + str(j + 1) + ', ' + str(self.stackup[j].get("angle")) + '°', fontsize=10)
                    ax[1].text(x_coord1, y_coord, 'ply ' + str(j + 1) + ', ' + str(self.stackup[j].get("angle")) + '°', fontsize=10)

        if calc_shear is True:
            ## VARIABLES INITIALIZATION
            # number of points in the ply: so far only the top and bottom values of stress in the ply where computed
            # to have a better plot, the stress at some internal points is calculated
            # the number of such points is length 
            length = length
            zed = np.zeros((length * len(self.thickness_vector)))
            Qbar_star = np.array(self.Qbars_star) ## new
            tau = np.zeros((length * len(self.thickness_vector), 2))
            gamma_out = np.zeros((length * len(self.thickness_vector), 2))
            
            for i in range(2) :

                for j in range(len(self.stackup)) :
                    
                    # overall z vector (abscissa of the graph) - created only the first time beacause the geometry does not change
                    # generated only once
                    if i == 0:
                        zed[length * j:length * j + length] = np.linspace(0, self.thickness_vector[j], length) + np.sum(self.thickness_vector[:j])   
                    
                    z2 = zed[length * j:length * j + length] - self.thickness / 2
                    
                    tau[length * j:length * j + length, i] = (Qbar_star[j].dot(self.gamma_oop_avrg))[i] * (5 / 4) * (1 - (z2 / (self.thickness / 2)) ** 2) 


                # laminate point gamma
                z3 = zed - self.thickness / 2 
                gamma = self.gamma_oop_avrg[i] * (5 / 4) * (1 - (z3 / (self.thickness / 2)) ** 2) 
                gamma_precent = gamma * 100
                gamma_out[:, i] = gamma


            self.tauxz = tau[:, 0]
            self.zed_ool = zed
            self.tauyz = tau[:, 1]
            self.tau_out = tau
            self.gamma_out = gamma_out
            self.gamma_out_percent = gamma_out * 100


        if print_shear is True: ## new
            # font size in the printed graphs
            font1 = font_title
            font2 = font_ax

            # loop over the three directions (tau1, tau2, and shear)
            for i in range(2) :
    
            # graph definition (2 graphs: one for the gamma and another for the stresses) 
                fig, ax = plt.subplots(1, 2)
                if title != 'None':
                    fig.suptitle(title, fontsize = font_title)
                ax[0].axhline(0, color = 'red')
                ax[1].axhline(0, color = 'red')

                # loop over the plies  
                for j in range(len(self.stackup)) :
                                        
                    # top abscissa of each ply
                    ply_top_coord = sum(self.thickness_vector[:j + 1]) 
                    # horizontal lines to indicate the plies  
                    ax[0].axhline(ply_top_coord, color = 'red')
                    ax[1].axhline(ply_top_coord, color = 'red')

                gamma_precent = gamma_out[:, i] * 100

                # plots 
                if i == 0:
                    ax[1].set_title('tau xz', fontsize = font1) 
                    ax[0].set_title('gamma xz', fontsize = font1)
                elif i == 1:
                    ax[1].set_title('tau yz', fontsize = font1) 
                    ax[0].set_title('gamma yz', fontsize = font1)
                ax[1].plot(tau[:, i], zed)
                ax[0].plot(gamma_precent, zed)
                ax[1].set_xlabel('Stress [MPa]', fontsize = font2)
                #ax[1].set_ylabel('z [mm]')
                ax[0].set_xlabel('Deformation [%]', fontsize = font2)
                ax[0].set_ylabel('z [mm]', fontsize = font2)
                ax[0].tick_params(axis = 'x', labelsize = font2)
                ax[0].tick_params(axis = 'y', which = 'both', labelsize = font2)
                ax[1].tick_params(axis = 'x', labelsize = font2)
                ax[1].tick_params(axis = 'y', which = 'both', labelsize = font2)
                fig.set_figwidth(figwidth)
                fig.set_figheight(figheight)

                delta_def = max(gamma_precent) - min(gamma_precent)
                mean_def = (max(gamma_precent) + min(gamma_precent)) / 2
                delta_tau =  max(tau[:, i]) - min(tau[:, i])
                mean_tau = (max(tau[:, i]) + min(tau[:, i])) / 2
                ax[0].axvline(mean_def, color = 'black', ls = '--', lw = .5)
                ax[1].axvline(mean_tau, color = 'black', ls = '--', lw = .5)

                # Definition of boundaries for plots
                # Check if there is an appereciable difference between max and min
                if  delta_def <= .0004:
                    # if not set the boundaries as +/- .04 about the max value (which is equal or almost equal to the min value)
                    boundary0_min = -.04 + max(gamma_precent) 
                    boundary0_max = .04 + max(gamma_precent)
                else: 
                    # if yes check the sign of the maximum
                    if max(gamma_precent) <= 0: 
                        # if it is negative take a 10% margin by multipling it by .9 
                        boundary0_max = max(gamma_precent) + max(abs(gamma_precent)) * .1
                    else:
                        # else multiply it by 1.1 to take the same margin
                        boundary0_max = max(gamma_precent) * 1.1
                    # to have a sym graph take the interval from before 
                    increment = abs(boundary0_max - max(gamma_precent))
                    boundary0_min = min(gamma_precent) - increment

                if delta_tau <= .1:
                    if max(tau[:, i]) >= 0:
                        boundary1_min = - 5 + max(tau[:, i]) 
                        boundary1_max = 5 + max(tau[:, i]) 
                    # else:
                    #     boundary1_min = - 5 + max(tau) 
                    #     boundary1_max = 5 + max(tau) 
                else: 
                    boundary1_max = max(tau[:, i]) + max(abs(tau[:, i])) * .1
                    increment2 = abs(boundary1_max - max(tau[:, i]))
                    boundary1_min = min(tau[:, i]) - increment2 

                max_x_0 = boundary0_max 
                min_x_0 = boundary0_min
                max_x_1 = boundary1_max
                min_x_1 = boundary1_min

                ax[0].set_xlim(min_x_0, max_x_0)
                ax[1].set_xlim(min_x_1, max_x_1)
                x_coord0 = max_x_0 + (max_x_0 - min_x_0) * .01
                #print('x_coord0 ' + str(x_coord0))
                x_coord1 = max_x_1 + (max_x_1 - min_x_1) * .01 
                #print('x_coord1 ' + str(x_coord1))
                for j in range(len(self.stackup)) :
                    ply_top_coord = sum(self.thickness_vector[:j + 1]) 
                    y_coord = ply_top_coord - self.thickness_vector[j] / 2
                    ax[0].text(x_coord0, y_coord, 'ply ' + str(j + 1) + ', ' + str(self.stackup[j].get("angle")) + '°', fontsize=10)
                    ax[1].text(x_coord1, y_coord, 'ply ' + str(j + 1) + ', ' + str(self.stackup[j].get("angle")) + '°', fontsize=10)

                        # accessible only if print is "True"
    # perform the FPF structural verification of the laminate
    def FPF(self, *varargs:list[list or np.ndarray], N_in:list or np.array, M_in:list or np.array, \
        print_stress_state:bool = False, V_in:str or np.array='None', criteria:str='TsaiWu', F12:float or int="-.5/FxtFxcFytFyc_max**2", \
            print_margins:bool = False, cntrl_external_actions:int = 0, T_in:int or float = 0, m_in:int or float = 0):
        '''
# DESCRIPTION:
this method is used verify the stackup and to identify the first failing ply in the laminate. 

    # INPUTS:

    ## Required
        -   *varargs: n block each containing the strenght values of the laminae. The number of block n is the number of different types of plies 
        in the laminate: how many different types of plies are present indepentenly on how many times such plies are repeated. 
        e.g. for a laminate with 10 plies of a given material and 2 plies of another one the code requires 2 blocks 

    ### ATTENTION: the blocks containing the strenght properties should be oredered as the referred materials (in order of apperance from bottom to top). 
        
        -   N_in: input forces per unit lenght. Three elements vector or list [Nx, Ny, Nxy].
        -   M_in: input moments per unit lenght. Three elements vector or list [Mx, My, Mxy]. 

    ## Optional
        -   V_in: input out-of-plane shear forces per unit lenght. Three elements vector or list [Vx, Vy]. 
        -   print_stress_state: print the stress distribution in the thickness.
        -   criteria: criteria used to verify the plies. Either 'TsaiWu', 'MaxStress' or 'MaxStrain'. Default is 'TsaiWu'.
        -   F12: coefficient of the Tsai Wu formula. Default is 0. 
        -   print_margins: whether to print the margin for each ply (if True prints). When equal or less than 0 the ply fails. Default is False. 
        -   cntrl_external_actions: wheather the thermal and/or the humidity effects should be included in the calculation. 
        (0 only mech, 1 mech and thermal, 2 mech and hygro, 3 mech and hygrothermal)

    ### ATTENTION: the "print margins" feature is available only for the TsaiWu verification

    # OUTPUTS: 
        -   None

    # ATTRIBUTES:
        self.margin_top         # margin at the top of each ply
        self.reserve_factor_top         # reserve factor at the top of each ply  
        self.inv_reserve_factor_top         # inverse reserve factor at the top of each ply 
        self.margin_bot         # margin at the bottom of each ply 
        self.reserve_factor_bot         # reserve factor at the bottom of each ply  
        self.inv_reserve_factor_bot         # inverse reserve factor at the bottom of each ply  
        self.margin_average         # margin in the middle of each ply
        self.average_reserve_factor         # reserve factor in the middle of each ply  
        self.average_inv_reserve_factor         # inverse reserve factor in the middle of each ply 

    # Example:
        import PyComp as comp
        ply_name = 'ANSYS Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
        N = [230, .10, -2.5]
        M = [-160, .012, -0.3]
        V = [.0005, 3.5]
        strght = [414, 414 * .5, 414, 414 * .5, 81.41, 40, 40]
        strain = [1/100, .5/100, 1/100, .5/100, 5/100, 5/100, 5/100]
        laminate.FPF([strght], N_in=N, M_in=M, V_in=V, criteria='TsaiWu')
    '''         


        ## CHECKS
        if criteria != 'TsaiWu' and criteria != 'MaxStress' and criteria != 'MaxStrain':
            raise ValueError('The variable criteria must be either equal to "MaxStress" or to "TsaiWu"')
        if isinstance(N_in, list) is True:
            for i in N_in:
                if isinstance(i, float) is False and isinstance(i, int) is False:
                    raise Exception('Error. Some elements in N_in are neither floats nor integers.')
        elif isinstance(N_in, np.ndarray) is True:
            for i in N_in:
                if i.dtype != 'int32' and i.dtype != 'float64':
                    raise Exception('Error. N_in is a numpy array but its elements are neither "int32" nor "float64".')
        else:
            raise Exception('Error. The input N_in must be either a list or a numpy array.')
        if isinstance(M_in, list) is True:
            for i in M_in:
                if isinstance(i, float) is False and isinstance(i, int) is False:
                    raise Exception('Error. Some elements in M_in are neither floats nor integers.')
        elif isinstance(M_in, np.ndarray) is True:
            for i in M_in:
                if i.dtype != 'int32' and i.dtype != 'float64':
                    raise Exception('Error. M_in is a numpy array but its elements are neither "int32" nor "float64".')
        else:
            raise Exception('Error. The input M_in must be either a list or a numpy array.')
        if len(N_in) != 3:
                raise Exception('Error. N must have 3 elements')
        if len(M_in) != 3:
                raise Exception('Error. M must have 3 elements')      
        if isinstance(V_in, str) is False:
            if isinstance(V_in, list) is True:
                for i in V_in:
                    if isinstance(i, float) is False and isinstance(i, int) is False:
                        raise Exception('Error. Some elements in V_in are neither floats nor integers.')
            elif isinstance(V_in, np.ndarray) is True:
                for i in V_in:
                    if i.dtype != 'int32' and i.dtype != 'float64':
                        raise Exception('Error. V_in is a numpy array but its elements are neither "int32" nor "float64".')
            else:
                raise Exception('Error. The input V_in must be either a list or a numpy array.')       
            
            if len(V_in) != 2:
                    raise Exception('Error. V must have 2 elements')
        else:
            if V_in != 'None':
                raise Exception('Error. The input V_in must be either a list or a numpy array.')       
        if F12 != "-.5/FxtFxcFytFyc_max**2":
            if isinstance(F12, float) is False and isinstance(F12, int) is False:
                raise Exception('Error. The input F12 must be either a float or an integer.')
        if isinstance(print_margins, bool) is False:
            raise Exception('Error. The input print_margins must be boolean.')
        if isinstance(cntrl_external_actions, int) is False:
            raise Exception('Error. The variable cntrl_external_actions must be an integer.')
        if cntrl_external_actions != 0 and cntrl_external_actions != 1 and cntrl_external_actions != 2 and cntrl_external_actions != 3:
            raise Exception('Error. The variable "cntrl_external_actions" must be either equal to 0, 1, 2 or 3.')
        if isinstance(T_in, float) is False and isinstance(T_in, int) is False:
            raise Exception('Error. The variable T is neither a float nor an integer.')
        if isinstance(m_in, float) is False and isinstance(m_in, int) is False:
            raise Exception('Error. The variable m is neither a float nor an integer.')

        ## Stresses calculation 
        if isinstance(V_in, str) is False: 
            print_stress_state2 = print_stress_state
            cntrl_shear = True
        else:
            print_stress_state2 = False
            cntrl_shear = False
        
        if criteria == 'TsaiWu' or criteria == 'MaxStress':
            # internal method which check, extract and store the plies' max strenght data
            self.__add_max_strenght_properties(varargs[0], cntrl_shear=cntrl_shear)
            # user warnings
            print('\033[35m', 'Input strenghts are assumed to be provided in MPa', '\033[37m',' ')
            print('\033[35m', 'N_in is assumed to be provided in N/mm', '\033[37m',' ')
            print('\033[35m', 'M_in is assumed to be provided in N', '\033[37m',' ')
        else:
            # internal method which check, extract and store the plies' max strain data
            self.__add_max_strain_properties(varargs[0], cntrl_shear=cntrl_shear)
            # user warnings
            print('\033[35m', 'Input strains are assumed to be provided in mm/mm', '\033[37m',' ')
            print('\033[35m', 'N_in is assumed to be provided in N/mm', '\033[37m',' ')
            print('\033[35m', 'M_in is assumed to be provided in N', '\033[37m',' ')

        self.criteria = criteria

        self.calculate_stress_state(N_in, M_in, V_in=V_in, print=print_stress_state, print_shear=print_stress_state2, \
            cntrl_external_actions=cntrl_external_actions, T_in=T_in, m_in=m_in)

        sigma_lcl_bot = self.sigma_bot_ply_ref
        sigma_lcl_top = self.sigma_top_ply_ref
        sigma_lcl_avrg = self.sigma_avrg_ply_ref
        def_lcl_bot = self.def_bot_ply_ref
        def_lcl_top = self.def_top_ply_ref
        def_lcl_avrg = self.def_avrg_ply_ref
        tau_lcl_bot = self.tau_oop_bot_ply_ref
        tau_lcl_top = self.tau_oop_top_ply_ref
        tau_lcl_avrg = self.tau_oop_avrg_ply_ref

        ## PLY-BY-PLY verification 
        ###########################################################
        # Maximum stress criteria
        if criteria == 'MaxStress':
            # margin vectors initialization 
            ## ATTENTION: since this criteria DOES NOT consider the multiaxial state of stress there is MORE THAN one margin value per ply
            margin_top = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            margin_bot = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            inv_reserve_factor_top = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            inv_reserve_factor_bot = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            reserve_factor_top = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            reserve_factor_bot = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))

            # for each ply
            for i in range (len(self.stackup)):
                # exctract the strenght properties
                Xt = self.stackup_long[i].get('Xt')
                Xc = self.stackup_long[i].get('Xc')
                Yt = self.stackup_long[i].get('Yt')
                Yc = self.stackup_long[i].get('Yc')
                S = self.stackup_long[i].get('S12')
                
                # it has to be verified if the ply is either in tension or compression 
                # to use the correct check value
                if sigma_lcl_top[i, 0] > 0: 
                    X = Xt    
                else:
                    X = Xc
                # same sign check 
                if sigma_lcl_top[i, 1] > 0:
                    Y = Yt
                else: 
                    Y = Yc

                # shear failure is independent on the sign therfore the values are turned positive before proceeding 
                if sigma_lcl_top[i, 2] < 0:
                    sigma_lcl_top[i, 2] = - sigma_lcl_top[i, 2] 
                
                #margin_top[i, :] = abs(np.array([X, Y, S])) - abs(sigma_lcl_top[i, :])
                # the margin is defined as the ratio of the stress and the maximum stress. It stats from 1 and goes down.
                # If the margin is less or equal than 0 the ply fails.

                inv_reserve_factor_top[i, :] = abs(sigma_lcl_top[i, :]) / abs(np.array([X, Y, S]))

                if sigma_lcl_bot[i, 0] > 0: 
                    X = Xt    
                else :
                    X = Xc
                if sigma_lcl_bot[i, 1] > 0:
                    Y = Yt
                else : 
                    Y = Yc

                # shear failure is independent on the sign therfore the values are turned positive before proceeding 
                if sigma_lcl_bot[i, 2] < 0:
                    sigma_lcl_bot[i, 2] = - sigma_lcl_bot[i, 2] 

                inv_reserve_factor_bot[i, :] = abs(sigma_lcl_bot[i, :]) / abs(np.array([X, Y, S]))
                
            reserve_factor_bot = inv_reserve_factor_bot ** - 1   
            reserve_factor_top = inv_reserve_factor_top ** - 1    
            # store the margins as class variables
            self.reserve_factor_bot = reserve_factor_bot
            self.reserve_factor_top = reserve_factor_top
            self.inv_reserve_factor_bot = inv_reserve_factor_bot
            self.inv_reserve_factor_top = inv_reserve_factor_top
            self.margin_top = reserve_factor_top - 1
            self.margin_bot = reserve_factor_bot - 1

            # perform the check on the margin bot
            # create a direction vector for the possible warning
            direction = ['x', 'y', 'xy']
            # check if "margin_bot" has any value below or at zero
            check_bot = np.where(self.margin_bot <= 0)
            # same operation for the top
            check_top = np.where(self.margin_top <= 0)
            # if there are such values
            print(' ')
            if len(np.where(self.margin_bot < 0)[0]) != 0 or len(np.where(self.margin_top < 0)[0]) != 0:
                # for each of those plies display the following warning
                for i in range (int(len(check_bot[0]))):
                    print('\033[31m', 'Failure expected on the bottom side of ply no. ' + str(check_bot[0][i] + 1) + \
                        ' along the ' + direction[check_bot[1][i]] + ' direction ' + 'according to the Maximum stress criteria', '\033[37m',' ')          
                for i in range (int(len(check_top[0]))):
                    print('\033[31m', 'Failure expected on the top side of ply no. ' + str(check_top[0][i] + 1) + \
                        ' along the ' + direction[check_top[1][i]] + ' direction ' + 'according to the Maximum stress criteria', '\033[37m',' ')
            else: 
                print('\033[32m', 'No failure expected in the analysed laminate', '\033[37m', ' ')

            if print_margins is True:
                print(' ')
                print('\033[33m', 'The option "print_margins" is available only for the TsaiWu FPF.')
        ###########################################################
        # Maximum strain criteria
        if criteria == 'MaxStrain':
            # margin vectors initialization 
            ## ATTENTION: since this criteria DOES NOT consider the multiaxial state of stress there is MORE THAN one margin value per ply
            margin_top = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            margin_bot = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            inv_reserve_factor_top = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            inv_reserve_factor_bot = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            reserve_factor_top = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))
            reserve_factor_bot = np.zeros((sigma_lcl_top.shape[0], sigma_lcl_top.shape[1]))

            # for each ply
            for i in range (len(self.stackup)):
                # exctract the strenght properties
                epsx_t = self.stackup_long[i].get('epsx_t')
                epsx_c = self.stackup_long[i].get('epsx_c')
                epsy_t = self.stackup_long[i].get('epsy_t')
                epsy_c = self.stackup_long[i].get('epsy_c')
                gamma_12 = self.stackup_long[i].get('gamma_12')
            
                # it has to be verified if the ply is either in tension or compression 
                # to use the correct check value
                if def_lcl_top[i, 0] > 0: 
                    epsx = epsx_t    
                else :
                    epsx = epsx_c
                # same sign check 
                if def_lcl_top[i, 1] > 0:
                    epsy = epsy_t    
                else : 
                    epsy = epsy_c  

                # shear failure is independent on the sign therfore the values are turned positive before proceeding 
                if def_lcl_top[i, 2] < 0:
                    def_lcl_top[i, 2] = - def_lcl_top[i, 2] 
                
                #margin_top[i, :] = abs(np.array([X, Y, S])) - abs(sigma_lcl_top[i, :])
                # the margin is defined as the ratio of the stress and the maximum stress. It stats from 1 and goes down.
                # If the margin is less or equal than 0 the ply fails.
                inv_reserve_factor_top[i, :] = abs(def_lcl_top[i, :]) / abs(np.array([epsx, epsy, gamma_12]))

                if def_lcl_bot[i, 0] > 0: 
                    epsx = epsx_t    
                else :
                    epsx = epsx_c
                if def_lcl_bot[i, 1] > 0:
                    epsy = epsy_t    
                else : 
                    epsy = epsy_c  

                # shear failure is independent on the sign therfore the values are turned positive before proceeding 
                if def_lcl_bot[i, 2] < 0:
                    def_lcl_bot[i, 2] = - def_lcl_bot[i, 2] 

                inv_reserve_factor_bot[i, :] = abs(def_lcl_bot[i, :]) / abs(np.array([epsx, epsy, gamma_12]))
            
            # store the margins as class variables
            self.inv_reserve_factor_top = inv_reserve_factor_top
            self.inv_reserve_factor_bot = inv_reserve_factor_bot
            self.reserve_factor_top = inv_reserve_factor_top ** -1
            self.reserve_factor_bot = inv_reserve_factor_bot ** -1
            self.margin_top = inv_reserve_factor_top ** -1 - 1
            self.margin_bot = inv_reserve_factor_bot ** -1 - 1

            # perform the check on the margin bot
            # create a direction vector for the possible warning
            direction = ['x', 'y', 'xy']
            # check if "margin_bot" has any value below or at zero
            check_bot = np.where(self.margin_bot <= 0)
            # same operation for the top
            check_top = np.where(self.margin_top <= 0)
            # if there are such values
            print(' ')
            if len(np.where(self.margin_bot < 0)[0]) != 0 or len(np.where(self.margin_top < 0)[0]) != 0:
                # for each of those plies display the following warning
                for i in range (int(len(check_bot[0]))):
                    print('\033[31m', 'Failure expected on the bottom side of ply no. ' + str(check_bot[0][i] + 1) + \
                        ' along the ' + direction[check_bot[1][i]] + ' direction ' + 'according to the Maximum strain criteria', '\033[37m',' ')          
                for i in range (int(len(check_top[0]))):
                    print('\033[31m', 'Failure expected on the top side of ply no. ' + str(check_top[0][i] + 1) + \
                        ' along the ' + direction[check_top[1][i]] + ' direction ' + 'according to the Maximum strain criteria', '\033[37m',' ')
            else: 
                print('\033[32m', 'No failure expected in the analysed laminate', '\033[37m', ' ')

            if print_margins is True:
                print(' ')
                print('\033[33m', 'The option "print_margins" is available only for the TsaiWu FPF.')   
        ###########################################################
        # Tsai-Wu criteria
        if criteria == 'TsaiWu':
            # initialization of the margin variables
            ## ATTENTION: since this criteria considers the multiaxial state of stress there is only one margin value per ply
            margin_top = np.zeros((sigma_lcl_top.shape[0]))
            reserve_factor_top = np.zeros((sigma_lcl_top.shape[0]))
            margin_bot = np.zeros((sigma_lcl_bot.shape[0]))
            reserve_factor_bot = np.zeros((sigma_lcl_bot.shape[0]))
            average_reserve_factor_temp = np.zeros((len(self.stackup)))
            average_strength_index_temp = np.zeros((len(self.stackup)))
            
            # loop on the plies 
            for i in range (len(self.stackup)):
                # exctraction of strenghts values
                Xt = self.stackup_long[i].get('Xt')
                Xc = self.stackup_long[i].get('Xc')
                if Xc < 0: 
                    print('\033[33m', "The provided Xc is negative. The sign was autmatically changed.", '\033[37m')
                    Xc = - Xc
                Yt = self.stackup_long[i].get('Yt')
                Yc = self.stackup_long[i].get('Yc')
                if Yc < 0: 
                    print('\033[33m', "The provided Yc is negative. The sign was autmatically changed.", '\033[37m')
                    Yc = - Yc
                S = self.stackup_long[i].get('S12')
                if cntrl_shear is True: 
                    S13 = self.stackup_long[i].get('S13')
                    S23 = self.stackup_long[i].get('S23')

                ## Calculation of Tsai coefficients
                F1 = 1 / Xt - 1 / Xc
                F11 = 1 / (Xt * Xc)
                F2 = 1 / Yt - 1 / Yc
                F22 = 1 / (Yt * Yc)
                F66 = 1 / (S ** 2)
                
                if cntrl_shear is True:
                    F44 = 1 / (S23) ** 2
                    F55 = 1 / (S13) ** 2
                if F12 != "-.5/FxtFxcFytFyc_max**2":
                    F12 = F12
                else:
                    # F12 = -.5 / Xt ** 2
                    F12 = - .5 / (Xt * Xc *  Yt * Yc) ** .5 

                ### TOP ###
                # exctration of stresses (x, y and shear)
                sigma1 = sigma_lcl_top[i, 0]
                sigma2 = sigma_lcl_top[i, 1]
                sigma3 = sigma_lcl_top[i, 2]
                if cntrl_shear is True:
                    tau1 = tau_lcl_top[i, 0]
                    tau2 = tau_lcl_top[i, 1]
                
                # calculation of the margin at the top of the ply
                a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2 + 2 * F12 * sigma1 * sigma2
                if cntrl_shear is True: 
                    a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2  + F44 * tau1 ** 2  + F55 * tau2 ** 2 + 2 * F12 * sigma1 * sigma2
                b = F1 * sigma1 + F2 * sigma2
                reserve_factor_top[i] = (- b / (2 * a) + ((b / (2 * a)) ** 2 + 1 / a) ** .5)
                # reserve_factor_top[i] = 1 / a * (- b + (b ** 2 + a) ** .5)
                margin_top[i] = reserve_factor_top[i] - 1
                del a 
                del b

                ### BOT ###
                sigma1 = sigma_lcl_bot[i, 0]
                sigma2 = sigma_lcl_bot[i, 1]
                sigma3 = sigma_lcl_bot[i, 2]
                if cntrl_shear is True:
                    tau1 = tau_lcl_bot[i, 0]
                    tau2 = tau_lcl_bot[i, 1]

                # calculation of the margin at the bottom of the ply
                # a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2 + 2 * F12 * sigma1 * sigma2
                # if cntrl_shear is True: 
                #     a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2  + F44 * tau1 ** 2  + F55 * tau2 ** 2 + 2 * F12 * sigma1 * sigma2
                a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2 + 2 * F12 * sigma1 * sigma2
                if cntrl_shear is True: 
                    a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2  + F44 * tau1 ** 2  + F55 * tau2 ** 2 + 2 * F12 * sigma1 * sigma2
                b = F1 * sigma1 + F2 * sigma2
                reserve_factor_bot[i] = (- b / (2 * a) + ((b / (2 * a)) ** 2 + 1 / a) ** .5)
                margin_bot[i] = reserve_factor_bot[i] - 1
                del a 
                del b

                # average reserve factor
                ### AVERAGE ###
                sigma1 = sigma_lcl_avrg[i, 0]
                sigma2 = sigma_lcl_avrg[i, 1]
                sigma3 = sigma_lcl_avrg[i, 2]
                if cntrl_shear is True:
                    tau1 = tau_lcl_avrg[i, 0]
                    tau2 = tau_lcl_avrg[i, 1]
                # calculation of the margin at the bottom of the ply
                # a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2 + 2 * F12 * sigma1 * sigma2
                # if cntrl_shear is True: 
                #     a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2  + F44 * tau1 ** 2  + F55 * tau2 ** 2 + 2 * F12 * sigma1 * sigma2
                a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2 + 2 * F12 * sigma1 * sigma2
                if cntrl_shear is True: 
                    a = F11 * sigma1 ** 2 + F22 * sigma2 ** 2 + F66 * sigma3 ** 2  + F44 * tau1 ** 2  + F55 * tau2 ** 2 + 2 * F12 * sigma1 * sigma2
                b = F1 * sigma1 + F2 * sigma2
                average_reserve_factor_temp[i] = (- b / (2 * a) + ((b / (2 * a)) ** 2 + 1 / a) ** .5)
                # margin_bot[i] = reserve_factor_bot[i] - 1
                average_strength_index_temp[i] = a + b
                del a 
                del b

            # definition of the margins as ply properties
            self.margin_top = margin_top
            self.reserve_factor_top = reserve_factor_top 
            self.inv_reserve_factor_top = (reserve_factor_top) ** -1 
            #
            self.margin_bot = margin_bot
            self.reserve_factor_bot = reserve_factor_bot 
            self.inv_reserve_factor_bot = (reserve_factor_bot) ** -1 
            #
            self.average_reserve_factor = average_reserve_factor_temp
            self.average_inv_reserve_factor = self.average_reserve_factor ** -1
            self.average_strenght_index = average_strength_index_temp
            self.margin_average = self.average_reserve_factor - 1

            # top and bottom check
            check_top = np.where(self.margin_top < 0)
            check_bot = np.where(self.margin_bot < 0)

            # if any of the top plies fails display a message  
            print(' ')
            if len(np.where(self.margin_top < 0)[0]) != 0 or len(np.where(self.margin_bot < 0)[0]) != 0:
                for i in range (int(len(check_top[0]))):
                    print('\033[31m', 'Failure expected on the top side of ply no. ' + str(check_top[0][i] + 1) + ' according to the Tsai-Wu criteria', '\033[37m',' ')
            # if any of the bot plies fails display a message  
                for i in range (int(len(check_bot[0]))):
                    print('\033[31m', 'Failure expected on the bottom side of ply no. ' + str(check_bot[0][i] + 1) + ' according to the Tsai-Wu criteria', '\033[37m',' ')
            # if none of the plies fails display a message  
            else: 
                print('\033[32m', 'No failure expected in the analysed laminate', '\033[37m', ' ')

            # print margins or not depending on the optional control variable
            if print_margins is True:
                self.print_margins_TsaiWu()
    # print the margins from a Tsai Wu verification (FPF required)
    def print_margins(self):
        '''
# DESCRIPTION:
Internal method to print the margins computed with the TsaiWu method

# INPUTS:

## Required
    -   None

# OUTPUTS:
    -   None

# Example:
    import PyComp as comp
    ply_name = 'ANSYS Epoxy Carbon Woven (230 GPa) Prepreg'
    ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
    ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

    laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
    N = [230, .10, -2.5]
    M = [-160, .012, -0.3]
    V = [.0005, 3.5]
    strght = [414, 414 * .5, 414, 414 * .5, 81.41, 40, 40]
    laminate.FPF([strght], N_in=N, M_in=M, V_in=V, criteria='TsaiWu')
    laminate.print_margins()

    or
   
    import PyComp as comp
    ply_name = 'ANSYS Epoxy Carbon Woven (230 GPa) Prepreg'
    ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
    ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

    laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
    N = [230, .10, -2.5]
    M = [-160, .012, -0.3]
    V = [.0005, 3.5]
    strght = [414, 414 * .5, 414, 414 * .5, 81.41, 40, 40]
    strain = [1/100, .5/100, 1/100, .5/100, 5/100, 5/100, 5/100]
    laminate.FPF([strght], N_in=N, M_in=M, V_in=V, criteria='MaxStress')
    laminate.print_margins()
    '''    
        error = 0  
        skip = 0
        try: 
            self.margin_top.shape[0] > 1
            self.margin_bot.shape[0] > 1
        # in case the attributes are note defined set the cntrl variable "error" to 1
        except AttributeError:
            error = 1  
        finally:
            # if error is 1 and so the attributes are not defined raise exception
            if error == 1:
                raise Exception('Either one between the class properties "margin_top" and "margin_bot" is not defined. Such properties are attached to the class by the FPF method. Please run the FPF method before running this one.')

        # internal variable 
        lgt = len(self.stackup)
        
        margin_top = self.margin_top
        margin_bot = self.margin_bot

        if self.criteria != 'TsaiWu':
            margin_top = np.min(self.margin_top, 1)            
            margin_bot = np.min(self.margin_bot, 1)            

        # loop aver the plies and print the ply name and then the ply margin 
        print('  ')
        for i in range(lgt):
            # if one of the values is equal or less than 0 print it in red, otherwise in the default white
            if margin_top[lgt - i - 1] <= 0:
                print('\033[31m', 'ply ' + str(lgt - i) + ' top margin: ' + str(margin_top[lgt - i - 1]), '\033[37m',' ')
            elif margin_top[lgt - i - 1] <= .5 and margin_top[lgt - i - 1] > 0:
                print('\033[33m', 'ply ' + str(lgt - i) + ' top margin: ' + str(margin_top[lgt - i - 1]), '\033[37m',' ')
            else:
                print('ply ' + str(lgt - i) + ' top margin: ' + str(margin_top[lgt - i - 1]))
            
            if margin_bot[lgt - i - 1] <= 0:
                print('\033[31m', 'ply ' + str(lgt - i) + ' bottom margin: ' + str(margin_bot[lgt - i - 1]), '\033[37m',' ')   
            elif margin_bot[lgt - i - 1] <= .5 and margin_bot[lgt - i - 1] > 0:
                print('\033[33m', 'ply ' + str(lgt - i) + ' bottom margin: ' + str(margin_bot[lgt - i - 1]), '\033[37m',' ')   
            else:
                print('ply ' + str(lgt - i) + ' bottom margin: ' + str(margin_bot[lgt - i - 1]))   
        print('  ')

