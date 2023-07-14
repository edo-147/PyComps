import unittest
import numpy as np
import sys, os 
my_path = os.path.dirname(os.path.abspath("PyComp"))
sys.path.append(my_path)
import PyComp as comp


class TestLamCalc(unittest.TestCase):
    def test_init_inputs(self):
        ply_name = 'Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

        # mech_props_units is neither GPa nor MPa 
        with self.assertRaises(Exception):
            comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='Pa', hide_text=True)
        # ply_block has too many arguments
        with self.assertRaises(Exception):
            comp.Laminate([ply_name, ply_mech_props, ply_stkup, ply_stkup], mech_prop_units='GPa', hide_text=True)
        # ply_name is not a str
        with self.assertRaises(Exception):
            comp.Laminate([ply_mech_props, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)
        # ply_mech_props is not a list nor a numpy array
        with self.assertRaises(Exception):
            comp.Laminate([ply_name, ply_name, ply_stkup], mech_prop_units='GPa', hide_text=True)
        # ply_stkup is not a list
        with self.assertRaises(Exception):
            comp.Laminate([ply_name, ply_mech_props, ply_name], mech_prop_units='GPa', hide_text=True)
        # ply_mech_props is too long
        with self.assertRaises(Exception):
            ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, 1420, .275]
            comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)
        # one of ply_mech_props is not a float nor a integer
        with self.assertRaises(Exception):
            ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, '2.7', 1420, .275]
            comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)
    
    def test_init_outputs(self):
        ply_name = 'Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]
        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)

        ref_Q = np.array([[61438.3, 2457.53, 0], [2457.53, 61438.3, 0], [0, 0, 3300]])
        np.testing.assert_array_equal(np.round(laminate.Q[0], 2) - ref_Q, np.zeros((ref_Q.shape)))

        ref_A = np.array([[106354.84, 34215.994, 0], [34215.994, 106354.84, 0], [0, 0, 36069.423]])
        np.testing.assert_array_equal(laminate.A - ref_A, np.zeros((ref_A.shape)))

        ref_B = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
        np.testing.assert_array_equal(laminate.B - ref_B, np.zeros((ref_B.shape)))

        ref_D = np.array([[47253.877,  9443.026,     0.   ], [ 9443.026, 47253.877,     0.   ], [    0.   ,     0.   , 10190.575]])
        np.testing.assert_array_equal(laminate.D - ref_D, np.zeros((ref_D.shape)))

        ref_H = np.array([[4950, 0], [0, 4950]])
        np.testing.assert_array_equal(laminate.H - ref_H, np.zeros((ref_H.shape)))

        ref_stiff = np.array([[106354.84 , 34215.994, 0, 0, 0, 0], \
                              [34215.994, 106354.84 , 0, 0, 0, 0], \
                                [0, 0, 36069.423, 0, 0, 0], [0, 0, 0, 47253.877, 9443.026, 0], \
                                    [0, 0, 0, 9443.026, 47253.877, 0], [0, 0, 0, 0, 0, 10190.575]])
        np.testing.assert_array_equal(laminate.stiff - ref_stiff, np.zeros((ref_stiff.shape)))

        ref_inv_stiff = np.array([[ 1.05e-05, -3.37e-06,  0.00e+00,  0.00e+00,  0.00e+00,  0.00e+00], \
                [-3.37e-06,  1.05e-05,  0.00e+00,  0.00e+00,  0.00e+00,  0.00e+00], \
                [ 0.00e+00,  0.00e+00,  2.77e-05,  0.00e+00,  0.00e+00,  0.00e+00], \
                [ 0.00e+00,  0.00e+00,  0.00e+00,  2.20e-05, -4.40e-06,  0.00e+00], \
                [ 0.00e+00,  0.00e+00,  0.00e+00, -4.40e-06,  2.20e-05,  0.00e+00], \
                [ 0.00e+00,  0.00e+00,  0.00e+00,  0.00e+00,  0.00e+00,  9.81e-05]])
        np.testing.assert_array_equal(laminate.inv_stiff - ref_inv_stiff, np.zeros((ref_inv_stiff.shape)))

    def test_calc_eq_props(self):
        ply_name = 'Toray T300 - Epoxy 8552'
        ply_mech_props = [133.15, 16.931, 16.931, .264, .4361, .264, 5.8944, 5.7868, 5.8944, 1556.9, .275]
        ply_stkup = [0, 90, 45, -45, 45, -45, 45, 90, 0, 45]
        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)
        
        #### inputs ####

        # print_cntrl is not boolean 
        with self.assertRaises(Exception):
            laminate.calc_equivalent_properties(print_cntrl='Barbero', disp_waring=False)
        # method is not a string 
        with self.assertRaises(Exception):
            laminate.calc_equivalent_properties(method=1, disp_waring=False)
        # method is neither 'Barbero' nor 'ANSYS' 
        with self.assertRaises(Exception):
            laminate.calc_equivalent_properties(method='Wrong_input', disp_waring=False)

        #### outputs ####
        laminate.calc_equivalent_properties(method='Barbero', disp_waring=False)
        self.assertEqual(laminate.G_flex_eq, 17301.683)
        self.assertEqual(laminate.Ex_flex_eq, 69501.344)
        self.assertEqual(laminate.Ey_flex_eq, 51037.942)
        self.assertEqual(laminate.ni_flex_eq, .292)
        self.assertEqual(laminate.G_eq, 23718.279)
        self.assertEqual(laminate.Ex_eq, 49271.059)
        self.assertEqual(laminate.Ey_eq, 49271.059)
        self.assertEqual(laminate.ni_eq, .386)
        self.assertEqual(laminate.rm, 0.0297)
        self.assertEqual(laminate.rb, 0.00644)
        self.assertEqual(laminate.rn, 0.0103)

    def test_calc_sss_state_inputs(self):
        ply_name = 'Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]
        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)
        N = [230, .10, -2.5]
        M = [-160, .012, -0.3]
        V = [.0005, 3.5]       


        # N is not a list nor a numpy array
        with self.assertRaises(Exception):
            laminate.calculate_stress_state('N', M, V, print=True, print_shear=True)
        # M is not a list nor a numpy array
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, 'M', V, print=True, print_shear=True)
        # V is not a list nor a numpy array
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, M, 'V', print=True, print_shear=True)
        
        # N is a list and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.calculate_stress_state([230, '.10', -2.5], M, V, print=True, print_shear=True)
        # M is a list and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, [-160, '.012', -0.3], V, print=True, print_shear=True)
        # V is a list and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, M, [.0005, '3.5'], print=True, print_shear=True)
        
        # N is an array and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(np.array([230, '.10', -2.5]), M, V, print=True, print_shear=True)
        # M is an array and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, np.array([-160, '.012', -0.3]), V, print=True, print_shear=True)
        # V is an array and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, M, np.array([.0005, '3.5']), print=True, print_shear=True)

        # N has the wrong length
        with self.assertRaises(Exception):
            laminate.calculate_stress_state([230, -2.5], M, V, print=True, print_shear=True)
        # M has the wrong length
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, [-160, -0.3], V, print=True, print_shear=True)
        # V has the wrong length
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, M, [.0005], print=True, print_shear=True)
        
        # cntr_external_actions is not an int
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, M, V, print=True, print_shear=True, cntrl_external_actions='2')
        # cntr_external_actions is not 0, 1, 2 or 3
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, M, V, print=True, print_shear=True, cntrl_external_actions=4)
        # T_in is nor a float an int
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, M, V, print=True, print_shear=True, T_in='1')
        # m_in is nor a float an int
        with self.assertRaises(Exception):
            laminate.calculate_stress_state(N, M, V, print=True, print_shear=True, m_in='0')
        
    def test_calc_sss_state_outputs(self):
        ply_name = 'Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]
        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)
        N = [230, .10, -2.5]
        M = [-160, .012, -0.3]
        V = [.0005, 3.5]      
        laminate.calculate_stress_state(N, M, V, print=False)    

        sigma_bot_0_ref = np.array([ 3.824e+02, -7.970e+01, -1.000e-01])
        np.testing.assert_array_equal(np.round(laminate.sigma_bot[0], 1), sigma_bot_0_ref)
        sigma_top_3_ply_ref_ref = np.array([ 50.4,  54.5, -10.5])
        np.testing.assert_array_equal(np.round(laminate.sigma_top_ply_ref[3], 1), sigma_top_3_ply_ref_ref)
        tau_oop_bot_6_ref = np.array([0. , 1.2])
        np.testing.assert_array_equal(np.round(laminate.tau_oop_bot[6], 1), tau_oop_bot_6_ref)

        def_bot_4_ref = np.array([2.41e-03, -7.70e-04, -7.00e-05])
        np.testing.assert_array_equal(np.round(laminate.def_bot[4], 5), def_bot_4_ref)
        def_top_2_ref = np.array([ 3.38e-03, -9.70e-04, -6.00e-05])
        np.testing.assert_array_equal(np.round(laminate.def_top[2], 5), def_top_2_ref)
        gamma_out_2_ref = np.array([0, 0.00039])
        np.testing.assert_array_equal(np.round(laminate.gamma_out[2], 5), gamma_out_2_ref)

    def test_FPF_inputs(self):
        ply_name = 'Epoxy Carbon Woven (230 GPa) Prepreg'
        ply_mech_props = [61.34, 61.34, 6.9, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1420, .275]
        ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]
        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)
        N = [230, .10, -2.5]
        M = [-160, .012, -0.3]
        V = [.0005, 3.5]       
        strght = [414, 414 * .5, 414, 414 * .5, 81.41, 40, 40]
        strain = [1/100, .5/100, 1/100, .5/100, 5/100, 5/100, 5/100]


        # one of the strength vector elements is not a int nor a float
        with self.assertRaises(Exception):
            laminate.FPF([414, 414 * .5, 414, 414 * .5, 81.41, 40, '40'], N, M, V, criteria='TsaiWu')
        # strength vector length is wrong
        with self.assertRaises(Exception):
            laminate.FPF([414, 414 * .5, 414, 414 * .5, 81.41, 40], N, M, V, criteria='TsaiWu')
        # strength vector length is wrong
        with self.assertRaises(Exception):
            laminate.FPF([414, 414 * .5, 414, 414 * .5], N, M, criteria='TsaiWu')

        # two strength vectors are expected but one is provided
        with self.assertRaises(Exception):
            laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], ['ply_name', ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)
            laminate.FPF(strght, N, M, criteria='TsaiWu')
        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)

        # one of the strain vector elements is not a int nor a float
        with self.assertRaises(Exception):
            laminate.FPF([1/100, .5/100, 1/100, .5/100, 5/100, 5/100, '5/100'], N, M, V, criteria='MaxStrain')
        # strain vector length is wrong
        with self.assertRaises(Exception):
            laminate.FPF([1/100, .5/100, 1/100, .5/100, 5/100, 5/100], N, M, V, criteria='MaxStrain')
        # strain vector length is wrong
        with self.assertRaises(Exception):
            laminate.FPF([1/100, .5/100, 1/100, .5/100], N, M, criteria='MaxStrain')

        # two strength vectors are expected but one is provided
        with self.assertRaises(Exception):
            laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], ['ply_name', ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)
            laminate.FPF(strain, N, M, criteria='MaxStrain')
        laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa', hide_text=True)

        # N is not a list nor a numpy array
        with self.assertRaises(Exception):
            laminate.FPF(strght, 'N', M, V, criteria='TsaiWu')
        # M is not a list nor a numpy array
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, 'M', V, criteria='TsaiWu')
        # V is not a list nor a numpy array
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, M, 'V', criteria='TsaiWu')
        
        # N is a list and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.FPF(strght, [230, '.10', -2.5], M, V, criteria='TsaiWu')
        # M is a list and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, [-160, '.012', -0.3], V, criteria='TsaiWu')
        # V is a list and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, M, [.0005, '3.5'], criteria='TsaiWu')
        
        # N is an array and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.FPF(strght, np.array([230, '.10', -2.5]), M, V, criteria='TsaiWu')
        # M is an array and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, np.array([-160, '.012', -0.3]), V, criteria='TsaiWu')
        # V is an array and its elemnts are neither float nor int
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, M, np.array([.0005, '3.5']), criteria='TsaiWu')

        # N has the wrong length
        with self.assertRaises(Exception):
            laminate.FPF(strght, [230, -2.5], M, V, criteria='TsaiWu')
        # M has the wrong length
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, [-160, -0.3], V, criteria='TsaiWu')
        # V has the wrong length
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, M, [.0005], criteria='TsaiWu')

        # criteria is not one of the available ('MaxStress', 'MaxStrain' or 'TsaiWu')
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, M, V, criteria='wrongCriteria', cntrl_external_actions='2')

        # cntr_external_actions is not an int
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, M, V, criteria='TsaiWu', cntrl_external_actions='2')
        # cntr_external_actions is not 0, 1, 2 or 3
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, M, V, criteria='TsaiWu', cntrl_external_actions=4)
        # T_in is nor a float an int
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, M, V, criteria='TsaiWu', T_in='1')
        # m_in is nor a float an int
        with self.assertRaises(Exception):
            laminate.FPF(strght, N, M, V, criteria='TsaiWu', m_in='0')
    



if __name__ == '__main__':
    unittest.main()