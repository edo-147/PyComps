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

if __name__ == '__main__':
    unittest.main()