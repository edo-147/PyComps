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


if __name__ == '__main__':
    unittest.main()