import unittest
import numpy as np
import sys, os 
my_path = os.path.dirname(os.path.abspath("PyComp"))
sys.path.append(my_path)
import PyComp as comp


class TestMicromech(unittest.TestCase):

    def test_init_inputs(self):
        # fiber_properties is neither a string nor a numpy array 
        with self.assertRaises(Exception):
            comp.PlyDef('[513, .3, 1910]', 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, grams_per_square_meter=160, compute_cured_thickness=True)
        # fiber_name is not a string
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 5, [5, .3, 1156], 'Ex-1515', fiber_frac=.69, grams_per_square_meter=160, compute_cured_thickness=True)
        # matrix_properties is neither a string nor a numpy array 
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', '[5, .3, 1156]', 'Ex-1515', fiber_frac=.69, grams_per_square_meter=160, compute_cured_thickness=True)
        # matrix_name is not a string
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 5, fiber_frac=.69, grams_per_square_meter=160, compute_cured_thickness=True)
        # fiber_frac is not float nor zero
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac='.69', grams_per_square_meter=160, compute_cured_thickness=True)
        # matrix_frac is not float nor zero
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', matrix_frac='.69', grams_per_square_meter=160, compute_cured_thickness=True)
        # draping is not a string
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, draping=5, grams_per_square_meter=160, compute_cured_thickness=True)
        # mass_or_vol_frac is not a string
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac=5, draping=5, grams_per_square_meter=160, compute_cured_thickness=True)
        # mass_or_vol_frac is neither "vol" nor "wgt"
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='test', draping=5, grams_per_square_meter=160, compute_cured_thickness=True)
        # grams_per_square_meter is neither float nor integer 
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, grams_per_square_meter='160', compute_cured_thickness=True)
        # compute_cured_thickness is not boolean
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, grams_per_square_meter=160, compute_cured_thickness=0)
        # fiber_frac and matrix_frac are both zero
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=0, matrix_frac=0, grams_per_square_meter=160, compute_cured_thickness=0)
        # fiber_frac is negative
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=-.25, grams_per_square_meter=160, compute_cured_thickness=0)
        # martix_frac is negative
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', matrix_frac=-.25, grams_per_square_meter=160, compute_cured_thickness=0)
        # fiber_frac is larger than 1
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=1.25, grams_per_square_meter=160, compute_cured_thickness=0)
        # martix_frac is larger than 1
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', matrix_frac=1.25, grams_per_square_meter=160, compute_cured_thickness=0)
        # fiber_frac + martix_frac is not 1
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.37, matrix_frac=.25, grams_per_square_meter=160, compute_cured_thickness=0)
        # draping is not "UD"
        with self.assertRaises(Exception):
            comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, draping='test', grams_per_square_meter=160, compute_cured_thickness=0)
    
    def test_init_outputs(self):
        # attibutes calculation    
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        example_ply_1.cured_thickness = np.round(example_ply_1.cured_thickness, 5)
        example_ply_1.fiber_G = np.round(example_ply_1.fiber_G, 3)
        example_ply_1.matrix_G = np.round(example_ply_1.matrix_G, 3)
        example_ply_1.rho = np.round(example_ply_1.rho, 3)

        self.assertEqual(example_ply_1.cured_thickness, 0.14595)
        self.assertEqual(example_ply_1.fiber_G, 197.308)
        self.assertEqual(example_ply_1.matrix_G, 1.923)
        self.assertEqual(example_ply_1.rho, 1588.758)

    def test_ROM_input(self):
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        # print_cntrl is not boolean 
        with self.assertRaises(Exception):
            example_ply_1.ROM(print_cntrl=7)
        with self.assertRaises(Exception):
            example_ply_1.ROM(print_cntrl='7')

    def test_ROM_outputs(self):
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        example_ply_1.ROM()

        self.assertEqual(np.round(example_ply_1.E1, 4), 296.5661)
        self.assertEqual(np.round(example_ply_1.E2, 4), 11.5836)
        self.assertEqual(np.round(example_ply_1.E3, 4), 11.5836)
        self.assertEqual(np.round(example_ply_1.ni12, 4), .3)
        self.assertEqual(np.round(example_ply_1.ni13, 4), .3)
        self.assertEqual(example_ply_1.ni23, 'NA')
        self.assertEqual(np.round(example_ply_1.G12, 4), 4.4552)
        self.assertEqual(np.round(example_ply_1.G23, 4), 1.9231)
        self.assertEqual(np.round(example_ply_1.G13, 4), 4.4552)
    
    def test_Halphin_Tsai_input(self):
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        # print_cntrl is not boolean 
        with self.assertRaises(Exception):
            example_ply_1.Halphin_Tsai(print_cntrl=7)
        with self.assertRaises(Exception):
            example_ply_1.Halphin_Tsai(print_cntrl='7')
        # csi_E is neither float nor integer
        with self.assertRaises(Exception):
            example_ply_1.Halphin_Tsai(csi_E='7')
        # csi_G is neither float nor integer nor 'def'
        with self.assertRaises(Exception):
            example_ply_1.Halphin_Tsai(csi_G='7')
        # csi_E is negative
        with self.assertRaises(Exception):
            example_ply_1.Halphin_Tsai(csi_E=-7)
        # csi_G is negative
        with self.assertRaises(Exception):
            example_ply_1.Halphin_Tsai(csi_G=-7)

    def test_Halphin_Tsai_outputs(self):
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        example_ply_1.Halphin_Tsai()

        self.assertEqual(np.round(example_ply_1.E1, 4), 296.5661)
        self.assertEqual(np.round(example_ply_1.E2, 4), 23.8974)
        self.assertEqual(np.round(example_ply_1.E3, 4), 23.8974)
        self.assertEqual(np.round(example_ply_1.ni12, 4), .3)
        self.assertEqual(np.round(example_ply_1.ni13, 4), .3)
        self.assertEqual(np.round(example_ply_1.ni23, 4), .3)
        self.assertEqual(np.round(example_ply_1.G12, 4), 7.1673)
        self.assertEqual(np.round(example_ply_1.G23, 4), 5.8132)
        self.assertEqual(np.round(example_ply_1.G13, 4), 7.1673)

    def test_PMM_input(self):
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        # print_cntrl is not boolean 
        with self.assertRaises(Exception):
            example_ply_1.PMM(print_cntrl=7)
        with self.assertRaises(Exception):
            example_ply_1.PMM(print_cntrl='7')

    def test_PMM_outputs(self):
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        example_ply_1.PMM()

        self.assertEqual(np.round(example_ply_1.E1, 4), 296.5661)
        self.assertEqual(np.round(example_ply_1.E2, 4), 17.4726)
        self.assertEqual(np.round(example_ply_1.E3, 4), 17.4726)
        self.assertEqual(np.round(example_ply_1.ni12, 4), .3)
        self.assertEqual(np.round(example_ply_1.ni13, 4), .3)
        self.assertEqual(np.round(example_ply_1.ni23, 4), 0.3524)
        self.assertEqual(np.round(example_ply_1.G12, 4), 6.9391)
        self.assertEqual(np.round(example_ply_1.G23, 4), 6.46)
        self.assertEqual(np.round(example_ply_1.G13, 4), 6.9391)

    def test_print_properties_outputs(self):
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        with self.assertRaises(Exception):
            example_ply_1.print_properties()        
        
    def test_error_percent_inputs(self):
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        example_ply_1.PMM()

        with self.assertRaises(Exception):
            data_in = 'data'
            example_ply_1.error_percent(data_in)
        with self.assertRaises(Exception):
            data_in = [295, 17, 17, .25, .25, .27, 6.8, 6.2, 6.8, 1575, 1575]
            example_ply_1.error_percent(data_in)
        with self.assertRaises(Exception):
            data_in = [295, 17, 17, .25, .25, .27, 6.8, 6.2, 6.8, '1575']
            example_ply_1.error_percent(data_in)
        with self.assertRaises(Exception):
            data_in = [295, 17, 17, .25, .25, .27, 6.8, 6.2, 6.8, 1575]
            example_ply_1.error_percent(data_in, print_cntrl='False')

    def test_error_percent_outputs(self):
        example_ply_1 = comp.PlyDef([513, .3, 1910], 'M55J/Toray', [5, .3, 1156], 'Ex-1515', fiber_frac=.69, mass_or_vol_frac='wgt', grams_per_square_meter=160, compute_cured_thickness=True)
        example_ply_1.PMM()
        data_in = [295, 17, 17, .25, .25, .27, 6.8, 6.2, 6.8, 1575]
        example_ply_1.error_percent(data_in)
        self.assertEqual(np.round(example_ply_1.error_percent_E1, 4), .5309)
        self.assertEqual(np.round(example_ply_1.error_percent_E2, 4), 2.7801)
        self.assertEqual(np.round(example_ply_1.error_percent_E3, 4), 2.7801)
        self.assertEqual(np.round(example_ply_1.error_percent_ni12, 4), 20)
        self.assertEqual(np.round(example_ply_1.error_percent_ni13, 4), 20)
        self.assertEqual(np.round(example_ply_1.error_percent_ni23, 4), 0.0121)
        self.assertEqual(np.round(example_ply_1.error_percent_G12, 4), 2.046)
        self.assertEqual(np.round(example_ply_1.error_percent_G23, 4), 4.1931)
        self.assertEqual(np.round(example_ply_1.error_percent_G13, 4), 2.046)
        self.assertEqual(np.round(example_ply_1.error_percent_rho, 4), .8735)

if __name__ == '__main__':
    unittest.main()
