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

        self.assertEqual(example_ply_1.cured_thickness, 0.14595)
        self.assertEqual(example_ply_1.fiber_G, 197.308)
        self.assertEqual(example_ply_1.matrix_G, 1.923)


if __name__ == '__main__':
    unittest.main()
