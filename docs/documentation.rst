PyComps - easy calculation of composite properties
---------------------------------------------

- calculation of composite stiffness matrix and equivalent stiffness properties
- calculation of composite ply-by-ply stresses and deformation
- FPF verification according to the max stress, max strain and, Tsai-Wu criteria

Package usage: 
-----------------------------------------------------------

Download the package folder, append the link to your "sys" and import PyComp:

.. code-block:: console

    import sys
    sys.path.append(<folder_path>)
    import PyComps as comp

LaminateCalculator.py
---------------------------------------------

**Example_1: composite generation**

.. code-block:: python

    import PyComps as comp
    ply_name = 'Epoxy Carbon Woven (230 GPa) Prepreg'
    ply_mech_props = [60, 60, 7, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1450, .275]
    ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

    laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')

**Example_2: composite equivalent properties calculation**

.. code-block:: python

    import PyComps as comp
    ply_name = 'Epoxy Carbon Woven'
    ply_mech_props = [60, 60, 7, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1450, .275]
    ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

    laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
    laminate.calc_equivalent_properties(print_cntrl=True, method='Barbero')

**Example_3: calculate ply-by-ply stressess following the FSDT**

.. code-block:: python

    import PyComps as comp
    ply_name = 'Epoxy Carbon Woven'
    ply_mech_props = [60, 60, 7, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1450, .275]
    ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

    laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
    N = [230, .10, -2.5]
    M = [-160, .012, -0.3]
    V = [.0005, 3.5]
    laminate.calculate_stress_state(N, M, V, print=True, print_shear=True)

**Example_4: perform a Fist-Ply-Failure verification following Tsai-Wu criteria**

.. code-block:: python

    import PyComps as comp
    ply_name = 'Epoxy Carbon Woven'
    ply_mech_props = [60, 60, 7, 0.04, 0.3, 0.3, 3.3, 2.7, 2.7, 1450, .275]
    ply_stkup = [0, 45, 0, 45, 45, 0, 45, 0]

    laminate = comp.Laminate([ply_name, ply_mech_props, ply_stkup], mech_prop_units='GPa')
    N = [230, .10, -2.5]
    M = [-160, .012, -0.3]

    strght = [414, 414 * .5, 414, 414 * .5, 81.41]
    strain = [1/100, .5/100, 1/100, .5/100, 5/100]
    laminate.FPF([strght], N_in=N, M_in=M,  criteria='TsaiWu')

