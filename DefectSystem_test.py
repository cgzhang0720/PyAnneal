if __name__ == '__main__' :

    '''
    '1. plotting elastic energy ( utotal: eV ) among defectObjects from MD cascade simulation of bulk tungsten, 
    when a ILoop111 is moved with the vector of its burgers_vector_unit * j * burgers.'

    #------------------------------------- anneal.json BEG ----------------------------------------

    {
    "_comment": "individual cascade or muti-cascade, units: s, ion.nm^-2.s^-1, ion.nm^-2, 6.25E-4, fulence=1.0E-2",

    "individual_cascade":  true,
    "time":                1.0E0,

    "_comment": "cascades from files",
    "cascade_from_file":   true,
    "sys_defects_cfg":     "/home/cgzhang/anneal/MD-cascade/bulk/defects-cfg",
    "defects_cfg_prefix":  "defect",
    "systems":             "/home/cgzhang/anneal/program-now/program/examples/cascade",
    "cascade_prefix":      "cascade-IV",
    "cluster":             true,
    "IC_VC_rcut":          [ 0.3165, 0.3165 ],
    "VLoop_lower":         14,
    "ICluster_upper":      3,
    "probability":         0.2,
    "burgers_vector_unit": [ [  1.0,          0.0,          0.0         ],
                             [  0.0,          1.0,          0.0         ],
                             [  0.0,          0.0,          1.0         ],
                             [  0.577350269,  0.577350269,  0.577350269 ],
                             [ -0.577350269,  0.577350269,  0.577350269 ],
                             [  0.577350269, -0.577350269,  0.577350269 ],
                             [  0.577350269,  0.577350269, -0.577350269 ]  ],
    "burgers_vector_scale_length":  [ 1.0, 1.0, 1.0, 0.866025404, 0.866025404, 0.866025404, 0.866025404 ],


    "initial_recombine":  true,

    "_comment": "model parameters, units: s, K, nm",
    "temperature":        300.0,
    "seed":               10,

    "_comment": "reset parameters of materials and defects",
    "lattice":            "bcc",
    "mass":               183.84,
    "symbol":             "W",
    "alatt" :             0.3165,


    "_comment": "trap parameters, units: Volume ~ (nm)^-3 or Atomic ~ ppm, nm, eV",
    "trap":               true,
    "trap_density":       [ [ "Atomic", 30.0 ] ],
    "trap_radius":        [   0.3165           ],
    "trap_energy":        [   2.0              ],


    "_comment": "simulation box, units: nm",
    "box":                [ 100,   100,   100   ],
    "linkcell_space":     [ 5.0,   5.0,   5.0   ],
    "periodic":           [ false, false, false ],
    "max_box_image_int":  [ 1000,  1000,  1000  ],


    "_comment": "elastic interaction distance of loops, must .lt. linkcell_space, units: nm",
    "rcutoff_elastic":      4.9,


    "max_numb_in_27_cells": 5000,


    "_comment": "defect properties",
    "defecttype": ["I",                     "V",           "ICluster",                     "VCluster",           "ILoop100",           "VLoop100",           "ILoop111",            "VLoop111"],
    "events":     [["migrate_1D","rotate"], ["migrate_3D"],["migrate_1D","emit","rotate"], ["migrate_3D","emit"],["emit","transform"], ["emit","transform"], ["migrate_1D","emit"], ["emit","transform"]],


    "_comment": "output and frequence",
    "out_dir":            "/home/cgzhang/anneal/program-now/program/examples/",
    "out_style_freq":     [ "step_uniform", 1 ],

    "out_cfg_dir_name":            [ "cfg", "okmc"                      ],
    "out_defectStr_dir_name":      [ "defectStr", "defectStr"           ],
    "out_volume_density_dir_name": [ "density_volume", "density_volume" ],
    "out_area_density_dir_name":   [ "density_area", "density_area"     ],

    "out_point_defects_dir_name":  [ "pointdefect", "cascade-IV"        ],
    "out_steptime_dir_name":       [ "steptime", "steptime.dat"         ],

    "out_size_distri_dir_name":    "size_distri",
    "out_size_distri_defecttype":  [ "ICluster", "VCluster", "ILoop100", "VLoop100", "ILoop111", "VLoop111", "VAll", "IAll", "VLoopAll", "ILoopAll", "All" ],


    "_comment": "that's all"
    }

    # ------------------------------------ anneal.json END -------------------------------------

    parser = argparse.ArgumentParser( description = '--- Class DefectSystem detecting ---' )
    parser.add_argument( 'INPUT', help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    cascade = Cascade( jdata )

    tmp_ConstNumber = ConstNumber( jdata )
    delta = 1.0e-7
    for i in range( cascade.return_cascade_numb( ) ) :
        'build a defectSys of Class DefectSystem for okmc simulation.'
        model = DefectSystem( jdata, cascade.return_defectStr_list( jdata, i ) )
        'model.print_properties( )'

        'jumping i.'
        if i in [ 4, 11, 15, 20, 26, 28, 33 ] : continue
        index = -1
        max_size = 0
        for it in range( model.numb_defectObject ) :
            defect = model.defectObject_list[ it ]
            'defect.print_properties( )'

            'select index according to defect.nsize( maximum ) and defect.defecttype.'
            if max_size < defect.nsize and defect.defecttype == 'ILoop111' :
                max_size = defect.nsize
                index = it

        if index == -1 :
            raise RuntimeError( 'No defect is selected, wrong !' )
      
        defect = model.defectObject_list[ index ]
        print( '# i, defecttype, nsize, burgers_vector_unit:', i, defect.defecttype, defect.nsize, defect.burgers_vector_unit )

        old_defectObject = copy.deepcopy( defect )

        outputFile = OutputFile( model, jdata, i )

        dataStr_list = [ ]

        for j in range( -50, 50 ) :

            defectObject_virtual_migrate = copy.deepcopy( old_defectObject )
            move_vector = j * tmp_ConstNumber.first * old_defectObject.burgers_vector_unit
            defectObject_virtual_migrate.move_defectObject( jdata, move_vector )

            model.substitute_defectObject( jdata, index, defectObject_virtual_migrate )

            if j == -50 : 
                utotal_inf_1 = model.utotal

            if j == 50 : 
                utotal_inf_2 = model.utotal
                judge_equal( utotal_inf_1, utotal_inf_2, delta )

            dataStr = str( j ) + ' ' + str( model.utotal - utotal_inf_1 ) + '\n'
            dataStr_list.append( dataStr )

            'config output.'
            if outputFile.is_output_now( model, jdata ) :
                outputFile.output_results( model, jdata )

            model.reset_iter( )

        filenameStr = './utotal/ILoop111/utotal-' + str( i ) + '.dat'
        output_data_into_file( dataStr_list, filenameStr )
    '''



    '''
    '2. test self.cascade_from_file and self.cluster == True'

    'first step'
    numb_defectObject = 200
    seed = 30
    defectObject_list = create_random_defectObject_list( jdata, numb_defectObject, seed )

    model = DefectSystem( jdata, defectObject_list )

    model.output_cfg( jdata )
    model.output_defectStr( jdata )

    'second step'
    cascade = Cascade( jdata )

    for i in range( cascade.return_cascade_numb( ) ) :

        print( '# cascade', i )
        defectStr_list = cascade.return_defectStr_list( jdata, i )
        'cascade.print_properties( )'
    
    model = DefectSystem( jdata, defectStr_list )

    model.output_cfg( jdata )
    model.output_defectStr( jdata )

    '''


    '''
    '3. Test function: substitute_defectObject( self, jdata, index, tmp_defectObject )'

    # ------------------------------ anneal.json ----------------------------------------
    {
    "_comment": "individual cascade or muti-cascade, units: s, ion.nm^-2.s^-1, ion.nm^-2, 6.25E-4, fulence=1.0E-2",

    "individual_cascade":               true,
    "time":                             1.0E0,

    "_comment": "cascades from files",
    "cascade_from_file":   true,
    "sys_defects_cfg":     "/home/cgzhang/anneal/MD-cascade/bulk/defects-cfg",
    "defects_cfg_prefix":  "defect",
    "systems":             "/home/cgzhang/anneal/program-now/program/examples/cascade",
    "cascade_prefix":      "cascade-IV",
    "cluster":             true,
    "IC_VC_rcut":          [ 0.3165, 0.3165 ],
    "VLoop_lower":         14,
    "ICluster_upper":      3,
    "probability":         0.2,
    "burgers_vector_unit": [ [  1.0,          0.0,          0.0         ],
                             [  0.0,          1.0,          0.0         ],
                             [  0.0,          0.0,          1.0         ],
                             [  0.577350269,  0.577350269,  0.577350269 ],
                             [ -0.577350269,  0.577350269,  0.577350269 ],
                             [  0.577350269, -0.577350269,  0.577350269 ],
                             [  0.577350269,  0.577350269, -0.577350269 ]  ],
    "burgers_vector_scale_length": [ 1.0, 1.0, 1.0, 0.866025404, 0.866025404, 0.866025404, 0.866025404 ],


    "initial_recombine":  false,

    "_comment": "model parameters, units: s, K, nm",
    "temperature":      300.0,
    "seed":             10,

    "_comment": "reset parameters of materials and defects",
    "lattice":          "bcc",
    "mass":             183.84,
    "symbol":           "W",
    "alatt" :           0.3165,


    "_comment": "trap parameters, units: Volume ~ (nm)^-3 or Atomic ~ ppm, nm, eV",
    "trap":             true,
    "trap_density":     [ [ "Atomic", 30.0 ] ],
    "trap_radius":      [   0.3165           ],
    "trap_energy":      [   2.0              ],


    "_comment": "simulation box, units: nm",
    "box":               [ 100,  100,  100  ],
    "linkcell_space":    [ 5.0,  5.0,  5.0  ],
    "periodic":          [ true, true, true ],
    "max_box_image_int": [ 1000, 1000, 1000 ],


    "_comment": "elastic interaction distance of loops, must .lt. linkcell_space, units: nm",
    "rcutoff_elastic":      4.9,

    "max_numb_in_27_cells": 5000,

    "_comment": "elastic interaction distance of loops, must .lt. linkcell_space, units: nm",
    "rcutoff_elastic":      4.9,

    "max_numb_in_27_cells": 5000,

    "_comment": "defect properties",
    "defecttype": ["I",                     "V",           "ICluster",                     "VCluster",           "ILoop100",           "VLoop100",           "ILoop111",            "VLoop111"],
    "events":     [["migrate_1D","rotate"], ["migrate_3D"],["migrate_1D","emit","rotate"], ["migrate_3D","emit"],["emit","transform"], ["emit","transform"], ["migrate_1D","emit"], ["emit","transform"]],

    "_comment": "that's all"
    }
    # -------------------------------------- anneal.json End ----------------------------------------------------------

    parser = argparse.ArgumentParser( description = '--- Class DefectSystem detecting ---' )
    parser.add_argument( 'INPUT', help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    tmp_ConstNumber = ConstNumber( jdata )
    numb_defectObject = 1000
    delta = 1.0e-7
    for i in range( 1000 ) :

        np.random.seed( )
        'build defectObject_list via create_defectObject_list( jdata, numb_defectObject )'
        defectObject_list = create_random_defectObject_list( jdata, numb_defectObject )
        old_defectObject_list = copy.deepcopy( defectObject_list )

        index = np.random.randint( 1000 )

        print( '# substituded index:', index )

        'A. real value'
        tmp_defectObject = create_random_defectObject( jdata )
        copy_tmp_defectObject = copy.deepcopy( tmp_defectObject )

        defectObject_list[ index ] = tmp_defectObject
        model1 = DefectSystem( jdata, defectObject_list )
       
        print( '# --------------------------- model1 -----------------------' )
        model1.print_properties( 'utotal' )

        'B. substitude'
        model2 = DefectSystem( jdata, old_defectObject_list )
        model2.substitute_defectObject( jdata, index, copy_tmp_defectObject )

        print( '# --------------------------- model2 -----------------------' )
        model2.print_properties( 'utotal' )

        model2.set_which_defect_action( )

        if not model1.equal( model2, delta ) :
            raise RuntimeError( '# substitude model1 not equal model2, wrong:', i )

        print( '# iter: %d is OK' % i )
    '''





    '''
    '4. Test function: delete_defectObject( self, jdata, index_list )'

    # ------------------------------ anneal.json ----------------------------------------
    {
    "_comment": "individual cascade or muti-cascade, units: s, ion.nm^-2.s^-1, ion.nm^-2, 6.25E-4, fulence=1.0E-2",

    "individual_cascade":               true,
    "time":                             1.0E0,

    "_comment": "cascades from files",
    "cascade_from_file": true,
    "sys_defects_cfg":    "/home/cgzhang/anneal/MD-cascade/bulk/defects-cfg",
    "defects_cfg_prefix": "defect",
    "systems":           "/home/cgzhang/anneal/program-now/program/examples/cascade",
    "cascade_prefix":    "cascade-IV",
    "cluster":           true,
    "IC_VC_rcut":        [ 0.3165, 0.3165 ],
    "VLoop_lower":       14,
    "ICluster_upper":    3,
    "probability":         0.2,
    "burgers_vector_unit": [ [  1.0,          0.0,          0.0         ],
                             [  0.0,          1.0,          0.0         ],
                             [  0.0,          0.0,          1.0         ],
                             [  0.577350269,  0.577350269,  0.577350269 ],
                             [ -0.577350269,  0.577350269,  0.577350269 ],
                             [  0.577350269, -0.577350269,  0.577350269 ],
                             [  0.577350269,  0.577350269, -0.577350269 ]  ],
    "burgers_vector_scale_length": [ 1.0, 1.0, 1.0, 0.866025404, 0.866025404, 0.866025404, 0.866025404 ],


    "initial_recombine":  false,

    "_comment": "model parameters, units: s, K, nm",
    "temperature":      300.0,
    "seed":             10,

    "_comment": "reset parameters of materials and defects",
    "lattice":          "bcc",
    "mass":             183.84,
    "symbol":           "W",
    "alatt" :           0.3165,


    "_comment": "trap parameters, units: Volume ~ (nm)^-3 or Atomic ~ ppm, nm, eV",
    "trap":             true,
    "trap_density":     [ [ "Atomic", 30.0 ] ],
    "trap_radius":      [   0.3165           ],
    "trap_energy":      [   2.0              ],


    "_comment": "simulation box, units: nm",
    "box":               [ 100,  100,  100  ],
    "linkcell_space":    [ 5.0,  5.0,  5.0  ],
    "periodic":          [ true, true, true ],
    "max_box_image_int": [ 1000, 1000, 1000 ],


    "_comment": "elastic interaction distance of loops, must .lt. linkcell_space, units: nm",
    "rcutoff_elastic":      4.9,

    "max_numb_in_27_cells": 5000,

    "_comment": "elastic interaction distance of loops, must .lt. linkcell_space, units: nm",
    "rcutoff_elastic":      4.9,

    "max_numb_in_27_cells": 5000,

    "_comment": "defect properties",
    "defecttype": ["I",                     "V",           "ICluster",                     "VCluster",           "ILoop100",           "VLoop100",           "ILoop111",            "VLoop111"],
    "events":     [["migrate_1D","rotate"], ["migrate_3D"],["migrate_1D","emit","rotate"], ["migrate_3D","emit"],["emit","transform"], ["emit","transform"], ["migrate_1D","emit"], ["emit","transform"]],

    "_comment": "that's all"
    }
    # -------------------------------------- anneal.json End ----------------------------------------------------------

    parser = argparse.ArgumentParser( description = '--- Class DefectSystem detecting ---' )
    parser.add_argument( 'INPUT', help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    tmp_ConstNumber = ConstNumber( jdata )
    numb_defectObject = 1000
    delta = 1.0e-7
    for i in range( 1000 ) :

        np.random.seed( )
        'build defectObject_list via create_defectObject_list( jdata, numb_defectObject )'
        defectObject_list = create_random_defectObject_list( jdata, numb_defectObject )
        old_defectObject_list = copy.deepcopy( defectObject_list )

        'defect_index_list = [ 4, 7, 0, 40, numb_defectObject - 1, 507, 7 ]'
        defect_index_list = np.random.randint( numb_defectObject, size = ( 5 ) )
        old_defect_index_list = copy.deepcopy( defect_index_list )

        'A. real value'
        defect_index_list = list( set( defect_index_list ) )
        defect_index_list.sort( reverse = True )

        for item in defect_index_list :
            del defectObject_list[ item ]
        model1 = DefectSystem( jdata, defectObject_list )

        print( '# --------------------------- model1 -----------------------' )
        model1.print_properties( 'utotal' )

        'B. delete_defectObject'
        model2 = DefectSystem( jdata, old_defectObject_list )
        model2.delete_defectObject( jdata, old_defect_index_list )

        model2.set_which_defect_action( )

        print( '# --------------------------- model2 -----------------------' )
        model2.print_properties( 'utotal' )

        if not model1.equal( model2, delta ) :
            raise RuntimeError( '# delete_defectObject model1 not equal model2, wrong' )

        print( '# iter: %d is OK' % i )
 
    '''




    '''
    '5. Test function add_defectObject( self, jdata, tmp_defectObject_list, position_random=False )'

    # ------------------------------ anneal.json ----------------------------------------
    {
    "_comment": "individual cascade or muti-cascade, units: s, ion.nm^-2.s^-1, ion.nm^-2, 6.25E-4, fulence=1.0E-2",

    "individual_cascade":               true,
    "time":                             1.0E0,

    "_comment": "cascades from files",
    "cascade_from_file": true,
    "sys_defects_cfg":    "/home/cgzhang/anneal/MD-cascade/bulk/defects-cfg",
    "defects_cfg_prefix": "defect",
    "systems":           "/home/cgzhang/anneal/program-now/program/examples/cascade",
    "cascade_prefix":    "cascade-IV",
    "cluster":           true,
    "IC_VC_rcut":        [ 0.3165, 0.3165 ],
    "VLoop_lower":       14,
    "ICluster_upper":    3,
    "probability":         0.2,
    "burgers_vector_unit": [ [  1.0,          0.0,          0.0         ],
                             [  0.0,          1.0,          0.0         ],
                             [  0.0,          0.0,          1.0         ],
                             [  0.577350269,  0.577350269,  0.577350269 ],
                             [ -0.577350269,  0.577350269,  0.577350269 ],
                             [  0.577350269, -0.577350269,  0.577350269 ],
                             [  0.577350269,  0.577350269, -0.577350269 ]  ],
    "burgers_vector_scale_length": [ 1.0, 1.0, 1.0, 0.866025404, 0.866025404, 0.866025404, 0.866025404 ],


    "initial_recombine":  false,

    "_comment": "model parameters, units: s, K, nm",
    "temperature":      300.0,
    "seed":             10,

    "_comment": "reset parameters of materials and defects",
    "lattice":          "bcc",
    "mass":             183.84,
    "symbol":           "W",
    "alatt" :           0.3165,


    "_comment": "trap parameters, units: Volume ~ (nm)^-3 or Atomic ~ ppm, nm, eV",
    "trap":             true,
    "trap_density":     [ [ "Atomic", 30.0 ] ],
    "trap_radius":      [   0.3165           ],
    "trap_energy":      [   2.0              ],


    "_comment": "simulation box, units: nm",
    "box":               [ 100,  100,  100  ],
    "linkcell_space":    [ 5.0,  5.0,  5.0  ],
    "periodic":          [ true, true, true ],
    "max_box_image_int": [ 1000, 1000, 1000 ],


    "_comment": "elastic interaction distance of loops, must .lt. linkcell_space, units: nm",
    "rcutoff_elastic":      4.9,

    "max_numb_in_27_cells": 5000,

    "_comment": "elastic interaction distance of loops, must .lt. linkcell_space, units: nm",
    "rcutoff_elastic":      4.9,

    "max_numb_in_27_cells": 5000,

    "_comment": "defect properties",
    "defecttype": ["I",                     "V",           "ICluster",                     "VCluster",           "ILoop100",           "VLoop100",           "ILoop111",            "VLoop111"],
    "events":     [["migrate_1D","rotate"], ["migrate_3D"],["migrate_1D","emit","rotate"], ["migrate_3D","emit"],["emit","transform"], ["emit","transform"], ["migrate_1D","emit"], ["emit","transform"]],

    "_comment": "that's all"
    }
    # -------------------------------------- anneal.json End ----------------------------------------------------------

    parser = argparse.ArgumentParser( description = '--- Class DefectSystem detecting ---' )
    parser.add_argument( 'INPUT', help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    tmp_ConstNumber = ConstNumber( jdata )
    numb_defectObject = 1000
    delta = 1.0e-7
    for i in range( 1000 ) :

        np.random.seed( )
        'build defectObject_list via create_defectObject_list( jdata, numb_defectObject )'
        defectObject_list = create_random_defectObject_list( jdata, numb_defectObject )
        old_defectObject_list = copy.deepcopy( defectObject_list )

        np.random.seed( )
        add_defectObject_list = create_random_defectObject_list( jdata, 7 )
        old_add_defectObject_list = copy.deepcopy( add_defectObject_list )

        'A. real'
        defectObject_list.extend( add_defectObject_list )
        model1 = DefectSystem( jdata, defectObject_list )

        print( '# --------------------------- model1 -----------------------' )
        model1.print_properties( 'utotal' )

        'B. add_defectObject'
        model2 = DefectSystem( jdata, old_defectObject_list )
        model2.add_defectObject( jdata, old_add_defectObject_list )

        model2.set_which_defect_action( )

        print( '# --------------------------- model2 -----------------------' )
        model2.print_properties( 'utotal' )

        if not model1.equal( model2, delta ) :
            raise RuntimeError( '# add_defectObject model1 not equal model2, wrong' )

        print( '# iter: %d is OK' % i )

    '''
