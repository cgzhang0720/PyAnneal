#!/usr/bin/env python3

import numpy as np
import copy
import math
import argparse
import json

from collections import Iterable

from ConstNumber import ConstNumber

from Auxiliary import j_must_have, j_have, judge_equal, to_str, create_random_unit_100_direction,         \
                      create_random_unit_111_direction, second_elliptic_integral, first_elliptic_integral, \
                      inner_product, my_str_match, create_random_unit_vector_array_vertical,                \
                      create_random_unit_vector_array, create_unit_111_direction,                            \
                      function_reshape, solve_function

import tungsten.Create
import MyGlobalParameter

class DefectObject( object ) :

    def __init__( self, jdata, defectStr ) :
        ''' 
        defectStr format (1) or (2):
        (1): ILoop100 13 10.0 10.0 10.0
        (2): ILoop111 13 10.0 10.0 10.0 0.5 0.5 0.5

        Class DefectObject properties:

        a. set_defecttype_properties( self, jdata, tmp_defecttype )
          1. self.actions_list
          2. self.numb_actions, 
 
        b. self.set_position_properties( self, jdata, tmp_ConstNumber, tmp_position )
          3. self.position
             self.set_position_scale_image_int( self, jdata )

             4. self.position_image_int
             5. self.position_scale
                self._set_position_next_list( self, tmp_ConstNumber )
                  6. self.position_next_list
 
        c. 7. self.nsize
              set_nsize_properties( self, jdata, tmp_ConstNumber, tmp_nsize )

        d. 8. set self.trap 
        '''

        NumbConst = ConstNumber( jdata )

        defectStr = defectStr.lstrip( ).rstrip( )
        defectStr_list = defectStr.split( )

        tmp_defecttype = defectStr_list[ 0 ]

        'a. set defecttype properties'
        self.set_defecttype_properties( jdata, tmp_defecttype )

        tmp_nsize = int( defectStr_list[ 1 ] )
        tmp_position = np.array( list( map( float, defectStr_list[ 2 : 5 ] ) ) )

        if len( defectStr_list ) == 5 :
            burgers_vector_scale = None
        elif len( defectStr_list ) == 8 :
            burgers_vector_scale = np.array( list( map( float, defectStr_list[ 5 : 8 ] ) ) )
        else :
            raise RuntimeError( 'defectStr: %s, format is wrong' % defectStr )

        'b. 1. set self.burgers_vector_scale, 2. self.burgers_vector_unit, 3. self.burgers'
        self.set_burgers_vector( NumbConst, burgers_vector_scale )
        
        'c. set position properties'
        self.set_position_properties( jdata, NumbConst, tmp_position )

        'd. set nsize properties'
        self.set_nsize_properties( jdata, NumbConst, tmp_nsize )

        'e. default self.trap'
        self.trap = ( False, None )

    def set_defecttype_properties( self, 
                                   jdata,
                                   tmp_defecttype ) :
        'Set 1. self.actions_list, 2. self.numb_actions'

        defecttype_list = j_must_have( jdata, 'defecttype' )
        events_list = j_must_have( jdata, 'events' )

        if len( defecttype_list ) != len( events_list ) :
            raise RuntimeError( 'len( defecttype ) != len( events ), wrong!' )

        '1. self.actions_list, 2. self.numb_actions'
        for i in range( len( defecttype_list ) ) :
            if tmp_defecttype == defecttype_list[ i ] :

                self.actions_list = events_list[ i ]

                "reset self.actions_list if 'migrate_3D' or 'migrate_1D' in self.actions_list"
                if 'migrate_1D' in self.actions_list :

                    index = self.actions_list.index( 'migrate_1D' )
                    del self.actions_list[ index ]
                    self.actions_list.insert( index, 'migrate_1D_0' )
                    self.actions_list.insert( index + 1, 'migrate_1D_1' )

                elif 'migrate_3D' in self.actions_list :

                    index = self.actions_list.index( 'migrate_3D' )
                    del self.actions_list[ index ]

                    for j in range( 8 ) :
                        action = 'migrate_3D_' + str( j )
                        self.actions_list.insert( index + j, action )

                self.numb_actions = len( self.actions_list )
                break

    def set_burgers_vector( self, 
                            tmp_ConstNumber, 
                            tmp_burgers_vector_scale = None ) :
        '''
        Set 1. self.burgers_vector_scale, 2. self.burgers_vector_unit and 3. self.burgers
        for Iterable tmp_burgers_vector_scale or for defecttype 'Loop100', 'Loop111', 'I' or 'ICluster'.

        for defecttype of 'VCluster', has no burgers properties.
        '''

        if isinstance( tmp_burgers_vector_scale, Iterable ) :
            self.burgers_vector_scale = tmp_burgers_vector_scale
            self.burgers = np.linalg.norm( tmp_burgers_vector_scale ) * tmp_ConstNumber.alatt
            self.burgers_vector_unit = tmp_burgers_vector_scale * tmp_ConstNumber.alatt / self.burgers

            '''
            'Test_cgzhang, for detecting burgers properties.'
            print( ' defecttype:', self.defecttype )
            print( ' burgers', self.burgers )
            print( ' burgers_vector_scale', self.burgers_vector_scale )
            print( ' burgers_vector_unit', self.burgers_vector_unit )
            '''

        else :
            'set 1. self.burgers_vector_unit, 2. self.burgers, 3. self.burgers_vector_scale'
            if 'Loop100' in self.defecttype :

                self.burgers_vector_unit = create_random_unit_100_direction( )
                self.burgers = tmp_ConstNumber.alatt
                self.burgers_vector_scale = self.burgers_vector_unit

            elif 'Loop111' in self.defecttype or self.defecttype in [ 'I', 'ICluster' ]:

                self.burgers_vector_unit = create_random_unit_111_direction( )
                self.burgers = tmp_ConstNumber.first
                self.burgers_vector_scale = self.burgers_vector_unit * self.burgers / tmp_ConstNumber.alatt

                '''
                'Test_cgzhang, for detecting burgers properties.'
                print( ' defecttype:', self.defecttype )
                print( ' burgers', self.burgers )
                print( ' burgers_vector_scale', self.burgers_vector_scale )
                print( ' burgers_vector_unit', self.burgers_vector_unit )
                '''

    def _set_position_virtual_migrate( self,
                                       jdata,
                                       tmp_position ) :
        ''' 
        This fuction is used when virtural migration happens. In this situation, self.position_next_list should not be calculated.
        set self.position related properties: 1. self.position, 2. self.position_image_int, 3. self.position_scale 
        '''

        self.position = np.array( [ x for x in tmp_position ] )
        self.set_position_scale_image_int( jdata )

    def set_position_properties( self, 
                                 jdata,
                                 tmp_ConstNumber,
                                 tmp_position     ) :
        ''' 
        Set self.position related properties: 
        1. self.position, 2. self.position_image_int, 3. self.position_scale, 4. self.position_next_list. 
        '''

        self.position = np.array( [ x for x in tmp_position ] ) 
        self.set_position_scale_image_int( jdata )

        self._set_position_next_list( tmp_ConstNumber )

    def set_position_scale_image_int( self, jdata ) :
        'Set 1. self.position_image_int, 2. self.position_scale'

        tmp_box = np.array( j_must_have( jdata, 'box' ) )

        self.position_image_int = np.floor_divide( self.position, tmp_box )
        self.position_scale = np.divide( self.position, tmp_box ) - self.position_image_int

    def _set_position_next_list( self, tmp_ConstNumber ) :
        ''' 
        self.position_next_list is a list of next positions to migrate.
        1. 'migrate_1D' in self.actions_list, len( self.position_next_list ) eq 2, forwards and afterwards.
        2. 'migrate_3D' in self.actions_list, len( self.position_next_list ) eq 8, 8 nearest neighbour directions.
        3. others do not move, set self.position_next_list = [ ]. 
        '''

        self.position_next_list = [ ]

        if 'migrate_1D_0' in self.actions_list :

            step = tmp_ConstNumber.first * self.burgers_vector_unit
            next_position = [ self.position + step, self.position - step ]
            self.position_next_list.extend( next_position )

        elif 'migrate_3D_0' in self.actions_list :

            for i in range( 4 ) :
                step = create_unit_111_direction( i ) * tmp_ConstNumber.first
                next_position = [ self.position + step, self.position - step ]
                self.position_next_list.extend( next_position )

    def set_nsize_properties( self, 
                              jdata, 
                              NumbConst, 
                              tmp_nsize  ) :
        ''' 
        nsize related properties:
        1.  self.nsize
        2.  self.migration_energy
        3.  self.pref_migration
        4.  self.emission_energy
        5.  self.pref_emission
        6.  self.rotation_energy
        7.  self.pref_rotation
        8.  self.transform_to
        9.  self.transform_energy
        10. self.pref_transform
        11. self.rho
        12. self.delta_rho
        13. self.radius
        14. self.relax_volume
        15. self.area          
        '''

        self.nsize = tmp_nsize

        'first set shape properties.'
        self._set_shape_properties( jdata, NumbConst )

        if my_str_match( 'migrate', self.actions_list ) : self._set_migration( )

        if 'emit' in self.actions_list : self._set_emission( NumbConst )

        if 'rotate' in self.actions_list : self._set_rotation( )

        if 'transform' in self.actions_list : self._set_transformation( jdata )

    '1. need reshape'
    def _set_migration( self ) :
        'Set self.migration_energy and self.pref_migration for migrated DefectObject.'

        pass

    '2. need reshape'
    def _set_emission( self, tmp_NumbConst ) :
        'Set self.emission_energy and self.pref_emission for emiting DefectObject.'

        pass

    '3. need reshape'
    def _set_rotation( self ) :
        'Set self.rotation_energy and self.pref_rotation for rotating DefectObject.'

        pass

    '4. need reshape'
    def _reset_b_pn_rotation( self, tmp_ConstNumber ) :
        'Reset 1. self.burgers_vector_unit, 2. self.burgers_vector_scale and 3. self.position_next_list after rotation.'
        
        pass

    '5. need reshape'
    def _set_transformation( self, jdata ) :
        'Set self.transform_to, self.transform_energy and self.pref_transform for transformating DefectObject.'

        pass

    '6. need reshape'
    def _set_shape_properties( self, jdata, tmp_ConstNumber ) :
        ''' 
        Set 1. self.rho, 2. self.delta_rho, 3. self.radius, 4. self.relax_volume, 5. self.area 
        according to Class.shape : 1. particle, 2. sphere, 3. disk.
        1. self.rho, 2. self.delta_rho and 3. self.radius should be set for all defecttype.
        self.radius = self.rho + self.delta_rho
        4. self.relax_volume should be set for defecttypes of sphere VCluster and particle Vacancy.
        5. self.area should be set for defecttypes of disk Loops, ICluster and particle Interstitial. 
        '''

        pass

    def return_rate( self,
                     jdata,
                     temperature,
                     u0,
                     defectObject_list_27_cells ) :
        'Return a list of rate corresponding to self.actions_list'

        'default value is set True.'
        having_elastic_interaction = True
        if j_have( jdata, 'having_elastic_interaction' ) : having_elastic_interaction = jdata[ 'having_elastic_interaction' ]

        tmp_ConstNumber = ConstNumber( jdata )

        rate_list = [ ]

        for action in self.actions_list :

            if 'migrate' in action :
                rate = self.return_rate_migrate( jdata, tmp_ConstNumber, temperature, u0, defectObject_list_27_cells, action, having_elastic_interaction )
            elif action == 'emit' :
                rate = self.return_rate_emit( jdata, tmp_ConstNumber, temperature, u0, defectObject_list_27_cells, having_elastic_interaction )
            elif action == 'rotate' :
                rate = self.return_rate_rotate( jdata, tmp_ConstNumber, temperature, u0, defectObject_list_27_cells, having_elastic_interaction )
            elif action == 'transform' :
                rate = self.return_rate_transform( jdata, tmp_ConstNumber, temperature, u0, defectObject_list_27_cells, having_elastic_interaction )

            rate_list.append( rate )

        return rate_list

    def return_rate_test( self,
                          jdata,
                          temperature,
                          defectObject_list_27_cells ) :
        'Return a list of rate corresponding to self.actions_list for testing.'

        tmp_ConstNumber = ConstNumber( jdata )

        rate_list = [ ]

        for action in self.actions_list :

            if 'migrate' in action :
                rate = 10.0
            elif action == 'emit' :
                rate = 10.0
            elif action == 'rotate' :
                rate = 10.0
            elif action == 'transform' :
                rate = 10.0

            rate_list.append( rate )

        return rate_list

    def return_rate_migrate( self,
                             jdata,
                             tmp_ConstNumber,
                             temperature,
                             u0,
                             defectObject_list_27_cells,
                             action,
                             having_elastic_interaction  ) :
        ''' 
        Return rate of migration:
        u0: elastic interaction energy for self, i.e. ui = 0.5 * sigma_( j != i )^N U_elastic_energy, 
        j are neighbours of i, i.e. defectObject_list_27_cells.
        u1: elastic interaction energy for defectObject_virtual_migrate,
        this is true only when defectObject_list_27_cells is not changed for the migrating defectObject.
        '''

        activation = self.migration_energy
        if self.trap[ 0 ] : activation = activation + self.trap[ 1 ].energy

        if not having_elastic_interaction :
            return self.pref_migration * math.exp( -activation / ( tmp_ConstNumber.kb * temperature ) )

        defectObject_virtual_migrate = copy.deepcopy( self )
 
        '''
        'Test_cgzhang, for detecting virtual migration.'
        print( ' before migrate !' )
        defectObject_virtual_migrate.print_properties( )
        '''

        '''
        only reset properties: position, position_scale and position_image_int.
        not reset property position_next_list, because it has no effect on u1.
        '''
        defectObject_virtual_migrate._set_position_virtual_migrate( jdata, self.position_next_list[ int( action[ -1 ] ) ] )

        '''
        'Test_cgzhang, for detecting virtual migration.'
        print( ' after migrate !' )
        defectObject_virtual_migrate.print_properties( )
        '''

        'label_distance = V for virtual migration.'
        u1 = 0.5 * defectObject_virtual_migrate.return_u_neighbour( jdata, tmp_ConstNumber, defectObject_list_27_cells, 'V' )

        '''
        'Test_cgzhang, for detecting defectObject_list_27_cells'
        print( 'len(defectObject_list_27_cells)', len( defectObject_list_27_cells ) )
        for item in defectObject_list_27_cells :
            print( ' properties of defectObject_list_27_cells !' )
            item.print_properties( )
        '''

        '''
        attension ! activation_energy = 0.5 * delta_elastic_energy + migration_energy, 
        here delta_elastic_energy = ( u1 - u0 ) * 2.0, so activation = u1 - u0 + migration_energy. 
        '''
        activation = activation + ( u1 - u0 )

        '''
        'test_cgzhang.'
        if activation > 10 : 
            print( '# activation_migrate:', activation, 'defecttype:', self.defecttype )
            activation = 10.0
        if activation < -10 : 
            print( '# activation_migrate:', activation, 'defecttype:', self.defecttype )
            activation = -10.0
        '''

        '''
        'Test_cgzhang, for detecting activation energy.'
        print( '# u0: %0.10f, u1: %0.10f' % ( u0, u1 ) )
        if self.defecttype == 'ILoop111' : print( '# delta_utotal:', 2.0 * ( u1 - u0 ) )
        '''

        if activation > 10 : 
            activation = 10.0
        if activation < -10 : 
            activation = -10.0

        rate = self.pref_migration * math.exp( -activation / ( tmp_ConstNumber.kb * temperature ) )

        "print( '# rate_migrate per ps:', rate / 1.0E12 )"

        return rate

    def return_u_neighbour( self,
                            jdata,
                            tmp_ConstNumber,
                            defectObject_list,
                            label_distance = 'A' ) :
        '''
        Return elastic interaction energy between self and DefectObject_list, unit: eV.
        1. the default value of label_distance is set 'A', meaning actual distance, so delta depends periodic nearest images. 
        2. the label_distance is inputted with 'V' when getting the u_neighbour for virtual migration or emitting, 
           delta is nearest distance no mater periodic is true or false.
        '''

        u = 0.0

        for defectObject in defectObject_list :
            judge = self.return_judge_elastic_interaction( jdata, defectObject, label_distance )
            if judge[ 0 ] : u += self.return_u12( tmp_ConstNumber, defectObject, judge[ 1 ], judge[ 2 ] )

        return u

    def return_judge_elastic_interaction( self,
                                          jdata,
                                          tmp_defectObject,
                                          label_distance = 'A' ) :
        'Return_judge_elastic_interaction between two DefectObjects, return ( judge, delta, R_square ).'

        delta = self.return_delta( jdata, tmp_defectObject, label_distance )

        rcutoff_elastic = j_must_have( jdata, 'rcutoff_elastic' )

        for item in delta :
            if item > rcutoff_elastic :
                return ( False, )

        R_square = np.inner( delta, delta )

        return ( R_square < rcutoff_elastic * rcutoff_elastic, delta, R_square )

    def return_effect_or_not( self,
                              jdata,
                              tmp_defectObject ) :
        'return effect or not according to distance between self and tmp_defectObject.'

        delta = self.return_delta( jdata, tmp_defectObject )

        rcutoff_elastic = j_must_have( jdata, 'rcutoff_elastic' ) + 0.01

        for item in delta :
            if item > rcutoff_elastic :
                return False

        R_square = np.inner( delta, delta )

        return R_square < rcutoff_elastic * rcutoff_elastic

    '7. need reshape'
    def return_u12( self,
                    tmp_ConstNumber,
                    defectObject,
                    R_vector,
                    R_square         ) :
        'Return elastic interaction energy between two DefectObjects, unit: eV.'

        pass

    def _return_u_disk_disk( self,
                             tmp_ConstNumber,
                             tmp_defectObject,
                             R_vector,
                             R_square          ) :
        '''
        Return the elastic interaction energy between two disk DefectObjects.
        The derivations assume the defects are separated by distances greater than their size.
        Ref-1: D.R. Mason et al., J. Phys.: Condens. Matter 26 375701. P5.
        Ref-2: S.L. Dudarev and A.P. Sutton. Elastic interactions between nano-scale defects in irradiated materials. Acta Materialia 125 (2017) 425-430, P427.

        the length of burgers vector, b1 = self.burgers, b2 = tmp_defectObject.burgers. the area of loops, a1 = self.area, a2 = tmp_defectObject.area.
        the radius of loops, rho1 = self.rho, rho2 = tmp_defectObject.rho.
        '''

        R = math.sqrt( R_square )
        R_tri_vers = 1.0 / ( R * R_square + ( self.rho + tmp_defectObject.rho )**3.0 )

        'test_cgzhang, should considering zero R.'
        R_unit = R_vector / ( R + 0.000001 )

        'the unit of burgers vector.'
        n1 = self.burgers_vector_unit
        n2 = tmp_defectObject.burgers_vector_unit

        n1_dot_n2 = np.inner( n1, n2 )

        n1_dot_R_unit = np.inner( n1, R_unit )
        n2_dot_R_unit = np.inner( n2, R_unit )

        n1_dot_R_unit_square = n1_dot_R_unit * n1_dot_R_unit
        n2_dot_R_unit_square = n2_dot_R_unit * n2_dot_R_unit

        '''
        #Ref: D.R. Mason et al., J. Phys.: Condens. Matter 26 375701. P5.
        n1_cross_n2 = np.cross( n1, n2 )
        n1_cross_n2_dot_R_unit = np.inner( n1_cross_n2, R_unit )
        u_disk_disk = tmp_ConstNumber.u_disk_disk_prefactor *                                                           \
                      self.burgers * tmp_defectObject.burgers * self.area * tmp_defectObject.area * R_tri_vers *         \
                      (   6.0 * ( 1.0 - tmp_ConstNumber.poisson_ratio ) * n1_cross_n2_dot_R_unit * n1_cross_n2_dot_R_unit \
                        - 2.0 * ( 1.0 - tmp_ConstNumber.poisson_ratio ) * np.inner( n1_cross_n2, n1_cross_n2 )             \
                        - 3.0                                                                                               \
                        + 2.0 * n1_dot_n2 * n1_dot_n2                                                                        \
                        + 3.0 * n1_dot_R_unit_square                                                                          \
                        + 3.0 * n2_dot_R_unit_square                                                                           \
                        + 15.0 * n1_dot_R_unit_square * n2_dot_R_unit_square                                                    \
                        - 12.0 * n1_dot_R_unit * n2_dot_R_unit * n1_dot_n2    )
        '''

        '''
        #Ref: S.L. Dudarev and A.P. Sutton. Elastic interactions between nano-scale defects in irradiated materials. Acta Materialia 125 (2017) 425-430.
        '''
        u_disk_disk = tmp_ConstNumber.u_disk_disk_prefactor *                                                   \
                      self.burgers * tmp_defectObject.burgers * self.area * tmp_defectObject.area * R_tri_vers * \
                      ( 15.0 * n1_dot_R_unit_square * n2_dot_R_unit_square                                        \
                       -( 4 * tmp_ConstNumber.poisson_ratio - 1.0 )                                                \
                       -12.0 * tmp_ConstNumber.poisson_ratio * n1_dot_R_unit * n2_dot_R_unit * n1_dot_n2            \
                       -( 1.0 - 2.0 * tmp_ConstNumber.poisson_ratio ) * ( 3.0 * n1_dot_R_unit_square + 3.0 * n2_dot_R_unit_square + 2.0 * n1_dot_n2 * n1_dot_n2 ) )

        return u_disk_disk

    '''
    def _return_u_sphere_disk( self,
                               tmp_ConstNumber,
                               tmp_defectObject,
                               R_vector,
                               R_square          ) :
        'Return the elastic interaction energy between a voids( self ) and a disk( tmp_defectObject ). Ref: D.R. Mason et al., J. Phys.: Condens. Matter 26 375701. P5.'

        'the unit of burgers vector of loop tmp_defectObject.'
        n_unit = tmp_defectObject.burgers_vector_unit

        'z is the projection of R_vector on n_unit.'
        z = np.inner( R_vector, n_unit )
        z_square = z * z

        'the length of burgers vector of loop tmp_defectObject.'
        b = tmp_defectObject.burgers

        'the radius of loop tmp_defectObject.'
        rho = tmp_defectObject.rho
        r_square = abs( R_square - z_square )

        'Test_cgzhang, considering r_square .lt. 0, when R_square and z_square are very closer, so we add abs().'
        'the projection of R_vector on the plane of loop tmp_defectObject.'
        r = math.sqrt( r_square )

        rho_plus_r = rho + r

        rho_plus_r_square = rho_plus_r * rho_plus_r

        rho_substract_r = rho - r

        tmp_1 = 4.0 * r * rho

        tmp_2 = 1.0 / ( z_square + rho_plus_r_square )

        tmp = tmp_1 * tmp_2

        u_sphere_disk = tmp_ConstNumber.u_sphere_disk_prefactor * self.relax_volume * b * math.sqrt( tmp_2 )  \
                       * (    ( rho * rho - R_square ) / ( rho_substract_r * rho_substract_r + z_square )      \
                            * second_elliptic_integral( tmp )                                                   \
                          + first_elliptic_integral( tmp )                                                ) 

        return u_sphere_disk
    '''

    def _return_u_sphere_disk( self,
                               tmp_ConstNumber,
                               tmp_defectObject,
                               R_vector,
                               R_square          ) :
        'Return the elastic interaction energy between a voids( self ) and a disk( tmp_defectObject ). Ref: S.L. Dudarev and A.P. Sutton. Elastic interactions between nano-scale defects in irradiated materials. Acta Materialia 125 (2017) 425-430, P427.'

        'the unit of burgers vector of loop tmp_defectObject.'
        n_unit = tmp_defectObject.burgers_vector_unit

        'the length of burgers vector and area of loop tmp_defectObject.'
        b = tmp_defectObject.burgers
        A = tmp_defectObject.area

        R = math.sqrt( R_square )
        R_unit = R_vector / R

        R_tri = R * R_square

        R_unit_dot_n_unit = np.inner( R_unit, n_unit )

        p2 = 0.5 * ( 3.0 * R_unit_dot_n_unit * R_unit_dot_n_unit - 1.0 )

        u_sphere_disk = tmp_ConstNumber.u_sphere_disk_prefactor * A * b * self.relax_volume * p2 / R_tri

        return u_sphere_disk

    '''
    def _return_u_sphere_disk( self,
                               tmp_ConstNumber,
                               tmp_defectObject,
                               R_vector,
                               R_square          ) :
        #Return the elastic interaction energy between a voids( self ) and a disk( tmp_defectObject ).
        #Ref: D.R. Mason et al., J. Phys.: Condens. Matter 26 375701. P5.

        #original projections: z_original, r_original, z_square_original and r_square_original; 
        #reshape projections: z_reshape, r_reshape, z_square_reshape and r_square_reshape.

        #1. detect z_reshape, r_reshape. ????

        # 0. test_cgzhang
        #parameter used in function_reshape, data_
        #print( 'r0 in return_u12:', data_ )
        data_ = MyGlobalParameter.get_value( )

        'the unit of burgers vector of loop tmp_defectObject.'
        n_unit = tmp_defectObject.burgers_vector_unit

        'z is the projection of R_vector on n_unit, z .lt. 0 or z .gt. 0.'
        z_original = np.inner( R_vector, n_unit )

        # 1. test_cgzhang, added for reshape function f(r).
        z_original = abs( z_original )
        #z_reshape = solve_function( z_original, data_, delta = 1.0E-3 )
        z_reshape = z_original * function_reshape( z_original, data_ )

        z_square_original = z_original * z_original
        z_square_reshape = z_reshape * z_reshape

        'the length of burgers vector of loop tmp_defectObject.'
        b = tmp_defectObject.burgers

        'the radius of loop tmp_defectObject.'
        rho = tmp_defectObject.rho
        r_square_original = abs( R_square - z_square_original )

        'Test_cgzhang, considering r_square_original .lt. 0, when R_square and z_square_original are very closer, so we add abs().'
        #r_square_original = R_square - z_square_original

        'the projection of R_vector on the plane of loop tmp_defectObject.'
        r_original = math.sqrt( r_square_original )

        # 2. test_cgzhang, added for reshape function f(r).
        r_reshape = r_original * function_reshape( r_original, data_ )
        #r_reshape = solve_function( r_original, data_, delta = 1.0E-3 )

        r_square_reshape = r_reshape * r_reshape

        rho_plus_r = rho + r_reshape

        rho_plus_r_square = rho_plus_r * rho_plus_r

        rho_substract_r = rho - r_reshape

        tmp_1 = 4.0 * r_reshape * rho

        tmp_2 = 1.0 / ( z_square_reshape + rho_plus_r_square )

        tmp = tmp_1 * tmp_2

        u_sphere_disk = tmp_ConstNumber.u_sphere_disk_prefactor * self.relax_volume * b * math.sqrt( tmp_2 )                               \
                       * ( ( rho * rho - r_square_reshape - z_square_reshape ) / ( rho_substract_r * rho_substract_r + z_square_reshape )   \
                            * second_elliptic_integral( tmp ) + first_elliptic_integral( tmp )                                            ) 

        # test_cgzhang, for detecting original and reshape value.
        print( 'n_unit:', n_unit, 'R_square:', R_square, 'R_vector:', R_vector )
        print( 'z_original:', z_original, 'z_square_original:', z_square_original, 'z_reshape:', z_reshape, 'z_square_reshape:', z_square_reshape )
        print( 'r_original:', r_original, 'r_square_original:', r_square_original, 'r_reshape:', r_reshape, 'r_square_reshape:', r_square_reshape )

        return u_sphere_disk
    '''

    def return_distance( self, 
                         jdata, 
                         tmp_defectObject,
                         label_distance = 'A' ) :
        '''
        Return the distance between two defectObject according to label_distance and periodic nearest images.
        1. the default value of label_distance is set 'A', meaning actual distance, so delta depends periodic nearest images. 
        2. the label_distance is inputted with 'V' when getting the u_neighbour for virtual migration or emitting, 
           delta is nearest distance no mater periodic is true or false.
        '''

        return math.sqrt( self.return_distance_square( jdata, tmp_defectObject, label_distance ) )

    def return_distance_square( self,
                                jdata,
                                tmp_defectObject,
                                label_distance = 'A' ) :
        'Return the distance_square between two defectObject according to label_distance and periodic nearest images.'

        delta = self.return_delta( jdata, tmp_defectObject, label_distance )

        return np.inner( delta, delta )

    def return_delta( self,
                      jdata,
                      tmp_defectObject,
                      label_distance = 'A' ) :
        '''
        Return delta = ( dx, dy, dz ) > 0 between two DefectObjects according to label_distance and periodic nearest images.
        1. the default value of label_distance is set 'A', meaning actual distance, so delta depends periodic nearest images. 
        2. the label_distance is inputted with 'V' when getting the u_neighbour for virtual migration or emitting, 
           delta is nearest distance no mater periodic is true or false.
        '''

        tmp_periodic = np.array( j_must_have( jdata, 'periodic' ) )
        tmp_box = np.array( j_must_have( jdata, 'box' ) )

        p1 = self.return_position_in_box( jdata )
        p2 = tmp_defectObject.return_position_in_box( jdata )

        tmp_nearest = [ ]

        if label_distance == 'A' :
            'actual distance, so take care of periodic.'
            for i in range( 3 ) :
                if tmp_periodic[ i ] :
                    tmp_ =  min( [ abs( p1[ i ] - p2[ i ] + j * tmp_box[ i ] ) for j in range( -1, 2 ) ] )
                else :
                    tmp_ = p1[ i ] - p2[ i ]
                tmp_nearest.append( tmp_ )
        else :
            'virtual distance due to virtual migration or emitted, so do not take care of periodic, only getting the nearest distance.'
            for i in range( 3 ) :
                tmp_ =  min( [ abs( p1[ i ] - p2[ i ] + j * tmp_box[ i ] ) for j in range( -1, 2 ) ] )
                tmp_nearest.append( tmp_ )

        return np.array( tmp_nearest )

    def return_rate_emit( self,
                          jdata,
                          tmp_ConstNumber,
                          temperature,
                          u0,
                          defectObject_list_27_cells,
                          having_elastic_interaction  ) :
        ''' 
        Return rate of emit:
        u0: elastic interaction energy for self, i.e. ui = 0.5 * sigma_( j != i )^N U_elastic_energy, 
        j are neighbours of i, i.e. defectObject_list_27_cells.

        u1: elastic interaction energy for defectObject_virtual_emitting,
        u2: elastic interaction energy for defectObject_virtual_emitted,
        this is true only when defectObject_list_27_cells is not changed for defectObject_virtual_emitting and defectObject_virtual_emitted.
        '''

        activation = self.emission_energy

        if not having_elastic_interaction :
            return self.pref_emission * math.exp( -activation / ( tmp_ConstNumber.kb * temperature ) )

        defectObject_virtual_emitting, defectObject_virtual_emitted = self.take_emission( jdata, tmp_ConstNumber )

        u1 = 0.5 * defectObject_virtual_emitting.return_u_neighbour( jdata, tmp_ConstNumber, defectObject_list_27_cells )
        u2 = 0.5 * defectObject_virtual_emitted.return_u_neighbour( jdata, tmp_ConstNumber, defectObject_list_27_cells, 'V' )

        '''
        attension ! activation_energy = 0.5 * delta_elastic_energy + self.emission_energy, 
        here delta_elastic_energy = ( u2 + u1 - u0 ) * 2.0, so activation = u2 + u1 - u0 + self.emission_energy.
        '''
        activation = activation + ( u2 + u1 - u0 )

        'Test_cgzhang'
        if activation > 10 :
            #print( '# activation_emit:', activation )
            activation = 10.0
        if activation < -10:
            #print( '# activation_emit:', activation )
            activation = -10.0

        return self.pref_emission * math.exp( -activation / ( tmp_ConstNumber.kb * temperature ) )

    def return_rate_rotate( self,
                            jdata,
                            tmp_ConstNumber,
                            temperature,
                            u0,
                            defectObject_list_27_cells,
                            having_elastic_interaction  ) :
        ''' 
        Return rate of rotation :
        u0: elastic interaction energy for self, i.e. ui = 0.5 * sigma_( j != i )^N U_elastic_energy, 
        j are neighbours of i, i.e. defectObject_list_27_cells.
        u1: elastic interaction energy for defectObject_virtual_rotate.
        '''

        activation = self.rotation_energy

        if not having_elastic_interaction :
            return self.pref_rotation * math.exp( -activation / ( tmp_ConstNumber.kb * temperature ) )

        defectObject_virtual_rotate = self.take_rotation( tmp_ConstNumber )

        u1 = 0.5 * defectObject_virtual_rotate[0].return_u_neighbour( jdata, tmp_ConstNumber, defectObject_list_27_cells )

        '''
        attension ! activation_energy = 0.5 * delta_elastic_energy + rotation_energy, 
        here delta_elastic_energy = ( u1 - u0 ) * 2.0, so activation = u1 - u0 + rotation_energy. 
        '''
        activation = activation + ( u1 - u0 )

        'Test_cgzhang'
        if activation > 10 :
            #print( '# activation_rotation:', activation )
            activation = 10.0
        if activation < -10:
            #print( '# activation_rotation:', activation )
            activation = -10.0

        return self.pref_rotation * math.exp( -activation / ( tmp_ConstNumber.kb * temperature ) )

    def return_rate_transform( self,
                               jdata,
                               tmp_ConstNumber,
                               temperature,
                               u0,
                               defectObject_list_27_cells,
                               having_elastic_interaction  ) :
        ''' 
        Return rate of transformation:
        u0: elastic interaction energy for self, i.e. ui = 0.5 * sigma_( j != i )^N U_elastic_energy, 
        j are neighbours of i, i.e. defectObject_list_27_cells.
        u1: elastic interaction energy for defectObject_virtual_transform.
        '''

        activation = self.transform_energy

        if not having_elastic_interaction :
            return self.pref_transform * math.exp( -activation / ( tmp_ConstNumber.kb * temperature ) )

        defectObject_virtual_transform = self.take_transformation( jdata )

        u1 = 0.5 * defectObject_virtual_transform[0].return_u_neighbour( jdata, tmp_ConstNumber, defectObject_list_27_cells )

        '''
        attension ! activation_energy = 0.5 * delta_elastic_energy + transform_energy, 
        here delta_elastic_energy = ( u1 - u0 ) * 2.0, so activation = u1 - u0 + transform_energy. 
        '''
        activation = activation + ( u1 - u0 )

        'Test_cgzhang'
        if activation > 10 :
            #print( '# activation_transform:', activation )
            activation = 10.0
        if activation < -10:
            #print( '# activation_transform:', activation )
            activation = -10.0

        return self.pref_transform * math.exp( -activation / ( tmp_ConstNumber.kb * temperature ) )

    def defect_take_action( self,
                            jdata,
                            action_i ) :
        'Taking action according to self.actions_list[ action_i ].'

        tmp_ConstNumber = ConstNumber( jdata )

        action = self.actions_list[ action_i ]

        if 'migrate' in action :
            return self.take_migration( jdata, tmp_ConstNumber, action )

        if action == 'emit' :
            return self.take_emission( jdata, tmp_ConstNumber )

        if action == 'rotate' :
            return self.take_rotation( tmp_ConstNumber )

        if action == 'transform' :
            return self.take_transformation( jdata )

    def take_migration( self,
                        jdata,
                        tmp_ConstNumber,
                        action           ) :
        'Taking migration according to action.'

        i_migrate = int( action[ -1 ] )

        new_defectObject = copy.deepcopy( self )
        new_defectObject.set_position_properties( jdata, tmp_ConstNumber, self.position_next_list[ i_migrate ] )

        return ( new_defectObject, )

    def take_emission( self, jdata, tmp_ConstNumber ) :
        ''' 
        Taking emission.
        a. if tmp_defecttyle having burgers_vector_scale 
        then set the burgers_vector_scale of new_defectObject to self.burgers_vector_scale. 
        else randomly set burgers_vector_scale.
        b. emitting point defect apart from body defect with the distance of 2.0 * tmp_ConstNumber.alatt.
        c. random position selected according to unit_vector_array and unit_vector_array_vertical to self.burgers_vector_unit for hody defect having self.burgers_vector_scale.
        d. random position selected according to unit_vector_array for body defect having no self.burgers_vector_scale.
        '''

        '''
        'Test_cgzhang, for creating fixed-emitted point defect, when detecting sub-, add-, delete- defectObject.'
        seed = 101
        np.random.seed( seed )
        '''

        tmp_defecttype = self.transite_defecttype_decrease_nsize( self.defecttype, self.nsize - 1  )

        if self.have_burgers_vector_scale( tmp_defecttype ) :
            new_defectObject = tungsten.Create.create_defectObject( jdata, tmp_defecttype, self.nsize - 1, self.position, self.burgers_vector_scale )

            'for body defect having burgers_vector_scale, such as loops.'
            unit_vector_array_vertical = create_random_unit_vector_array_vertical( self.burgers_vector_unit )

            'random unit_vector_array.'
            'unit_vector_array = create_random_unit_vector_array( )'

            'random 111 unit_vector_array.'
            unit_vector_array = create_random_unit_111_direction( )

            position = self.position + self.rho * unit_vector_array_vertical + 2.0 * tmp_ConstNumber.alatt * unit_vector_array

        else :
            new_defectObject = tungsten.Create.create_defectObject( jdata, tmp_defecttype, self.nsize - 1, self.position )

            'for no burgers_vector_scale, such as vacancy clluster.'
            unit_vector_array = create_random_unit_vector_array( )
            position = self.position + ( self.rho + 2.0 * tmp_ConstNumber.alatt ) * unit_vector_array

        'the burgers_vector_scale of emitted I is randomly set.'
        'emit_defectObject = tungsten.Create.create_defectObject( jdata, self.defecttype[ 0 ], 1, position )'

        'the burgers_vector_scale of emitted I is the same as body defect.'
        point = self.defecttype[ 0 ]

        if point == 'I' :
            emit_defectObject = tungsten.Create.create_defectObject( jdata, point, 1, position, self.burgers_vector_scale )
        else :
            emit_defectObject = tungsten.Create.create_defectObject( jdata, point, 1, position )

        return ( new_defectObject, emit_defectObject )

    def take_rotation( self, tmp_ConstNumber ) :
        'Taking rotation.'

        '''
        'Test_cgzhang, for creating fixed-rotated defectObject, when detecting sub-, add-, delete- defectObject.'
        seed = 211
        np.random.seed( seed )
        '''

        new_defectObject = copy.deepcopy( self )
        new_defectObject._reset_b_pn_rotation( tmp_ConstNumber )

        return ( new_defectObject, )

    def take_transformation( self, jdata ) :
        'Taking transformation, the burgers_vector_scale of transform_to is randomly set.'

        '''
        'Test_cgzhang, for creating fixed-transformed defectObject, when detecting sub-, add-, delete- defectObject.'
        seed = 985
        np.random.seed( seed )
        '''

        return ( tungsten.Create.create_defectObject( jdata, self.transform_to, self.nsize, self.position ), )

    '8. need reshape'
    def transite_defecttype_decrease_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize       ) :
        'Tansite from type tmp_defecttype to type transite_type because decrease tmp_nsize get to the number range.'
   
        pass

    '9. need reshape'
    def transite_defecttype_increase_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize       ) :
        'Tansite from type tmp_defecttype to type transite_type because increase tmp_nsize >= item[ 2 ].'

        pass

    def reset_trap_defectObject( self, tmp_trap ) :
        'Reset trap according to tuple tmp_trap.'

        self.trap = tmp_trap
   
    def move_defectObject( self, 
                           jdata, 
                           tmp_move_vector, 
                           tmp_ConstNumber  ) :
        '''
        Move the defectObject with the vector tmp_move_vector, 
        so reset all related properties, such as self.position, self.position_image_int, self.position_scale and self.position_next_list.
        '''

        self.position += tmp_move_vector
        self.set_position_scale_image_int( jdata )
        self._set_position_next_list( tmp_ConstNumber )

    def return_judge_in_box( self, 
                             jdata,
                             tmp_max_box_int ) :
        'Judge the defectObject in box or not according to self.position_image_int.'

        tmp_box = np.array( j_must_have( jdata, 'box' ) )

        in_box = abs( self.position_image_int ) <= tmp_max_box_int
        delta_l = np.array( [ max( x, y ) for ( x, y ) in np.vstack( ( self.position - tmp_box * ( tmp_max_box_int + 1 ), -self.position - tmp_box * tmp_max_box_int ) ).transpose( 1, 0 ) ] )

        return ( in_box, delta_l )

    def _return_judge_recombine_sphere_sphere( self,
                                               jdata,
                                               tmp_defectObject ) :
        'Judge two sphere defectObject recombine or not.'

        rcutoff = self.radius + tmp_defectObject.radius

        delta = self.return_delta( jdata, tmp_defectObject )

        if any( delta > rcutoff ) :
            return False
        else :
            rcutoff_square = rcutoff * rcutoff

            delta_square = self.return_distance_square( jdata, tmp_defectObject )

            recombine = delta_square < rcutoff_square

            if recombine :
                print( '# ---------- void and void recombine ------------' )
                self.print_properties( )
                tmp_defectObject.print_properties( )
        
                print( '# delta_square:', delta_square )
                print( '# rcutoff_square:', rcutoff_square )

            return recombine

    def _return_judge_recombine_sphere_disk_ring( self,
                                                  jdata,
                                                  tmp_defectObject ) :
        '''
        Judge two defectObject recombine or not, self.shape == sphere and tmp_defectObject.shape == disk.
        recombine condition sphere interact with the ring of the disk.
        '''

        # Test_cgzhang
        #print( '# recombine ring !' )

        rcutoff = self.radius + tmp_defectObject.delta_rho
 
        rcutoff_square = rcutoff * rcutoff

        delta = self.return_delta( jdata, tmp_defectObject )

        'projection along the unit vector of disk plane.'
        z = np.inner( tmp_defectObject.burgers_vector_unit, delta )
  
        z_square = z * z

        delta_square = np.inner( delta, delta )

        'projection on the disk plane.'
        d =  math.sqrt( abs( delta_square - z_square ) )

        y = d - tmp_defectObject.rho

        r_square = y * y + z_square
      
        recombine = r_square < rcutoff_square

        if recombine :
            print( '# ---------- void and the ring of dislocation loop recombine ------------' )
            self.print_properties( )
            tmp_defectObject.print_properties( )

            print( '# d<rho?:', d < tmp_defectObject.rho )
            print( '# z=<R>.<n>:', z )
            print( '# d=sqrt(R^2-z^2):', d )
            print( '# r_square:', r_square )
            print( '# rcutoff_square:', rcutoff_square )

        return recombine

    def _return_judge_recombine_sphere_disk_all( self,
                                                 jdata,
                                                 tmp_defectObject ) :
        '''
        Judge two sphere defectObject recombine or not, self.shape == sphere and tmp_defectObject.shape == disk.
        sphere interact with all parts of the disk.
        '''

        # Test_cgzhang
        #print( '# recombine all !' )
       
        rcutoff = self.radius + tmp_defectObject.delta_rho

        rcutoff_square = rcutoff * rcutoff

        delta = self.return_delta( jdata, tmp_defectObject )

        'projection along the unit vector of disk plane.'
        z = abs( np.inner( tmp_defectObject.burgers_vector_unit, delta ) )

        z_square = z * z

        delta_square = np.inner( delta, delta )

        'projection on the disk plane.'
        d =  math.sqrt( abs( delta_square - z_square ) )

        y = d - tmp_defectObject.rho

        if y < 0 : 
            recombine = z < rcutoff
            if recombine :
                print( '# ----------void and all plane of dislocation loop recombine -----------' )
                self.print_properties( )
                tmp_defectObject.print_properties( )

                print( '# d<rho?:', d < tmp_defectObject.rho )
                print( '# z=<R>.<n>:', z )
                print( '# rcutoff:', rcutoff )

            return recombine

        else :

            r_square = y * y + z_square
        
            recombine = r_square < rcutoff_square
        
            if recombine :
                print( '# ----------void and all plane of dislocation loop recombine -----------' )
                self.print_properties( )
                tmp_defectObject.print_properties( )
        
                print( '# d<rho?:', d < tmp_defectObject.rho )
                print( '# z=<R>.<n>:', z )
                print( '# d=sqrt(R^2-z^2):', d )
                print( '# r_square:', r_square )
                print( '# rcutoff_square:', rcutoff_square )

            return recombine

    def _return_judge_recombine_sphere_disk( self,
                                             jdata,
                                             tmp_defectObject ) :

        'choose function _return_judge_recombine_sphere_disk_X depending on void_disk_recombine_all_part.'

        'default value True'
        void_disk_recombine_all_part = True
        if j_have( jdata, 'void_disk_recombine_all_part' ) : void_disk_recombine_all_part = jdata[ 'void_disk_recombine_all_part' ]

        if void_disk_recombine_all_part : 
            return self._return_judge_recombine_sphere_disk_all( jdata, tmp_defectObject )
        else :
            return self._return_judge_recombine_sphere_disk_ring( jdata, tmp_defectObject )

    def _return_judge_recombine_disk_disk( self,
                                           jdata,
                                           tmp_defectObject ) :
        'Judge two disk defectObject recombine or not.'

        if self._return_judge_recombine_sphere_disk( jdata, tmp_defectObject ) :
            return True
        else :
            return tmp_defectObject._return_judge_recombine_sphere_disk( jdata, self )

    def return_judge_recombine( self,
                                jdata,
                                tmp_defectObject ) :
        'Judge recombine or not for self.shape == disk, others will reset this function.'

        if tmp_defectObject.shape == 'sphere' or tmp_defectObject.defecttype == 'V' :
            return tmp_defectObject._return_judge_recombine_sphere_disk( jdata, self )
        else :
            return self._return_judge_recombine_disk_disk( jdata, tmp_defectObject )

    def return_recombine_results( self,
                                  jdata,
                                  tmp_defectObject ) :
        ''' 
        Three conditions:
        1. judge == False, return ( False, None ), nothing to do
        2. judge == True and two defectObjects disapear, return ( True, None )
        3. judge == True and new_defectObject is created, return ( True, new_defectObject ), 

        judge self.defecttype[ 0 ] == tmp_defectObject.defecttype[ 0 ] ? 
        set nsize, defecttype and position.

        if tmp_defecttype having burgers_vector_scale, then set new_defectObject to obj.burgers_scale_scale.
        else randomly set burgers_vector_scale of new_defectObject.
        '''

        judge = self.return_judge_recombine( jdata, tmp_defectObject )

        if not judge : return ( False, None )

        if self.nsize >= tmp_defectObject.nsize :
            obj = self
        else :
            obj = tmp_defectObject

        if self.defecttype[ 0 ] == tmp_defectObject.defecttype[ 0 ] :
            tmp_nsize = self.nsize + tmp_defectObject.nsize
            tmp_defecttype = obj.transite_defecttype_increase_nsize( obj.defecttype, tmp_nsize )
            tmp_position = ( self.position * self.nsize + tmp_defectObject.position * tmp_defectObject.nsize ) / ( self.nsize + tmp_defectObject.nsize )

            '''
            # Test_cgzhang
            print( '# increasing !' )
            '''
        else :
            tmp_nsize = abs( self.nsize - tmp_defectObject.nsize )
            if tmp_nsize == 0 : return ( True, None )

            '''
            # Test_cgzhang
            print( '# decrease!' )
            '''

            tmp_defecttype = obj.transite_defecttype_decrease_nsize( obj.defecttype, tmp_nsize )
            tmp_position = obj.position

        if self.have_burgers_vector_scale( tmp_defecttype ) :
            new_defectObject = tungsten.Create.create_defectObject( jdata, tmp_defecttype, tmp_nsize, tmp_position, obj.burgers_vector_scale )
        else :
            new_defectObject = tungsten.Create.create_defectObject( jdata, tmp_defecttype, tmp_nsize, tmp_position )

        return ( True, new_defectObject )

    def convert_point_defects( self, jdata ) :
        'Convert into point defects to output cfg file, format: position_scale + IV'

        alatt = j_must_have( jdata, 'alatt' )
        itmax = int( self.radius / alatt ) + 1

        IV_dict = { 'I' : 1, 'V' : 0 }

        IV = IV_dict[ self.defecttype[ 0 ] ]

        convert_point_defects_str = 'self.convert_point_defects_' + self.shape + '(jdata,IV,itmax)'

        return eval( convert_point_defects_str )

    def convert_point_defects_particle( self, 
                                        jdata, 
                                        IV,
                                        itmax  ) :
        ''' 
        Convert defectObject with the shape of particle into point defect str list.
        Format: xs ys zs IV. 
        '''

        return [ to_str( self.position_scale ) + str( IV ) + '\n' ]
 
    def convert_point_defects_sphere( self,
                                      jdata,
                                      IV,
                                      itmax  ) :
        ''' 
        Convert defectObject with the shape of sphere into point defect str list.
        Format: xs ys zs IV. 
        '''

        box = np.array( j_must_have( jdata, 'box' ) )
        lattice = j_must_have( jdata, 'lattice' )
        alatt = j_must_have( jdata, 'alatt' )

        if lattice == 'bcc' :
            basis = [ [ 0.0, 0.0, 0.0 ],
                      [ 0.5, 0.5, 0.5 ] ]

        elif lattice == 'fcc' :
            basis = [ [ 0.0, 0.0, 0.0 ],
                      [ 0.5, 0.5, 0.0 ],
                      [ 0.5, 0.0, 0.5 ],
                      [ 0.0, 0.5, 0.5 ] ]

        elif lattice == 'hcp' :
            print( '# need change for hcp' )

        position_scale_list = [ ]

        for i in range( -itmax, itmax ) :
            for j in range( -itmax, itmax ) :
                for k in range( -itmax, itmax ) :
                    for m in range( len( basis ) ) :
                        relative = np.array( [ i, j, k ] ) + np.array( basis[ m ] )
                        position_scale_list.append( relative )

        position_scale_list.sort( key = inner_product, reverse = False )

        if len( position_scale_list ) < self.nsize :
            raise RuntimeError( 'itmax is too small, wrong !' )

        point_defect_position_list = position_scale_list[ 0 : self.nsize ]

        result = [ ]

        for item in point_defect_position_list :
            tmp_position = self.position + item * alatt

            image_int = np.floor_divide( tmp_position, box )
            scale = np.divide( tmp_position, box ) - image_int

            point_defect_str = to_str( scale ) + str( IV ) + '\n'

            result.append( point_defect_str )

        return result
                                                                     
    def convert_point_defects_disk( self,
                                    jdata,
                                    IV,
                                    itmax  ) :
        ''' 
        Convert defectObject with the shape of disk into point defect str list.
        Format: xs ys zs IV. 
        '''

        box = np.array( j_must_have( jdata, 'box' ) )
        lattice = j_must_have( jdata, 'lattice' )
        alatt = j_must_have( jdata, 'alatt' )

        if lattice == 'bcc' :
            basis = [ [ 0.0, 0.0, 0.0 ],
                      [ 0.5, 0.5, 0.5 ]  ]

        elif lattice == 'fcc' :
            basis = [ [ 0.0, 0.0, 0.0 ],
                      [ 0.5, 0.5, 0.0 ],
                      [ 0.5, 0.0, 0.5 ],
                      [ 0.0, 0.5, 0.5 ]  ]

        elif lattice == 'hcp' :
            print( '# need change for hcp' )

        position_scale_list = [ ]

        scale_cutoff = self.burgers / alatt * 0.5

        for i in range( -itmax, itmax ) :
            for j in range( -itmax, itmax ) :
                for k in range( -itmax, itmax ) :
                    for m in range( len( basis ) ) :
                        relative = np.array( [ i, j, k ] ) + np.array( basis[ m ] )
                        if abs( np.inner( relative, self.burgers_vector_unit ) - 0.01 ) < scale_cutoff :
                            position_scale_list.append( relative )

        position_scale_list.sort( key = inner_product, reverse = False )

        if len( position_scale_list ) < self.nsize :
            raise RuntimeError( 'itmax is too small, wrong !' )

        point_defect_position_list = position_scale_list[ 0 : self.nsize ]

        result = [ ]

        for item in point_defect_position_list :
            tmp_position = self.position + item * alatt

            image_int = np.floor_divide( tmp_position, box )
            scale = np.divide( tmp_position, box ) - image_int

            point_defect_str = to_str( scale ) + str( IV ) + '\n'

            result.append( point_defect_str )

        return result

    def have_burgers_vector_scale( self, tmp_defecttype ) :
        'judge defecttype have burgers_vector_scale, return True or False.'
        defecttype_have_burgers_vector_scale = [ 'I', 'ICluster', 'ILoop100', 'ILoop111', 'VLoop100', 'VLoop111' ]

        return tmp_defecttype in defecttype_have_burgers_vector_scale

    def return_numb_actions( self ) :
        return self.numb_actions
 
    def return_actions_list( self ) :
        return self.actions_list

    def return_position( self ) :
        return self.position
 
    def return_position_image_int( self ) :
        return self.position_image_int

    def return_position_scale( self ) :
        return self.position_scale

    def return_position_next_list( self ) :
        return self.position_next_list

    def return_position_in_box( self, jdata ) :
        tmp_box = np.array( j_must_have( jdata, 'box' ) )
        return self.position_scale * tmp_box

    '10. need reshape'
    def _return_properties( self ) :

        properties_list = [ [ 'actions_list',         self.actions_list         ],
                            [ 'numb_actions',         self.numb_actions         ],
                            [ 'position',             self.position             ],
                            [ 'position_image_int',   self.position_image_int   ],
                            [ 'position_scale',       self.position_scale       ],
                            [ 'position_next_list',   self.position_next_list   ],
                            [ 'nsize',                self.nsize                ],
                            [ 'trap',                 self.trap                 ]  ]

        return properties_list

    def equal( self,
               tmp_defectObject,
               delta             ) :
        'Judge equal between two objects of Class DefectObject'

        self_properties = self._return_properties( )
        tmp_properties = tmp_defectObject._return_properties( )

        for i in range( len( self_properties ) ) :
            if not judge_equal( self_properties[ i ][ 1 ], tmp_properties[ i ][ 1 ], delta ) :
                print( '# ' + self_properties[ i ][ 0 ] + ' not equal !' )
                return False

        return True

    def print_properties( self ) :

        print( '# ' )
        print( '# ------------- properties for DefectObject --------------' )
        print( '\n'.join( ( '# %s : %s' % ( item[ 0 ], item[ 1 ] ) for item in self._return_properties( ) ) ) )
        print( '# ' )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class DefectObject detecting ---' )

    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
