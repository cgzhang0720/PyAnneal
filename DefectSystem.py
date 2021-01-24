#!/usr/bin/env python3
'This is an DefectSystem module.'

import os
import numpy as np
import argparse
import json
import math
import copy
import time

from collections import Iterable

from Cascade import Cascade
from ConstNumber import ConstNumber
from BoxLinkcell import BoxLinkcell
from tungsten.DefectObject import DefectObject
from Auxiliary import j_must_have, j_have, judge_equal, output_data_into_file, save_float_several_digits
from tungsten.Create import create_defectObject_of_subclass, create_random_defectObject_list, create_random_defectObject, create_defectObject

from OutputFile import OutputFile
import MyGlobalParameter

class DefectSystem( object ) :

    def __init__( self, 
                  jdata, 
                  defectObjectStr_list ) :
        '''
        Properties of DefectSystem :
         0. self.box
         1. self.c_step
         2. self.c_time
         3. self.temperature_evolution
         4. self.individual_cascade
         5. self.multi_cascade
         6. self.time_evolution
         7. self.step_evolution
         8. self.having_elastic_interaction
         9. self.max_box_int
         10. self.numb_defectObject
         11. self.defectObject_list
         12. self.linkcell
         13. self.initial_recombine
         14. self.total_defects
         15. self.sum_numb_actions_list
         16. self.ui_array
         17. self.utotal
         18. self.rate_array
         19. self.size_rate_array
         20. self.sum_rate_array
         21. self.sum_rate
         22. self.delta_time
         23. self.which_defect
         24. self.what_action
        '''

        tmp_ConstNumber = ConstNumber( jdata )

        'set 0.self.box, 1.self.c_step, 2.self.c_time, 3.self.temperature_evolution, 4.self.individual_cascade, 5.self.multi_cascade. 6.self.time_evolution, 7.self.step_evolution.'

        self.box = np.array( j_must_have( jdata, 'box' ) )

        self.c_step = 0
        self.c_time = 0.0
        self.temperature_evolution = j_must_have( jdata, 'temperature' )

        self.individual_cascade = False
        if j_have( jdata, 'individual_cascade' ) : self.individual_cascade = jdata[ 'individual_cascade' ]

        self.multi_cascade = False
        if j_have( jdata, 'multi_cascade' ) : self.multi_cascade = jdata[ 'multi_cascade' ]

        if self.individual_cascade : 
            self.time_evolution = j_must_have( jdata, 'time' )
            self.step_evolution = math.inf
            if j_have( jdata, 'step' ) : self.step_evolution = jdata[ 'step' ]

        if self.multi_cascade : self.time_evolution = 1.0 / ( j_must_have( jdata, 'fluence_rate' ) * self.box[ 0 ] * self.box[ 1 ] )

        'set 8.self.having_elastic_interaction.'
        self.having_elastic_interaction = True
        if j_have( jdata, 'having_elastic_interaction' ) : self.having_elastic_interaction = jdata[ 'having_elastic_interaction' ]

        'set 9.self.max_box_int.'
        tmp_periodic = j_must_have( jdata, 'periodic' )

        if any( tmp_periodic ) : 
            self.max_box_int = np.array( j_must_have( jdata, 'max_box_image_int') )
        else :
            self.max_box_int = np.array( [ 0, 0, 0 ] )

        self.max_box_int = self.max_box_int * tmp_periodic

        'set 10.self.numb_defectObject.'
        self.numb_defectObject = len( defectObjectStr_list )

        if self.numb_defectObject == 0 :
            self.defectObject_list = [ ]
        else :
            'set 11.self.defectObject_list.'
            if isinstance( defectObjectStr_list[ 0 ], str ) :
                self.defectObject_list = [ create_defectObject_of_subclass( jdata, item ) for item in defectObjectStr_list ]
            else :
      
                '''
                'Test_cgzhang, for detecting variable assignment operation. [ item for item in list ] is the same as copy.copy( ).'
                self.defectObject_list = [ item for item in defectObjectStr_list ]
                '''
                self.defectObject_list = defectObjectStr_list

        ( in_box, delta_l ) = self.return_judge_in_box( jdata, self.max_box_int )

        for j in range( 3 ) :
            if not in_box[ j ] and not tmp_periodic[ j ]:
                raise RuntimeError( 'defectSystem not in box, please increase box length in %d direction with %2.3f nm' % ( j, delta_l[ j ] ) )

        self.remap_defectSystem( jdata )

        'very important ! set 12.self.linkcell after ( self.remap_defectSystem: position, position_scale, position_image_int ) has been set.'
        self.linkcell = BoxLinkcell( jdata, self.defectObject_list )

        '''
        'Test_cgzhang, the default value is set to True. when detecting delete-, add-, sub- defectObject, set False.'
        self.initial_recombine = False
        '''

        'set 13.self.initial_recombine.'
        self.initial_recombine = True
        if j_have( jdata, 'initial_recombine' ) : self.initial_recombine = jdata[ 'initial_recombine' ]

        'set related properties according to initial recombine or not.'
        self.set_initial_recombine_and_related_properties( jdata, tmp_ConstNumber )

    def set_initial_recombine_and_related_properties( self, 
                                                      jdata, 
                                                      tmp_ConstNumber ) :
        '''
        set related properties, such as :
        14. self.total_defects, 15. self.sum_numb_actions_list,
        reset 12. self.linkcell,
        16. self.ui_array, 17. utotal,
        18. self.rate_array, 19. self.size_rate_array, 
        20. self.sum_rate_array, 21. self.sum_rate, 
        22. self.delta_time, 23. self.which_defect, 24. self.what_action.
        '''

        if self.initial_recombine : 
            numb_recombine_results = self.numb_defectObject
            while numb_recombine_results != 0 :
                numb_defectObject_back = self.numb_defectObject
                self.set_initial_recombine( jdata )
                numb_recombine_results = numb_defectObject_back - self.numb_defectObject

                '''
                'Test_cgzhang, for detecting recombine results.'
                print( '# numb_recombine_results:', numb_recombine_results )
                '''

        'set 14.self.total_defects.'
        self.set_total_defects( )

        'set 15.self.sum_numb_actions_list.'
        self.set_sum_numb_actions_list( )

        '''
        'Test_cgzhang, set geometric center of defect system to the center of tmp_box for self.individual_cascade==True.
        no move_defectSystem when detecting add-, delete- or sub- defectObject.'
        '''
        if self.individual_cascade :
            tmp_move_vector = self.box * 0.5 - self.return_center_defectSystem( self.defectObject_list )
            self.move_defectSystem( jdata, self.defectObject_list, tmp_move_vector, tmp_ConstNumber )

        ( in_box, delta_l ) = self.return_judge_in_box( jdata, self.max_box_int )

        for j in range( 3 ) :
            if not in_box[ j ] and not tmp_periodic[ j ]:
                raise RuntimeError( 'DefectSystem not in box, please increase box length in %d direction with %2.3f nm' % ( j, delta_l[ j ] ) )

        self.remap_defectSystem( jdata )

        'very important ! reset 12.self.linkcell after ( self.remap_defectSystem: position, position_scale, position_image_int ) has been set.'
        self.linkcell.linked( self.defectObject_list )

        'set 16.ui_array for every defectObject and 17.self.utotal for defectSystem, including elastic energy and others.'
        self.set_ui_and_utotal( jdata, tmp_ConstNumber )

        'set 18.self.rate_array, 19.self.size_rate_array.'
        self.set_rate_array( jdata )

        'set 20.self.sum_rate_array, 21.self.sum_rate.'
        self.sum_rate_array = np.zeros( self.size_rate_array )
        self.set_sum_rate_array( )

        'set 22.self.delta_time, 23.self.which_defect, 24.self.what_action.'
        self.set_which_defect_action( )

    def set_initial_recombine( self, jdata ) :
        '''
        Recombine defectObjects initially. set related properties, such as 
        8. self.defectObject_list, 9. self.numb_defectObject
        reset 10. self.linkcell.
        '''

        indicate = np.zeros( self.numb_defectObject, dtype = np.int32 )
        recombine_results_list = [ ]

        for i in range( self.numb_defectObject ) :
            if indicate[ i ] == 1 : continue
            index_list = self.return_27_cells_index_list( jdata, i )
        
            for j in index_list :
                if indicate[ j ] == 1 : continue
                tmp_judge = self.defectObject_list[ i ].return_recombine_results( jdata, self.defectObject_list[ j ] )
                if tmp_judge[ 0 ] :
        
                    indicate[ i ] = 1
                    indicate[ j ] = 1
        
                    if tmp_judge[ 1 ] == None :
                        pass
                    else :
                        recombine_results_list.append( tmp_judge[ 1 ] )

                        '''
                        print( '# disappear two defectObjects, and create a new !' )
        
                        print( '# ----- i -------' )
                        self.defectObject_list[ i ].print_properties( )
        
                        print( '# ----- j -------' )
                        self.defectObject_list[ j ].print_properties( )
        
                        print( '# ----- R -------' )
                        tmp_judge[ 1 ].print_properties( )
                        '''
        
                    break
        
        for i in range( self.numb_defectObject ) :
            if indicate[ i ] == 0 :
                recombine_results_list.append( self.defectObject_list[ i ] )
        
        'reset 8.self.defectObject_list.'
        self.defectObject_list = recombine_results_list
        
        'reset 9.self.numb_defectObject.'
        self.numb_defectObject = len( self.defectObject_list )
        
        self.remap_defectSystem( jdata )
        
        'very important ! reset 10.self.linkcell after ( self.remap_defectSystem: position, position_scale, position_image_int ) has been set.'
        self.linkcell.linked( self.defectObject_list )

    def set_ui_and_utotal( self, 
                           jdata, 
                           tmp_ConstNumber ) :
        '''
        set ui for every defectObject and utotal for defectSystem, including elastic energy.
        here ui = 0.5 * sigma_( j != i )^N U_elastic_( ij ), j is the neighbours of i.
        utotal = sigma_(i=1)^N ui.
        '''

        self.utotal = 0.0
        self.ui_array = np.zeros( self.numb_defectObject, dtype = np.float64 )

        if not self.having_elastic_interaction : return

        for index in range( self.numb_defectObject ) :
            ui = self.return_ui_index( jdata, tmp_ConstNumber, index )
            self.ui_array[ index ] = ui
            self.utotal = self.utotal + ui

    def set_total_defects( self ) :
        'compute total point defects.'

        self.total_defects = { 'I' : 0, 'V': 0 }

        for item in self.defectObject_list :
            IV = item.defecttype[ 0 ]
            self.total_defects[ IV ] = self.total_defects[ IV ] + item.nsize

    def set_sum_numb_actions_list( self ) :
        'set self.sum_numb_actions_list.'

        self.sum_numb_actions_list = [ ]

        tmp_sum = 0

        for item in self.defectObject_list :
            tmp_sum += item.return_numb_actions( )
            self.sum_numb_actions_list.append( tmp_sum )

    def _return_m_n( self, index ) :
        'return ( m, n ), m : which_defect, n : what_action, index in rate_array.'

        if index < 0 or index >= self.size_rate_array :
            raise RuntimeError( 'defectObject index: %d is not in range [ 0, %d ]' %( index, self.size_rate_array ) )

        if self.sum_numb_actions_list[ 0 ] > index : return ( 0, index )

        for i in range( 1, self.numb_defectObject ) :
            if self.sum_numb_actions_list[ i ] > index :
                return ( i, index - self.sum_numb_actions_list[ i - 1 ] )

    def _return_begin_end_index( self, index_defectObject ) :
        'return begin index and end index in rate array, index_defectObject in self.defectObject_list.'

        if index_defectObject < 0 or index_defectObject >= self.numb_defectObject :
            raise RuntimeError( 'defectObject index is out of range [ 0, %d ]' %( self.numb_defectObject ) )

        if index_defectObject == 0 : return ( 0, self.sum_numb_actions_list[ 0 ] )

        return ( self.sum_numb_actions_list[ index_defectObject - 1 ], self.sum_numb_actions_list[ index_defectObject ] )

    def return_rate_defectObject( self, 
                                  jdata, 
                                  index  ) :
        'return rate of defectObject with index in self.defectObject_list.'

        index_list = self.return_27_cells_index_list( jdata, index )
        defectObject_list_27_cells = self.return_defectObject_list_from_index( index_list )

        rate_defectObject = self.defectObject_list[ index ].return_rate( jdata, self.temperature_evolution, self.ui_array[ index ], defectObject_list_27_cells )

        return rate_defectObject

    def return_rate_defectObject_test( self, 
                                       jdata, 
                                       index  ) :
        'return rate of defectObject with index in self.defectObject_list.'

        index_list = self.return_27_cells_index_list( jdata, index )
        defectObject_list_27_cells = self.return_defectObject_list_from_index( index_list )

        rate_defectObject = self.defectObject_list[ index ].return_rate_test( jdata, self.temperature_evolution, defectObject_list_27_cells )

        return rate_defectObject

    def set_rate_defectObject_index( self, 
                                     jdata, 
                                     index  ) :
        'reset rate of defectObject with index in self.defectObject_list.'

        ( begin_index, end_index ) = self._return_begin_end_index( index )

        rate_array_defectObject = np.array( self.return_rate_defectObject( jdata, index ) )

        self.rate_array[ begin_index : end_index ] = rate_array_defectObject
   
    def set_rate_array( self, jdata ) :
        'set rate array for all defectObject in self.defectObject_list.'

        rate_list = [ ]

        for i in range( self.numb_defectObject ) :
            rate_defectObject = self.return_rate_defectObject( jdata, i )
            rate_list.extend( rate_defectObject )

        self.rate_array = np.array( rate_list )
        self.size_rate_array = self.rate_array.size

        ''' 
        'Test_cgzhang, for detecting self.rate_array.'
        print( '# set_rate_array' )
        print( '# ( self.rate_array, self.size_rate_array ): ', ( self.rate_array, self.size_rate_array ) )   
        '''

    def set_sum_rate_array( self, 
                            begin_index_of_defectObject = 0 ) :
        '''
        set self.sum_rate_array and self.sum_rate.
        a. if self.numb_defectObject == 0, self.sum_rate_array=np.array([]), self.sum_rate = 0.0.
        b. if self.numb_defectObject != 0, set these two properties from begin_index_of_sum_rate_array to the end index of sum_rate_array.
        '''

        if self.numb_defectObject == 0 :
            self.sum_rate_array = np.array( [ ] )
            self.sum_rate = 0.0
            return

        if begin_index_of_defectObject == 0 : 
            self.sum_rate_array[ -1 ] = 0.0
            begin_index_of_sum_rate_array = 0
        else :
            begin_index_of_sum_rate_array = self.sum_numb_actions_list[ begin_index_of_defectObject - 1 ]

        for i in range( begin_index_of_sum_rate_array, self.size_rate_array ) :
            self.sum_rate_array[ i ] = self.sum_rate_array[ i - 1 ] + self.rate_array[ i ]

        self.sum_rate = self.sum_rate_array[ -1 ]

        '''
        'Test_cgzhang, for detecting self.sum_rate_array.'
        print( '# efficiency:', begin_index_of_sum_rate_array / self.size_rate_array )
        print( '# set_sum_rate_array' )
        print( '# ( self.sum_rate_array, self.sum_rate ):', ( self.sum_rate_array, self.sum_rate ) )   
        '''

    def set_which_defect_action( self ) :
        '''
        set self.which_defect, self.what_action and self.delta_time.
        a. if self.numb_defectObject == 0, self.which_defect = -1, self.what_action = -1, self.delta_time = 0.0.
        b. if self.numb_defectObject != 0, set these three properties as following algorithm.
        '''

        '''
        'Test_cgzhang, for detecting self.delta_time, self.which_defect and self.what_action.'
        seed = 777
        np.random.seed( seed )
        '''

        if self.numb_defectObject == 0 :
            self.which_defect = -1
            self.what_action = -1
            self.delta_time = +math.inf
            return

        kesi_1 = np.random.rand( ) * self.sum_rate

        kesi_2 = np.random.rand( )
        self.delta_time = -math.log( kesi_2, math.e ) / self.sum_rate

        '''
        'Another Algorithm :'

        for i in range( self.size_rate_array ) :
            if self.sum_rate_array[ i ] >= kesi_1 :
                ( self.which_defect, self.what_action ) = self._return_m_n( i )

                'Test_cgzhang, for detecting this algorithm.'
                which_defect_1, what_action_1 = self._return_m_n( i )
                #print( '# set_which_defect_action real:' )
                #print( '# end_i real:', i )
                #print( '# ( kesi_1, self.delta_time ): ', ( kesi_1, self.delta_time ) )
                #print( '# ( self.which_defect, self.what_action ): ', ( self.which_defect, self.what_action ) )

                break
        '''        

        if kesi_1 <= self.sum_rate_array[ 0 ] :
            self.which_defect = 0
            self.what_action = 0
        else :
            begin_i = 0
            end_i = self.size_rate_array - 1
       
            while True :
                half_i = int( 0.5 * ( begin_i + end_i ) )
                half = self.sum_rate_array[ half_i ]
       
                if begin_i == half_i :
                    ( self.which_defect, self.what_action ) = self._return_m_n( end_i )

                    'Test_cgzhang, for detecting this algorithm.'
                    #which_defect_2, what_action_2 = self._return_m_n( end_i )

                    #if not judge_equal( which_defect_1, which_defect_2, 1.0E-3 ) : print( '# which_defect not equal!' )
                    #if not judge_equal( what_action_1, what_action_2, 1.0E-3 ) : print( '# what_action not equal!' )

                    #print( '# end_i:', end_i )
                    #print( '# ( kesi_1, self.delta_time ): ', ( kesi_1, self.delta_time ) )
                    #print( '# ( self.which_defect, self.what_action ): ', ( self.which_defect, self.what_action ) )
                   
                    'return is important!' 
                    return
       
                if kesi_1 >  half : begin_i = half_i
                if kesi_1 <= half : end_i = half_i

    def sys_take_action( self, 
                         jdata, 
                         tmp_trapSys = None ) :
        '''
        return new_defectObject_list after defect_take_action( ), 
        new_defectObject_list[ 0 ] substitute index: self.which_defect.
        new_defectObject_list[ 1 ] emitted and be added, when len( new_defectObject_list ) == 2.

        return None or first_delete, which is used to the function sys_take_recombine().
        '''

        'return 0 when self.numb_defectObject eq 0.'
        if self.numb_defectObject == 0 : return 0

        new_defectObject_list = self.defectObject_list[ self.which_defect ].defect_take_action( jdata, self.what_action )

        self.substitute_defectObject( jdata, self.which_defect, new_defectObject_list[ 0 ] )

        'judge the first defectObject out of box and deleted, and then return the deleted defectObject.'
        first_delete = self.deleted_due_2_out_max_box( jdata, self.which_defect )

        if tmp_trapSys != None and first_delete == None : 
            'set self.trap of the DefectObject that substituded.'
            self.defectObject_list[ self.which_defect ].reset_trap_defectObject( tmp_trapSys.return_judge_trap( new_defectObject_list[ 0 ] ) )

        if len( new_defectObject_list ) == 2 :
            'len() == 2, take emitting, then add the emited defectObject.'
            self.add_defectObject( jdata, new_defectObject_list[ 1 ] )

            index_second_delete = self.numb_defectObject - 1
            second_delete = self.deleted_due_2_out_max_box( jdata, index_second_delete )

            if tmp_trapSys != None and second_delete == None :
                'set self.trap of the DefectObject that added.'
                self.defectObject_list[ index_second_delete ].reset_trap_defectObject( tmp_trapSys.return_judge_trap( new_defectObject_list[ 1 ] ) )

        return first_delete

    def sys_take_recombine( self, jdata ) :
        '''
        take recombine:
        1. tmp_judge[ 0 ] == False, nothing to do
        2. tmp_judge[ 0 ] == True and tmp_judge[ 1 ] == None, two defectObjects disapear, 
                                                              delete defectObject for index: [ self.which_defect, item ]
        3. tmp_judge[ 0 ] == True and tmp_judge[ 1 ] == new_defectObject, new_defectObject is created, 
                                                              substitude new_defectObject for index: self.which_defect 
                                                              and delete defectObject for index: item
        '''

        'self.which_defect is not in the range [ 0 : self.numb_defectObject ].'
        if self.which_defect >= self.numb_defectObject or self.which_defect < 0 : 
            print( '# self.which_defect %d not in the range [0:%d]' %( self.which_defect, self.numb_defectObject ) )

        index_list = self.return_27_cells_index_list( jdata, self.which_defect )

        for item in index_list :

            tmp_judge = self.defectObject_list[ self.which_defect ].return_recombine_results( jdata, self.defectObject_list[ item ] )

            if tmp_judge[ 0 ] :

                # test_cgzhang, for recombine detecting
                print( '#--------------------------- Take recombination ---------------------------' )

                print( '# react-1:' )
                self.defectObject_list[ self.which_defect ].print_properties( )

                print( '# react-2:' )
                self.defectObject_list[ item ].print_properties( )

                if tmp_judge[ 1 ] != None :
                    print( '# product:' )
                    tmp_judge[ 1 ].print_properties( )

                if tmp_judge[ 1 ] == None :
                    self.delete_defectObject( jdata, [ self.which_defect, item ] )
                else :
                    self.substitute_defectObject( jdata, self.which_defect, tmp_judge[ 1 ] )
                    self.delete_defectObject( jdata, item )

                break

    def add_defectObject( self, 
                          jdata, 
                          tmp_defectObject_list, 
                          position_random = None ) :
        '''
        Add defectObject with the type of input parameter : 
         1. str
         2. DefectObject
         3. list of str
         4. list of DefectObject
        position_random == None, add defectObject with real position.
        position_random == 'xyz', the center of tmp_defectObject_list with random positions in tmp_box.
        position_random == 'xy', the center of tmp_defectObject_list with random positions in xy direction in tmp_box, and z is the real position.
        '''

        tmp_ConstNumber = ConstNumber( jdata )
 
        if isinstance( tmp_defectObject_list, str ) :

            tmp_defectObject_list = [ create_defectObject_of_subclass( jdata, tmp_defectObject_list ) ]
       
        elif isinstance( tmp_defectObject_list, DefectObject ) :

            tmp_defectObject_list = [ tmp_defectObject_list ]
      
        elif isinstance( tmp_defectObject_list, list ) and isinstance( tmp_defectObject_list[ 0 ], str ) :

            tmp_defectObject_list = [ create_defectObject_of_subclass( jdata, tmp_defectObject ) for tmp_defectObject in tmp_defectObject_list ]
      
        if position_random != None :

            tmp_move_vector = np.random.rand( 3 ) * self.box - self.return_center_defectSystem( tmp_defectObject_list )
            if position_random == 'xy' : tmp_move_vector[ 2 ] = 0.0
 
            self.move_defectSystem( jdata, tmp_defectObject_list, tmp_move_vector, tmp_ConstNumber )
      
        self._add_defectObject_list( jdata, tmp_defectObject_list )

    def _add_defectObject_list( self, 
                                jdata,
                                tmp_defectObject_list ) :
        '''
        Add a list of defectObject.
        No changed roperties of DefectSystem :
         1. self.c_step
         2. self.c_time
         3. self.time_evolution
         4. self.temperature_evolution
         5. self.individual_cascade
         6. self.multi_cascade
         7. self.max_box_int
         11. self.initial_recombine

        Changed roperties of DefectSystem :
         8. self.numb_defectObject
         9. self.defectObject_list
         10. self.linkcell in self.linkcell.add_defectObject_list( tmp_defectObject_list ).
         12. self.total_defects
         13. self.sum_numb_actions_list
         14. self.ui_array
         15. self.utotal

         16. self.rate_array in self.set_rate_defectObject_index( self, jdata, index ) with effect_index_list.
         17. self.size_rate_array

         18. self.sum_rate_array
         19. self.sum_rate
         18-19 in self.set_sum_rate_array( self, begin_index_of_defectObject ) not in this function.

         20. self.delta_time
         21. self.which_defect
         22. self.what_action
         20-22 in self.set_which_defect_action( self ) not in this function.
        '''

        if len( tmp_defectObject_list ) == 0 :
            raise RuntimeError( 'Add no defectObject when trying to add defectObject.' )

        tmp_ConstNumber = ConstNumber( jdata )

        'reset 10.self.linkcell.'
        self.linkcell.add_defectObject_list( tmp_defectObject_list )

        'reset 9.self.defectObject_list.'
        self.defectObject_list.extend( tmp_defectObject_list )

        old_numb_defectObject = self.numb_defectObject

        'reset 8.self.numb_defectObject.'
        self.numb_defectObject = len( self.defectObject_list )

        'reset 12.self.total_defects.'
        for item in tmp_defectObject_list :
            IV = item.defecttype[ 0 ]
            self.total_defects[ IV ] = self.total_defects[ IV ] + item.nsize

        tmp_sum = self.sum_numb_actions_list[ -1 ]
        effect_index_list = [ ]

        for i in range( len(  tmp_defectObject_list ) ) :

            real_index = i + old_numb_defectObject

            effect_index_list.extend( self.return_effect_index_list( jdata, real_index ) )

            tmp_sum += self.defectObject_list[ real_index ].return_numb_actions( )

            'reset 13.self.sum_numb_actions_list.'
            self.sum_numb_actions_list.append( tmp_sum )

        numb_actions_add = self.sum_numb_actions_list[ -1 ] - self.sum_numb_actions_list[ old_numb_defectObject - 1 ]

        'add ui_array space with length of ( self.numb_defectObject - old_numb_defectObject ).'
        self.ui_array = np.concatenate( ( self.ui_array, np.zeros( ( self.numb_defectObject - old_numb_defectObject ), dtype = np.float64 ) ), axis = 0 )

        'add numb_actions_add space for self.rate_array and self.sum_rate_array.'
        self.rate_array = np.concatenate( ( self.rate_array, np.zeros( ( numb_actions_add ), dtype = np.float64 ) ), axis = 0 )
        self.sum_rate_array = np.concatenate( ( self.sum_rate_array, np.zeros( ( numb_actions_add ), dtype = np.float64 ) ), axis = 0 )

        'reset 17.self.size_rate_array.'
        self.size_rate_array = self.size_rate_array + numb_actions_add

        effect_index_list = list( set( effect_index_list ) )
        effect_index_list.extend( list( range( old_numb_defectObject, self.numb_defectObject ) ) )
        'for set_sum_rate_array, minimum element of effect_index_list.'
        begin_index_of_defectObject = min( effect_index_list ) 

        'reset 14.self.ui_array and 15.self.utotal.'
        for item in effect_index_list :
            ui = self.return_ui_index( jdata, tmp_ConstNumber, item )
            self.utotal = self.utotal + ui - self.ui_array[ item ]
            self.ui_array[ item ] = ui

            'reset 16.self.rate_array.'
            self.set_rate_defectObject_index( jdata, item )

        'reset 18.self.sum_rate_array and 19.self.sum_rate.'
        self.set_sum_rate_array( begin_index_of_defectObject )

    def return_ui_index( self, 
                         jdata, 
                         tmp_ConstNumber, 
                         tmp_index        ) :
        '''
        Here ui = 0.5 * sigma_( j != i )^N U_elastic_( ij ), j is the neighbours of i.
        return ui of defectObject with tmp_index of self.defectObject_list.
        '''

        if not self.having_elastic_interaction : return 0.0

        effect_index_list = self.return_effect_index_list( jdata, tmp_index  )
        effect_defectObject_list = self.return_defectObject_list_from_index( effect_index_list )

        return 0.5 * self.defectObject_list[ tmp_index ].return_u_neighbour( jdata, tmp_ConstNumber, effect_defectObject_list )

    def delete_defectObject( self, 
                             jdata, 
                             index_list ) :
        '''
        Delete defectObject, the type of input parameter :
         1. int
         2. list of int 
        '''

        if self.numb_defectObject == 0 :
            raise RuntimeError( 'Try to delete defectObject when numb_defectObject eq 0 in the defectSystem !' )
        if not isinstance( index_list, Iterable ) : index_list = [ index_list ]

        return self._delete_defectObject_list( jdata, index_list )

    def _delete_defectObject_list( self, 
                                   jdata, 
                                   index_list ) :
        '''
        Delete defectObject_list, 
        No changed properties of DefectSystem :
         1. self.c_step
         2. self.c_time
         3. self.time_evolution
         4. self.temperature_evolution
         5. self.individual_cascade
         6. self.multi_cascade
         7. self.max_box_int
         11. self.initial_recombine

        Changed properties of DefectSystem :
         8. self.numb_defectObject
         9. self.defectObject_list
         10. self.linkcell in self.linkcell.delete_defectObject_list( tmp_delete_defectObject_list, index_list ).
         12. self.total_defects
         13. self.sum_numb_actions_list in self.delete_sum_numb_actions_list_index( item ).
         14. self.ui_array
         15. self.utotal

         16. self.rate_array in self.set_rate_defectObject_index( self, jdata, index ).
         17. self.size_rate_array

         18. self.sum_rate_array
         19. self.sum_rate
         18-19 in self.set_sum_rate_array( self, begin_index_of_defectObject ).

         20. self.delta_time
         21. self.which_defect
         22. self.what_action
         20-22 in self.set_which_defect_action( self ) not in this function.
        '''

        if len( index_list ) < 1 :
            raise RuntimeError( 'Delete no defectObject in delete_defectObject_list( ) of Class DefectSystem.' )

        tmp_ConstNumber = ConstNumber( jdata )

        'delete repeated elements.'
        index_list = list( set( index_list ) )

        'reset 12.self.total_defects.'
        for item in index_list :
            IV = self.defectObject_list[ item ].defecttype[ 0 ]
            self.total_defects[ IV ] = self.total_defects[ IV ] - self.defectObject_list[ item ].nsize

        'get effect_index_list in old self.defectObject_list index.'
        old_effect_index_list = [ ]

        for item in index_list :
            old_effect_index_list.extend( self.return_effect_index_list( jdata, item ) )

        old_effect_index_list = list( set( old_effect_index_list ) )

        'old_effect_index delete index in index_list.'
        old_effect_index_list = [ item for item in old_effect_index_list if item not in index_list ]

        tmp_delete_defectObject_list = [ ]

        for i in index_list :
            if i < 0 or i > self.numb_defectObject - 1 :
                raise RuntimeError( 'Index %d out of range [ 0 : %d ] when trying to delete a defectObject' % ( i, self.numb_defectObject ) )

            tmp_delete_defectObject_list.append( self.defectObject_list[ i ] )

        'reset 10.self.linkcell.'
        self.linkcell.delete_defectObject_list( tmp_delete_defectObject_list, index_list )

        'decreasing index_list[ index ] according value.'
        index_list.sort( reverse = True )
 
        for item in index_list :

            ( begin_index, end_index ) = self._return_begin_end_index( item )

            'delete space of self.rate_array.'
            self.rate_array = np.delete( self.rate_array, slice( begin_index, end_index ) )

            'delete self.ui_array[ item ] and 15.reset self.utotal.'
            self.utotal = self.utotal - self.ui_array[ item ]
            self.ui_array = np.delete( self.ui_array, item )

            'delete space of self.sum_numb_actions_list, reset 13.self.sum_numb_actions_list.'
            self.delete_sum_numb_actions_list_index( item )

            'delete self.defectObject_list, reset 9.self.defectObject_list.'
            del self.defectObject_list[ item ]

            'set 8.self.numb_defectObject.'
            self.numb_defectObject -= 1

        'delete ( self.size_rate_array - self.rate_array.size ) spaces for self.sum_rate_array.'
        self.sum_rate_array = np.delete( self.sum_rate_array, slice( self.rate_array.size, self.size_rate_array ) )

        'set 17.self.size_rate_array.'
        self.size_rate_array = self.rate_array.size

        if len( old_effect_index_list ) == 0 :
            begin_index_of_defectObject = index_list[ -1 ]

        else :
            new_effect_index_list = self._return_new_index_list( old_effect_index_list, index_list )
            'for set_sum_rate_array, minimum element of new_effect_index_list.'
            begin_index_of_defectObject = min( min( new_effect_index_list ), index_list[ -1 ] )
        
            'reset 14.self.ui_array and 15.self.utotal.'
            for item in new_effect_index_list :
                ui = self.return_ui_index( jdata, tmp_ConstNumber, item )
                self.utotal = self.utotal + ui - self.ui_array[ item ]
                self.ui_array[ item ] = ui

                'reset 16.self.rate_array.'
                self.set_rate_defectObject_index( jdata, item )

        '''
        'Test_cgzhang, detecting rate array.'
        print( '# new_effect_index_list:', new_effect_index_list )
        print( '# begin_index_of_defectObject:', begin_index_of_defectObject )
        'self.set_rate_array( jdata )'
        '''

        'reset 18.self.sum_rate_array and 19.self.sum_rate.'
        self.set_sum_rate_array( begin_index_of_defectObject )

        "print( '# self.rate_array.size: %d, self.size_rate_array: %d' % ( self.rate_array.size, self.size_rate_array ) )"

        return tmp_delete_defectObject_list

    def delete_sum_numb_actions_list_index( self, index ) :
        'delete and reset sum_numb_actions_list after delete defectObject with index.'

        for i in range( self.numb_defectObject ) :
           if i > index : self.sum_numb_actions_list[ i ] -= self.defectObject_list[ index ].return_numb_actions( )

        del self.sum_numb_actions_list[ index ]

    def _return_new_index( self,
                           old_index,
                           delete_index_list ) :
        'return new_index after delete defectObjects with delete_index_list.'

        if len( delete_index_list ) == 0 : 
            return old_index

        'increasing according to value.'
        delete_index_list = sorted( delete_index_list, reverse = False )

        for i in range( len( delete_index_list ) ) :

            if old_index < delete_index_list[ i ] : 
                return old_index - i
            elif old_index == delete_index_list[ i ] :
                raise RuntimeError( 'Old_index in delete_index_list, wrong !' )
        
        return old_index - len( delete_index_list )

    def _return_new_index_list( self,
                                old_index_list, 
                                delete_index_list ) :
        'return_new_index_list after delete defectObjects with delete_index_list.'

        new_index_list = [ ]

        for item in old_index_list :
            new_index_list.append( self._return_new_index( item, delete_index_list ) )
       
        return new_index_list
   
    def substitute_defectObject( self, 
                                 jdata, 
                                 index, 
                                 tmp_defectObject ) :
        '''
        Substitute defectObject.
        The type of input parameter :
         1. str
         2. DefectObject

        No changed properties of DefectSystem :
         1. self.c_step
         2. self.c_time
         3. self.time_evolution
         4. self.temperature_evolution
         5. self.individual_cascade
         6. self.multi_cascade
         7. self.max_box_int
         8. self.numb_defectObject
         11. self.initial_recombine

        Changed properties of DefectSystem :
         9. self.defectObject_list
         10. self.linkcell in self.linkcell.substitute_defectObject( index, old_defectObject, tmp_defectObject ).
         12. self.total_defects
         13. self.sum_numb_actions_list
         14. self.ui_array
         15. self.utotal

         16. self.rate_array in self.set_rate_defectObject_index( self, jdata, index ).
         17. self.size_rate_array

         18. self.sum_rate_array
         19. self.sum_rate
         18-19 in self.set_sum_rate_array( self, begin_index_of_defectObject ).

         20. self.delta_time
         21. self.which_defect
         22. self.what_action
         20-22 in self.set_which_defect_action( self ) not in this function.
        '''

        if index > self.numb_defectObject - 1 or index < 0 :
            raise RuntimeError( 'Index %d out of range [ 0 : %d ] when trying to substitute a defectObject' % ( index, self.numb_defectObject ) )

        if isinstance( tmp_defectObject, str ) :
           tmp_defectObject = create_defectObject_of_subclass( jdata, tmp_defectObject )

        if not ( isinstance( tmp_defectObject, DefectObject ) ) :
            raise RuntimeError( 'Wrong data type !' )

        tmp_ConstNumber = ConstNumber( jdata )

        'reset 12.self.total_defects.'
        IV = self.defectObject_list[ index ].defecttype[ 0 ]
        self.total_defects[ IV ] = self.total_defects[ IV ] - self.defectObject_list[ index ].nsize

        IV = tmp_defectObject.defecttype[ 0 ]
        self.total_defects[ IV ] = self.total_defects[ IV ] + tmp_defectObject.nsize

        'get effect_index_list with old self.linkcell.'
        effect_index_list = [ ]
        effect_index_list.extend( self.return_effect_index_list( jdata, index ) )

        ''' 
        'Test_cgzhang, for detecting variable assignment operation.'
        old_defectObject = copy.copy( self.defectObject_list[ index ] )
        old_defectObject = copy.deepcopy( self.defectObject_list[ index ] )  
        '''

        old_defectObject = self.defectObject_list[ index ]

        'set 9.self.defectObject_list.'
        self.defectObject_list[ index ] = tmp_defectObject

        '''
        'Test_cgzhang, for detecting variable assignment operation.'
        print( 'id(old_defectObject)=', id( old_defectObject ) )
        print( 'id(self.defectObject_list(index))=', id( self.defectObject_list[index] ) )
        print( 'id(tmp_defectObject)=', id( tmp_defectObject ) )
        old_defectObject.print_properties( )
        self.defectObject_list[index].print_properties( )
        tmp_defectObject.print_properties( )
        '''

        delta_numb_actions = tmp_defectObject.return_numb_actions( ) - old_defectObject.return_numb_actions( )

        'get the begin_index and end_index in self.rate_array using old self.sum_numb_actions_list.'
        ( begin_index, end_index ) = self._return_begin_end_index( index )

        'set space and size of self.rate_array.'
        if delta_numb_actions > 0 : 

            'insert delta_numb_actions element in index: begin_index.'
            self.rate_array = np.insert( self.rate_array, np.ones( delta_numb_actions, dtype = np.int ) * begin_index, 0.0 )
            self.size_rate_array = self.size_rate_array + delta_numb_actions 

            'add delta_numb_actions space for self.sum_rate_array.'
            self.sum_rate_array = np.concatenate( ( self.sum_rate_array, np.zeros( ( delta_numb_actions ), dtype = np.float64 ) ), axis = 0 )

        elif delta_numb_actions < 0 :

            'delete abs( delta_numb_actions ) element in from index : begin_index to begin_index + abs( delta_numb_actions ).'
            self.rate_array = np.delete( self.rate_array, slice( begin_index, begin_index + abs( delta_numb_actions ) ) )
            'reset 17.self.size_rate_array.'
            self.size_rate_array = self.size_rate_array - abs( delta_numb_actions )
         
            'delete -delta_numb_actions space for self.sum_rate_array.'
            self.sum_rate_array = np.delete( self.sum_rate_array, slice( self.size_rate_array + delta_numb_actions, self.size_rate_array ) )

        'reset 13.self.sum_numb_actions_list.'
        for i in range( self.numb_defectObject ) :
            if i >= index : 
                self.sum_numb_actions_list[ i ] = self.sum_numb_actions_list[ i ] + delta_numb_actions

        'reset 10.self.linkcell.'
        self.linkcell.substitute_defectObject( index, old_defectObject, tmp_defectObject )

        'add effect_index_list due to new linkcell, and then add index into effect_index_list.'
        effect_index_list.extend( self.return_effect_index_list( jdata, index ) )
        effect_index_list.append( index )

        effect_index_list = list( set( effect_index_list ) )
        'for set_sum_rate_array, minimum element of effect_index_list.'
        begin_index_of_defectObject = min( effect_index_list )

        '''
        'Test_cgzhang, for detecting self.ui_array and self.utotal.'
        'a.'
        'self.set_ui_and_utotal( jdata, tmp_ConstNumber )'

        'b.'
        self.utotal = 0.0
        for item in range( self.numb_defectObject ) :
            ui = self.return_ui_index( jdata, tmp_ConstNumber, item )
            self.utotal = self.utotal + ui
            self.ui_array[ item ] = ui
        '''

        'reset 14.self.ui_array and 15.self.utotal.'
        for item in effect_index_list :
            ui = self.return_ui_index( jdata, tmp_ConstNumber, item )
            self.utotal = self.utotal + ui - self.ui_array[ item ]

            '''
            'Test_cgzhang, for detecting self.ui_array and self.utotal.'
            print( '# delta_ui', ui - self.ui_array[ item ] )
            '''

            self.ui_array[ item ] = ui
            'reset 16.self.rate_array.'
            self.set_rate_defectObject_index( jdata, item )

        '''
        'Test_cgzhang, for detecting rate array.'
        self.set_rate_array( jdata )
        print( '# effect_index_list:', effect_index_list )
        '''

        'reset 18.self.sum_rate_array and 19.self.sum_rate.'
        self.set_sum_rate_array( begin_index_of_defectObject )

    def return_defectObject_list_from_index( self, index_list ) :
        'return defectObject_list from index_list.'

        cells_defectObject_list = [ ]

        for i in index_list :
            cells_defectObject_list.append( self.defectObject_list[ i ] )

        return cells_defectObject_list

    def return_effect_index_list( self, 
                                  jdata, 
                                  index  ) :
        'return index list taking effect, including the rate of emitting point defect, migration, rotation and transformation.'

        cells_defectObject_index_list = self.return_27_cells_index_list( jdata, index )

        return cells_defectObject_index_list

        '''
        'Sth. may be wrong when getting the rate of emitting point defect 
        because the emitted point defect having interaction with the changed defectObject( with index ) 
        while the emitting(body) defect having no interaction with changed defectObject.'

        effect_index_list = [ ]

        for item in cells_defectObject_index_list :
            if self.defectObject_list[ index ].return_effect_or_not( jdata, self.defectObject_list[ item ] ) :
                effect_index_list.append( item )

        return effect_index_list
        '''

    def return_27_cells_index_list( self, 
                                    jdata, 
                                    index  ) :
        'return 27_cells_defectObject_index_list in 27 cells of the index: index, no including its self index.'

        cell_index = self.linkcell.return_cell_index_defectObject( self.defectObject_list [ index ] )
        ip = self.linkcell.return_ip_cell_index( cell_index )

        i = self.linkcell.ltop[ ip ]
   
        cells_defectObject_index_list = [ ]

        'central cell.'
        if i != -1 :
            while True :
                if i != index :
                    cells_defectObject_index_list.append( i )
                i = self.linkcell.linkmp[ i ]
                if i == -1 : break

        for kc in range( 1, 27 ) :

            jx = cell_index[ 0 ] + self.linkcell.nix[ kc ]
            jy = cell_index[ 1 ] + self.linkcell.niy[ kc ]
            jz = cell_index[ 2 ] + self.linkcell.niz[ kc ]
           
            if cell_index[ 0 ] == self.linkcell.nlc_vector[ 0 ] - 1 and jx > cell_index[ 0 ] :
                jx = 0
                if not self.linkcell.periodic[ 0 ] or jx == cell_index[ 0 ] : break
            elif cell_index[ 0 ] == 0 and jx < cell_index[ 0 ] :
                jx = self.linkcell.nlc_vector[ 0 ] - 1
                if not self.linkcell.periodic[ 0 ] or jx == cell_index[ 0 ] : break
           
            if cell_index[ 1 ] == self.linkcell.nlc_vector[ 1 ] - 1 and jy > cell_index[ 1 ] :
                jy = 0
                if not self.linkcell.periodic[ 1 ] or jy == cell_index[ 1 ] : break
            elif cell_index[ 1 ] == 0 and jy < cell_index[ 1 ] :
                jy = self.linkcell.nlc_vector[ 1 ] - 1
                if not self.linkcell.periodic[ 1 ] or jy == cell_index[ 1 ] : break
           
            if cell_index[ 2 ] == self.linkcell.nlc_vector[ 2 ] - 1 and jz > cell_index[ 2 ] :
                jz = 0
                if not self.linkcell.periodic[ 2 ] or jz == cell_index[ 2 ] : break
            elif cell_index[ 2 ] == 0 and jz < cell_index[ 2 ] :
                jz = self.linkcell.nlc_vector[ 2 ] - 1
                if not self.linkcell.periodic[ 2 ] or jz == cell_index[ 2 ] : break

            jc = self.linkcell.return_ip_cell_index( [ jx, jy, jz ] )

            j = self.linkcell.ltop[ jc ]

            'Bypass this neighbouring cell if it is empty.'
            if j != -1 :
                while True :
                    cells_defectObject_index_list.append( j )
                    j = self.linkcell.linkmp[ j ]
                    if j == -1 : break

        if len( cells_defectObject_index_list ) > self.linkcell.nnbrs :
            raise RuntimeError( 'Max1: %d .gt. nnbrs: %d something wrong' % ( len( cells_defectObject_index_list ), self.linkcell.nnbrs ) )

        return cells_defectObject_index_list

    def return_center_defectSystem( self, tmp_defectObject_list ) :
        'return the center position of tmp_defectObject_list. when len(tmp_defectObject_list) eq 0, return np.zeros( ( 3 ), dtype = np.float ).'

        tmp_center_defectSystem = np.zeros( ( 3 ), dtype = np.float )
        tmp_numb_defectObject = len( tmp_defectObject_list )

        if tmp_numb_defectObject == 0 : 
            return tmp_center_defectSystem 

        for tmp_defectObject in tmp_defectObject_list :
            tmp_center_defectSystem += tmp_defectObject.return_position( )

        tmp_center_defectSystem /= tmp_numb_defectObject

        return tmp_center_defectSystem

    def move_defectSystem( self, 
                           jdata,
                           tmp_defectObject_list,
                           tmp_move_vector,
                           tmp_ConstNumber         ) :
        'move defectObjects of tmp_defectObject_list with the vector of tmp_move_vector.'

        for tmp_defectObject in tmp_defectObject_list :
            tmp_defectObject.move_defectObject( jdata, tmp_move_vector, tmp_ConstNumber )

    def return_judge_in_box( self,
                             jdata,
                             tmp_max_box_int ) :
        'return judge in box or not, and corresponding delta_l.'
 
        'add this line to solve the problem of no defectObject in defectSystem.' 
        delta_l = np.zeros( 3 )

        for item in self.defectObject_list :
            in_box, delta_l = item.return_judge_in_box( jdata, tmp_max_box_int )
            if not all( in_box ) : return ( in_box, delta_l )

        return ( [ True, True, True ], delta_l )

    def remap_defectSystem( self, jdata ) :
        'remap the defectSystem, i.e. set_position_scale_image_int.'

        for item in self.defectObject_list :
            item.set_position_scale_image_int( jdata )

    def reset_trap_defectSystem( self, 
                                 jdata, 
                                 tmp_trapSys ) :
       'reset trap of defectSystem.'

       for item in self.defectObject_list :
          item.reset_trap_defectObject( tmp_trapSys.return_judge_trap( jdata, item ) )
       
    def deleted_due_2_out_max_box( self, 
                                   jdata, 
                                   index  ) :
        'delete self.defectObject_list[ index ] because it migrates out of maximum simulation box, return the deleted defectObject.'
        'when self.defectObject_list[ index ] out of box delete it, and then return the deleted defectObject. if in the box, then return None.'

        delete_ = self.defectObject_list[ index ]

        if not all( delete_.return_judge_in_box( jdata, self.max_box_int )[ 0 ] ) :

            print( '# ---------- delete defectObject because it migrates out of simulation box ----------' )
            delete_.print_properties( )

            self.delete_defectObject( jdata, index )
            
            return delete_

        return None

    def reset_iter( self ) :
        'reset self.c_step and self.c_time.'

        self.c_step += 1
        self.c_time += self.delta_time

    def set_zero_c_time( self ) :
        self.c_time = 0.0

    def _return_properties( self ) :
        '''
        Properties of DefectSystem :
         0. self.box
         1. self.c_step
         2. self.c_time
         3. self.temperature_evolution
         4. self.individual_cascade
         5. self.multi_cascade
         6. self.time_evolution
         7. self.step_evolution
         8. self.having_elastic_interaction
         9. self.max_box_int
         10. self.numb_defectObject
         11. self.defectObject_list
         12. self.linkcell
         13. self.initial_recombine
         14. self.total_defects
         15. self.sum_numb_actions_list
         16. self.ui_array
         17. self.utotal
         18. self.rate_array
         19. self.size_rate_array
         20. self.sum_rate_array
         21. self.sum_rate
         22. self.delta_time
         23. self.which_defect
         24. self.what_action
        '''

        properties_list = [ [ 'box'                       , self.box                        ],
                            [ 'c_step'                    , self.c_step                     ],  
                            [ 'c_time'                    , self.c_time                     ],
                            [ 'temperature_evolution'     , self.temperature_evolution      ],
                            [ 'individual_cascade'        , self.individual_cascade         ],
                            [ 'multi_cascade'             , self.multi_cascade              ],
                            [ 'time_evolution'            , self.time_evolution             ],
                            [ 'step_evolution'            , self.step_evolution             ],
                            [ 'having_elastic_interaction', self.having_elastic_interaction ],
                            [ 'max_box_int'               , self.max_box_int                ],
                            [ 'numb_defectObject'         , self.numb_defectObject          ],
                            [ 'initial_recombine'         , self.initial_recombine          ],
                            [ 'total_defects'             , self.total_defects              ],
                            [ 'sum_numb_actions_list'     , self.sum_numb_actions_list      ],
                            [ 'ui_array'                  , self.ui_array                   ],
                            [ 'utotal'                    , self.utotal                     ],
                            [ 'rate_array'                , self.rate_array                 ],
                            [ 'size_rate_array'           , self.size_rate_array            ],
                            [ 'sum_rate_array'            , self.sum_rate_array             ],
                            [ 'sum_rate'                  , self.sum_rate                   ],
                            [ 'delta_time'                , self.delta_time                 ],
                            [ 'which_defect'              , self.which_defect               ],
                            [ 'what_action'               , self.what_action                ]  ]

        return properties_list

    def equal( self, 
               tmp_defectSystem, 
               delta             ) :
        'judge equal between two objects of Class DefectSystem.'

        if not self.linkcell.equal( tmp_defectSystem.linkcell, delta ) :
            print( '# linkcell not equal !' )
            return False

        self_properties = self._return_properties( )
        tmp_properties = tmp_defectSystem._return_properties( )

        for i in range( len( self_properties ) ) :
            if not judge_equal( self_properties[ i ][ 1 ], tmp_properties[ i ][ 1 ], delta ) :
                print( '# ' + self_properties[ i ][ 0 ] + ' not equal !' )
                return False

        return True

    def print_properties( self, which_item = None ) :

        if which_item == None :
            print( '# ' )
            print( '# ------------- properties for DefectSystem --------------' )
            self.linkcell.print_properties( )
            print( '\n'.join( ( '# %s : %s' % ( item[ 0 ], item[ 1 ] ) for item in self._return_properties( ) ) ) )
            print( '# ' )
        else :
            for item in self._return_properties( ) :
                if item[ 0 ] == which_item : print( '# ', which_item, ':', item[ 1 ] )
 
    def evolution( self,
                   jdata,
                   outputFile,
                   tmp_trapSys = None ) :
        'the center of okmc simulation.'

        while self.c_time < self.time_evolution and self.c_step < self.step_evolution :

            if outputFile.is_output_now( self, jdata ) : 
                outputFile.output_results( self, jdata )

            'if first_delete == None, i.e. the defectObject with index of self.which_defect is in the box, therefore the function sys_take_recombine() must take place.'
            first_delete = self.sys_take_action( jdata, tmp_trapSys = None )
            if first_delete == None : self.sys_take_recombine( jdata )

            self.reset_iter( )

            self.set_which_defect_action( )

        if outputFile.is_output_now( self, jdata ) :
            outputFile.output_results( self, jdata )

if __name__ == '__main__' :

    '''
    1. plotting elastic energy ( utotal: eV ) among defectObjects from defectStr_list, 
    when a defectObject is moved with a vector, such as :
    1. void and ILoop111
    2. void and ILoop100
    3. void and VLoop111
    4. void and VLoop100
    5. V and ILoop111
    6. V and ILoop100
    7. V and VLoop111
    8. V and VLoop100

    #------------------------------------- anneal.json BEG ----------------------------------------

    # ------------------------------------ anneal.json END -------------------------------------
    '''
    parser = argparse.ArgumentParser( description = '--- Class DefectSystem detecting ---' )
    parser.add_argument( 'INPUT', help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    tmp_ConstNumber = ConstNumber( jdata )
    delta = 1.0e-2
    several_digits = 4

    #parameter_list = [ 0.1, 0.3, 0.5, 0.8, 1.0, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0 ]
    #parameter_list = [ 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0 ]
    #parameter_list = [ 7.0, 8.0, 9.0, 11.0, 12.0, 13.0, 14.0 ]
    parameter_list = [ 7.0 ]

    for r0 in parameter_list :

        'change global data_ value.'

        MyGlobalParameter.set_value( r0 )
        data_ = MyGlobalParameter.get_value( )
        #print( ' r0: %0.4f begin!' % data_ )
        #MyGlobalParameter.print_parameters( )
 
        dataStr_list = [ ]
 
        # size of Vn
        Vn_list = [ 9, 59, 137 ]
 
        for i in Vn_list :
 
            #tmp_i = ( i + 8 ) % 10
            #if tmp_i != 0 : continue
 
            t1 = time.time( )
 
            'defectObject_1, VCluster, size from 2 to 1000.'
            d1_defecttype = 'VCluster'
            d1_size = i
            'd1 position not in the boundary of simulation box.'
            d1_position = np.array( [ 0.1, 0.1, 0.1 ] )
            d1_burgers_vector_scale = None
            d1 = create_defectObject( jdata, d1_defecttype, d1_size, d1_position, d1_burgers_vector_scale )
 
            print( '# i, defecttype, nsize:', i, d1.defecttype, d1.nsize )
            print( 'properties of d1' )
            d1.print_properties( )

            'the index must not be changed.' 
            for k in range( 4, 1000 ) :
  
                'defectObject_2, ILoop111, size from 4 to 1000.'
                d2_defecttype = 'ILoop111'
                d2_size = k
                d2_position = np.array( [ 0.1, 0.1, 0.1 ] )
                'pay attention to this: burgers_vector_scale of Loop111 [ 0.5, 0.5, 0.5 ].'
                d2_burgers_vector_scale = np.array( [ 0.5, 0.5, 0.5 ] )
 
                d2 = create_defectObject( jdata, d2_defecttype, d2_size, d2_position, d2_burgers_vector_scale )
 
                #print( '# k, defecttype, nsize, burgers_vector_unit:', k, d2.defecttype, d2.nsize, d2.burgers_vector_unit )

                if k in [ 19, 43, 85 ] : d2.print_properties( )
  
                for j in range( 1, 500 ) :
  
                    defectObject_virtual_migrate = copy.deepcopy( d2 )
                    move_vector = j * tmp_ConstNumber.first * d2.burgers_vector_unit
                    defectObject_virtual_migrate.move_defectObject( jdata, move_vector, tmp_ConstNumber )
 
                    #print( 'j:', j, 'move_vector:', move_vector )
                    #defectObject_virtual_migrate.print_properties( )
 
                    defectObject_list = [ d1, defectObject_virtual_migrate ]
                    model = DefectSystem( jdata, defectObject_list )
 
                    #model.print_properties( 'numb_defectObject' )
 
                    #reset model.initial_recombine=True.
                    #model.initial_recombine = True
                    #model.set_initial_recombine_and_related_properties( jdata, tmp_ConstNumber )
 
                    'this means these two defectObjects recombine, therefore jump this data of distance and utotal.'
                    if model.numb_defectObject <= 1 : 
                        #print( '# recombine happens !' )
                        continue
 
                    'the distance between these two defectObjects.'
                    distance = d1.return_distance( jdata, defectObject_virtual_migrate )
 
                    ''' 
                    if j == 300 : 
                        if not judge_equal( 0.0, model.utotal, delta ) : print( '# utotal: %.7f when at the distance: %.7f nm.' % ( model.utotal, distance ) )
                    '''
 
                    if judge_equal( 0.0, model.utotal, delta ) : 
                        #print( '# utotal: %.7f when at the distance: %.7f nm, i: %d, k: %d, j: %d' % ( model.utotal, distance, i, k, j ) )
                        break
 
                    tmp_distance = save_float_several_digits( str( distance ), several_digits ) 
                    tmp_utotal = save_float_several_digits( str( model.utotal ), several_digits )
 
                    #tmp_distance = str( distance )
                    #tmp_utotal = str( model.utotal )
 
                    dataStr = tmp_distance + ' ' + tmp_utotal + ' ' + '\n'
                    dataStr_list.append( dataStr )
  
                    model.reset_iter( )
 
                    #print( '# move %d is OK !' % j )
 
                dataStr_list.append( '  \n' )
  
            #filenameStr = './utotal-ch/void-ILoop111/'+'V' + str( i ) + '/utotal-ch-' + str( data_ ) + '/utotal-V-' + str( d1_size ) + '.dat'
            filenameStr = './utotal-new/void-ILoop111/'+'V' + str( i ) + '/utotal' + '/utotal-V-' + str( d1_size ) + '.dat'
            output_data_into_file( dataStr_list, filenameStr )
 
            t2 = time.time( )
 
            print( '# r0: %0.4f, i: %d using time: %0.6f' % ( data_, i, t2 - t1 ) )
