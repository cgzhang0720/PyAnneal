#!/usr/bin/env python3

import numpy as np
import math
import argparse
import json
import time
import os

from Auxiliary import to_str

'''
'for detecting Class Cluster and ClusterSystem.'
import DefectSystem
from Trap import TrapSystem
import Cascade
from OutputFile import OutputFile
from Auxiliary import j_must_have, j_have
'''

class Cluster( object ) :

    def __init__( self,
                  tmp_alatt,
                  tmp_ICluster_upper,
                  tmp_VLoop_lower,
                  tmp_IV,
                  tmp_rcut,
                  tmp_position,
                  burgers_vector_unit_array, 
                  burgers_vector_scale_length_array, 
                  d_probability                      ) :
        '''
        Class Cluster set defectStr and related properties. 
        Properties:
         1. self.position
         2. self.rcut
         3. self.IV
         4. self.defecttype
         5. self.atoms_numb
         6. self.position_center
         7. self.burgers_vector_scale
         8. self.defectStr

        Methods:
         self.set_defectStr( self, alatt, ICluster_upper, VLoop_lower )
             1. self.return_sum_probability_in_burgers_range( self, tmp_alatt, burgers_vector_unit, burgers_vector_scale_length )
                 a. self.return_projection_atoms( self, tmp_alatt, burgers_vector_unit )
             2. self.find_index_of_burgers_vector( self, d_probability, probability_distri_array )
        Length units: nm.
        '''

        self.position = tmp_position
        self.rcut = tmp_rcut
        self.IV = tmp_IV
        self.atoms_numb = len( self.position )

        self.position_center = np.mean( self.position, axis = 0 )

        self.set_defectStr( tmp_alatt, tmp_ICluster_upper, tmp_VLoop_lower, burgers_vector_unit_array, burgers_vector_scale_length_array, d_probability )

    def set_defectStr( self,
                       alatt,
                       ICluster_upper,
                       VLoop_lower,
                       burgers_vector_unit_array,
                       burgers_vector_scale_length_array,
                       d_probability                        ) :
        '''
        set self.defecttype and self.burgers_vector_scale.
        a. burgers_vector_unit_array, b. burgers_vector_scale_length_array, c. d_probability.
        '''

        if burgers_vector_unit_array.shape[ 0 ] != burgers_vector_scale_length_array.size :
            raise RuntimeError( 'the size of burgers_vector_unit_array and burgers_vector_scale_length_array is not equal, wrong !' )

        '1. get sum probability of atoms belong to burgers range.'
        probability_distri = [ ]
        for i in range( burgers_vector_unit_array.shape[ 0 ] ) :
            probability_distri.append( self.return_sum_probability_in_burgers_range( alatt, burgers_vector_unit_array[ i ], burgers_vector_scale_length_array[ i ], i ) )

        probability_distri_array = np.array( probability_distri )

        '''
        # test_cgzhang
        print( '# self.atoms_numb:', self.atoms_numb )
        print( '# self.position:', self.position )
        print( '# self.position_center:', self.position_center )
        print( '# self.IV:', self.IV )
        print( '# self.rcut:', self.rcut )
        print( '# probability_distri_array:', probability_distri_array )
        '''

        probability_distri_array = probability_distri_array / np.sum( probability_distri_array )

        '2. get find_index_list which is relate to selected burgers vector unit.'
        find_index_list = self.find_index_of_burgers_vector( d_probability, probability_distri_array )
        length_find_index_list = len( find_index_list )

        '3. set self.defecttype and self.burgers_vector_scale according to length_find_index_list.'
        if length_find_index_list == 1 :
            if find_index_list[ 0 ] >= 0 and find_index_list[ 0 ] < 3 : 
                self.defecttype = self.IV + 'Loop100'
            else :
                'find_index_list[ 0 ] >= 3 and find_index_list[ 0 ] < probability_distri_array.size'
                self.defecttype = self.IV + 'Loop111'

            self.burgers_vector_scale = burgers_vector_unit_array[ find_index_list[ 0 ] ] * burgers_vector_scale_length_array[ find_index_list[ 0 ] ]
            burgers_set = True
        else :
            'length_find_index_list == 0'
            self.defecttype = self.IV + 'Cluster'
            burgers_set = False

        '4. reset self.defecttype according to self.IV, self.atoms_numb, ICluster_upper and VLoop_lower.'
        if self.IV == 'I' :

            if self.atoms_numb == 1 :
                self.defecttype = self.IV
            elif self.atoms_numb > 1 and self.atoms_numb <= ICluster_upper :
                self.defecttype = 'ICluster'
                burgers_set = False
            else :
                'transite ICluster to ILoop111 or ILoop100, because self.atoms_numb > ICluster_upper.'
                if burgers_set :
                    pass
                else :
                    self.defecttype = 'ILoop111'

        if self.IV == 'V' :

            if self.atoms_numb == 1 : 
                self.defecttype = self.IV
            elif self.atoms_numb > 1 and self.atoms_numb < VLoop_lower :
                self.defecttype = 'VCluster'
                burgers_set = False

        '5. set self.defectStr according to burgers_set.'
        if burgers_set :
            self.defectStr = self.defecttype + ' ' + str( self.atoms_numb ) + ' ' + to_str( self.position_center ) + ' ' + to_str( self.burgers_vector_scale ) + '\n'
        else :
            self.defectStr = self.defecttype + ' ' + str( self.atoms_numb ) + ' ' + to_str( self.position_center ) + '\n'

        '''
        # Test_cgzhang
        if length_find_index_list == 1 : 
            raise RuntimeError( '# set defecttype wrong !' )

        if self.atoms_numb != 1 and self.defecttype != 'VCluster' :
            raise RuntimeError( '# set defecttype wrong !' )

        print( '#-------------- cluster begin ! ------------------------------' )
        print( '# atoms_numb:', self.atoms_numb )
        print( '# probability_distri_array:', probability_distri_array )
        print( '# find_index_list:', find_index_list )
        print( '# length_find_index_list:', length_find_index_list )
        print( '# burgers_set:', burgers_set )
        print( '# defecttype:', self.defecttype )
        print( '# defectStr:', self.defectStr )
        '''

    def return_sum_probability_in_burgers_range( self,
                                                 tmp_alatt,
                                                 burgers_vector_unit,
                                                 burgers_vector_scale_length,
                                                 ii                            ) :
        '''
        Return sum probability of atoms in the range of burgers width, 
        with the unit vector of burgers: burgers_vector_unit 
        and the lattice length of burgers: burgers_vector_scale_length.
        units: lattice.
        '''

        max_projection, min_projection, projection_atoms = self.return_projection_atoms( tmp_alatt, burgers_vector_unit )

        '''
        # test_cgzhang
        print( '# projection_atoms:', projection_atoms, 'i:', ii )
        '''

        'delta_projection = max_projection - min_projection'

        if self.atoms_numb == 1 :
            '''
            sum_numb_atoms = np.array( [ 1.0 ] )
            projection_distribution_outputstr_list = [ '0.0 1.0 \n', '   \n' ]
            '''

            prob = 1.0

        else :

            numb_atoms_in_burgers_range = 0

            for projection in projection_atoms :
                if abs( projection ) < burgers_vector_scale_length * 0.5 :
                    numb_atoms_in_burgers_range = numb_atoms_in_burgers_range + 1

            prob = numb_atoms_in_burgers_range / self.atoms_numb


            '''
            '1nn * 0.5 = 3.0**0.5 * 0.5 * 0.5' 
            delta_projection_dr = 0.866025404 * 0.5
            sum_numb_atoms_size = int( delta_projection / delta_projection_dr ) + 1

            sum_numb_atoms = np.zeros( sum_numb_atoms_size, dtype = np.int32 )

            for projection in projection_atoms :
                index_sum_numb_atoms = int( ( projection - min_projection ) / delta_projection_dr )

                if index_sum_numb_atoms < 0 or index_sum_numb_atoms >= sum_numb_atoms_size :
                    raise RuntimeError( 'index %d out of in range [ %d, %d ]' %( index_sum_numb_atoms, 0, sum_numb_atoms_size ) )
                sum_numb_atoms[ index_sum_numb_atoms ] = sum_numb_atoms[ index_sum_numb_atoms ] + 1
        
            sum_numb_atoms = sum_numb_atoms / self.atoms_numb
       
            # Test_cgzhang, for getting output data str list
            projection = min_projection + delta_projection_dr * 0.5
            projection_distribution_outputstr_list = [ ]
         
            for item in sum_numb_atoms :
         
                projection_distribution_outputstr_list.append( str( projection ) + ' ' + str( item ) + '\n' )
                projection = projection + delta_projection_dr
         
            projection_distribution_outputstr_list.append( '    \n' )
            '''

        '''
        # Test_cgzhang, for projection distribution
        fileout_projection = './projection/projection_distri.' + str( ii ) + '.dat'

        if os.path.exists( fileout_projection ) :
            fileout_projection_distribution = open( fileout_projection, 'a' )
        else :
            fileout_projection_distribution = open( fileout_projection, 'w' )
        
        fileout_projection_distribution.writelines( projection_distribution_outputstr_list )
        
        fileout_projection_distribution.close( )
        '''

        return prob

    def return_projection_atoms( self, 
                                 tmp_alatt, 
                                 burgers_vector_unit ) :
        '''
        Return max_projection, min_projection, projection_atoms 
        between relative positions to self.position_center and burgers_vector_unit.
        units: lattice.
        '''

        max_projection = -math.inf
        min_projection = +math.inf
        projection_atoms = [ ]

        for item in self.position :

            relative_position = item - self.position_center
            projection = np.inner( relative_position, burgers_vector_unit ) / tmp_alatt

            if projection > max_projection : max_projection = projection
            if projection < min_projection : min_projection = projection

            projection_atoms.append( projection )

        return ( max_projection, min_projection, projection_atoms )

    def find_index_of_burgers_vector( self, 
                                      d_probability, 
                                      probability_distri_array ) :
        'Return find_index_list according to item > d_probability, item in probability_distri_array.'

        find_index_list = [ ]

        for i in range( probability_distri_array.size ) :
            if probability_distri_array[ i ] > d_probability :
                find_index_list.append( i )

        return find_index_list

    def return_position( self ) :
        return self.position

    def return_rcut( self ) :
        return self.rcut

    def return_IV( self ) :
        return self.IV

    def return_defecttype( self ) :
        return self.defecttype

    def return_atoms_numb( self ) :
        return self.atoms_numb

    def return_position_center( self ) :
        return self.position_center

    def return_burgers_unit_scale( self ) :
        return self.burgers_vector_scale

    def return_defectStr( self ) :
        return self.defectStr

    def _return_properties( self ) :

        properties_list = [ [ 'position',             self.position             ],
                            [ 'rcut',                 self.rcut                 ],
                            [ 'IV',                   self.IV                   ],
                            [ 'defecttype',           self.defecttype           ],
                            [ 'atoms_numb',           self.atoms_numb           ],
                            [ 'position_center',      self.position_center      ],
                            [ 'burgers_vector_scale', self.burgers_vector_scale ],
                            [ 'defectStr',            self.defectStr            ]  ]

        return properties_list

    def print_properties( self ) :

        print( '# ' )
        print( '# ------------- properties for Cluster --------------' )
        print( '\n'.join( ( '# %s : %s' % ( item[ 0 ], item[ 1 ] ) for item in self._return_properties( ) ) ) )

        print( '# ' )

class ClusterSystem( object ) :

    def __init__( self, 
                  tmp_alatt,
                  tmp_ICluster_upper,
                  tmp_VLoop_lower,
                  tmp_IV,
                  tmp_rcut,
                  tmp_position,
                  tmp_box,
                  tmp_burgers_vector_unit_array,
                  tmp_burgers_vector_scale_length_array,
                  tmp_d_probability                      ) :
        '''
        Class ClusterSystem properties:
         1. self.cluster_list
         2. self.cluster_numb
         3. self.burgers_vector_unit_array
         4. self.burgers_vector_scale_length_array
         5. d_probability
        '''

        self.burgers_vector_unit_array = tmp_burgers_vector_unit_array
        self.burgers_vector_scale_length_array = tmp_burgers_vector_scale_length_array
        self.d_probability = tmp_d_probability

        if tmp_position.size == 0 :
            self.cluster_list = [ ]
            self.cluster_numb = 0
        else :
            self.cluster_list = self.identify_clusters( tmp_alatt, tmp_ICluster_upper, tmp_VLoop_lower, tmp_IV, tmp_rcut, tmp_position, tmp_box )
            self.cluster_numb = len( self.cluster_list )

    def identify_clusters( self,
                           tmp_alatt,
                           tmp_ICluster_upper,
                           tmp_VLoop_lower,
                           tmp_IV,
                           tmp_rcut,
                           tmp_position,
                           tmp_box         ) :
        'Note that the distance is the nearest distance between images. distance unit: nm.'

        rcut_square = tmp_rcut * tmp_rcut

        nn1 = len( tmp_position )

        ltop = [ i for i in range( nn1 ) ]

        cluster_head = [ ]

        for i in range( nn1 - 1 ) :
       
            if i != ltop[ i ] : continue
       
            cluster_head.append( i )
       
            j = i

            while True :
       
                for k in range( i + 1, nn1 ) :
           
                    if k != ltop[ k ] : continue

                    delta = tmp_position[ j ] - tmp_position[ k ]

                    tmp_nearest = [ ]
                    for m in range( 3 ) :
                        tmp_ = min( [ abs( delta[ m ] + jj * tmp_box[ m ] ) for jj in range( -1, 2 ) ] )
                        tmp_nearest.append( tmp_ )
                    delta = np.array( tmp_nearest )

                    if any( delta > tmp_rcut ) : continue
           
                    R_square = np.inner( delta, delta )
                    if R_square > rcut_square : continue
           
                    tmp = ltop[ j ]
                    ltop[ j ] = ltop[ k ]
                    ltop[ k ] = tmp
           
                j = ltop[ j ]
                if j == i : break

        if ltop[ nn1 - 1 ] == nn1 - 1 : cluster_head.append( nn1 - 1 )

        cluster_list = [ ]

        for i in range( len( cluster_head ) ) :

            'very important.'
            position = tmp_position[ cluster_head[ i ] ].reshape( 1, 3 )
            j = cluster_head[ i ]

            while True :

                j = ltop[ j ]
                if j != cluster_head[ i ] :
                    position = np.vstack( ( position, tmp_position[ j ] ) )
                else :
                    break

            cluster_list.append( Cluster( tmp_alatt, tmp_ICluster_upper, tmp_VLoop_lower, tmp_IV, tmp_rcut, position, self.burgers_vector_unit_array, self.burgers_vector_scale_length_array, self.d_probability ) )

        return cluster_list

    def add_clusters( self, 
                      tmp_alatt,
                      tmp_ICluster_upper,
                      tmp_VLoop_lower,
                      tmp_IV,
                      tmp_rcut,
                      tmp_position,
                      tmp_box          ) :

       if tmp_position.size == 0 :
           pass
       else :
           self.cluster_list.extend( self.identify_clusters( tmp_alatt, tmp_ICluster_upper, tmp_VLoop_lower, tmp_IV, tmp_rcut, tmp_position, tmp_box ) )
           self.cluster_numb = len( self.cluster_list )

    def return_cluster_numb( self ) :
        return self.cluster_numb

    def return_defectStr_list( self ) :

        defectStr_list = [ ]

        for item in self.cluster_list :
            defectStr_list.append( item.return_defectStr( ) )

        return defectStr_list

    def _return_properties( self ) :

        properties_list = [ [ 'cluster_numb',   self.cluster_numb   ] ]

        return properties_list

    def print_properties( self ) :

        print( '# ' )
        print( '# ------------- properties for ClusterSystem --------------' )
        print( '\n'.join( ( '# %s : %s' % ( item[ 0 ], item[ 1 ] ) for item in self._return_properties( ) ) ) )

        print( '# ' )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class Cluster and ClusterSystem detecting ---' )
    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    seed = None
    if j_have( jdata, 'seed' ) :
        seed = jdata[ 'seed' ]
        seed = seed % ( 2 **32 )
    np.random.seed( seed )

    cascade = Cascade.Cascade( jdata )

    trap = False
    trapsys = None

    if j_have( jdata, 'trap' ) : trap = jdata[ 'trap' ]
    if trap : trapsys = TrapSystem( jdata )

    for i in range( cascade.return_cascade_numb( ) ) :

        start_time = time.time( )

        'build a defectSys of Class DefectSystem for okmc simulation.'
        defectsys = DefectSystem.DefectSystem( jdata, cascade.return_defectStr_list( jdata, i ) )

        outputFile = OutputFile( defectsys, jdata, i )

        outputFile.reset_steptime( defectsys, jdata )

        outputFile.output_results( defectsys, jdata )

        'reset trap for all defectObject in model of Class DefectSystem.'
        if trap : defectsys.reset_trap_defectSystem( jdata, trapsys )

        end_time = time.time( )

        print( '# cascade %d running time: %.3f s' % ( i, end_time - start_time ) )
        print( '#    ' )

    print( '# All cascades annealing done !' )
