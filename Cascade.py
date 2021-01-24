#!/usr/bin/env python3

import glob
import numpy as np
import os.path
import shutil
import math
import argparse
import json

from scipy import stats
from collections import Iterable

from Cluster import ClusterSystem
from Auxiliary import j_have, j_must_have, R_convert_to_position, to_str, return_serial_numb_in_str

class Cascade( object ) :

    def __init__( self, jdata ) :
        '''
        Class Cascade properties:

         1. self.defectStr_list
            format:

            (1): ILoop100 13 10.0 10.0 10.0
            (2): ILoop111 13 10.0 10.0 10.0 0.5 0.5 0.5

        mode :  
         1. cascade from formula

         2. cascade from file
            a. self.cluster == True, defects of I and V convert into clusters
            b. self.cluster == False, defects of I and V does not convert into clusters

         3. cascade I and V randomly and uniformly distributed in simulation box, for electron irradiation
        '''

        self.cascade_from_formula = False

        if j_have( jdata, 'cascade_from_formula' ) : self.cascade_from_formula = jdata[ 'cascade_from_formula' ]

        if self.cascade_from_formula : 

            self.cascade_numb = j_must_have( jdata, 'cascade_numb' )
            self.I_total = j_must_have( jdata, 'I_total' )
            self.SI = j_must_have( jdata, 'SI' )
            self.mu = j_must_have( jdata, 'mu' )
            self.sigma = j_must_have( jdata, 'sigma' )
            self.SV = j_must_have( jdata, 'SV' )

            self.ICluster_upper = j_must_have( jdata, 'ICluster_upper' )
            self.ratio_ILoop111_ILoop100 = j_must_have( jdata, 'ratio_ILoop111_ILoop100' )

            self.VLoop_lower = j_must_have( jdata, 'VLoop_lower' )
            self.ratio_VCluster_VLoop100_VLoop111 = j_must_have( jdata, 'ratio_VCluster_VLoop100_VLoop111' )
       
        self.cascade_from_file = False

        if j_have( jdata, 'cascade_from_file' ) : self.cascade_from_file = jdata[ 'cascade_from_file' ]

        if self.cascade_from_file : 

            sys_path = j_must_have( jdata, 'systems' )
            cascade_prefix = j_must_have( jdata, 'cascade_prefix' )

            if j_have( jdata, 'sys_defects_cfg' ) :
                sys_defects_cfg = jdata[ 'sys_defects_cfg' ]
                defects_cfg_prefix = j_must_have( jdata, 'defects_cfg_prefix' )
                self.defects_cfg_file = glob.glob( os.path.join( sys_defects_cfg, defects_cfg_prefix + '.*' ) )
                self.defects_cfg_file.sort( key = return_serial_numb_in_str )

                if os.path.exists( sys_path ) : shutil.rmtree( sys_path )
                os.mkdir( sys_path )

                self.cascade_cfg_2_points_defectStr( sys_path, cascade_prefix ) 

            self.cascade_file = glob.glob( os.path.join( sys_path, cascade_prefix + '.*' ) )
            self.cascade_file.sort( key = return_serial_numb_in_str )
            self.cascade_numb = len( self.cascade_file )

            self.cluster = j_must_have( jdata, 'cluster' )

            if self.cluster : 
                self.IC_VC_rcut = j_must_have( jdata, 'IC_VC_rcut' )
                self.VLoop_lower = j_must_have( jdata, 'VLoop_lower' )
                self.ICluster_upper = j_must_have( jdata, 'ICluster_upper' )
 
        self.point_defects_uniform_distribution = False

        if j_have( jdata, 'point_defects_uniform_distribution' ) : self.point_defects_uniform_distribution = jdata[ 'point_defects_uniform_distribution' ]

        if self.point_defects_uniform_distribution :

            self.cascade_numb = j_must_have( jdata, 'cascade_numb' )
            self.I_total = j_must_have( jdata, 'I_total' )

    def create_defects_formula_distribution( self ) :
        '''
        create defects from distribution formula.
        the appendix A. of Ref: D.R. Mason et al., Elastic trapping of dislocation loops in cascades in ion-irradiated tungsten foils, 
        J. Phys.: Condense. Matter 26, 375701( 2014 ).
        '''
 
        nmax = int( self.I_total / 10 )
     
        I_numb = 0
        self.int_defect_size_R = [ ]
     
        while True :
 
            tmp_defect_size_R = self.create_interstitial_clusters( self.SI, self.mu, self.sigma, 1 )[ 0 ]
            I_numb = I_numb + tmp_defect_size_R[ 0 ]
     
            if I_numb == self.I_total :
                self.int_defect_size_R.append( tmp_defect_size_R )
                break
            elif I_numb < self.I_total :
                self.int_defect_size_R.append( tmp_defect_size_R )
            elif I_numb > self.I_total :
                I_numb = I_numb - tmp_defect_size_R[ 0 ]
 
        V_numb = 0
        self.vac_defect_size_R = [ ]
     
        while True :
 
            tmp_defect_size_R = self.create_vacancy_clusters( self.SV, self.sigma, nmax, 1 )[ 0 ]
            V_numb = V_numb + tmp_defect_size_R[ 0 ]
     
            if V_numb == self.I_total :
                self.vac_defect_size_R.append( tmp_defect_size_R )
                break
            elif V_numb < self.I_total :
                self.vac_defect_size_R.append( tmp_defect_size_R )
            elif V_numb > self.I_total :
                V_numb = V_numb - tmp_defect_size_R[ 0 ]

    def create_interstitial_clusters( self, 
                                      SI, 
                                      mu, 
                                      sigma, 
                                      nsize, 
                                      seed = None ) :
        '''
        interstitial clusters(loops) are initially distributed with position ( truncnormal ) and defect_size ( powerlaw ) uncorrelated.
        i.e., ( 1 / N^SI ) * exp( -( R - mu )^2 /(  2 * sigma^2 ) ), R_lower = 0.0, R_upper = np.inf
        '''
     
        lower = 0.0
     
        defect_size = self.create_powerlaw_random_numb_2( SI - 1, nsize )
        defect_R = self.create_truncnormal_random_numb( mu, sigma, nsize, lower )
     
        defect_size_R = [ ]
     
        for i in range( nsize ) :
            defect_size_R.append( ( int( defect_size[ i ] ), defect_R[ i ] ) )
     
        return defect_size_R
     
    def create_vacancy_clusters( self,
                                 SV, 
                                 sigma, 
                                 nmax, 
                                 nsize ) :
        '''
        nmax and Rmax trunc parameter
     
        vacancy clusters( loops ) are initially distributed with position ( truncnormal ) and defect_size ( powerlaw ) uncorrelated.
        i.e., exp( -1/2 * ( R / sigma )^SV * N ) * exp( -R^2 / ( 2 * sigma ^2 ) ), defect_size in range[ 2, nmax ], defect_R in range [ 0.0, Rmax ]
        '''
     
        defect_size_lower = 1
        defect_size_upper = nmax
        defect_size_upper_lower = defect_size_upper - defect_size_lower
     
        defect_R_lower = 0.0
        defect_R_upper = 5.0 * sigma
        defect_R_upper_lower = defect_R_upper - defect_R_lower
     
        m = 0
        defect_size_R = [ ]
     
        while m < nsize :
     
            kesi = np.random.random( 3 )
     
            kesi_1 = kesi[ 0 ] * defect_size_upper_lower + defect_size_lower
            kesi_2 = kesi[ 1 ] * defect_R_upper_lower + defect_R_lower

            kesi_2_vers_sigma = kesi_2 / sigma
 
            '''
            f_kesi_1_2 = exp( -0.5 * kesi_2_vers_sigma**SV * kesi_1 ) * exp( -0.5 * kesi_2_vers_sigma * kesi_2_vers_sigma ) / Lmax
            g_kesi_1_2 = math.log( f_kesi_1_2 )
            '''
 
            g_kesi_1_2 = -0.5 * ( kesi_2_vers_sigma**SV * kesi_1 + kesi_2_vers_sigma * kesi_2_vers_sigma )
 
            if math.log( kesi[ 2 ] ) < g_kesi_1_2 :
                m = m + 1
                defect_size_R.append( ( int( kesi_1 ), kesi_2 ) )

        return defect_size_R

    def create_uniform_random_numb( self,
                                    a, 
                                    b, 
                                    nsize, 
                                    seed = None ) :
        'create a np.array with size of nsize, the data is uniformly and randomly distributed in [ a, b ).'
     
        if seed != None : np.random.seed( seed )
     
        return np.random.random( nsize ) * ( b - a ) + a
     
    def create_normal_random_numb( self,
                                   mu, 
                                   sigma, 
                                   nsize, 
                                   seed = None ) :
        'normal distribution: f( x ) = 1.0 / ( sqrt( 2 * pi ) * sigma ) * exp( -( x - mu )^2 / ( 2 * sigma^2 ) ).'
  
        if seed != None : np.random.seed( seed )
     
        return stats.norm.rvs( mu, sigma, nsize )
     
    def create_truncnormal_random_numb( self,
                                        mu, 
                                        sigma, 
                                        nsize, 
                                        lower, 
                                        upper = np.inf, 
                                        seed = None     ) :
        '''
        truncnormal distribution: 
        f( x ) = 1.0 / ( sqrt( 2 * pi ) * sigma ) * exp( -( x - mu )^2 / ( 2 * sigma^2 ) ) in [ lower, upper ]
        a, b = (myclip_a - my_mean) / my_std, (myclip_b - my_mean) / my_std
        '''
     
        if seed != None : np.random.seed( seed )
     
        lower = ( lower - mu ) / sigma
        if upper != np.inf : upper = ( upper - mu ) / sigma
     
        return stats.truncnorm.rvs( lower, upper, mu, sigma, nsize )
     
    def create_powerlaw_random_numb_1( self,
                                       a, 
                                       nsize, 
                                       seed = None ) :
        '''
        power law distribution: f( x ) = a * x^( a - 1 ), 0 <= x <= 1, a > 0
        create a np.array with size of nsize, the data is power law and randomly distributed with ( a )
        '''
     
        if seed != None : np.random.seed( seed )
 
        return stats.powerlaw.rvs( a, size = nsize )
 
    def create_powerlaw_random_numb_2( self,
                                       b, 
                                       nsize, 
                                       seed = None ) :
        '''
        power law distribution: f( x ) = b / x^( b + 1 ), x >= 1, b > 0
        create a np.array with size of nsize, the data is power law and randomly distributed with ( b )
        '''
     
        if seed != None : np.random.seed( seed )
     
        return stats.pareto.rvs( b, size = nsize )

    def load_cascade_from_formula( self,
                                   jdata,
                                   seed = None ) :
        '''
        self.defectStr_list, defectStr format: 
         'I        1  10.0 10.0 10.0\n'
         'ICluster 7  10.0 10.0 10.0\n'
         'ILoop111 13 10.0 10.0 10.0\n'
         'ILoop100 13 10.0 10.0 10.0\n'

         'V        1  10.0 10.0 10.0\n'
         'VCluster 7  10.0 10.0 10.0\n'
         'VLoop111 13 10.0 10.0 10.0\n'
         'VLoop100 13 10.0 10.0 10.0\n'
        '''

        if seed != None : np.random.seed( seed )

        tmp_box = np.array( j_must_have( jdata, 'box' ) )

        self.create_defects_formula_distribution( )

        self.defectStr_list = [ ]

        for item in self.int_defect_size_R :

            if item[ 0 ] == 1 :
                defecttype = 'I'
            elif item[ 0 ] >= 2 and item[ 0 ] <= self.ICluster_upper :
                defecttype = 'ICluster'
            elif item[ 0 ] > self.ICluster_upper :
                if np.random.rand( ) < self.ratio_ILoop111_ILoop100[ 0 ] :
                    defecttype = 'ILoop111'
                else :
                    defecttype = 'ILoop100'

            position = 0.5 * tmp_box + R_convert_to_position( item[ 1 ] )
            defectStr = defecttype + ' ' + str( item[ 0 ] ) + ' ' + to_str( position ) + '\n'

            self.defectStr_list.append( defectStr )

        for item in self.vac_defect_size_R :

            if item[ 0 ] == 1 :
                defecttype = 'V'
            elif item[ 0 ] >= 2 and item[ 0 ] < self.VLoop_lower :
                defecttype = 'VCluster'
            elif item[ 0 ] >= self.VLoop_lower :

                kesi = np.random.rand( )

                if kesi < self.ratio_VCluster_VLoop100_VLoop111[ 0 ] :
                    defecttype = 'VCluster'
                elif kesi < self.ratio_VCluster_VLoop100_VLoop111[ 0 ] + self.ratio_VCluster_VLoop100_VLoop111[ 1 ] :
                    defecttype = 'VLoop100'
                else :
                    defecttype = 'VLoop111'

            position = 0.5 * tmp_box + R_convert_to_position( item[ 1 ] )
            defectStr = defecttype + ' ' + str( item[ 0 ] ) + ' ' + to_str( position ) + '\n'

            self.defectStr_list.append( defectStr )

    def load_cascade_from_file( self, 
                                jdata,
                                serial_numb ) :
        '''
        1. self.cluster == True, cascade file format:
         100.0 100.0 100.0 # b011, b022, b033, units: nm
         I 1 10.0 10.0 10.0
         V 1 10.0 10.0 10.0
         ...        

         reading defects and convert clusters

        2. self.cluster == False, cascade file format:
         100.0 100.0 100.0 # b011, b022, b033, units: nm
         (1): ILoop100 13 10.0 10.0 10.0
         (2): ILoop111 13 10.0 10.0 10.0 0.5 0.5 0.5
         ...

         reading defects directly
        '''

        if serial_numb > self.cascade_numb :
            raise RuntimeError( 'serial_numb %d not in range [ 0, %d ]' % ( serial_numb, self.cascade_numb ) )

        file_in = open( self.cascade_file[ serial_numb ], 'r' )

        alatt = j_must_have( jdata, 'alatt' )

        burgers_vector_unit_array = np.array( j_must_have( jdata, 'burgers_vector_unit' ) )

        burgers_vector_scale_length_array = np.array( j_must_have( jdata, 'burgers_vector_scale_length' ) )

        'set d_probability = 1.0 / burgers_vector_scale_length_array.size + 0.037142857 formatly.'
        d_probability = 0.037142857 + 1.0 / burgers_vector_scale_length_array.size
        if j_have( jdata, 'probability' ) : d_probability = jdata[ 'probability' ]

        if self.cluster :

            # test_cgzhang
            print( 'defect cluster algorithm is done!' )

            I_position_list = [ ]
            V_position_list = [ ]

            s = file_in.readline( )
            defectStr = s.lstrip( ).rstrip( )
            defectStr_list = defectStr.split( )
            tmp_box = np.array( [ np.float64( item ) for item in defectStr_list[ 0 : 3 ] ] )

            while True :
                s = file_in.readline( )
                if s == '':
                    break

                defectStr = s.lstrip( ).rstrip( )
                defectStr_list = defectStr.split( )

                position = [ np.float64( item ) for item in defectStr_list[ 2 : 5 ] ]

                if defectStr_list[ 0 ] == 'I' :
                    I_position_list.append( position )
                else :
                    V_position_list.append( position ) 

            I_position_array = np.array( I_position_list )
            V_position_array = np.array( V_position_list )

            clusterSys = ClusterSystem( alatt, self.ICluster_upper, self.VLoop_lower, 'I', self.IC_VC_rcut[ 0 ], I_position_array, tmp_box, burgers_vector_unit_array, burgers_vector_scale_length_array, d_probability )
            clusterSys.add_clusters( alatt, self.ICluster_upper, self.VLoop_lower, 'V', self.IC_VC_rcut[ 1 ], V_position_array, tmp_box )

            self.defectStr_list = clusterSys.return_defectStr_list( ) 

        else :

            self.defectStr_list = [ ]
 
            s = file_in.readline( )
            while True :
                s = file_in.readline( )
                if s == '':
                    break
                self.defectStr_list.append( s )

    def load_point_defects_uniform_distribution( self,
                                                 jdata,
                                                 seed = None ) :
        '''
        self.defectStr_list, format: 
         'I 1 10.0 10.0 10.0\n'
         'V 1 10.0 10.0 10.0\n'

         len( ) == 5
        '''

        if seed != None : np.random.seed( seed )

        tmp_box = np.array( j_must_have( jdata, 'box' ) )

        self.defectStr_list = [ ]

        for i in range( self.I_total ) :
            position = np.random.random( 3 ) * tmp_box
            defectStr = 'I 1 ' + to_str( position ) + '\n'
            self.defectStr_list.append( defectStr )

        for i in range( self.I_total ) :
            position = np.random.random( 3 ) * tmp_box
            defectStr = 'V 1 ' + to_str( position ) + '\n'
            self.defectStr_list.append( defectStr )
     
    def load_cascade( self,
                      jdata,
                      serial_numb = 0,
                      seed = None      ) :

        if self.cascade_from_formula : self.load_cascade_from_formula( jdata, seed )

        if self.cascade_from_file : self.load_cascade_from_file( jdata, serial_numb )

        if self.point_defects_uniform_distribution : self.load_point_defects_uniform_distribution( jdata, seed )


    def cascade_cfg_2_points_defectStr( self,
                                        sys_path,
                                        cascade_prefix ) :
        '''
        cfg format:
        Number of particles =    564               0   4
        A = 1.0 Angstrom (basic length-scale)      1
        H0(1,1) =    481.11039925 A                2   2   
        H0(1,2) =      0.00000000 A
        H0(1,3) =      0.00000000 A
        H0(2,1) =      0.00000000 A
        H0(2,2) =    481.11039925 A                6   2
        H0(2,3) =      0.00000000 A
        H0(3,1) =      0.00000000 A
        H0(3,2) =      0.00000000 A
        H0(3,3) =    664.69199896 A               10   2
        .NO_VELOCITY.
        entry_count =  4
        auxiliary[0] = IV                         13
        183.84                                    14                    
        W                                         15
          0.45574900000000002        4.4982300000000003E-002  0.27725538095238089                1
        183.84
        W
         0.47206700000000001       0.10702200000000001       0.24135238095238087                1

        points_defectStr format:
        I 1 10.0 10.0 10.0
        V 1 10.0 10.0 10.0
        '''

        for i in range( len( self.defects_cfg_file ) ) :

            file_in = open( self.defects_cfg_file[ i ], 'r' )
            file_points_defectStr = os.path.join( sys_path, cascade_prefix + '.' + str( i ) )
            file_out = open( file_points_defectStr, 'w' )
 
            lines = [ ]
            while True :
                s = file_in.readline( )
                if s == '':
                    break
                lines.append( s )
            file_in.close( )
 
            for i in range( 11 ) :
 
                defectStr = lines[ i ].lstrip( ).rstrip( )
                defectStr_list = defectStr.split( )
 
                if i == 0 : 
                    defects_numb = np.int64( defectStr_list[ 4 ] )
                elif i == 2  : 
                    b011 = np.float64( defectStr_list[ 2 ] )
                elif i == 6  : 
                    b022 = np.float64( defectStr_list[ 2 ] )
                elif i == 10 : 
                    b033 = np.float64( defectStr_list[ 2 ] )
                    tmp_box = np.array( [ b011, b022, b033 ] ) / 10.0

            points_defectStr_list = [ to_str( tmp_box ) + ' # b011, b022, b033\n' ]

            for i in range( 14, len( lines ) ) :
                less = ( i - 13 ) % 3
                defectStr = lines[ i ].lstrip( ).rstrip( )
                defectStr_list = defectStr.split( )
                if less == 1 : 
                    mass = np.float64( defectStr_list[ 0 ] )
                elif less == 2 : 
                    symbol = defectStr_list[ 0 ]
                elif less == 0 : 
                    position_scale = np.array( [ np.float64( item ) for item in defectStr_list[ 0 : 3 ] ] )
                    IV = defectStr_list[ 3 ]
            
                    if IV == '0' : 
                        points_defectStr_list.append( 'V 1 ' + to_str( position_scale * tmp_box ) + '\n' )
                    elif IV == '1' :
                        points_defectStr_list.append( 'I 1 ' + to_str( position_scale * tmp_box ) + '\n' )
 
            file_out.writelines( points_defectStr_list )
 
            file_out.close( )

    def return_cascade_numb( self ) :
        return self.cascade_numb

    def return_defectStr_list( self,
                               jdata,
                               serial_numb,
                               seed = None  ) :

        self.load_cascade( jdata, serial_numb, seed )

        return self.defectStr_list

    def _return_properties( self ) :

        if self.cascade_from_formula :

            properties_list = [ [ 'cascade_from_formula',             self.cascade_from_formula             ],
                                [ 'cascade_numb',                     self.cascade_numb                     ],
                                [ 'I_total',                          self.I_total                          ],
                                [ 'SI',                               self.SI                               ],
                                [ 'mu',                               self.mu                               ],
                                [ 'sigma',                            self.sigma                            ],
                                [ 'SV',                               self.SV                               ],
                                [ 'ICluster_upper',                   self.ICluster_upper                   ],
                                [ 'ratio_ILoop111_ILoop100',          self.ratio_ILoop111_ILoop100          ],
                                [ 'VLoop_lower',                      self.VLoop_lower                      ],
                                [ 'ratio_VCluster_VLoop100_VLoop111', self.ratio_VCluster_VLoop100_VLoop111 ]  ]

        if self.cascade_from_file :

            if self.cluster :
                properties_list = [ [ 'cascade_from_file', self.cascade_from_file ],
                                    [ 'cascade_numb',      self.cascade_numb      ],
                                    [ 'cascade_file',      self.cascade_file      ],
                                    [ 'cluster',           self.cluster           ],
                                    [ 'IC_VC_rcut',        self.IC_VC_rcut        ],
                                    [ 'VLoop_lower',       self.VLoop_lower       ],
                                    [ 'ICluster_upper',    self.ICluster_upper    ]  ]
            else :
                properties_list = [ [ 'cascade_from_file', self.cascade_from_file ],
                                    [ 'cascade_numb',      self.cascade_numb      ],
                                    [ 'cascade_file',      self.cascade_file      ],
                                    [ 'cluster',           self.cluster           ]  ]

        if self.point_defects_uniform_distribution :

            properties_list = [ [ 'point_defects_uniform_distribution', self.point_defects_uniform_distribution ],
                                [ 'cascade_numb',                       self.cascade_numb                       ],
                                [ 'I_total',                            self.I_total                            ]  ]

        return properties_list

    def print_properties( self ) :

        print( '# ' )
        print( '# ------------- properties for Cascade --------------' )
        print( '\n'.join( ( '# %s : %s' % ( item[ 0 ], item[ 1 ] ) for item in self._return_properties( ) ) ) )

        print( '# defectStr_list:', self.defectStr_list )
        print( '# ' )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class Cascade detecting ---' )
    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    '''
    cascade = Cascade( jdata )

    for i in range( cascade.return_cascade_numb( ) ) :

        print( '# cascade', i )
        cascade.return_defectStr_list( jdata, i )
        cascade.print_properties( )
    '''
