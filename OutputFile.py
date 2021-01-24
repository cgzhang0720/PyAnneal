#!/usr/bin/env python3

import math
import os
import shutil
import numpy as np

from Auxiliary import j_must_have, j_have, to_str

class OutputFile( object ) :

    def __init__( self,
                  defectsys,
                  jdata,
                  index = None ) :
        'property'

        self.out_dir = os.path.join( j_must_have( jdata, 'out_dir' ), 'output' )

        if index == None or index == 0 :
            if os.path.exists( self.out_dir ) : shutil.rmtree( self.out_dir )
            os.mkdir( self.out_dir )

        if index != None : 
            self.out_dir = os.path.join( self.out_dir, str( index ) )
            if os.path.exists( self.out_dir ) : shutil.rmtree( self.out_dir )
            os.mkdir( self.out_dir )

        self.step = defectsys.c_step
        self.time = defectsys.c_time

        self.out_style_freq = j_must_have( jdata, 'out_style_freq' )

        self.out_cfg = False

        if j_have( jdata, 'out_cfg_dir_name' ) : 
            '0. self.out_cfg, 1. self.cfg_path, 2. self.cfg_name, 3. self.target_cfg'

            self.out_cfg = True
            out_cfg_dir_name = jdata[ 'out_cfg_dir_name' ]
            self.cfg_path = os.path.join( self.out_dir, out_cfg_dir_name[ 0 ] )
            os.mkdir( self.cfg_path )

            self.cfg_name = out_cfg_dir_name[ 1 ]
            self.target_cfg = os.path.join( self.cfg_path, self.cfg_name  + '.' + str( self.step ) +'.cfg' )

        self.out_defectStr = False

        if j_have( jdata, 'out_defectStr_dir_name' ) : 
            '0. self.out_defectStr, 1. self.defectStr_path, 2. self.defectStr_name, 3. self.target_defectStr'

            self.out_defectStr = True

            out_defectStr_dir_name = jdata[ 'out_defectStr_dir_name' ]

            self.defectStr_path = os.path.join( self.out_dir, out_defectStr_dir_name[ 0 ] )
            os.mkdir( self.defectStr_path )

            self.defectStr_name = out_defectStr_dir_name[ 1 ]
            self.target_defectStr = os.path.join( self.defectStr_path, self.defectStr_name  + '.' + str( self.step ) )

        self.out_volume_density = False

        if j_have( jdata, 'out_volume_density_dir_name' ) : 
            '0. self.out_volume_density, 1. self.volume_density_path, 2. self.volume_density_name, 3. self.target_volume_density'

            self.out_volume_density = True

            out_volume_density_dir_name = jdata[ 'out_volume_density_dir_name' ]

            self.volume_density_path = os.path.join( self.out_dir, out_volume_density_dir_name[ 0 ] )
            os.mkdir( self.volume_density_path )

            self.volume_density_name = out_volume_density_dir_name[ 1 ]
            self.target_volume_density = os.path.join( self.volume_density_path, self.volume_density_name  + '.' + str( self.step ) )

        self.out_area_density = False

        if j_have( jdata, 'out_area_density_dir_name' ) :
            '0. self.out_area_density, 1. self.area_density_path, 2. self.area_density_name, 3. self.target_area_density'

            self.out_area_density = True

            out_area_density_dir_name = jdata[ 'out_area_density_dir_name' ]

            self.area_density_path = os.path.join( self.out_dir, out_area_density_dir_name[ 0 ] )
            os.mkdir( self.area_density_path )

            self.area_density_name = out_area_density_dir_name[ 1 ]
            self.target_area_density = os.path.join( self.area_density_path, self.area_density_name  + '.' + str( self.step ) )

        self.out_points = False

        if j_have( jdata, 'out_point_defects_dir_name' ) :
            '0. self.out_points, 1. self.points_path, 2. self.points_name, 3. self.target_points'

            self.out_points = True

            out_points_dir_name = jdata[ 'out_point_defects_dir_name' ]

            self.points_path = os.path.join( self.out_dir, out_points_dir_name[ 0 ] )
            os.mkdir( self.points_path )

            self.points_name = out_points_dir_name[ 1 ]
            self.target_points = os.path.join( self.points_path, self.points_name  + '.' + str( self.step ) )

        self.out_size_distri = False

        if j_have( jdata, 'out_size_distri_dir_name' ) : 

            '0. self.out_size_distri, 1. self.size_distri_path, 2. self.size_distri_name, 3. self.target_size_distri.'

            self.out_size_distri = True

            out_size_distri_dir_name = jdata[ 'out_size_distri_dir_name' ]
            self.size_distri_path = os.path.join( self.out_dir, out_size_distri_dir_name )
            os.mkdir( self.size_distri_path )

            self.out_size_distri_defecttype = jdata[ 'out_size_distri_defecttype' ]

            for item in self.out_size_distri_defecttype :
                os.mkdir( os.path.join( self.size_distri_path, item ) )

        '1. self.target_steptime.'

        out_steptime_dir_name = j_must_have( jdata, 'out_steptime_dir_name' )

        steptime_path = os.path.join( self.out_dir, out_steptime_dir_name[ 0 ] )
        os.mkdir( steptime_path )

        self.target_steptime = os.path.join( steptime_path, out_steptime_dir_name[ 1 ] )
        self.file_out_steptime = open( self.target_steptime, 'w' )
        self.file_out_steptime.writelines( [ '# step, c_time(s), total_defects_I, total_defects_V, numb_defectObject' + '\n' ] )
        self.file_out_steptime.close( )

        self.set_min_max_nsize_dict( defectsys, jdata )
        self.set_sum_numb_nsize_dict( defectsys, jdata )

    def is_output_now( self,
                       defectsys,
                       jdata      ) :
        '''
        time to output, return True; else return False.
        "out_style_freq":     [ "step_uniform", 100     ],
        "out_style_freq":     [ "time_uniform", 1.0E-2  ],
        "out_style_freq":     [ "step_exp",     10      ],
        "out_style_freq":     [ "time_exp",     10.0    ],
        '''

        judge = False

        c_step = defectsys.c_step
        if c_step == 0 : judge = True

        c_time = defectsys.c_time

        if self.out_style_freq[ 0 ] == 'step_uniform' :
            if c_step % self.out_style_freq[ 1 ] == 0 : judge = True

        elif self.out_style_freq[ 0 ] == 'time_uniform' :
            if c_time >= self.time + self.out_style_freq[ 1 ] : judge = True
        
        elif self.out_style_freq[ 0 ] == 'step_exp' :
            if c_step >= self.step *  self.out_style_freq[ 1 ] : judge = True

        elif self.out_style_freq[ 0 ] == 'time_exp' :
            if c_time >= self.time * self.out_style_freq[ 1 ] : judge = True

        if judge : self.reset_steptime( defectsys, jdata )

        return judge

    def reset_steptime( self,
                        defectsys,
                        jdata      ) :
        '''
        reset 1. self.step, 2. self.time, 3. self.target_cfg, 4. self.defectStr, 
              5. self.target_volume_density, 5-1. self.target_area_density, 
              6. self.target_points, 7. self.target_size_distri
        open files
        '''

        self.step = defectsys.c_step
        self.time = defectsys.c_time

        if self.out_cfg :         
            self.target_cfg = os.path.join( self.cfg_path, self.cfg_name  + '.' + str( self.step ) + '.cfg' )
            self.file_out_cfg = open( self.target_cfg, 'w' )

        if self.out_defectStr :   
            self.target_defectStr = os.path.join( self.defectStr_path, self.defectStr_name  + '.' + str( self.step ) )
            self.file_out_defectStr = open( self.target_defectStr, 'w' )

        if self.out_volume_density :     
            self.target_volume_density = os.path.join( self.volume_density_path, self.volume_density_name  + '.' + str( self.step ) )
            self.file_out_volume_density = open( self.target_volume_density, 'w' )

        if self.out_area_density :
            self.target_area_density = os.path.join( self.area_density_path, self.area_density_name  + '.' + str( self.step ) )
            self.file_out_area_density = open( self.target_area_density, 'w' )

        if self.out_points :
            self.target_points = os.path.join( self.points_path, self.points_name  + '.' + str( self.step ) )
            self.file_out_points = open( self.target_points, 'w' )

        if self.out_size_distri :

            self.file_out_size_distri = { }

            for item in self.out_size_distri_defecttype :
                target_size_distri = os.path.join( self.size_distri_path, item, 'size_distri.' + str( self.step ) )
                self.file_out_size_distri[ item ] = open( target_size_distri, 'w' )

        self.file_out_steptime = open( self.target_steptime, 'a' )

        self.set_min_max_nsize_dict( defectsys, jdata )
        self.set_sum_numb_nsize_dict( defectsys, jdata )

    def _return_points( self,
                        defectsys,
                        jdata      ) :
        '''
        format:
        0.1 0.1 0.1 0
        0.2 0.2 0.2 1
        '''

        total_point_defects_str_list = [ ]

        for item in defectsys.defectObject_list :
            total_point_defects_str_list.extend( item.convert_point_defects( jdata ) )

        return total_point_defects_str_list

    def output_cfg( self,
                    defectsys,
                    jdata      ) :
        '''
        output configuration file with format *.cfg:

        Number of particles = 10000
        A = 1 Angstrom (basic length-scale)
        H0(1,1) = 760.164 A
        H0(1,2) = 0 A 
        H0(1,3) = 0 A 
        H0(2,1) = 0 A 
        H0(2,2) = 760.164 A
        H0(2,3) = 0 A 
        H0(3,1) = 0 A 
        H0(3,2) = 0 A 
        H0(3,3) = 760.164 A
        .NO_VELOCITY.
        entry_count = 4
        auxiliary[0] = IV
        183.839996 
        W 
        0 0 0 1 
        183.839996 
        W 
        0.00208333 0.00208333 0.00208 0
        '''

        box = j_must_have( jdata, 'box' )
        mass = j_must_have( jdata, 'mass' )
        symbol = j_must_have( jdata, 'symbol' )

        total_point_defects_str_list = self._return_points( defectsys, jdata )
        numb_point_defects = len( total_point_defects_str_list )

        head_str_list = [ ]

        head_str_list.append( 'Number of particles = ' + str( numb_point_defects ) + '\n' )
        head_str_list.append( 'A = 1 Angstrom (basic length-scale)'                + '\n' )
        head_str_list.append( 'H0(1,1) = ' + str( box[ 0 ] * 10.0 ) + ' ' + 'A'    + '\n' )
        head_str_list.append( 'H0(1,2) = 0 A'                                      + '\n' )
        head_str_list.append( 'H0(1,3) = 0 A'                                      + '\n' )
        head_str_list.append( 'H0(2,1) = 0 A'                                      + '\n' )
        head_str_list.append( 'H0(2,2) = ' + str( box[ 1 ] * 10.0 ) + ' ' + 'A'    + '\n' )
        head_str_list.append( 'H0(2,3) = 0 A'                                      + '\n' )
        head_str_list.append( 'H0(3,1) = 0 A'                                      + '\n' )
        head_str_list.append( 'H0(3,2) = 0 A'                                      + '\n' )
        head_str_list.append( 'H0(3,3) = ' + str( box[ 2 ] * 10.0 ) + ' ' + 'A'    + '\n' )
        head_str_list.append( '.NO_VELOCITY.'                                      + '\n' )
        head_str_list.append( 'entry_count =  4'                                   + '\n' )
        head_str_list.append( 'auxiliary[0] = IV'                                  + '\n' )

        self.file_out_cfg.writelines( head_str_list )

        new_total_point_defects_str_list = [ ]
        for item in total_point_defects_str_list :
            new_total_point_defects_str_list.append( str( mass ) + '\n' )
            new_total_point_defects_str_list.append( symbol + '\n' )
            new_total_point_defects_str_list.append( item )

        self.file_out_cfg.writelines( new_total_point_defects_str_list )

    def output_defectStr( self, defectsys ) :
        '''
        output defectStr for next cascade annealing input
        format : ILoop111 13 10.0 10.0 10.0 0.5 0.5 0.5
        '''

        defectStr_list = [ ]

        for item in defectsys.defectObject_list :

            if hasattr( item, 'burgers_vector_scale' ) :
                defectStr = item.defecttype + ' '                    \
                          + str( item.nsize ) + ' '                   \
                          + to_str( item.position ) + ' '              \
                          + to_str( item.burgers_vector_scale ) + '\n'
            else :
                defectStr = item.defecttype + ' '          \
                          + str( item.nsize ) + ' '         \
                          + to_str( item.position ) + '\n'

            defectStr_list.append( defectStr )

        self.file_out_defectStr.writelines( defectStr_list )

    def set_min_max_nsize_dict( self,
                                defectsys,
                                jdata       ) :
        'set min nsize and max nsize for every defecttype.'

        ICluster_lower = 2
        ICluster_upper = j_must_have( jdata, 'ICluster_upper' )       

        VCluster_lower = 2
        VCluster_upper = VCluster_lower #

        ILoop100_lower = ICluster_upper + 1
        ILoop100_upper = ILoop100_lower #

        ILoop111_lower = ICluster_upper + 1
        ILoop111_upper = ILoop111_lower #

        VLoop100_lower = j_must_have( jdata, 'VLoop_lower' )
        VLoop100_upper = VLoop100_lower #
 
        VLoop111_lower = VLoop100_lower
        VLoop111_upper = VLoop111_lower #

        self.min_nsize_dict = { 'I'        : 1,              'V'        : 1,
                                'ICluster' : ICluster_lower, 'VCluster' : VCluster_lower,
                                'ILoop100' : ILoop100_lower, 'VLoop100' : VLoop100_lower,
                                'ILoop111' : ILoop111_lower, 'VLoop111' : VLoop111_lower  }

        self.max_nsize_dict = { 'I'        : 1,              'V'        : 1,
                                'ICluster' : ICluster_upper, 'VCluster' : VCluster_upper,
                                'ILoop100' : ILoop100_upper, 'VLoop100' : VLoop100_upper,
                                'ILoop111' : ILoop111_upper, 'VLoop111' : VLoop111_upper  }

        for item in defectsys.defectObject_list :
            tmp_defecttype = item.defecttype
            tmp_nsize = item.nsize
            if self.max_nsize_dict[ tmp_defecttype ] < tmp_nsize : self.max_nsize_dict[ tmp_defecttype ] = tmp_nsize

        self.max_nsize = max( item for item in self.max_nsize_dict.values( ) )

    def set_sum_numb_nsize_dict( self,
                                 defectsys,
                                 jdata      ) :

        defecttypes = j_must_have( jdata, 'defecttype' )
        self.sum_numb_nsize_dict = { }

        for item in defecttypes :
            self.sum_numb_nsize_dict[ item ] = np.zeros( self.max_nsize + 1 )

        for item in defectsys.defectObject_list :
            tmp_defecttype = item.defecttype
            tmp_nsize = item.nsize
            self.sum_numb_nsize_dict[ tmp_defecttype ][ tmp_nsize ] += 1.0

    def output_volume_density( self, jdata ) :
        ''' format:
        I 1 1.0
        V 1 1.0
        ICluster nsize 1.0
        VCluster nsize 1.0
        ILoop100 nsize 1.0
        ILoop111 nsize 1.0
        VLoop100 nsize 1.0
        VLoop111 nsize 1.0
        '''

        box = j_must_have( jdata, 'box' )
        v = box[ 0 ] * box[ 1 ] * box[ 2 ]

        density_str_list = [ '#--(1)Dtype--(2) nsize--(3) density( nm^-3 )---\n' ]
        for item in self.sum_numb_nsize_dict.keys( ) :

            min_size = self.min_nsize_dict[ item ]
            max_size = self.max_nsize_dict[ item ]
            for i in range( self.sum_numb_nsize_dict[ item ].size ) :
                if i >= min_size and i <= max_size :
                    density_str = item + ' ' + str( i ) + ' ' + str( self.sum_numb_nsize_dict[ item ][ i ] / v ) + '\n'
                    density_str_list.append( density_str )

        self.file_out_volume_density.writelines( density_str_list )

    def output_area_density( self, jdata ) :
        ''' format:
        I        1     1.0
        V        1     1.0
        ICluster nsize 1.0
        VCluster nsize 1.0
        ILoop100 nsize 1.0
        ILoop111 nsize 1.0
        VLoop100 nsize 1.0
        VLoop111 nsize 1.0
        '''

        box = j_must_have( jdata, 'box' )
        area = box[ 0 ] * box[ 1 ]

        density_str_list = [ '#--(1) Dtype--(2) nsize--(3) density( nm^-2 )---\n' ]
        for item in self.sum_numb_nsize_dict.keys( ) :

            min_size = self.min_nsize_dict[ item ]
            max_size = self.max_nsize_dict[ item ]
            for i in range( self.sum_numb_nsize_dict[ item ].size ) :
                if i >= min_size and i <= max_size :
                    density_str = item + ' ' + str( i ) + ' ' + str( self.sum_numb_nsize_dict[ item ][ i ] / area ) + '\n'
                    density_str_list.append( density_str )

        self.file_out_area_density.writelines( density_str_list )

    def output_points( self,
                       defectsys,
                       jdata      ) :
        ''' format:
        I 1 10.0 10.0 10.0
        V 1 10.0 10.0 10.0
        '''

        box = np.array( j_must_have( jdata, 'box' ) )
        total_point_defects_str_list = self._return_points( defectsys, jdata )

        new_total_point_defects_str_list = [ to_str( box ) + ' # b011, b022, b033\n' ]

        for item in total_point_defects_str_list :

            defectStr = item.lstrip( ).rstrip( )
            defectStr_list = defectStr.split( )
            position = np.array( [ np.float64( item2 ) for item2 in defectStr_list[ 0 : 3 ] ] ) * box

            if defectStr_list[ 3 ] == '0' :
                new_defectStr = 'V ' + ' 1 ' + to_str( position ) + '\n'
            else :
                new_defectStr = 'I ' + ' 1 ' + to_str( position ) + '\n'

            new_total_point_defects_str_list.append( new_defectStr )

        self.file_out_points.writelines( new_total_point_defects_str_list )

    def output_size_distri( self ) :

        for item in self.out_size_distri_defecttype :

            if item in self.sum_numb_nsize_dict.keys( ) :

                min_size = self.min_nsize_dict[ item ]
                max_size = self.max_nsize_dict[ item ]
                self._out_cluster_size_distri( item, min_size, max_size, self.sum_numb_nsize_dict[ item ] )

            else :

                min_size, max_size, sum_nsize = self._return_size_distri( item )
                self._out_cluster_size_distri( item, min_size, max_size, sum_nsize )

    def _return_size_distri( self, cluster_name ) :

        if cluster_name == 'VAll' :

            min_size = min( [ self.min_nsize_dict[ 'VCluster' ], self.min_nsize_dict[ 'VLoop100' ], self.min_nsize_dict[ 'VLoop111' ] ] )
            max_size = min( [ self.max_nsize_dict[ 'VCluster' ], self.max_nsize_dict[ 'VLoop100' ], self.max_nsize_dict[ 'VLoop111' ] ] )

            sum_nsize = self.sum_numb_nsize_dict[ 'VCluster' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'VLoop100' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'VLoop111' ][ min_size : max_size + 1 ] 

        elif cluster_name == 'IAll' :

            min_size = min( [ self.min_nsize_dict[ 'ICluster' ], self.min_nsize_dict[ 'ILoop100' ], self.min_nsize_dict[ 'ILoop111' ] ] )
            max_size = min( [ self.max_nsize_dict[ 'ICluster' ], self.max_nsize_dict[ 'ILoop100' ], self.max_nsize_dict[ 'ILoop111' ] ] )

            sum_nsize = self.sum_numb_nsize_dict[ 'ICluster' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'ILoop100' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'ILoop111' ][ min_size : max_size + 1 ] 

        elif cluster_name == 'ILoopAll' :

            min_size = min( [ self.min_nsize_dict[ 'ILoop100' ], self.min_nsize_dict[ 'ILoop111' ] ] )
            max_size = min( [ self.max_nsize_dict[ 'ILoop100' ], self.max_nsize_dict[ 'ILoop111' ] ] )

            sum_nsize = self.sum_numb_nsize_dict[ 'ILoop100' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'ILoop111' ][ min_size : max_size + 1 ] 

        elif cluster_name == 'VLoopAll' :

            min_size = min( [ self.min_nsize_dict[ 'VLoop100' ], self.min_nsize_dict[ 'VLoop111' ] ] )
            max_size = min( [ self.max_nsize_dict[ 'VLoop100' ], self.max_nsize_dict[ 'VLoop111' ] ] )

            sum_nsize = self.sum_numb_nsize_dict[ 'VLoop100' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'VLoop111' ][ min_size : max_size + 1 ] 

        elif cluster_name == 'LoopAll' :

            min_size = min( [ self.min_nsize_dict[ 'ILoop100' ], self.min_nsize_dict[ 'ILoop111' ], self.min_nsize_dict[ 'VLoop100' ], self.min_nsize_dict[ 'VLoop111' ] ] )
            max_size = min( [ self.max_nsize_dict[ 'ILoop100' ], self.max_nsize_dict[ 'ILoop111' ], self.max_nsize_dict[ 'VLoop100' ], self.max_nsize_dict[ 'VLoop111' ] ] )

            sum_nsize = self.sum_numb_nsize_dict[ 'ILoop100' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'ILoop111' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'VLoop100' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'VLoop111' ][ min_size : max_size + 1 ]

        elif cluster_name == 'All' :

            min_size = 1
            max_size = self.max_nsize

            sum_nsize = self.sum_numb_nsize_dict[ 'ILoop100' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'ILoop111' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'VLoop100' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'VLoop111' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'ICluster' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'VCluster' ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'I'        ][ min_size : max_size + 1 ] \
                      + self.sum_numb_nsize_dict[ 'V'        ][ min_size : max_size + 1 ]

        return ( min_size, max_size, sum_nsize )

    def _out_cluster_size_distri( self,
                                  cluster_name,
                                  min_size,     
                                  max_size,     
                                  sum_nsize     ) :
         density_str_list = [ ]
         for i in range( sum_nsize.size ) :
             if i >= min_size and i <= max_size :
                 density_str = str( i ) + ' ' + str( sum_nsize[ i ] ) + '\n'
                 density_str_list.append( density_str )

         self.file_out_size_distri[ cluster_name ].writelines( density_str_list )

    def output_steptime( self, defectsys ) :

        steptime_str = str( self.step ) + ' ' + str( self.time ) + ' '                                     \
                      + str( defectsys.total_defects[ 'I' ] ) + ' ' + str( defectsys.total_defects[ 'V' ] ) \
                      + ' ' + str( defectsys.numb_defectObject ) + '\n'

        self.file_out_steptime.writelines( [ steptime_str ] )

    def analyse_results( self,
                         defectsys,
                         jdata      ) :
        'self.sum_numb_nsize_dict = { }, visible defects or others corresponding to Ref.'
        return None

    def close_files( self ) :
        'close all result files.'

        if self.out_cfg : self.file_out_cfg.close( )

        if self.out_defectStr : self.file_out_defectStr.close( )

        if self.out_volume_density : self.file_out_volume_density.close( )
        if self.out_area_density : self.file_out_area_density.close( )

        if self.out_points : self.file_out_points.close( )

        if self.out_size_distri :
            for item in self.file_out_size_distri.values( ) : item.close( )

        self.file_out_steptime.close( )

    def output_results( self,
                        defectsys,
                        jdata      ) :

        if self.out_cfg : self.output_cfg( defectsys, jdata )

        if self.out_defectStr : self.output_defectStr( defectsys )

        if self.out_volume_density : self.output_volume_density( jdata )

        if self.out_area_density : self.output_area_density( jdata )

        if self.out_points : self.output_points( defectsys, jdata )

        if self.out_size_distri : self.output_size_distri( )

        self.output_steptime( defectsys )

        self.close_files( )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class OutputFile detecting ---' )
    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
