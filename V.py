#!/usr/bin/env python3

import numpy as np
import math
import argparse
import json

from tungsten.DefectObject import DefectObject

from Auxiliary import j_have

class V( DefectObject ) :
    'Class V properties: 1. defecttype, 2. shape'

    defecttype = 'V'
    shape = 'particle'

    '1. reshape'
    def _set_migration( self ) :
        'Set self.migration_energy ( unit: eV ) and self.pref_migration( unit: Hz ) for migrated DefectObject.'

        self.migration_energy = 1.66
        self.pref_migration = 6.0E12

    '2. reshape'
    def _set_shape_properties( self, jdata, tmp_ConstNumber ) :
        ''' Set 1. self.rho( unit: nm ), 2. self.delta_rho( unit: nm ), 3. self.radius( unit: nm ), 4. self.relax_volume( nm^-3 ) 
        according to Class.shape : 1. particle, and defecttype.
        Ref: D.R. Mason et al., J. Phys.: Condens. Matter 26(2014)375701.   '''

        self.rho = tmp_ConstNumber.alatt * ( 3.0 / ( 8.0 * math.pi ) )**0.3333333

        'default value: 0.5*tmp_ConstNumber.alatt'
        self.delta_rho = 0.5 * tmp_ConstNumber.alatt

        if j_have( jdata, 'V_delta_rho' ) : self.delta_rho = jdata[ 'V_delta_rho' ]

        self.radius = self.rho + self.delta_rho

        self.relax_volume = -0.37 * tmp_ConstNumber.volume

    '3. reshape'
    def return_u12( self,
                    tmp_ConstNumber,
                    defectObject,
                    R_vector,
                    R_square         ) :
        'Return elastic interaction energy between two DefectObjects, unit: eV.'

        if defectObject.shape == 'disk' or defectObject.defecttype == 'I' :
            return self._return_u_sphere_disk( tmp_ConstNumber, defectObject, R_vector, R_square )
        else :
            return 0.0

    '4. reshape'
    def transite_defecttype_increase_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize       ) :
        'Tansite from type tmp_defecttype to type transite_type because increase tmp_nsize >= item[ 0 ].'

        transite_increase_nsize = ( 2, 'VCluster' )

        if tmp_nsize >= transite_increase_nsize[ 0 ] :
            return transite_increase_nsize[ 1 ]
        else :
            return tmp_defecttype

    '5. reshape'
    def _return_properties( self ) :

        properties_list = super( )._return_properties( )

        'Auxiliary properties.'
        properties_list.extend( [ [ 'migration_energy', self.migration_energy ],
                                  [ 'pref_migration',   self.pref_migration   ],
                                  [ 'rho',              self.rho              ],
                                  [ 'delta_rho',        self.delta_rho        ],
                                  [ 'radius',           self.radius           ],
                                  [ 'relax_volume',     self.relax_volume     ],
                                  [ 'defecttype',       self.defecttype       ]  ] )

        return properties_list

    def return_judge_recombine( self,
                                jdata,
                                tmp_defectObject ) :
        'Judge recombine or not for self.shape == sphere.'

        if tmp_defectObject.shape == 'sphere' or tmp_defectObject.defecttype == 'V' :
            return self._return_judge_recombine_sphere_sphere( jdata, tmp_defectObject )
        else :
            return self._return_judge_recombine_sphere_disk( jdata, tmp_defectObject )

    def _set_emission( self, tmp_ConstNumber ) :
        'Set self.emission_energy( unit: eV ) and self.pref_emission( unit: Hz ) for emiting DefectObject.'
        raise RuntimeError( 'emission_energy and pref_emission have not been set for Class ' + self.defecttype )

    def _reset_b_pn_rotation( self, tmp_ConstNumber ) :
        'Reset 1. self.burgers_vector_unit, 2. self.burgers_vector_scale and 3. self.position_next_list after rotation.'
        raise RuntimeError( 'function : _reset_b_pn_rotation has not been set for Class ' + self.defecttype )

    def _set_rotation( self ) :
        'Set self.rotation_energy( unit: eV ) and self.pref_rotation( unit: Hz ) for rotating DefectObject.'
        raise RuntimeError( 'rotation_energy and pref_rotation have not been set for Class ' + self.defecttype )

    def _set_transformation( self, jdata ) :
        'Set self.transform_to( defecttype ), self.transform_energy( unit: eV ) and self.pref_transform( unit: Hz ) for transformating DefectObject.'
        raise RuntimeError( 'transform_to, transform_energy and pref_transform have not been set for Class ' + self.defecttype )

    def transite_defecttype_decrease_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize       ) :
        'Tansite from type tmp_defecttype to type transite_type because decrease tmp_nsize get to the number range.'
        raise RuntimeError( 'function : transite_defecttype_decrease_nsize( ) has not been set for Class ' + self.defecttype )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class V detecting ---' )

    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
