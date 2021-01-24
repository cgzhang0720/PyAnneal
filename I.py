#!/usr/bin/env python3

import numpy as np
import math
import argparse
import json

from tungsten.DefectObject import DefectObject

from Auxiliary import j_have

class I( DefectObject ) :
    'Class I properties: 1. defecttype, 2. shape'

    defecttype = 'I'
    shape = 'particle'

    '1. reshape.'
    def _set_migration( self ) :
        'Set self.migration_energy ( unit: eV ) and self.pref_migration( unit: Hz ) for migrated DefectObject.'

        self.migration_energy = 0.013
        self.pref_migration = 6.0E12

    '2. reshape.'
    def _set_rotation( self ) :
        ''' Set self.rotation_energy( unit: eV ) and self.pref_rotation( unit: Hz ) for rotating DefectObject.
        Ref: C.S. Becquart et al., J. Nucl. Mater. 403(2010)75. '''

        self.rotation_energy = 0.38
        self.pref_rotation = 6.0E12

    '3. reshape.'
    def _reset_b_pn_rotation( self, tmp_ConstNumber ) :
        'Reset 1. self.burgers_vector_unit, 2. self.burgers_vector_scale and 3. self.position_next_list after rotation.'
        
        index = np.random.randint( 0, 3 )
        self.burgers_vector_unit[ index ] = -self.burgers_vector_unit[ index ]
        self.burgers_vector_scale[ index ] = -self.burgers_vector_scale[ index ]

        self._set_position_next_list( tmp_ConstNumber )

    '4. reshape.'
    def _set_shape_properties( self, jdata, tmp_ConstNumber ) :
        ''' Set 1. self.rho, 2. self.delta_rho, 3. self.radius, 4. self.area according to Class.shape : 1. particle, defecttype I.
        Ref: D.R. Mason et al., J. Phys.: Condens. Matter 26 (2014) 375701. '''

        self.rho = tmp_ConstNumber.alatt * math.sqrt( 1.0 / ( math.pi * math.sqrt( 3.0 ) ) )

        'default value: tmp_ConstNumber.alatt'
        self.delta_rho = tmp_ConstNumber.alatt

        if j_have( jdata, 'I_delta_rho' ) : self.delta_rho = jdata[ 'I_delta_rho' ]

        self.radius = self.rho + self.delta_rho

        self.area = math.pi * self.rho * self.rho

    '5. reshape.'
    def return_u12( self,
                    tmp_ConstNumber,
                    defectObject,
                    R_vector,
                    R_square         ) :
        'Return elastic interaction energy between two DefectObjects, unit: eV.'

        if defectObject.shape == 'disk' or defectObject.defecttype == 'I' :
            return self._return_u_disk_disk( tmp_ConstNumber, defectObject, R_vector, R_square )
        elif defectObject.shape == 'sphere' or defectObject.defecttype == 'V' :
            return defectObject._return_u_sphere_disk(tmp_ConstNumber, self, R_vector, R_square )
        else :
            return 0.0

    '6. reshape.'
    def transite_defecttype_increase_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize        ) :
        'Tansite from type tmp_defecttype to type transite_type because increase tmp_nsize >= item[ 0 ].'

        transite_increase_nsize = ( 2, 'ICluster' )

        if tmp_nsize >= transite_increase_nsize[ 0 ] :
            return transite_increase_nsize[ 1 ]
        else :
            return tmp_defecttype

    '7. reshape.'
    def _return_properties( self ) :

        properties_list = super( )._return_properties( )

        'Auxiliary properties.'
        properties_list.extend( [ [ 'migration_energy',     self.migration_energy     ],
                                  [ 'pref_migration',       self.pref_migration       ],
                                  [ 'rho',                  self.rho                  ],
                                  [ 'delta_rho',            self.delta_rho            ],
                                  [ 'radius',               self.radius               ],
                                  [ 'area',                 self.area                 ],
                                  [ 'burgers',              self.burgers              ],
                                  [ 'burgers_vector_scale', self.burgers_vector_scale ],
                                  [ 'burgers_vector_unit',  self.burgers_vector_unit  ],
                                  [ 'rotation_energy',      self.rotation_energy      ],
                                  [ 'pref_rotation',        self.pref_rotation        ],
                                  [ 'defecttype',           self.defecttype           ]  ] )

        return properties_list

    def _set_emission( self ) :
        'Set self.emission_energy( unit: eV ) and self.pref_emission( unit: Hz ) for emiting DefectObject.'
        raise RuntimeError( 'emission_energy and pref_emission have not been set for Class ' + self.defecttype )

    def _set_transformation( self, jdata ) :
        'Set self.transform_to, self.transform_energy and self.pref_transform for transformating DefectObject.'
        raise RuntimeError( 'transform_to, transform_energy and pref_transform have not been set for Class ' + self.defecttype )

    def transite_defecttype_decrease_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize       ) :
        'Tansite from type tmp_defecttype to type transite_type because decrease tmp_nsize get to the number range.'
        raise RuntimeError( 'function : transite_defecttype_decrease_nsize( ) has not been set for Class ' + self.defecttype )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class I detecting ---' )

    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
