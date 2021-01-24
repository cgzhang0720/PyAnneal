#!/usr/bin/env python3

import numpy as np
import math
import argparse
import json

from tungsten.DefectObject import DefectObject

from Auxiliary import j_have

class ICluster( DefectObject ) :
    'Class ICluster properties: 1. defecttype, 2. shape'

    defecttype = 'ICluster'
    shape = 'disk'

    '1. reshape.'
    def _set_migration( self ) :
        ''' Set self.migration_energy ( unit: eV ) and self.pref_migration( unit: Hz ) for migrated DefectObject. 
        Ref: C.S. Becquart et al., J. Nucl. Mater. 403(2010)75. '''

        self.migration_energy = 0.013
        v0 = 6.0E12
        s = 0.5
        self.pref_migration = v0 * self.nsize **( -s )

    '2. reshape.'
    def _set_rotation( self ) :
        ''' Set self.rotation_energy( unit: eV ) and self.pref_rotation( unit: Hz ) for rotating DefectObject.
        W.H. Zhou et al., J. Nucl. Mater. 453(2014)202–209. Potential: AT.'''

        if self.nsize == 2 : 
            self.rotation_energy = 1.307
            self.pref_rotation = 1.0E12
        elif self.nsize == 3 : 
            self.rotation_energy = 1.730
            self.pref_rotation = 1.0E12
        else :
            raise RuntimeError( 'The size of ICluster is not in the range [ 2, 3 ].' )

    '3. reshape.'
    def _reset_b_pn_rotation( self, tmp_ConstNumber ) :
        'Reset 1. self.burgers_vector_unit, 2. self.burgers_vector_scale and 3. self.position_next_list after rotation.'

        index = np.random.randint( 0, 3 )
        self.burgers_vector_unit[ index ] = -self.burgers_vector_unit[ index ]
        self.burgers_vector_scale[ index ] = -self.burgers_vector_scale[ index ]

        self._set_position_next_list( tmp_ConstNumber )

    '4. reshape.'
    def _set_emission( self, tmp_ConstNumber ) :
        ''' 
        Set self.emission_energy( unit: eV ) and self.pref_emission( unit: Hz ) for emiting DefectObject. 
        Binding energy of one SIA to a SIA cluster of size m according to the following reaction (m-1)I + I -> mI.
        emission_energy = 0.5 * binding_energy + migration_energy.
        pref_emission ~ nsize. Ref: C.S. Becquart et al., J. Nucl. Mater. 403(2010)75. 
        '''

        EmI = 0.013

        if self.nsize == 2 :
            binding_energy = 2.12
        elif self.nsize == 3 :
            binding_energy = 3.02
        else :
            raise RuntimeError( 'The size of ICluster: %d is not in the range [ 2, 3 ].' % self.nsize )

        'attention ! emission_energy = 0.5 * binding_energy + EmI for emission point defect.'
        self.emission_energy = 0.5 * binding_energy + EmI

        self.pref_emission = 6.0E12 * self.nsize

    '5. reshape.'
    def _set_shape_properties( self, jdata, tmp_ConstNumber ) :
        ''' Set 1. self.rho, 2. self.delta_rho, 3. self.radius, 4. self.area according to Class.shape : 1. particle, defecttype I.
        Ref: D.R. Mason et al., J. Phys.: Condens. Matter 26 (2014) 375701. '''

        self.rho = tmp_ConstNumber.alatt * math.sqrt( self.nsize / ( math.pi * math.sqrt( 3.0 ) ) )

        'default value: tmp_ConstNumber.alatt'
        self.delta_rho = tmp_ConstNumber.alatt

        if j_have( jdata, 'I_delta_rho' ) : self.delta_rho = jdata[ 'I_delta_rho' ]

        self.radius = self.rho + self.delta_rho

        self.area = math.pi * self.rho * self.rho

    '6. reshape.'
    def return_u12( self,
                    tmp_ConstNumber,
                    defectObject,
                    R_vector,
                    R_square         ) :
        'Return elastic interaction energy between two DefectObjects, unit: eV.'

        if defectObject.shape == 'disk' or defectObject.defecttype == 'I' :
            return self._return_u_disk_disk( tmp_ConstNumber, defectObject, R_vector, R_square )
        elif defectObject.shape == 'sphere' or defectObject.defecttype == 'V' :
            return defectObject._return_u_sphere_disk( tmp_ConstNumber, self, R_vector, R_square )
        else :
            return 0.0

    '7. reshape.'
    def transite_defecttype_increase_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize        ) :
        'Tansite from type tmp_defecttype to type transite_type because increase tmp_nsize >= item[ 0 ].'

        transite_increase_nsize = ( 4, 'ILoop111' )

        if tmp_nsize >= transite_increase_nsize[ 0 ] :
            return transite_increase_nsize[ 1 ]
        else :
            return tmp_defecttype

    '8. reshape.'
    def transite_defecttype_decrease_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize       ) :
        ''' Tansite from type tmp_defecttype to type transite_type because decrease tmp_nsize get to the number range, 
        i.e., tmp_nsize >= item[ 0 ] and tmp_nsize <= item[ 1 ], return item[ 2 ]. '''

        transite_decrease_nsize = [ ( 1, 1, 'I' ) ]

        for item in transite_decrease_nsize :
            if tmp_nsize >= item[ 0 ] and tmp_nsize <= item[ 1 ] : 
                return item[ 2 ]

        return tmp_defecttype

    '9. reshape'
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
                                  [ 'emission_energy',      self.emission_energy      ],
                                  [ 'pref_emission',        self.pref_emission        ],
                                  [ 'rotation_energy',      self.rotation_energy      ],
                                  [ 'pref_rotation',        self.pref_rotation        ],
                                  [ 'defecttype',           self.defecttype           ]  ] )

        return properties_list

    def _set_transformation( self, jdata ) :
        'Set self.transform_to, self.transform_energy and self.pref_transform for transformating DefectObject.'
        raise RuntimeError( 'transform_to, transform_energy and pref_transform have not been set for Class ' + self.defecttype )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class ICluster detecting ---' )

    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
