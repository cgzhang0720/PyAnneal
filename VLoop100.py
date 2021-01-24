#!/usr/bin/env python3

import numpy as np
import math
import argparse
import json

from tungsten.DefectObject import DefectObject

from Auxiliary import j_have

class VLoop100( DefectObject ) :
    'Class VLoop100 properties: 1. defecttype, 2. shape'

    defecttype = 'VLoop100'
    shape = 'disk'

    '1. reshape.'
    def _set_emission( self, tmp_ConstNumber ) :
        ''' 
        Set self.emission_energy( unit: eV ) and self.pref_emission( unit: Hz ) for emiting DefectObject. 
        Binding energy of one Vacancy to a Vacancy loop 100 of size nsize according to the following reaction 
        (nsize-1) V + V -> nsize V, nsize begins from 2. emission_energy = 0.5 * binding_energy + migration_energy. 
        pref_emission ~ nsize. 
        Ref-1: C.S. Becquart et al., J. Nucl. Mater. 403(2010)75. for nsize in [ 2 : 8 ].
        Ref-2: J. Fikar et al., Nuclear Materials and Energy 16(2018)60–65. for nsize > 8.
        Ref-3: D. R. Mason et al., J. Phys.: Condens. Matter 29(2017)505501.                 

        vac-100-loop:
        Rc_vac(x)=a0*(x/(2.0*3.1415926))**0.5-delta_Rc
        vacloop100(x)=Rc_vac(x)*b*b*(c1_vac+c2*log(Rc_vac(x)))

        Final set of parameters            Asymptotic Standard Error
        =======================            ==========================
        c1_vac          = 1863.07          +/- 8.838        (0.4744%)
        '''

        'm in [ 2 : 8 ]'
        binding_energy = [ -0.1, 0.04, 0.64, 0.72, 0.89, 0.72, 0.88 ]

        EfV = 3.23
        EmV = 1.66

        'nsize > 8.'
        c1 = 1863.07
        c2 = 22.41 / ( tmp_ConstNumber.volume * 2.0 )
        delta_Rc = tmp_ConstNumber.alatt / 4.0

        'Rc_nsize = tmp_ConstNumber.alatt * math.sqrt( self.nsize / ( 2.0 * math.pi ) ) - delta_Rc'
        Rc_nsize = self.rho - delta_Rc

        Rc_nsize_1 = tmp_ConstNumber.alatt * math.sqrt( ( self.nsize - 1.0 ) / ( 2.0 * math.pi ) ) - delta_Rc

        b2 = self.burgers * self.burgers

        loop_nsize = Rc_nsize * b2 * ( c1 + c2 * math.log( Rc_nsize ) )
        loop_nsize_1 = Rc_nsize_1 * b2 * ( c1 + c2 * math.log( Rc_nsize_1 ) )

        if self.nsize > 8 :
            Ebind = 3.727 + loop_nsize_1 - loop_nsize
        else :
            Ebind = binding_energy[ self.nsize - 2 ]

        'attention ! emission_energy = 0.5 * binding_energy + EmV for emission point defect.'
        self.emission_energy = Ebind + EmV
        self.pref_emission = 6.0E12 * self.nsize

    '2. reshape.'
    def _set_transformation( self, jdata ) :
        ''' Set self.transform_to( defecttype ), self.transform_energy( unit: eV ) and self.pref_transform( unit: Hz ) for transformating DefectObject. 
        Ref: this study. or other Refs ?      '''
        
        self.transform_to = 'VCluster'
        self.transform_energy = 2.0

        v0 = 6.0E12
        s = 0.5
        self.pref_transform = v0 * self.nsize **( -s )

    '3. reshape.'
    def _set_shape_properties( self, jdata, tmp_ConstNumber ) :
        ''' Set 1. self.rho, 2. self.delta_rho, 3. self.radius, 4. self.area according to Class.shape : 1. particle, defecttype I.
        Ref: our study. rho = alatt * sqrt( nsize / ( 2 * pi ) ), delta_rho = alatt.   '''

        self.rho = tmp_ConstNumber.alatt * math.sqrt( self.nsize / ( math.pi * 2.0 ) )

        'default value: 0.5*tmp_ConstNumber.alatt'
        self.delta_rho = 0.5 * tmp_ConstNumber.alatt

        if j_have( jdata, 'V_delta_rho' ) : self.delta_rho = jdata[ 'V_delta_rho' ]

        self.radius = self.rho + self.delta_rho

        self.area = math.pi * self.rho * self.rho

    '4. reshape.'
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

    '5. reshape.'
    def transite_defecttype_increase_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize       ) :
        'Tansite from type tmp_defecttype to type transite_type because increase tmp_nsize >= item[ 0 ].'
        return tmp_defecttype

    '6. reshape.'
    def transite_defecttype_decrease_nsize( self,
                                            tmp_defecttype,
                                            tmp_nsize       ) :
        ''' Tansite from type tmp_defecttype to type transite_type because decrease tmp_nsize get to the number range, 
        i.e., tmp_nsize >= item[ 0 ] and tmp_nsize <= item[ 1 ], return item[ 2 ]. '''

        transite_decrease_nsize = [ ( 1, 1,  'V'        ),
                                    ( 2, 13, 'VCluster' )  ]

        for item in transite_decrease_nsize :
            if tmp_nsize >= item[ 0 ] and tmp_nsize <= item[ 1 ] : 
                return item[ 2 ]

        return tmp_defecttype

    '7. reshape'
    def _return_properties( self ) :

        properties_list = super( )._return_properties( )

        'Auxiliary properties.'
        properties_list.extend( [ [ 'rho',                  self.rho                  ],
                                  [ 'delta_rho',            self.delta_rho            ],
                                  [ 'radius',               self.radius               ],
                                  [ 'area',                 self.area                 ],
                                  [ 'burgers',              self.burgers              ],
                                  [ 'burgers_vector_scale', self.burgers_vector_scale ],
                                  [ 'burgers_vector_unit',  self.burgers_vector_unit  ],
                                  [ 'emission_energy',      self.emission_energy      ],
                                  [ 'pref_emission',        self.pref_emission        ],
                                  [ 'transform_to',         self.transform_to         ],
                                  [ 'transform_energy',     self.transform_energy     ],
                                  [ 'pref_transform',       self.pref_transform       ],
                                  [ 'defecttype',           self.defecttype           ]  ] )

        return properties_list

    def _set_migration( self ) :
        'Set self.migration_energy ( unit: eV ) and self.pref_migration( unit: Hz ) for migrated DefectObject.'
        raise RuntimeError( 'migration_energy and pref_migration have not been set for Class ' + self.defecttype )

    def _set_rotation( self ) :
        'Set self.rotation_energy( unit: eV ) and self.pref_rotation( unit: Hz ) for rotating DefectObject.'
        raise RuntimeError( 'rotation_energy and pref_rotation have not been set for Class ' + self.defecttype )

    def _reset_b_pn_rotation( self, tmp_ConstNumber ) :
        'Reset 1. self.burgers_vector_unit, 2. self.burgers_vector_scale and 3. self.position_next_list after rotation.'
        raise RuntimeError( 'position_next_list have not been set after rotation for Class ' + self.defecttype )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class VLoop100 detecting ---' )

    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
