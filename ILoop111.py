#!/usr/bin/env python3

import numpy as np
import math
import argparse
import json

from tungsten.DefectObject import DefectObject

from Auxiliary import j_have

class ILoop111( DefectObject ) :
    'Class ILoop111 properties: 1. defecttype, 2. shape'

    defecttype = 'ILoop111'
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
    def _set_emission( self, tmp_ConstNumber ) :
        ''' 
        Set self.emission_energy( unit: eV ) and self.pref_emission( unit: Hz ) for emiting DefectObject. 
        Binding energy of one SIA to a SIA cluster of size nsize according to the following reaction 
        (nsize-1) I + I -> nsize I, nsize begins from 2. emission_energy = 0.5 * binding_energy + migration_energy. 
        pref_emission ~ nsize. 
        Ref-1: C.S. Becquart et al., J. Nucl. Mater. 403(2010)75. for nsize in [ 2 : 7 ].
        Ref-2: J. Fikar et al., Nuclear Materials and Energy 16(2018)60–65. for nsize > 7.
        Ref-3: D. R. Mason et al., J. Phys.: Condens. Matter 29(2017)505501.                 

        int-111-loop
        b=3.0**0.5/2.0*a0
        c2=22.41/(a0*a0*a0)
        delta_Rc=a0/(2.0*6.0**0.5)
        Rc_int(x)=a0*(x/(3.0**0.5*3.1415926))**0.5+delta_Rc
        formation energy of size x( eV ) : intloop111(x)=Rc_int(x)*b*b*(c1_int+c2*log(Rc_int(x)))

        Final set of parameters            Asymptotic Standard Error
        =======================            ==========================
        c1_int          = 1845.64          +/- 2.58         (0.1398%)
        '''
   
        if self.nsize < 4 : raise RuntimeError( 'The size of ILoop111 is not in the range [ 4 : ].' )

        'nsize in [ 2 : 7 ]'
        binding_energy = [ 2.12, 3.02, 3.60, 3.98, 4.27, 5.39 ]
        EfI = 9.96
        EmI = 0.013

        'nsize > 7.' 
        c1 = 1845.64 
        c2 = 22.41 / ( tmp_ConstNumber.volume * 2.0 )

        delta_Rc = tmp_ConstNumber.alatt / ( 2.0 * math.sqrt( 6.0 ) )

        'Rc_nsize = tmp_ConstNumber.alatt * math.sqrt( self.nsize / ( math.sqrt( 3.0 ) * math.pi ) ) + delta_Rc'
        Rc_nsize = self.rho + delta_Rc

        Rc_nsize_1 = tmp_ConstNumber.alatt * math.sqrt( ( self.nsize - 1.0 ) / ( math.sqrt( 3.0 ) * math.pi ) ) + delta_Rc

        b2 = self.burgers * self.burgers 

        loop_nsize = Rc_nsize * b2 * ( c1 + c2 * math.log( Rc_nsize ) )
        loop_nsize_1 = Rc_nsize_1 * b2 * ( c1 + c2 * math.log( Rc_nsize_1 ) )

        if self.nsize > 7 :
            Ebind = 9.31 + loop_nsize_1 - loop_nsize
        else :
            Ebind = binding_energy[ self.nsize - 2 ]

        'attention ! emission_energy = 0.5 * binding_energy + EmI for emission point defect.'
        self.emission_energy = Ebind + EmI
        self.pref_emission = 6.0E12 * self.nsize

    '3. reshape.'
    def _set_shape_properties( self, jdata, tmp_ConstNumber ) :
        ''' Set 1. self.rho, 2. self.delta_rho, 3. self.radius, 4. self.area according to Class.shape : 1. particle, defecttype I.
        Ref: D.R. Mason et al., J. Phys.: Condens. Matter 26(2014)375701. '''

        self.rho = tmp_ConstNumber.alatt * math.sqrt( self.nsize / ( math.pi * math.sqrt( 3.0 ) ) )

        'default value: tmp_ConstNumber.alatt'
        self.delta_rho = tmp_ConstNumber.alatt

        if j_have( jdata, 'I_delta_rho' ) : self.delta_rho = jdata[ 'I_delta_rho' ]

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

        transite_decrease_nsize = [ ( 1, 1, 'I'        ),
                                    ( 2, 3, 'ICluster' )  ]

        for item in transite_decrease_nsize :
            if tmp_nsize >= item[ 0 ] and tmp_nsize <= item[ 1 ] : 
                return item[ 2 ]

        return tmp_defecttype

    '7. reshape'
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
                                  [ 'defecttype',           self.defecttype           ]  ] )

        return properties_list

    def _set_rotation( self ) :
        'Set self.rotation_energy( unit: eV ) and self.pref_rotation( unit: Hz ) for rotating DefectObject.'
        raise RuntimeError( 'rotation_energy and pref_rotation have not been set for Class ' + self.defecttype )

    def _reset_b_pn_rotation( self, tmp_ConstNumber ) :
        'Reset 1. self.burgers_vector_unit, 2. self.burgers_vector_scale and 3. self.position_next_list after rotation.'
        raise RuntimeError( 'position_next_list have not been set after rotation for Class ' + self.defecttype )

    def _set_transformation( self, jdata ) :
        'Set self.transform_to, self.transform_energy and self.pref_transformation for transformating DefectObject.'
        raise RuntimeError( 'transform_to, transform_energy and pref_transformation have not been set for Class ' + self.defecttype )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class ILoop111 detecting ---' )

    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
