#!/usr/bin/env python3

import numpy as np
import math
import argparse
import json

from tungsten.DefectObject import DefectObject

from Auxiliary import j_have


class VCluster( DefectObject ) :
    'Class VCluster properties: 1. defecttype, 2. shape'

    defecttype = 'VCluster'
    shape = 'sphere'

    '1. reshape'
    def _set_migration( self ) :
        ''' Set self.migration_energy ( unit: eV ) and self.pref_migration( unit: Hz ) for migrated DefectObject. 
        nsize >= 2. Ref: C.S. Becquart et al., J. Nucl. Mater. 403(2010)75. '''

        self.migration_energy = 1.66

        v0 = 6.0E12
        q = 1000.0
        self.pref_migration = v0 * q ** ( 1 - self.nsize )

    '2. reshape.'
    def _set_emission( self, tmp_ConstNumber ) :
        ''' 
        Set self.emission_energy( unit: eV ) and self.pref_emission( unit: Hz ) for emiting DefectObject. 
        Binding energy of one Vacancy to a Vcluster of size m according to the following reaction (n-1)V + V -> nV, n begins from 2.
        The Capillary law is an empirical model assuming that defect cluster containing n atoms is a spherical object, 
        whose formatation energy is proportional to the surface of the cluster. 
        With this approximation, the binding energy Eb(n) corresponding to adding an atom to a cluster of size n-1 is given by 
        Eb( n -> ( n - 1 ) + 1 ) = Ef( n - 1 ) + Ef( 1 ) -Ef( n ) = Ef(1) + ( Eb(2) - Ef(1) ) / ( n^(2/3) -(n-1)^(2/3) ) / ( 2^(2/3) -1 ). 
        When n becomes large the asymptotic limit will be Ef(1).
        emission_energy = 0.5 * binding_energy + migration_energy. pref_emission ~ nsize. 
        Ref-1: C.S. Becquart et al., J. Nucl. Mater. 403(2010)75. for m in [ 2 : 8 ]
        Ref-2: J. Fikar et al., Nuclear Materials and Energy 16(2018)60–65. for nsize > 8.
        Ref-3: D. R. Mason et al., J. Phys.: Condens. Matter 29(2017)505501.                     

        a0=0.3165
        gamma_a=23.0
        b1=gamma_a*(9.0*3.1415926)**(1.0/3.0)*a0*a0
        void(x) = b1*x**(2.0/3.0)
        '''

        'm in [ 2 : 8 ]'
        binding_energy = [ -0.1, 0.04, 0.64, 0.72, 0.89, 0.72, 0.88 ]

        EfV = 3.23
        EmV = 1.66

        e0 = 2.0 / 3.0

        ''' Capillary approximation, Ref-1.
        if self.nsize > 8 :
            Ebind = EfV + ( binding_energy[ 0 ] - EfV ) * ( self.nsize **e0 - ( self.nsize - 1 ) **e0 ) / 0.587401
        else :
            Ebind = binding_energy[ self.nsize - 2 ]

        average surface energy, unit: eV/ nm^2, Ref-2.
        gamma_a = 23.0
        b1 = gamma_a * ( 9.0 * math.pi ) **( 1.0 / 3.0 ) * tmp_ConstNumber.alatt * tmp_ConstNumber.alatt  '''

        b1 = 70.068899133 * tmp_ConstNumber.alatt * tmp_ConstNumber.alatt
        En = b1 * self.nsize**e0
        En_1 = b1 * ( self.nsize - 1.0 )**e0
        Ebind = En_1 + 3.727 - En

        'attention ! emission_energy = 0.5 * binding_energy + EmV for emission point defect.'
        self.emission_energy = 0.5 * Ebind + EmV
        self.pref_emission = 6.0E12 * self.nsize

    '3. reshape'
    def _set_shape_properties( self, jdata, tmp_ConstNumber ) :
        ''' Set 1. self.rho( unit: nm ), 2. self.delta_rho( unit: nm ), 3. self.radius( unit: nm ), 4. self.relax_volume( nm^-3 ) 
        according to Class.shape : 1. particle, and defecttype.
        Ref: D.R. Mason et al., J. Phys.: Condens. Matter 26 (2014) 375701. P15 '''

        self.rho = tmp_ConstNumber.alatt * ( 3.0 * self.nsize / ( 8.0 * math.pi ) )**0.3333333

        'default value: 0.5*tmp_ConstNumber.alatt'
        self.delta_rho = 0.5 * tmp_ConstNumber.alatt

        if j_have( jdata, 'V_delta_rho' ) : self.delta_rho = jdata[ 'V_delta_rho' ]

        self.radius = self.rho + self.delta_rho

        self.relax_volume = -0.531 * self.nsize **( 2.0 / 3.0 ) * tmp_ConstNumber.volume

    '4. reshape'
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

    '5. reshape'
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

        transite_decrease_nsize = [ ( 1, 1, 'V' ) ]

        for item in transite_decrease_nsize :
            if tmp_nsize >= item[ 0 ] and tmp_nsize <= item[ 1 ] :
                return item[ 2 ]

        return tmp_defecttype

    '7. reshape'
    def _return_properties( self ) :

        properties_list = super( )._return_properties( )

        'Auxiliary properties.'
        properties_list.extend( [ [ 'migration_energy',  self.migration_energy  ],
                                  [ 'pref_migration',    self.pref_migration    ],
                                  [ 'rho',               self.rho               ],
                                  [ 'delta_rho',         self.delta_rho         ],
                                  [ 'radius',            self.radius            ],
                                  [ 'relax_volume',      self.relax_volume      ],
                                  [ 'emission_energy',   self.emission_energy   ],
                                  [ 'pref_emission',     self.pref_emission     ],
                                  [ 'defecttype',        self.defecttype        ]  ] )

        return properties_list

    def return_judge_recombine( self,
                                jdata,
                                tmp_defectObject ) :
        'Judge recombine or not for self.shape == sphere.'

        if tmp_defectObject.shape == 'sphere' or tmp_defectObject.defecttype == 'V' :
            return self._return_judge_recombine_sphere_sphere( jdata, tmp_defectObject )
        else :
            return self._return_judge_recombine_sphere_disk( jdata, tmp_defectObject )

    def _set_rotation( self ) :
        'Set self.rotation_energy( unit: eV ) and self.pref_rotation( unit: Hz ) for rotating DefectObject.'
        raise RuntimeError( 'rotation_energy and pref_rotation have not been set for Class ' + self.defecttype )

    def _reset_b_pn_rotation( self, tmp_ConstNumber ) :
        'Reset 1. self.burgers_vector_unit, 2. self.burgers_vector_scale and 3. self.position_next_list after rotation.'
        raise RuntimeError( 'function : _reset_b_pn_rotation has not been set for Class ' + self.defecttype )

    def _set_transformation( self, jdata ) :
        'Set self.transform_to( defecttype ), self.transform_energy( unit: eV ) and self.pref_transform( unit: Hz ) for transformating DefectObject.'
        raise RuntimeError( 'transform_to, transform_energy and pref_transform have not been set for Class ' + self.defecttype )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class VCluster detecting ---' )

    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
