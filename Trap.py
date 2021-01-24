#!/usr/bin/env python3

import numpy as np
import argparse
import json

from BoxLinkcell import BoxLinkcell
from ConstNumber import ConstNumber
from Auxiliary import j_must_have, j_have

class TrapObject( object ) :

    def __init__( self, 
                  jdata,
                  tmp_type, 
                  tmp_position,
                  tmp_radius,
                  tmp_energy   ) :
        '''
        Class TrapObject properties:
         1. self.traptype
         2. self.position
         3. self.position_scale
         4. self.radius
         5. self.radius_square
         6. self.energy
        '''

        self.traptype = tmp_type

        self.position = np.array( [ item for item in tmp_position ] )

        self.set_position_scale( jdata )

        self.radius = tmp_radius

        self.radius_square = self.radius * self.radius

        self.energy = tmp_energy

    def return_traptype( self ) :
        return self.traptype

    def return_position( self ) :
        'All traps are in box, so only have property self.position, no self.position_image_int.'
        return self.position

    def return_position_scale( self ) :
        return self.position_scale

    def return_radius( self ) :
        return self.radius

    def return_radius_square( self ) :
        return self.radius_square
    
    def return_energy( self ) :
        return self.energy

    def return_judge_trap( self,
                           jdata,
                           tmp_defectObject ) :
        'return_judge_trap between a TrapObject and a DefectObject.'

        delta = self.return_delta( jdata, tmp_defectObject )

        for item in delta :
            if item > self.radius :
                return False

        return np.inner( delta, delta ) < self.radius_square

    def return_delta( self,
                      jdata,
                      tmp_defectObject ) :
        'return delta = ( dx, dy, dz ) between a TrapObject and a DefectObject, according to periodic return nearest image.'

        tmp_periodic = np.array( j_must_have( jdata, 'periodic' ) )
        tmp_box = np.array( j_must_have( jdata, 'box' ) )

        p1 = self.return_position( )
        p2 = tmp_defectObject.return_position( )

        tmp_delta = p1 - p2

        tmp_nearest = [ ]
        for i in range( 3 ) :
            if tmp_periodic[ i ] :
                tmp_ =  min( [ abs( tmp_delta[ i ] + j * tmp_box[ i ] ) for j in range( -1, 2 ) ] )
            else :
                tmp_ = abs( tmp_delta[ i ] )
            tmp_nearest.append( tmp_ )
        delta = np.array( tmp_nearest )

        return delta

    def set_position_scale( self, jdata ) :
        'set the value of self.position_scale.'

        tmp_box = np.array( j_must_have( jdata, 'box' ) )
        self.position_scale = np.divide( self.position, tmp_box )

    def print_properties( self ) :
        '''
        Class TrapObject properties:
         1. self.traptype
         2. self.position
         3. self.position_scale
         4. self.radius
         5. self.radius_square
         6. self.energy
        '''

        print( '# traptype: ', self.traptype )

        print( '# position:', self.position )
        print( '# position_scale:', self.position_scale )

        print( '# radius: ', self.radius )
        print( '# radius_square: ', self.radius_square )
        print( '# energy: ', self.energy )

class TrapSystem( object ) :

    def __init__( self, jdata ) :
        '''
        Class TrapSystem properties:
         1. self.numb_trapObject
         2. self.trapObject_list
         3. self.linkcell
        '''

        ConstNumb = ConstNumber( jdata )

        density = j_must_have( jdata, 'trap_density' )
        radius = j_must_have( jdata, 'trap_radius' )
        energy = j_must_have( jdata, 'trap_energy' )

        if len( density ) != len( radius ) or len( density ) != len( energy ) :
            raise RuntimeError( 'the length of trap_density, trap_radius and trap_energy are not equal !' )

        tmp_box = np.array( j_must_have( jdata, 'box' ) )
        box_volume = tmp_box[ 0 ] * tmp_box[ 1 ] * tmp_box[ 2 ]

        self.trapObject_list = [ ]

        for i in range( len( density ) ) :

            if density[ i ][ 0 ] == 'Volume' :
                numb = int( density[ i ][ 1 ] * box_volume )
            elif density[ i ][ 0 ] == 'Atomic' :
                numb = int( density[ i ][ 1 ] * 1.0E-6 / ConstNumb.volume * box_volume )

            # Test
            #print( 'total traps:', numb )

            for j in range( numb ) :
                position = np.array( [ np.random.rand( ) * l for l in tmp_box ] )
                self.trapObject_list.append( TrapObject( jdata, i, position, radius[ i ], energy[ i ] ) )

        self.numb_trapObject = len( self.trapObject_list )

        if self.numb_trapObject == 0 :
            raise RuntimeError( 'self.numb_trapObject eq 0 while trap eq True' )

        self.linkcell = BoxLinkcell( jdata, self.trapObject_list )

    def return_judge_trap( self, 
                           jdata,
                           tmp_defectObject ) :

        if self.numb_trapObject < self.linkcell.nnbrs :
            select_trapObject_list = self.trapObject_list
        else :
            select_trapObject_list = self.return_cells_trapObject_list( jdata, tmp_defectObject )

        for item in select_trapObject_list :
           if item.return_judge_trap( jdata, tmp_defectObject ) : return ( True, item )

        return ( False, None )

    def return_cells_trapObject_list( self, 
                                      jdata, 
                                      tmp_defectObject ) :
        'return cells_trapObject_list in 27 cells of the tmp_defectObject.'

        cell_index = self.linkcell.return_cell_index_defectObject( tmp_defectObject )
        ip = self.linkcell.return_ip_cell_index( cell_index ) 

        i = self.linkcell.ltop[ ip ]
   
        cells_trapObject_list = [ ]

        if i != -1 :
            while True :
                cells_trapObject_list.append( self.trapObject_list[ i ] )
                i = self.linkcell.linkmp[ i ]
                if i == -1 : break

        for kc in range( 1, 27 ) :

            jx = cell_index[ 0 ] + self.linkcell.nix[ kc ]
            jy = cell_index[ 1 ] + self.linkcell.niy[ kc ]
            jz = cell_index[ 2 ] + self.linkcell.niz[ kc ]
           
            if cell_index[ 0 ] == self.linkcell.nlc_vector[ 0 ] - 1 and jx > cell_index[ 0 ] :
                jx = 0
                if not self.linkcell.periodic[ 0 ] : break
            elif cell_index[ 0 ] == 0 and jx < cell_index[ 0 ] :
                jx = self.linkcell.nlc_vector[ 0 ] - 1
                if not self.linkcell.periodic[ 0 ] : break
           
            if cell_index[ 1 ] == self.linkcell.nlc_vector[ 1 ] - 1 and jy > cell_index[ 1 ] :
                jy = 0
                if not self.linkcell.periodic[ 1 ] : break
            elif cell_index[ 1 ] == 0 and jy < cell_index[ 1 ] :
                jy = self.linkcell.nlc_vector[ 1 ] - 1
                if not self.linkcell.periodic[ 1 ] : break
           
            if cell_index[ 2 ] == self.linkcell.nlc_vector[ 2 ] - 1 and jz > cell_index[ 2 ] :
                jz = 0
                if not self.linkcell.periodic[ 2 ] : break
            elif cell_index[ 2 ] == 0 and jz < cell_index[ 2 ] :
                jz = self.linkcell.nlc_vector[ 2 ] - 1
                if not self.linkcell.periodic[ 2 ] : break

            jc = self.linkcell.return_ip_cell_index( [ jx, jy, jz ] )

            j = self.linkcell.ltop[ jc ]

            'Bypass this neighbouring cell if it is empty.'
            if j != -1 :
                while True :
                    cells_trapObject_list.append( self.trapObject_list[ j ] )
                    j = self.linkcell.linkmp[ j ]
                    if j == -1 : break

        if len( cells_trapObject_list ) > self.linkcell.nnbrs :
            raise RuntimeError( 'max1: %d .gt. nnbrs: %d something wrong' % ( len( cells_trapObject_list ), self.linkcell.nnbrs ) )

        return cells_trapObject_list

    def print_properties( self ) :
        print( '# numb_trapObject: %d' % ( self.numb_trapObject ), 'numb_defectObject: %d' % ( self.linkcell.numb_defectObject), \
               'ltop:', self.linkcell.ltop , 'linkmp:', self.linkcell.linkmp )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class TrapSystem detecting ---' )
    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
