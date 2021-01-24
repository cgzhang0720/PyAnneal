#!/usr/bin/env python3

import numpy as np
import argparse
import json

from Auxiliary import j_must_have, convert_int_array, judge_equal
from tungsten.Create import create_random_defectObject_list, create_defectObject_of_subclass

class BoxLinkcell( object ) :
    def __init__( self, 
                  jdata, 
                  tmp_defectObject_list ) :
        '''
        Properties of BoxLinkcell:
        No change:
         1.  self.box
         2.  self.linkcell_space
         3.  self.nlc_vector
         4.  self.fnlc_vector
         5.  self.nlc
         6.  self.periodic
         7.  self.nnbrs
         8.  self.hmeps
         9. self.nix, self.niy, self.niz
       
        Change :
         10. self.numb_defectObject
         11. self.ltop[ 0 : self.nlc ]
         12. self.linkmp[ 0 : self.numb_defectObject ]
        '''

        self.box = np.array( j_must_have( jdata, 'box' ) )

        self.linkcell_space = np.array( j_must_have( jdata, 'linkcell_space' ) )

        for item in self.linkcell_space > self.box :
            if item : raise RuntimeError( 'linkcell_space > box, wrong !' )

        rcutoff_elastic = j_must_have( jdata, 'rcutoff_elastic' )

        for item in self.linkcell_space < rcutoff_elastic :
            if item : raise RuntimeError( 'linkcell_space < rcutoff_elastic, wrong !' )

        'python3 np.divide, / and np.true_divide are same.'
        self.fnlc_vector = np.floor_divide( self.box, self.linkcell_space )

        self.nlc_vector = convert_int_array( self.fnlc_vector )

        self.nlc = self.nlc_vector[ 0 ] * self.nlc_vector[ 1 ] * self.nlc_vector[ 2 ]
       
        self.periodic = j_must_have( jdata, 'periodic' )
        
        self.nnbrs = j_must_have( jdata, 'max_numb_in_27_cells' )

        self.hmeps = -1e-9

        self.nix = [ 0, -1, -1, -1, 0, 0, -1, 1, -1,  0,  1, -1, 0, 1,  0, 1, 1,  1,  0,  0,  0,  1,  1,  1, -1, -1, -1 ]
        self.niy = [ 0,  0, -1,  1, 1, 0,  0, 0, -1, -1, -1,  1, 1, 1, -1, 0, 1, -1,  1,  0, -1,  0, -1,  1,  0, -1,  1 ]
        self.niz = [ 0,  0,  0,  0, 0, 1,  1, 1,  1,  1,  1,  1, 1, 1,  0, 0, 0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1 ]

        self.linked( tmp_defectObject_list )

    def linked( self, tmp_defectObject_list ) :
        'link self.defectObject_list'

        self.numb_defectObject = len( tmp_defectObject_list )

        self.ltop = -np.ones( ( self.nlc ), dtype = np.int )

        self.linkmp = -np.ones( ( self.numb_defectObject ), dtype = np.int )

        for i in range( self.numb_defectObject ) :
            ip = self.return_ip_defectObject( tmp_defectObject_list[ i ] )
           
            j = self.ltop[ ip ]
            self.ltop[ ip ] = i
            self.linkmp[ i ] = j

    def add_a_defectObject( self, tmp_defectObject ) :

        ip = self.return_ip_defectObject( tmp_defectObject )

        j = self.ltop[ ip ]
        self.ltop[ ip ] = self.numb_defectObject
        self.linkmp = np.append( self.linkmp, j )

        self.numb_defectObject += 1

    def add_defectObject_list( self, tmp_defectObject_list ) :

        for tmp_defectObject in tmp_defectObject_list :
           self.add_a_defectObject( tmp_defectObject )

    def substitute_defectObject( self,
                                 index,
                                 old_defectObject,
                                 tmp_defectObject  ) :

        'A. First delete index of old_defectObject in self.linkcell.'
        ip = self.return_ip_defectObject( old_defectObject )

        i = self.ltop[ ip ]

        if i == index :
            self.ltop[ ip ] = self.linkmp[ i ]
        else :
            while True :
                j = self.linkmp[ i ]
                if j == index :
                    self.linkmp[ i ] = self.linkmp[ j ]
                    break
                elif j == -1 :
                    raise RuntimeError( 'not find index %d in cell ip %d: ' % ( index, ip ) )
                i = j

        'B. Second add index of tmp_defectObject.'
        ip = self.return_ip_defectObject( tmp_defectObject )

        i = self.ltop[ ip ]

        if index > i :
            self.ltop[ ip ] = index
            self.linkmp[ index ] = i
        else :
            while True :
                j = self.linkmp[ i ]
                if index > j :
                    self.linkmp[ i ] = index
                    self.linkmp[ index ] = j
                    break
                i = j

    def delete_a_defectObject( self, 
                               tmp_defectObject, 
                               tmp_index         ) :

        if tmp_index < 0 or tmp_index > self.numb_defectObject - 1 :
            raise RuntimeError( 'tmp_index %d out of range [ 0, %d ] !' % ( tmp_index, self.numb_defectObject ) )

        ip = self.return_ip_defectObject( tmp_defectObject )

        i = self.ltop[ ip ]

        if i == tmp_index :
            self.ltop[ ip ] = self.linkmp[ i ]
        else :
            while True :
                j = self.linkmp[ i ]
                if j == tmp_index :
                    self.linkmp[ i ] = self.linkmp[ j ]
                    break
                i = j

        for i in range( self.nlc ) : 
            if self.ltop[ i ] > tmp_index :
                self.ltop[ i ] -= 1

        for i in range( self.numb_defectObject ) :
            self.linkmp[ i - int( i > tmp_index ) ] = self.linkmp[ i ] - int( self.linkmp[ i ] > tmp_index )

        self.linkmp = np.delete( self.linkmp, -1 )

        self.numb_defectObject -= 1

    def delete_defectObject_list( self,
                                  tmp_defectObject_list,
                                  tmp_index_list         ) :

       if len( tmp_defectObject_list ) != len( tmp_index_list ) :
           raise RuntimeError( 'The length of input parameters are not equal !' )

       tmp_add_list = [ [ tmp_index_list[ i ], tmp_defectObject_list[ i ] ] for i in range( len( tmp_defectObject_list ) ) ]

       'decreasing tmp_index_list[ index ] according index.'
       tmp_add_list.sort( reverse = True )

       for item in tmp_add_list :
           self.delete_a_defectObject( item[ 1 ], item[ 0 ] )

    def return_ip_defectObject( self, tmp_defectObject ) :

        cell_index = self.return_cell_index_defectObject( tmp_defectObject )
        ip = self.return_ip_cell_index( cell_index )

        return ip

    def return_cell_index_defectObject( self, tmp_defectObject ) :
        return convert_int_array( ( tmp_defectObject.return_position_scale( ) + self.hmeps ) * self.fnlc_vector )

    def return_ip_cell_index( self, cell_index ) :

        ip = cell_index[ 0 ] + self.nlc_vector[ 0 ] * ( cell_index[ 1 ] + self.nlc_vector[ 1 ] * cell_index[ 2 ] )

        if ip > self.nlc - 1 or ip < 0 :
            print( 'Error in ip:', ip, 'not in [ 0, %d ]' %d ( self.nlc ) )
            print( 'return_cell_index: ', cell_index )

        return ip

    def _return_properties( self ) :
        '''
        Properties of BoxLinkcell:
        No change :
         1.  self.box
         2.  self.linkcell_space
         3.  self.nlc_vector
         4.  self.fnlc_vector
         5.  self.nlc
         6.  self.periodic
         7.  self.nnbrs
         8.  self.hmeps
         9.  self.nix, self.niy, self.niz
       
        Change :
         10. self.numb_defectObject
         11. self.ltop[ 0 : self.nlc ]
         12. self.linkmp[ 0 : self.numb_defectObject ]
        '''

        properties_list = [ [ 'box',                 self.box               ],  
                            [ 'linkcell_space',      self.linkcell_space    ],
                            [ 'nlc_vector',          self.nlc_vector        ],
                            [ 'fnlc_vector',         self.fnlc_vector       ],
                            [ 'nlc',                 self.nlc               ],
                            [ 'periodic',            self.periodic          ],
                            [ 'nnbrs',               self.nnbrs             ],
                            [ 'hmeps',               self.hmeps             ],
                            [ 'nix',                 self.nix               ],
                            [ 'niy',                 self.niy               ],
                            [ 'niz',                 self.niz               ],
                            [ 'numb_defectObject',   self.numb_defectObject ],
                            [ 'ltop',                self.ltop              ],
                            [ 'linkmp',              self.linkmp            ] ]

        return properties_list
                                                                                 
    def equal( self,                                                              
               tmp_linkcell,
               delta         ) :
        'compare the value between tmp_linkcell and self.'

        self_properties = self._return_properties( )        
        tmp_properties = tmp_linkcell._return_properties( )        

        for i in range( len( self_properties ) ) :
            if not judge_equal( self_properties[ i ][ 1 ], tmp_properties[ i ][ 1 ], delta ) : 
                print( '# ' + self_properties[ i ][ 0 ] + ' not equal !' )
                return False

        return True

    def print_properties( self ) :

        print( '#   ' )
        print( '# ------------- properties for BoxLinkcell --------------' )

        print( '\n'.join( ( '# %s : %s' % ( item[ 0 ], item[ 1 ] ) for item in self._return_properties( ) ) ) )

        print( '#   ' )

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class BoxLinkcell detecting ---' )

    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    '''
    tmp_periodic = np.array( j_must_have( jdata, 'periodic' ) )

    delta = 1.0e-6
 
    numb_defectObject = 10
    defectObject_list = create_random_defectObject_list( jdata, numb_defectObject )
    old_defectObject_list = [ item for item in defectObject_list ]

    linkcell_1 = BoxLinkcell( jdata, defectObject_list )
    linkcell_2 = BoxLinkcell( jdata, old_defectObject_list ) 
   
    #1. add_a_defectObject( self, tmp_defectObject ) and linked( defectObject_list )
    # A. real value
    a_defectStr = 'ILoop111 13 20.0 20.0 20.0   0.57735 0.57735 0.57735'
    a_defectObject = create_defectObject_of_subclass( jdata, a_defectStr )
    defectObject_list.append( a_defectObject )
    linkcell_1.linked( defectObject_list )

    # B. add_a_defectObject
    linkcell_2.add_a_defectObject( a_defectObject )


    #2. delete_a_defectObject( self, tmp_defectObject, tmp_index )
    # A. real value
    del defectObject_list[ 3 ]
    linkcell_1.linked( defectObject_list )

    # B. delete_a_defectObject
    linkcell_2.delete_a_defectObject( old_defectObject_list[ 3 ], 3 )


    if not linkcell_1.equal( linkcell_2, delta ) :
        raise RuntimeError( '# linkcell_1 not equal linkcell_2, wrong !' )
    else :
        print( '# every is OK !' )
    '''
