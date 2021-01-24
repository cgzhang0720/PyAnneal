#!/usr/bin/env python3
'This is an Auxiliary module of functions.'

import numpy as np
import math
import argparse
import json
import os
import shutil

from scipy import stats
import matplotlib.pyplot as plt
from collections import Iterable

def function_reshape( x, x0 ) :
    'function reshape, default x0 = 2.0, units: nm.'

    return 1.0 - math.exp( -x / x0 )

def function_original( x, x_, x0 ) :

    return x_ * ( 1.0 - math.exp( -x / x0 ) )

def solve_function( x_, x0, delta = 1.0E-10 ) :

    if x_ < 0.0 or x0 < 0.0 : raise RuntimeError( 'x_ or x0 .lt. 0, wrong!' )

    'when x_ is very closer to 0.0, then return 0.'
    if x_ < delta : return 0.0

    x = x_
    while True :

        x_next = function_original( x, x_, x0 )
        print( 'x_next:', x_next )
        dx = abs( ( 1.0 - math.exp( -x_next / x0 ) ) - x_next / x_ )
        #dx = abs( x - x_next )

        if dx < delta : 
            print( 'x_:', x_, ', x:', x_next, ', dx:', dx, ', delta:', abs( x_next / ( 1.0 - math.exp( -x_next / x0 ) ) - x_ ) )
            #print( 'x_:', x_, 'x:', x, 'delta:', dx )
            return x_next
        else :
            x = x_next

def save_float_several_digits( data_float_str, several_digits ) :
    '''
    Save data_float_str with several_digits to save space when output to data files.
    data_float_str: the float data handed, several_digits: save the number of digits.
    '''

    "if several_digits < 0 : raise RuntimeError( '# Several_digits: %d is .lt. 0.' % several_digits )"

    point_index = data_float_str.find( '.' )

    'not find .'
    if point_index == -1 : return data_float_str

    end_index = min( len( data_float_str ), point_index + several_digits + 1 )

    return data_float_str[ : end_index ]

def output_data_into_file( dataStr_list, filenameStr ) :
    'Output dataStr_list into filenameStr, according to dataStr_list is a str or str_list.'

    'the index of the last /.'
    index_end = filenameStr.rfind( '/' )
    'detect the dir of filenameStr exists or not.'
    if not os.path.exists( filenameStr[ 0 : index_end ] ) :
        os.makedirs( filenameStr[ 0 : index_end ], exist_ok = True )

    file_out = open( filenameStr, 'w' )

    if type( dataStr_list ) == list :
        file_out.writelines( dataStr_list )
    elif type( dataStr_list ) == str :
        file_out.writeline( [ dataStr_list ] )
    else :
        raise RuntimeError( '# the type of dataStr_list is wrong, please choose str or str_list.' )

    file_out.close( )

def return_serial_numb_in_str( father_str, key_char = '.' ) :
    ''' 
    Take out all the numbers that after key_char in the string path, 
    connect them into a new string number and convert it into integer number.
    '''

    father_str_reverse = ''

    for i in range( len( father_str ) ) :
        father_str_reverse = father_str_reverse + father_str[ len( father_str ) - 1 - i ]

    'the index of last char / in string father_str.'
    index_ = len( father_str ) - 1 - father_str_reverse.find( '/' )

    index_begin = father_str.find( key_char, index_ )

    numb_str = ''
    for i in range( len( father_str ) ) :
        if father_str[ i ] >= '0' and father_str[ i ] <= '9' and i > index_begin :
            numb_str = numb_str + father_str[ i ]

    return int( numb_str )

def inner_product( position ) :
    'Return inner product of position.'
    return np.inner( position, position )

def my_str_match( str1, str_list ) :
    'Return str1 match element of str_list or not.'

    for item in str_list :
        if str1 in item :
            return True

    return False

def judge_equal( obj1, obj2, tol ) :
    '''
    judge two object equal or not
    1. number
    2. list
    3. np.array
    4. dict
    '''

    if isinstance( obj1, dict ) and isinstance( obj2, dict ) :
        return obj1 == obj2

    if not isinstance( obj1, np.ndarray ) :

        if isinstance( obj1, list ) :
            obj1 = np.array( obj1 )
        else :
            obj1 = np.array( [ obj1 ] )

    if not isinstance( obj2, np.ndarray ) :

        if isinstance( obj2, list ) :
            obj2 = np.array( obj2 )
        else :
            obj2 = np.array( [ obj2 ] )

    return _judge_equal_array( obj1, obj2, tol )

def _judge_equal_array( obj1, obj2, tol ) :
    'compare two array obj1 and obj2, tolerance < tol.'

    if isinstance( obj1[ 0 ], np.bool_ ) :
       obj1 = np.array( [ int( item ) for item in obj1 ] )

    if isinstance( obj2[ 0 ], np.bool_ ) :
       obj2 = np.array( [ int( item ) for item in obj2 ] )

    '''
    # Test_cgzhang, detecting the index of not equal data.
    obj = obj1 - obj2
    for i in range( obj.size ) :
        if abs( obj[ i ] ) > tol : print( 'obj[%d]: %f' %( i, obj[ i ] ) )
    '''

    return all( np.abs( ( obj1 - obj2 ) ) < tol )

def fun( x ) :
    'integral function'

    return x * x + math.sin( x )

def integral_legendre_Gauss( a, b, tol = 1.0e-10 ) :
    'legendre-Gauss integral method, integral f(x)dx from a to b, with tolerance tol.'

    t = [ -0.9061798459, -0.5384693101, 0.0,          0.5384693101,  0.9061798459 ]
    c = [  0.2369268851,  0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851  ]

    m = 1
    s = ( b - a ) * 0.001
    p = 0.0

    while True :

        h = ( b - a ) / m
        g = 0.0

        for i in range( m ) :

            aa = a + i * h
            bb = a + ( i + 1 ) * h
            w = 0.0

            for j in range( 5 ) :
                x = ( ( bb - aa ) * t[ j ] + ( bb + aa ) ) / 2.0
                w = w + fun( x ) * c[ j ]

            g = g + w

        g = g * h / 2.0
        q = abs( g - p ) / ( 1.0 + abs( g ) )

        if ( ( q >= tol ) and ( abs( h ) > abs( s ) ) ) :
            p = g
            m = m + 1
        else :
            break

    return g

def first_elliptic_integral( MM, eps = 1.0e-10 ) :
    '''
    K( k ), the elliptic integral of the first kind.
    K( k ) = integral from 0 to 0.5*pi for 1.0 / math.sqrt( 1.0 - k * k * math.sin( x ) * math.sin( x ) ) * dx
    k in [ 0, 1 ].
    K( k ) = first_elliptic_integral( MM ), here MM = k * k.
    '''
 
    t = [ -0.9061798459, -0.5384693101, 0.0,          0.5384693101,  0.9061798459 ]
    c = [  0.2369268851,  0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851  ]
 
    a = 0.0
    b = math.pi * 0.5
    m = 1
    s = ( b - a ) * 0.001
    p = 0.0 
    
    while True :

        h = ( b - a ) / m 
        g = 0.0
      
        for i in range( m ) :
      
            aa = a + i * h 
            bb = a + ( i + 1 ) * h 
            w = 0.0
      
            for j in range( 5 ) :
                x = ( ( bb - aa ) * t[ j ] + ( bb + aa ) ) / 2.0
                fk = 1.0 / math.sqrt( 1.0 - MM * math.sin( x ) * math.sin( x ) )
                w = w + fk * c[ j ]
       
            g = g + w 
      
        g = g * h / 2.0
        q = abs( g - p ) / ( 1.0 + abs( g ) )
     
        if ( ( q >= eps ) and ( abs( h ) > abs( s ) ) ) :
            p = g 
            m = m + 1
        else :
            break
 
    return g
        
def second_elliptic_integral( MM, eps = 1.0e-10 ) :
    '''
    E( k ), the elliptic integral of the second kind.
    E( k ) = integral from 0 to 0.5*pi for math.sqrt( 1.0 - k * k * math.sin( x ) * math.sin( x ) ) * dx
    k in [ 0, 1 ].
    E( k ) = second_elliptic_integral( MM ), here MM = k * k.
    '''

    t = [ -0.9061798459, -0.538469310, 0.0,          0.538469310,  0.9061798459 ]
    c = [ 0.2369268851,  0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851 ]

    a = 0.0
    b = math.pi * 0.5
    m = 1
    s = ( b - a ) * 0.02
    p = 0.0
    
    while True :
        h = ( b - a ) / m 
        g = 0.0

        for i in range( m ) :
            aa = a + i * h 
            bb = a + ( i + 1 ) * h 
            w = 0.0

            for j in range( 5 ) :
                x = ( ( bb - aa ) * t[ j ] + ( bb + aa ) ) / 2.0
                ek = math.sqrt( 1.0 - MM * math.sin( x ) * math.sin( x ) )
                w = w + ek * c[ j ]
        
            g = g + w 

        g = g * h / 2.0
        q = abs( g - p ) / ( 1.0 + abs( g ) )
    
        if ( ( q >= eps ) and ( abs( h ) > abs( s ) ) ) :
            p = g 
            m = m + 1
        else :
            break

    return g

def create_random_unit_vector_array( seed = None ) :
    'Return unit vector randomly.'

    if seed != None : np.random.seed( seed )

    theta = np.random.uniform( 0.0, math.pi )
    phi = np.random.uniform( 0.0, 2.0 * math.pi )

    v1 = math.sin( theta ) * math.cos( phi )
    v2 = math.sin( theta ) * math.sin( phi )
    v3 = math.cos(theta )

    return np.array( [ v1, v2, v3 ] )

def create_random_unit_vector_array_vertical( unit_vector_array, seed = None ) :
    'Randomly return unit vector array vertical to unit_vector_array.'

    unit_ = create_random_unit_vector_array( seed )

    return np.cross( unit_, unit_vector_array )

def create_random_unit_111_direction( seed = None ) :

    if seed != None : np.random.seed( seed )
    
    return create_unit_111_direction( np.random.randint( 0, 4 ) )

def create_random_unit_100_direction( seed = None ) :

    if seed != None : np.random.seed( seed )

    return create_unit_100_direction( np.random.randint( 0, 3 ) )

def create_unit_100_direction( index ) :

    choose = [ [ 1.0, 0.0, 0.0 ],
               [ 0.0, 1.0, 0.0 ],
               [ 0.0, 0.0, 1.0 ] ]

    return np.array( choose[ index ] )

def create_unit_111_direction( index ) :

    choose = [ [  0.57735,  0.57735,  0.57735 ],
               [ -0.57735,  0.57735,  0.57735 ],
               [  0.57735, -0.57735,  0.57735 ],
               [  0.57735,  0.57735, -0.57735 ]  ]

    return np.array( choose[ index ] )

def j_must_have( jdata, key ) :

    if not key in jdata.keys( ) :
        raise RuntimeError( 'json database must provide key ' + key )
    else :
        return jdata[ key ]

def j_have( jdata, key ) :
    return key in jdata.keys( )

def to_str( before ) :

    tmp_str = ''

    for item in before :
        tmp_str += ' ' + str( item ) + ' '

    return tmp_str

def convert_int_array( obj ) :
    return np.array( list( map( int, obj ) ) )

def convert_float_array( obj ) :
    return np.array( list( map( float, obj ) ) )

def R_convert_to_position( R ) :

    theta = np.random.rand( ) * math.pi
    phi = np.random.rand( ) * math.pi * 2.0

    x = R * math.sin( theta ) * math.cos( phi )
    y = R * math.sin( theta ) * math.sin( phi )
    z = R * math.cos( theta )

    return np.array( [ x, y, z ] )

def print_object( obj ) :
    print( '\n'.join( [ '%s:%s' % item for item in obj.__dict__.items( ) ] ) )

if __name__ == '__main__' :

    dx = 0.01
    x0 = 2.0

    for i in range( 1000 ) :
        x_ = i * dx
        x = solve_function( x_, x0, delta = 1.0E-10 )


    '''
    parser = argparse.ArgumentParser( description = '--- Moduler Auxiliary detecting ---' )
    parser.add_argument( 'INPUT', help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    tmp_periodic = np.array( j_must_have( jdata, 'periodic' ) )
    '''


    '''
    # 4. test create_normal_random_numb( mu, sigma, nsize ) and create_truncnormal_random_numb( mu, sigma, nsize, lower )
    fig, ax = plt.subplots( 1, 1 )

    nsize = 100000

    mu = 5.0
    sigma = 5.0

    lower = 0.0

    x = np.linspace( 0.0, 20.0, 10000 )

    ax.plot( x, stats.truncnorm.pdf( x, ( lower - mu ) / sigma, np.inf, mu, sigma ), 'r-', lw = 5, alpha = 0.6, label = 'truncnorm pdf' )
    ax.plot( x, stats.norm.pdf( x, mu, sigma ), 'b-', lw = 5, alpha = 0.6, label = 'norm pdf' )

    defect_R_trunc = create_truncnormal_random_numb( mu, sigma, nsize, lower )
    defect_R = create_normal_random_numb( mu, sigma, nsize )

    ax.hist( defect_R_trunc, 20, density = True, histtype = 'stepfilled', alpha = 0.2 )
    ax.hist( defect_R, 20, density = True, histtype = 'stepfilled', alpha = 0.2 )

    ax.legend( loc = 'best', frameon = False )

    plt.show( )

    '''


    '''
    # 3. test create_powerlaw_random_numb_2( b, nsize )
    fig, ax = plt.subplots( 1, 1 )

    nsize = 100000
    b = 0.63

    #nmax = 10
    #nmax = 20
    nmax = 30
    #nmax = 40

    x = np.linspace( 1.0, nmax, 10000 )
    ax.plot( x, stats.pareto.pdf( x, b ), 'r-', lw = 5, alpha = 0.6, label = 'pareto pdf' )
    
    defect_size =  create_powerlaw_random_numb_2( b, nsize )

    defect_size2 = [ int( item ) for item in defect_size if item <= nmax ]
    
    plt.xlim( 0.0, nmax )

    ax.hist( defect_size2, nmax, normed = True, histtype = 'stepfilled', alpha = 0.2 )
    ax.legend( loc = 'best', frameon = False )

    plt.show( )
    print( '# defect_size :', defect_size )

    '''


    '''
    # 2. test create_powerlaw_random_numb_1( a, nsize )

    fig, ax = plt.subplots( 1, 1 )

    nsize = 100000
    a = 1.66
    mean, var, skew, kurt = stats.powerlaw.stats( a, moments = 'mvsk' )
    
    x = np.linspace( stats.powerlaw.ppf( 0.01, a ), stats.powerlaw.ppf( 0.99, a ), 100 )
    ax.plot( x, stats.powerlaw.pdf( x, a ), 'r-', lw = 5, alpha = 0.6, label = 'powerlaw pdf' )
    
    rv = stats.powerlaw( a )
    ax.plot( x, rv.pdf( x ), 'k-', lw = 2, label = 'frozen pdf' )
    
    vals = stats.powerlaw.ppf( [ 0.001, 0.5, 0.999 ], a )
    np.allclose( [ 0.001, 0.5, 0.999 ], stats.powerlaw.cdf( vals, a ) )
    
    r = create_powerlaw_random_numb_1( a, nsize )
    
    ax.hist( r, density = True, histtype = 'stepfilled', alpha = 0.2 )
    ax.legend( loc = 'best', frameon = False )
    plt.show( )

    '''


    '''
    # 1. test elliptic integral
    for i in range( 2 ) :

        M = 1.0 * i / 2

        f = first_elliptic_integral( M )
        s = second_elliptic_integral( M )

        print( '# f(%2.6f): %2.6f' % ( M, f ) )
        print( '# s(%2.6f): %2.6f' % ( M, s ) )

    tol = 0.0001

    print( '# judge f:', judge_equal( first_elliptic_integral( 0.5 ), 1.8541, tol ) )
    print( '# judge s:', judge_equal( second_elliptic_integral( 0.5 ), 1.3506, tol ) )

    print( '# integral f(x) = x * x + sin( x ) from 2.5 to 8.4:', integral_legendre_Gauss( 2.5, 8.4, 1.0e-10 ) )

    '''
