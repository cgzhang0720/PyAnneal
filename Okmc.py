#!/usr/bin/env python3

import time
import numpy as np
import argparse
import json

from DefectSystem import DefectSystem
from Trap import TrapSystem
from Cascade import Cascade
from ConstNumber import ConstNumber
from OutputFile import OutputFile
from Auxiliary import j_must_have, j_have

def individual_cascade_annealing( jdata, cascade ) :

    trap = False
    trapSys = None

    if j_have( jdata, 'trap' ) : trap = jdata[ 'trap' ]
    if trap : trapSys = TrapSystem( jdata )

    for i in range( cascade.return_cascade_numb( ) ) :

        start_time = time.time( )

        'build a defectSys of Class DefectSystem for okmc simulation.'
        defectSys = DefectSystem( jdata, cascade.return_defectStr_list( jdata, i ) )

        outputFile = OutputFile( defectSys, jdata, i )

        'reset trap for all defectObject in model of Class DefectSystem.'
        if trap : defectSys.reset_trap_defectSystem( jdata, trapSys )

        defectSys.evolution( jdata, outputFile, trapSys )

        end_time = time.time( )

        print( '# cascade %d running time: %.3f s' % ( i, end_time - start_time ) )
        print( '#    ' )

    print( '# All cascades annealing done !' )
 
    return None

def multi_cascade_annealing( jdata, cascade ) :

    tmp_ConstNumber = ConstNumber( jdata )
    box = j_must_have( jdata, 'box' )

    trap = False
    trapSys = None

    if j_have( jdata, 'trap' ) : trap = jdata[ 'trap' ]
    if trap : trapSys = TrapSystem( jdata )

    total_cascade_numb = int( jdata[ 'fluence' ] * box[ 0 ] * box[ 1 ] )
    position_random = jdata[ 'position_random' ]

    cascade_numb = cascade.return_cascade_numb( )
    i_select = np.random.randint( 0, cascade_numb )

    'build a defectSys of Class DefectSystem for okmc simulation.'
    defectSys = DefectSystem( jdata, cascade.return_defectStr_list( jdata, i_select ) )

    outputFile = OutputFile( defectSys, jdata )

    print( '# %d cascades in database, total %d cascades annealing' %( cascade_numb, total_cascade_numb ) )

    for i in range( total_cascade_numb ) :

        start_time = time.time( )

        'reset trap for all defectObject in model of Class DefectSystem.'
        if trap : defectSys.reset_trap_defectSystem( jdata, trapSys )

        defectSys.evolution( jdata, outputFile, trapSys )

        i_select = np.random.randint( 0, cascade_numb )

        add_defectStr_list = cascade.return_defectStr_list( jdata, i_select )

        defectSys.add_defectObject( jdata, add_defectStr_list, position_random )
        defectSys.set_initial_recombine_and_related_properties( jdata, tmp_ConstNumber )

        defectSys.set_zero_c_time( )

        end_time = time.time( )

        print( '# cascade %d running time: %.3f s' % ( i, end_time - start_time ) )
        print( '#    ' )

    print( '# All cascades annealing done !' )

    return None

def _main( ) :

    parser = argparse.ArgumentParser( description = '--- cascade annealing ---' )

    parser.add_argument( 'INPUT', help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )

    seed = None
    if j_have( jdata, 'seed' ) : 
        seed = jdata[ 'seed' ]
        seed = seed % ( 2 **32 )
    np.random.seed( seed )

    individual_cascade = False
    if j_have( jdata, 'individual_cascade' ) : individual_cascade = jdata[ 'individual_cascade' ]

    multi_cascade = False
    if j_have( jdata, 'multi_cascade' ) : multi_cascade = jdata[ 'multi_cascade' ]

    if individual_cascade == False and multi_cascade == False :
        raise RuntimeError( 'individual_cascade and multi_cascade are all set False, wrong ! Please set true for one of them !' )

    cascade = Cascade( jdata )

    if individual_cascade : individual_cascade_annealing( jdata, cascade )
    if multi_cascade : multi_cascade_annealing( jdata, cascade )

if __name__ == '__main__' :
    _main( )
