#!/usr/bin/env python3

'''
Parameters used to detect Classes and functions.
At most 5 global parameters.
'''

global data_0_
global data_1_
global data_2_
global data_3_
global data_4_

def set_value( value, serial_numb = 0 ) :

    if serial_numb == 0 :
        global data_0_
        data_0_ = value
    elif serial_numb == 1 :
        global data_1_
        data_1_ = value
    elif serial_numb == 2 :
        global data_2_
        data_2_ = value
    elif serial_numb == 3 :
        global data_3_
        data_3_ = value
    elif serial_numb == 4 :
        global data_4_
        data_4_ = value
    else :
        raise RuntimeError( 'Over 5 global parameters, wrong!' )

def get_value( serial_numb = 0 ) :

    if serial_numb == 0 :
        global data_0_
        return data_0_
    elif serial_numb == 1 :
        global data_1_
        return data_1_
    elif serial_numb == 2 :
        global data_2_
        return data_2_
    elif serial_numb == 3 :
        global data_3_
        return data_3_
    elif serial_numb == 4 :
        global data_4_
        return data_4_
    else :
        raise RuntimeError( 'Over 5 global parameters, wrong!' )

def return_parameters( serial_numb = 0 ) :
    '''
    Parameters of this moduler:
    '''

    if serial_numb == 0 :
        global data_0_
        return [ 'data_0_', data_0_ ]
    elif serial_numb == 1 :
        global data_1_
        return [ 'data_1_', data_1_ ]
    elif serial_numb == 2 :
        global data_2_
        return [ 'data_2_', data_2_ ]
    elif serial_numb == 3 :
        global data_3_
        return [ 'data_3_', data_3_ ]
    elif serial_numb == 4 :
        global data_4_
        return [ 'data_4_', data_4_ ]
    else :
        raise RuntimeError( 'Over 5 global parameters, wrong!' )

def print_parameters( serial_numb = 0 ) :

    item = return_parameters( serial_numb )

    print( '# %s : %s' % ( item[ 0 ], item[ 1 ] ) )

if __name__ == '__main__' :

    set_value( 3, 4 )
    print_parameters( 4 )
