#!/usr/bin/env python3
'This is an Parameter_test moduler.'

import MyGlobalParameter

def fun( ) :
    global data_

    MyGlobalParameter.set_value( 1.0 )
    print( 'data_ in fun', MyGlobalParameter.get_value( ) )
    MyGlobalParameter.print_parameters()

if __name__ == '__main__' :

    fun( ) 
    MyGlobalParameter.set_value( 3.0 )
    MyGlobalParameter.print_parameters()
