#!/usr/bin/env python3

import math
from Auxiliary import j_must_have

class ConstNumber( object ) :

    def __init__( self, jdata ) :
        '''
        --- property    --  unit    -- 
        1.  Length      --   nm     --   self.alatt, self.first, self.second
        2.  Volume      --   nm^3   --   self.volume
        2.  Time        --   s      --
        3.  Temperature --   K      --
        4.  Energy      --   eV     --
        5.  Pressure    --   Gpa    --  

        1 eV = 1.60217663410E-19 J  ->  1 J = 1.0 / 1.60217663410 * 1.E+19 eV
        shear modulus:
        1.0 Gpa = 1.E9 Pa = 1.E9 N / m^2 = 1.e9 J/m^3 = 1.E9 * 1.0 / 1.60217663410 * 1.E+19 eV/(10^9 nm)^3 = 1.0 / 1.60217663410 * 10 eV/nm^3
        Boltzmann constant, unit: eV/ K
        kb = 1.380649 * 10 ** (-23) J / K = 1.380649 * 10 ** (-23) * 1.0 / 1.60217663410 * 10**( 19 ) eV/K = 1.380649 / 1.60217663410 * 1.E-4 eV/K
        '''

        'shear modulus, unit: eV/nm^3.'
        self.shear_modulus = 161.0 / 1.60217663410 * 10

        'poisson ratio'
        self.poisson_ratio = 0.28

        'lattice constant'
        self.alatt = j_must_have( jdata, 'alatt' )

        lattice = j_must_have( jdata, 'lattice' )

        if lattice == 'bcc' :

            'first nearest neighbour'
            self.first = math.sqrt( 3.0 ) * 0.5 * self.alatt
 
            'second nearest neighbour'
            self.second = self.alatt
 
            'atomic volume'
            self.volume = self.alatt * self.alatt * self.alatt * 0.5

        elif lattice == 'fcc' :
            print( '# should be change for lattice fcc !' )
            pass
        elif lattice == 'hcp' :
            print( '# should be change for lattice hcp !' )
            pass

        self.kb = 1.380649 / 1.60217663410 * 1.E-4

        'prefactor: Hz'
        self.prefactor = 10.0 * 10 ** 12

        self.u_disk_disk_prefactor = self.shear_modulus / ( 4.0 * math.pi * ( 1.0 - self.poisson_ratio ) )

        self.u_sphere_disk_prefactor = 1.0 / ( 3.0 * math.pi ) * ( 1.0 + self.poisson_ratio ) / ( 1.0 - self.poisson_ratio ) * self.shear_modulus

if __name__ == '__main__' :

    parser = argparse.ArgumentParser( description = '--- Class ConstNumber detecting ---' )
    parser.add_argument( 'INPUT',
                         help = 'the input json database' )

    args = parser.parse_args( )

    fp = open( args.INPUT, 'r' )
    jdata = json.load( fp )
