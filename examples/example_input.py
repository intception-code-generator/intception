### LICENSE INFORMATION
### This file is part of Intception, a molecular integral code generator
### Copyright (C) 2016 James C. Womack
### 
### Intception is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
### 
### Intception is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
### 
### You should have received a copy of the GNU General Public License
### along with Intception.  If not, see <http://www.gnu.org/licenses/>.
### 
### See the file named LICENSE for full details.
"""
EXAMPLE INPUT SCRIPT FOR INTCEPTION
James C. Womack, 01/2016

In order for the script to correctly import the intception modules, you 
must either run this script in the src/ directory containing the intception 
source code, or add this directory to your shell's PYTHONPATH environment 
variable, e.g.
  export PYTHONPATH=$PYTHONPATH:${REPODIR}/src
where ${REPODIR} is the directory into which the Intception git repository
has been cloned. 

To execute the script run it via the Python 3 interpreter:
  python3 example_input.py
Generated integral code will be output to the directory set in the 
'output_directory' argument of the generator_options object created at the
end of the script (by default, this is in a directory called "output" in the
current working directory).

The current version of this script generates C source code for evaluating 
two-index overlap and nuclear attraction integrals over primitive Cartesian
Gaussians, illustrating the key features and syntax used to specify integrals
and generate source code.
"""

# Import modules
from math import pi
from intception import *

# Define variables
ga = dsl_cartesian_gaussian('a', constant = True)
xa = ga.exp()
A  = ga.cen()
gb = dsl_cartesian_gaussian('b', constant = True)
xb = gb.exp()
B  = gb.cen()
gc = dsl_cartesian_gaussian('c', constant = True)
xc = gc.exp()
C  = gc.cen()

# Constants
pi3 = dsl_value(name = 'pi3', value = pi**3, vartype='double' )

# Derived variables
xp = dsl_scalar('xp', xa + xb )
o_o_xp = dsl_scalar('o_o_xp', 1.0/ xp )
xaxb_o_xp = dsl_scalar('xaxb_o_xp', xa * xb * o_o_xp )
o_o_xp3 = dsl_scalar('o_o_xp3', o_o_xp * o_o_xp * o_o_xp )
o_o_2xp = dsl_scalar('o_o_2xp', 0.5 * o_o_xp )
cAcB = dsl_position('cAcB', A - B )
cAcC = dsl_position('cAcC', A - C )
cBcC = dsl_position('cBcC', B - C )
PA = dsl_position('PA', xb * o_o_xp * -cAcB)
PB = dsl_position('PB', xa * o_o_xp * cAcB )
PC = dsl_position( 'PC', xa * o_o_xp * cAcC \
                     + xb * o_o_xp * cBcC )
RAB2 = dsl_scalar('RAB2', cAcB[0]*cAcB[0] + cAcB[1]*cAcB[1] + cAcB[2]*cAcB[2])
RPC2 = dsl_scalar('RPC2', PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2] )

### 2-idx overlap (overlap_2idx) ###
overlap_2idx = dsl_integral(ga,gb,name='overlap_2idx')
overlap_2idx.set_base( sqrt( pi3 * o_o_xp3 ) * exp( -xaxb_o_xp*RAB2) )
overlap_2idx.add_vrr('vrr1', ga, \
    PA * overlap_2idx.int(ga-1,gb) + o_o_2xp * (ga - 1) * overlap_2idx.int(ga-2,gb) \
                         + o_o_2xp * gb * overlap_2idx.int(ga-1,gb-1) \
  )
overlap_2idx.add_vrr('vrr2', gb, \
    PB * overlap_2idx.int(ga,gb-1) + o_o_2xp * (gb - 1) * overlap_2idx.int(ga,gb-2) \
                         + o_o_2xp * ga * overlap_2idx.int(ga-1,gb-1) \
  )

### 2-idx  nuclear attraction (nuclear_2idx) ###
m  = dsl_integer_index(name = 'm')
U = xp * RPC2
boys = boys_f( m, U )
nuclear_2idx = dsl_integral(ga,gb,m,C,name='nuclear_2idx')
nuclear_2idx.set_base( 2.0 * ( pi * o_o_xp ) * exp( -xaxb_o_xp*RAB2) * boys )
nuclear_2idx.add_vrr('vrr1', gb, \
    PB * nuclear_2idx.int(ga,gb-1,m) - PC * nuclear_2idx.int(ga,gb-1,m+1) \
  + o_o_2xp * ( gb - 1 ) \
    * ( nuclear_2idx.int(ga,gb-2,m) - nuclear_2idx.int(ga,gb-2,m+1) ) \
  + o_o_2xp * ga \
    * ( nuclear_2idx.int(ga-1,gb-1,m) - nuclear_2idx.int(ga-1,gb-1,m+1) ) 
)
nuclear_2idx.add_vrr('vrr2', ga, \
    PA * nuclear_2idx.int(ga-1,gb,m) - PC * nuclear_2idx.int(ga-1,gb,m+1) \
  + o_o_2xp * gb \
    * ( nuclear_2idx.int(ga-1,gb-1,m) - nuclear_2idx.int(ga-1,gb-1,m+1) ) \
  + o_o_2xp * ( ga - 1 ) \
    * ( nuclear_2idx.int(ga-2,gb,m) - nuclear_2idx.int(ga-2,gb,m+1) ) 
)


### Set generator options and create generator object ###
# Primitive, Cartesian
options = generator_options( output_directory = "./output", contracted = False )
gen = generator( overlap_2idx, nuclear_2idx, opt = options )
gen.out()
