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
Provides the boys_function_call, boys_function_prototype and boys_function_source 
classes which are derived from callable_* classes in intception.generator.wrappers.

The boys_function_source class contains the code necessary to generate a self-contained
function for calculation of the Boys function.
"""

import math
from intception.dsl import *
from intception.generator.wrappers import function_wrapper

class boys_function_call(function_wrapper.callable_function_call):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function call specific to the Boys function.
    """
    pass

class boys_function_prototype(function_wrapper.callable_prototype):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function prototypes specific to the Boys function.
    """
    pass

class boys_function_source(function_wrapper.callable_source):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate executable source code specific to the Boys function.
    """
    def __init__(self,f_dsl,boys_function_obj,\
                 local_variable_list,local_variable_dict,\
                 oof_data_array,smallx_data_array,mediumx_data_array,gen, const_dict = {}):
        """
        Setup a new instance of boys_function_source, extends 
        function_wrapper.callable_source.__init__.

        f_dsl:                  dsl_function object corresponding to Boys function
        boys_function_obj:      boys_function object (see boys.py)
        local_variable_list:    list of local variables (dsl_variable objects)
        local_variable_dict:    dictionary of local_variables (dsl_variable objects)
                                with dsl_variable.name() as keys.
        oof_data_array:         data array object with oof values (i.e. 1/N!).
        smallx_data_array:      data array containing tabulated values for small x
                                range for use in downward recursion.
        mediumx_data_array:     data array containing tabulated values for medium x
                                range for used in upward recursion.

        Built-in Boys function implementation based on Kaito Miyamoto's implementation
        in Fortran 2003, which is itself based on the scheme outlined in T. Helgaker, 
        P. Jorgensen, J. Olsen, "Molecular Electronic-Structure Theory", p.366.
        """
        function_wrapper.callable_source.__init__(self,f_dsl,gen, const_dict = const_dict)
        self._local_variable_list = local_variable_list
        self._local_variable_dict = local_variable_dict
        self._boys_function_obj   = boys_function_obj
        self._oof_data_array      = oof_data_array
        self._smallx_data_array   = smallx_data_array
        self._mediumx_data_array  = mediumx_data_array

    def __call__(self,p):
        """
        Output source code corresponding to built-in Boys function via printer object, p.

        For the associated generator object, if 
        generator_options().parallelism().vectorized() == True, output vectorized code,
        else output serial code.

        Implementation notes:
        * Only serial implementation currently supported.
        """
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            self.serial_out(p)

    def serial_out(self,p):
        """
        Output the serial (non-vectorized) version of the built-in Boys function 
        implementation via the printer object, p.
        """
        p.funcdef( self._f_dsl, self.const_dict() )
        # Local variable declarations
        for dvar in self._local_variable_list:
            p.declaration( dvar )
        # Local variable assignment
        for dvar in self._local_variable_list:
            if dvar.expr() != None: # assignment required
                if dvar.is_cartesian() == False:
                    p.out( dvar.assign( dvar.expr() ) )
                elif dvar.is_cartesian() == True:
                    for d in [0,1,2]:
                        p.converter.direction().set(d)
                        p.out( dvar.assign( dvar.expr() ) )

        # Give function arguments nice readable names
        f      = self._boys_function_obj.boys_output_pointer()
        ifskip = self._boys_function_obj.boys_ifskiparg()
        nmin   =  self._boys_function_obj.boys_nminarg()
        n      =  self._boys_function_obj.boys_narg()
        x      =  self._boys_function_obj.boys_xarg()

        # Assign local variables from local_variable_dict readable names
        oof = self._oof_data_array.pointer()
        i = self._local_variable_dict['i']
        ii = self._local_variable_dict['ii'] # counter
        xn = self._local_variable_dict['xn']
        dx = self._local_variable_dict['dx']
        expx = self._local_variable_dict['expx']
        ooxstep = self._local_variable_dict['ooxstep']
        oox = self._local_variable_dict['oox']
        # Set nice names for small x
        xmin = self._boys_function_obj.smallx_xmin() 
        xmax = self._boys_function_obj.smallx_xmax() 
        ntaylor = self._boys_function_obj.smallx_ntaylor() 
        xstep =  self._boys_function_obj.smallx_xstep() 
        xpoints = self._boys_function_obj.smallx_xpoints() 
        bfdata = self._smallx_data_array
        p.ifblock( dsl_binop(op_le, x, xmax ) )
        # small x Boys function (Taylor series + downward recursion)
        ooxstep_value = dsl_value(name = 'o_o_xstep',\
                value = eval( '1.0 /'+xstep.value() ) )
        p.out( ooxstep.assign( ooxstep_value ) )
        if( src( xmin ) == "0" ):
            p.out( i.assign( nint( x * ooxstep ) ) ) 
            p.out( dx.assign( dsl_unop(op_dbl,i) * xstep - x ) )
        else:
            p.out( i.assign( nint( ( x - xmin )* ooxstep ) ) ) 
            p.out( dx.assign( xmin + dsl_unop(op_dbl,i) * xstep - x ) )
        p.out( xn.assign( 1.0 ) )
        p.out( f[ n * ifskip ].assign( bfdata[ i ][ n ] ) )
        # Taylor expansion
        p.forloop( ii.assign(1), dsl_binop(op_lt,ii,ntaylor), ii.assign(ii+1) )
        p.out( xn.assign( xn * dx ) )
        p.out( f[ n * ifskip ].assign( f[ n * ifskip ] + bfdata[ i ][ n + ii ] * xn * oof[ ii ] ) )
        p.endforloop()
        p.out( expx.assign( exp( -x ) ) ) 
        # Downward recursion
        p.forloop( ii.assign(n-1), dsl_binop(op_ge,ii,nmin), ii.assign(ii-1) )
        p.out( f[ ii * ifskip ].assign( (2.0 * x * f[ ( ii + 1 )* ifskip ] + expx) / ( 2.0 * dsl_unop(op_dbl,ii) + 1.0 ) ) )
        p.endforloop()
        p.returnstmt()
        # Set nice names for medium x
        xmin = self._boys_function_obj.mediumx_xmin()
        xmax = self._boys_function_obj.mediumx_xmax()
        ntaylor = self._boys_function_obj.mediumx_ntaylor()
        xstep = self._boys_function_obj.mediumx_xstep()
        xpoints = self._boys_function_obj.mediumx_xpoints()
        bfdata = self._mediumx_data_array
        p.elifblock( dsl_binop(op_and,dsl_binop(op_gt,x,xmin),dsl_binop(op_lt,x,xmax) ) )
        # medium x Boys function (Taylor series + upward recursion)
        ooxstep_value = dsl_value(name = 'o_o_xstep',\
                value = eval( '1.0 /'+xstep.value() ) )
        p.out( ooxstep.assign( ooxstep_value ) )
        p.out( i.assign( nint( x * ooxstep ) ) )
        p.out( dx.assign( dsl_unop(op_dbl,i) * xstep - x ) )
        p.out( i.assign( i - int( eval( src(xmin)+' / '+xstep.value() ) ) ) )
        p.out( oox.assign( 1.0 / x ) )
        p.out( xn.assign( 1.0 ) )
        p.out( f[ 0 ].assign( bfdata[ i ][ 0 ] ) )
         # Taylor expansion
        p.forloop( ii.assign(1), dsl_binop(op_lt,ii,ntaylor), ii.assign(ii+1) )
        p.out( xn.assign( xn * dx ) )
        p.out( f[ 0 ].assign( f[ 0 ] + bfdata[ i ][ ii ] * xn * oof[ ii ] ) )
        p.endforloop()
        p.out( expx.assign( exp( -x ) ) ) 
        # Upward recursion
        p.forloop( ii.assign(0), dsl_binop(op_lt,ii,n), ii.assign(ii+1) )
        p.out( f[ ii*ifskip + ifskip ].assign( (dsl_unop(op_dbl,ii)+0.5 ) * oox * f[ ii*ifskip ] - expx * 0.5 * oox ) )
        p.endforloop()
        p.returnstmt()
        p.elseblock()
        # large x Boys function (approximate expression)
        p.out( oox.assign( 1.0 / x ) )
        # Approximate formula for large x
        sqrt = dsl_function('sqrt',[ math.pi * oox ])
        p.out( f[ 0 ].assign( 0.5 * sqrt ) )
        # Upward recursion (approximated exp(-x) = 0)
        p.forloop( ii.assign(0), dsl_binop(op_lt,ii,n), ii.assign(ii+1) )
        p.out( f[ ii*ifskip+ifskip ].assign( (dsl_unop(op_dbl,ii)+0.5 ) * oox * f[ ii*ifskip ] ) )
        p.endforloop()
        p.returnstmt()
        p.endifblock()
        p.endfuncdef()

    def vector_out(self,p):
        """
        Output the serial (non-vectorized) version of the built-in Boys function 
        implementation via the printer object, p.
        """
        #assert isinstance(p,printer), "p must be a printer object"
        f = self._f_dsl
        p.funcdef( f )
        p.out("NOTHING TO SEE HERE",endl='')
        p.endfuncdef( f )

