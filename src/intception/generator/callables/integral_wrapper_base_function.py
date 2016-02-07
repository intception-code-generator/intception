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
Provides the base_call, base_prototype and vbase_source 
classes which are derived from callable_* classes in intception.generator.wrappers.

The base_source class contains the code necessary to generate a self-contained
function for evaluating the zero-angular momentum case of an integral class.
"""

from intception.dsl import *
from intception.dsl_extensions import *
from intception.generator.wrappers import function_wrapper, boys_wrapper
from intception.boys import boys_function, boys_f

class base_call(function_wrapper.callable_function_call):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function call specific to the base function.
    """
    pass

class base_prototype(function_wrapper.callable_prototype):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function prototypes specific to the base function..
    """
    def __call__(self,p):
        """
        Extends callable_prototype.__call__ by adding suppression of
        expr_to_contract_convert for p.converter.
        """
        if self.generator().generator_options().contracted() == True:
            p.converter.set_contract_convert( False )
        function_wrapper.callable_prototype.__call__(self,p)
        if self.generator().generator_options().contracted() == True:
            p.converter.set_contract_convert( True )

class base_source(function_wrapper.callable_source):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate executable source code specific to the base function.

    Implementation notes:
    * The base class may only depend on a single auxiliary index.
    """
    def __init__(self,f_dsl,expr,wintegral,gen, const_dict = {},\
                 special_function_tuple_list=[],):
        function_wrapper.callable_source.__init__(self,f_dsl,gen, const_dict = const_dict)
        self._expr = expr
        self._wintegral = wintegral
        self._special_function_tuple_list = special_function_tuple_list
    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            if self.generator().generator_options().contracted() == True:
                # Suppress expr_to_contract_convert for code output (inside the base function
                # we treat primitive exponents as scalars, as the scalar value is passed in)
                p.converter.set_contract_convert( False )
            self.serial_out(p)
            if self.generator().generator_options().contracted() == True:
                # End suppression of expr_to_contract convert in p.converter
                p.converter.set_contract_convert( False )

    def serial_out(self,p):
        # [ Current implementation restrictions ]
        for t in self._special_function_tuple_list:
            spf = t[0]
            assert spf.m().name() in list( self._wintegral.windex_aux_dict().keys() ),\
                    "special function auxiliary index must be in integral_wrapper object "+\
                    "windex_aux_dict."
            windex = self._wintegral.windex_aux_dict()[ spf.m().name() ]
            assert isinstance(windex.index(),dsl_integer_index),\
                    "only dsl_integer_index objects supported currently"
            assert windex.index().is_cartesian() == False,\
                    "only non-Cartesian auxiliary index supported currently"
        f = self._f_dsl
        wintegral = self._wintegral
        work_array = wintegral.work_array()
        dintegral = wintegral.integral()
        expr = self._expr
        local_list = f.local_vars()[:]
        for windex in wintegral.windex_index_list():
            if windex.index().is_cartesian() == True:
                windex.index().set_angmom([0,0,0])
            else:
                windex.index().set_value(0)


        # Output top line of function definition
        p.funcdef( f, self.const_dict() )

        # Declare local variables
        for dvar in local_list:
            p.declaration( dvar )

        # Insert code for evaluating special functions (e.g. Boys function)
        work_array_index  = work_array.special_array_index_dict()['base']
        if len( self._special_function_tuple_list ) > 0:
            assert len( self._special_function_tuple_list ) <= 1,\
                    "Only one instance of special function currently supported"
            for spf, spf_wrapper in self._special_function_tuple_list: 
                assert isinstance( spf, boys_f ), \
                        "Only Boys function currently supported"
                assert isinstance( spf_wrapper, boys_wrapper ), \
                        "Only Boys function currently supported"
                aux_windex = wintegral.windex_aux_dict()[ spf.m().name() ]
                aux_min = aux_windex.length_dict()['loop min']
                aux_max = aux_windex.loop_length()
                aux_array_counter = aux_windex.array_index()
                # aux_array_skip always 1 in 2-layer reduced memory algorithm
                aux_array_skip = dsl_value( name='aux_array_skip', value='1' )
                x_arg = spf.x()
                # Insert call to special function which populates work array
                # with appropriate values
                spf_wrapper.special_function_call_out( p, \
                        array_pointer = dsl_unop(\
                            op_ref, work_array.pointer()[ work_array_index ] ),\
                        array_skip = aux_array_skip, 
                        aux_min = aux_min, \
                        aux_max = aux_max, \
                        x_arg = x_arg )
                # Modify base expression to replace call to boys_f with 
                # reference to the appropriate work_array elements.
                modified_expr = \
                    work_array[ work_array_index + aux_array_counter ].assign(\
                            expr_find_and_replace( expr, spf,\
                        work_array[ work_array_index + aux_array_counter ] ) ) 

                # Loop over auxiliary index (only 1 auxiliary index supported
                # currently)
                p.forloop( aux_array_counter.assign(0),\
                        dsl_binop(op_le,aux_array_counter,aux_max),\
                        aux_array_counter.assign(aux_array_counter+1) )
                p.out( modified_expr ) 
                p.endforloop()

        else: # no special function => no loop over auxiliary index
            modified_expr = work_array[ work_array_index ].assign( expr )
            p.out( modified_expr ) 

        # Close function definition block
        p.endfuncdef()


    def vector_out(self,p):
        f = self._f_dsl
        p.funcdef( f )
        p.out("NOTHING TO SEE HERE",endl='')
        p.endfuncdef()

