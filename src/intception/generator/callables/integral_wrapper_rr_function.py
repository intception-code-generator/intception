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
Provides the rr_call, rr_prototype and vrr_source 
classes which are derived from callable_* classes in intception.generator.wrappers.

The rr_source class contains the code necessary to generate a self-contained
VRR function containing an unrolled loop.
"""

from intception.dsl import *
from intception.generator.wrappers import function_wrapper

class rr_call(function_wrapper.callable_function_call):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function calls specific to a RR.
    """
    pass

class rr_prototype(function_wrapper.callable_prototype):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function prototypes specific to a RR.
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
class vrr_source(function_wrapper.callable_source):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate executable source code specific to a VRR.
    """
    def __init__(self,f_dsl,expr,windex,windex_src_max,wintegral,gen, const_dict = {}):
        function_wrapper.callable_source.__init__(self,f_dsl, gen, const_dict )
        self._expr = expr
        self._windex = windex
        self._windex_src_max = windex_src_max
        self._wintegral = wintegral

    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            if self.generator().generator_options().contracted() == True:
                # Suppress expr_to_contract_convert in p.converter object
                p.converter.set_contract_convert( False )
            self.serial_out(p)
            if self.generator().generator_options().contracted() == True:
                # End suppression of expr_to_contract_convert in p.converter object
                p.converter.set_contract_convert( True )

    def decide(self,a):
        """Return the index of the first element in an array of integers that is
        greater than 0. A dumb selector that does not optimize the route 
        through the RR.
        
        Having this as a method of rr_source allows this to be more flexible as the
        global method from the generator instance can be replaced by something else."""
        return self.generator().decide(a)


    def serial_out(self,p):
        # [ Current implementation restrictions ]
        assert isinstance(self._windex.index(),dsl_cartesian_gaussian),\
                "only dsl_cartesian_gaussian index objects supported currently"

        f = self._f_dsl
        dintegral = self._wintegral.integral()
        expr  = dintegral.int().assign( self._expr )
        index = self._windex.index()
        imax_var = self._windex.index_loop_max()
        imax_src = self._windex_src_max
        array = self._wintegral.work_array()
        array.set_current_index( index )
        direction = p.converter.direction()
        p.funcdef( f, self.const_dict() )
        # Set angmom of index being incremented to zero
        assert index.is_cartesian() == True,\
                "only Cartesian unrolled index supported currently"
        index.set_angmom([0,0,0])
        while index.angmom()[2] != imax_src:
            before_increment_total_angmom = sum( index.angmom() )
            index.increment()
            if sum( index.angmom() ) > before_increment_total_angmom:
                # total angular momentum incremented -- check if need to exit
                p.ifblock( dsl_binop(op_eq,imax_var,before_increment_total_angmom) )
                p.returnstmt()
                p.endifblock()
            direction.set( self.decide( index.angmom() ) )
            # Output RR expression
            p.out( expr )
        p.endfuncdef()

    def vector_out(self,p):
        f = self._f_dsl
        p.funcdef( f )
        p.out("NOTHING TO SEE HERE",endl='')
        p.endfuncdef()

