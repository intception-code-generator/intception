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
Provides the work_array_size_call, work_array_size_prototype and work_array_size_source 
classes which are derived from callable_* classes in intception.generator.wrappers.

The work_array_size_source class contains the code necessary to generate the work
array sizing function attached to each integral class.
"""

from intception.dsl import *
from intception.dsl_extensions import *
from intception.generator.wrappers import function_wrapper

class work_array_size_call(function_wrapper.callable_function_call):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate a function call for a function determining the size 
    of the work array required for a specific integral class.
    """
    pass

class work_array_size_prototype(function_wrapper.callable_prototype):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate a function prototype for a function determining the size
    of the work array required for a specific integral class.
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

class work_array_size_source(function_wrapper.callable_source):
    def __init__(self,f_dsl,wintegral,local_variable_list,gen, const_dict = {}):
        """
        Instances are callable objects carrying with them all the information
        necessary to generate executable source for a function determining the size
        of the work array required for a specific integral class.
        """
        function_wrapper.callable_source.__init__(self,f_dsl, gen, const_dict = const_dict )
        self._wintegral = wintegral
        self._local_variable_list = local_variable_list

    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            if self.generator().generator_options().contracted() == True:
                # Suppress expr_to_contract_convert in p.converter for output of function 
                # definition
                p.converter.set_contract_convert( False )
            self.serial_out(p)
            if self.generator().generator_options().contracted() == True:
                # End suppression of expr_to_contract_convert in p.converter
                p.converter.set_contract_convert( True )

    def serial_out(self,p):
        # Should be overloaded on a per-algorithm basis, see callables.main_source_algorithms
        pass
        

    def vector_out(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        f = self._f_dsl
        local_variable_list = self._local_variable_list
        p.funcdef( f )
        p.out( "PLACEHOLDER" )
        p.endfuncdef()

