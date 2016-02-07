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
Provides the main_call, main_prototype and main_source classes which are derived from
callable_* classes in intception.generator.wrappers.

Of particular importance is the main_source class, which contains code for outputting the
main function which is called when evaluating an integral class.
"""

from intception.dsl import *
from intception.generator.wrappers import function_wrapper
                                          
# Callable classes containing code necessary to generate prototype and source
# specific to main function
class main_call(function_wrapper.callable_function_call):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate a function call for the main integral function.
    """
    pass

class main_prototype(function_wrapper.callable_prototype):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate a function prototype for the main integral function.
    """
    pass

class main_source(function_wrapper.callable_source):
    def __init__(self,f_dsl,wintegral, gen, const_dict = {}):
        """
        Instances are callable objects carrying with them all the information
        necessary to generate executable source code for the main integral function.
        """
        function_wrapper.callable_source.__init__(self,f_dsl, gen, const_dict = const_dict )
        self._wintegral = wintegral

    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            self.serial_out(p)

    def serial_out(self,p):
        raise Exception('serial_out should only be called for a derived main_source class')
        
    def vector_out(self,p):
        raise Exception('serial_out should only be called for a derived main_source class')

