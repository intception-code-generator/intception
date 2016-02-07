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
Define set of modules that are imported when 
   import intception
is used in an input script

All identifiers are imported into the global namespace for ease of scripting.
To minimize the possibility of namespace collisions, the identifiers must 
be explicitly included in following from ... import statements.
"""
from .boys import boys_f, boys_function
from .dsl_extensions import dsl_rr, dsl_integral
from .dsl import dsl_base, dsl_binop, dsl_cartesian_gaussian, \
    dsl_function, dsl_integer_index, \
    dsl_position, dsl_scalar, dsl_unop, dsl_value, dsl_zero, \
    sqrt, exp, norm, nint # forms of these operators specific to DSL
from .generator.generator import generator, generator_options
