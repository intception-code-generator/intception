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
Provides parent classes with methods to be called by integral_wrapper setup_* routines, that set
the integral class specific variable lists. 

Each callable main_source_algorithm should provide an algorithm-specific set of classes which 
ensure that only the required variables are used in generated code.
"""

from intception import classtools

class algo_variable_setup(classtools.classtools):
    def __init__(self,wintegral,gen):
        self.wintegral = wintegral
        self.gen       = gen

    def __call__(self):
        """
        This is a placeholder, intended to be overidden in derived classes.
        """
        pass
