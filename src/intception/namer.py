#!/usr/bin/env python3
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
import sys
import intception.classtools

class namer:
    """Class which uniquely names DSL objects."""
    def __init__(self,separator='_',default_suffix=''):
        """A short-naming scheme is used where variables are simply numbered, e.g. 
        for integers
            i_1, i_2, i_3 etc.
        for double precision
            dp_1, dp_2, dp_3, etc.
        A dictionary of integer counters is updated for each named
        variable, with the keys set to the value of 
        str( variable.type() ) + str( variable.__class__ )
        which uniquely determines the type of the variable in 
        C source code and also the type of the object in the DSL
        (this distinguishes, for example, double precision scalars
        from double precision positions (3 index arrays)).
        Each DSL class which uses the namer class to name variables
        should supply a prefix() method which returns a short string
        dependent on the value of variable.type() and the object type
        (e.g. for dsl_vars with type() == 'int', prefix() == 'i')"""
        self._separator = separator
        self._default_suffix = default_suffix
        self._counter_dict = {}

    def __call__(self,obj):
        assert isinstance( obj.prefix(), str )
        assert isinstance( obj.type(), str )
        key = obj.type()+str( obj.__class__ )
        try:
            self._counter_dict[key] += 1
        except KeyError:
            # No dictionary entry yet
            self._counter_dict[key] = 1
        out = []
        out.append( obj.prefix() )
        out.append( self._separator )
        out.append( str( self._counter_dict[key] ) )
        out.append( self._default_suffix )
        return ''.join(out)

    def reset_counters(self):
        """Resets all counters for variable naming."""
        self._counter_dict = {}

