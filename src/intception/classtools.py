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
# Generic classes and methods for developing and debugging classes

class classtools:
    """
    Provides some generic methods which can help with debugging 
    classes and instances of classes.
    This class can be inherited by any other class and provide useful
    methods to it.
    """
    def debug_str(self,list_out=False):
        """
        Provides a nicely formatted string (with linebreaks) containing 
        information about the class, instance of the class
        and attributes of that instance to the standard output.
        If list_out = True, then return a list of lines.
        """
        assert isinstance(list_out,bool), 'list_out must be True or False'
        superclass_list = []
        def list_superclasses(c):
            nonlocal superclass_list
            if len( c.__bases__ ) > 0:
                superclass_list.append( c.__bases__ )
                for b in c.__bases__:
                    list_superclasses( b )
            else:
                pass
        out = []
        out.append( 'Class: '+ self.__class__.__name__ )
        out.append( 'Object: '+ self.__repr__() )
        if len( list(self.__dict__.items() ) ) > 0:
            out.append( 'Object attributes:' )
            for key, value in self.__dict__.items():
                out.append( str(key)+' => '+str(value) )
        out.append( 'Superclasses:' )
        list_superclasses(self.__class__)
        for c in reversed( superclass_list ):
            out.append( str( c ) )
        if list_out == False:
            return '\n'.join( out )
        elif list_out == True:
            return out


