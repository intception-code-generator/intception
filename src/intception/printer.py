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
from intception import classtools

# src( ) function, which allows the src() method of dsl objects to be called 
# in the same style as the __str__ method, e.g. src( arg ), rather than object.src()
def src(arg):
    """
    Calls the src() method of arg and returns the result (should be a string
    of valid C source code. If a src() method is unavailable (e.g. for non-dsl objects
    like integers or strings), then call the __str__() method.
    """
    try:
        return arg.src()
    except AttributeError:
        return str( arg )

class printer(classtools.classtools):
    """
    Basic output of source code from objects in the DSL language (dsl.py).
    The class provdes methods for outputting source code using the source code
    representations of DSL objects provided the dsl_object.src() methods provided
    with each DSL class that can be represented in source code.
    """
    def __init__(self,endl=';',output_file=sys.stdout,tab=2):
        """
        endl:         line end character, added to the end of each line from printer.out().
        output_file:  file object or file-like object to be used for text output
                      (must be open for writing).
        tab:          number of blank spaces for each indent/outdent step.
        """
        self._end = endl
        self._file = output_file
        self._tab  = tab
        self._indent = 0
        self._max_indent = 80

    def out(self,arg,endl=None,output_file=None):
        """
        printer.out() uses the built-in print function to output text to a file or 
        file-like object, surrounding the text representation of arg with the indentation
        and line end characters set for the printer object.

        arg:          object which is to be converted to a text (source code output)
                      representation; the src() function is called with arg as argument
        endl:         allows endl set in __init__() to be overidden for a single function call
        output_file:  allows output_file set in __init__() to be overidden for a single 
                      function call
        """
        if endl == None:
            endl = self._end
        if output_file == None:
            output_file = self._file
        print(self._indent*' '+src(arg), end=endl+'\n', file=output_file)

    def indent(self):
        """Increments the indentation level by self._tab spaces."""
        assert self._indent + self._tab <= self._max_indent,\
                'indent greater than _max_indent ('+str(self._max_indent)+')'
        self._indent += self._tab

    def outdent(self):
        """Decrements the indentation level by self._tab spaces."""
        assert self._indent - self._tab >= 0, 'negative indent not allowed'
        self._indent -= self._tab

    def blankline(self,output_file=None):
        """
        prints a blankline without a line ending character or indentation

        output_file:  allows output_file set in __init_() to be overidden
        """
        if output_file == None:
            output_file = self._file
        print( '', end='\n', file=output_file )
    
    def indent_value(self):
        """Returns the current indentation level set for the printer object."""
        return self._indent
        
