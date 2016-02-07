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
import unittest
import sys
from io import StringIO # mock file descriptor
from intception.printer import *

# Set up mock class for printing
class dsl_mock:
    """Has a src() method which should be preferred over
    __str__() when printer.out or src are called with
    a dsl_mock object as argument."""
    def __init__(self,string):
        self._string = string

    def __str__(self):
        return '__'+self._string+'__'

    def src(self):
        return self._string

class other_mock:
    """
    No src() method, so printer and src should fall back to
    __str__() when a other_mock object is passed as an 
    argument.
    """

    def __init__(self,string):
        self._string = string

    def __str__(self):
        return '__'+self._string+'__'

class TestPrinter(unittest.TestCase):
    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    def test_printer(self):
        # Test the printer class
        m1 = dsl_mock('m1')
        m2 = dsl_mock('m2')
        m3 = dsl_mock('m3')
        o1 = other_mock('o1')
        f  = StringIO()
        p = printer(endl=';',output_file=f,tab=2)
        p.out( m1 )
        self.assertEqual( f.getvalue(), 'm1;\n' )
        p.indent()
        p.out( 199 )
        self.assertEqual( f.getvalue(), 'm1;\n  199;\n' )
        p.indent()
        p.out( m2 )
        self.assertEqual( f.getvalue(), 'm1;\n  199;\n    m2;\n')
        p.outdent()
        p.outdent()
        p.out( m3 )
        self.assertEqual( f.getvalue(), 'm1;\n  199;\n    m2;\nm3;\n')
        p.out( o1 )
        self.assertEqual( f.getvalue(), 'm1;\n  199;\n    m2;\nm3;\n__o1__;\n')

        # Negative indent not allowed
        self.assertRaises(AssertionError, p.outdent)

    def test_src(self):
        # Test the src() subroutine
        i1  = 1
        dp1 = 2.945
        o1 = other_mock('o1')
        list_other = [ i1, dp1, o1 ]
        for item in list_other:
            # Default to __str__ when src() method not available
            self.assertEqual( src( item ), str( item ) )

        m1 = dsl_mock('m1')
        m2 = dsl_mock('m2')
        m3 = dsl_mock('m3')
        list_dsl = [ m1, m2, m3 ]
        list_srcstr = [ 'm1', 'm2', 'm3' ]
        for item, srcstr in zip( list_dsl, list_srcstr):
            # __str__() method should not be called for objects with src() method
            self.assertNotEqual( src( item ), str( item ) )
            self.assertEqual( src( item ), srcstr )
        
