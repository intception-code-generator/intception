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
from intception.supporting_functions import *

class TestSupporting_functions(unittest.TestCase):
    """Tests for supporting_functions.py classes and their methods."""
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_isinstance_list(self):
        # Test isinstance_list function
        # Create mock classes 
        class mock:
            pass
        class mock_child(mock):
            pass
        m = mock()
        mc = mock_child()
        self.assertTrue( isinstance( m, ( mock, mock_child, int, str, list ) ) ) 
        self.assertFalse( isinstance( mc, ( int, str, list ) ) )
        self.assertTrue( isinstance_list( m, [ mock, mock_child, int, str, list ] ) ) 
        self.assertFalse( isinstance_list( mc, [ int, str, list ] ) )
        # type_list must be of list type
        self.assertRaises(AssertionError, isinstance_list, m, ( mock, mock_child ) )
        self.assertRaises(AssertionError, isinstance_list, m, mock )

