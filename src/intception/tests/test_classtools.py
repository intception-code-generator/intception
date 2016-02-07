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
from intception import classtools

class TestClassTools(unittest.TestCase):
    """Tests for classtools.py classes and their methods."""
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_classtools(self):
        c = classtools.classtools()
        c.hw = 'hello, world!'
        # Test debug_print()
        self.assertTrue( isinstance( c.debug_str(), str ) )
        self.assertTrue( isinstance( c.debug_str(list_out = True), list ) )
        self.assertIn( 'hw => '+c.hw, c.debug_str( list_out = True ) )
        class test_class1(classtools.classtools):
            def __init__(self):
                self.a = 'apple'
                self.b = 100.0
                self.c = 42
        class test_class2(test_class1):
            def set_attributes(self):
                self.d = 'donut'
                self.e = 2.71828
                self.f = 101
        t1 = test_class1()
        self.assertIn( 'Class: '+str(t1.__class__.__name__), t1.debug_str( list_out = True ) )
        self.assertIn( 'Object: '+str(t1.__repr__()), t1.debug_str( list_out = True ) )
        t2 = test_class2()
        t2.set_attributes()
        for t in [ t1, t2 ]:
            # Check that class and object names are correctly output
            self.assertIn( 'Class: '+str(t.__class__.__name__), t.debug_str( list_out = True ) )
            self.assertIn( 'Object: '+str(t.__repr__()), t.debug_str( list_out = True ) )
            # Check all attributes are correctly output by debug_print
            for key, value in t.__dict__.items():
                self.assertIn( str(key) +' => '+str(value), t.debug_str( list_out = True ) )
        superclass_list = [ c.__class__, t1.__class__ ]
        # Check that list of superclasses is as expected
        for c in superclass_list:
            self.assertIn( str( (c,) ), t2.debug_str( list_out = True ) )

