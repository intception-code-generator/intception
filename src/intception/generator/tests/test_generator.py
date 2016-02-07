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
import math
from intception.boys import boys_f
from intception.dsl_extensions import dsl_integral
from intception.generator.generator import *

class TestGenerator(unittest.TestCase):
    """Tests for generator.py classes and their methods."""
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_generator_options(self):
        opt = generator_options()
        self.assertTrue( opt.contracted() )
        self.assertTrue( opt.spherical_transformed() )
        opt = generator_options( contracted = False, spherical_transformed = False )
        self.assertFalse( opt.contracted() )
        self.assertFalse( opt.spherical_transformed() )

    def test_generator(self):
        pass


