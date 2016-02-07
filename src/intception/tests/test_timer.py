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
import numbers
from intception.timer import timer

class TestTimer(unittest.TestCase):
    """Tests for generator/arrays.py classes and their methods."""
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_timer(self):
        with timer() as t:
            j = 1
            for i in range(1000):
                j *= i
        self.assertTrue( isinstance( t, timer ) )
        self.assertTrue( isinstance( t.start_time, numbers.Real ) )
        self.assertTrue( isinstance( t.end_time, numbers.Real ) )
        self.assertEqual( t.time_taken, t.end_time - t.start_time ) 
                
                





