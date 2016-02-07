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
from intception.namer import *

class TestNamer(unittest.TestCase):
    """Tests the namer class which creates objects
    for automatic naming of DSL objects."""
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_namer(self):
        """Test the namer class in isolation using mock objects."""
        class dsl_mock:
            """Mock object that can be processed by namer."""
            def __init__(self,prefix,vartype):
                self._prefix  = prefix
                self._vartype = vartype
            def prefix(self):
                return self._prefix
            def type(self):
                return self._vartype
        n = namer()
        i1 = dsl_mock(prefix = 'i', vartype = 'int')
        i2 = dsl_mock(prefix = 'i', vartype = 'int')
        i3 = dsl_mock(prefix = 'i', vartype = 'int')
        dp1 = dsl_mock(prefix = 'dp', vartype = 'double')
        dp2 = dsl_mock(prefix = 'dp', vartype = 'double')
        dp3 = dsl_mock(prefix = 'dp', vartype = 'double')
        self.assertEqual( 'i_1',n(i1) )
        self.assertEqual( 'i_2',n(i2) )
        self.assertEqual( 'i_3',n(i3) )
        self.assertEqual( 'dp_1',n(dp1) )
        self.assertEqual( 'dp_2',n(dp2) )
        self.assertEqual( 'dp_3',n(dp3) )
        n = namer(separator='X',default_suffix='000')
        for count in range(1,100):
            i = dsl_mock(prefix = 'i', vartype = 'int' )
            dp = dsl_mock(prefix = 'dp', vartype = 'double' )
            self.assertEqual('iX'+str(count)+'000',n(i) )
            self.assertEqual('dpX'+str(count)+'000',n(dp) )
        for count in range(100,200):
            i = dsl_mock(prefix = 'i', vartype = 'int' )
            dp = dsl_mock(prefix = 'dp', vartype = 'double' )
            self.assertEqual('iX'+str(count)+'000',n(i) )
            self.assertEqual('dpX'+str(count)+'000',n(dp) )

        # Test reset_counters()
        n.reset_counters()
        for count in range(1,100):
            i = dsl_mock(prefix = 'i', vartype = 'int' )
            dp = dsl_mock(prefix = 'dp', vartype = 'double' )
            self.assertEqual('iX'+str(count)+'000',n(i) )
            self.assertEqual('dpX'+str(count)+'000',n(dp) )

