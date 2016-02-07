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
from intception.read_json import *
from unittest.mock import mock_open # allows mocking of open() built-in function

class TestRead_json(unittest.TestCase):
    """Tests for read_json.py classes and their methods."""
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_read_json(self):
        # Test read_json() function
        json_text = """{ "a" : 1, "b" : { "c" : 2, "d" : "donut", "e": 42 }, "f" : "fjord" }"""
        data = read_json('filename',opener = mock_open(read_data=json_text) )
        self.assertEqual( data["a"], 1 )
        self.assertTrue( isinstance( data["b"], dict ) )
        self.assertEqual( data["b"]["c"], 2 )
        self.assertEqual( data["b"]["d"], "donut" )
        self.assertEqual( data["b"]["e"], 42 )
        self.assertEqual( data["f"], "fjord" )
        # FileNotFound
        self.assertRaises( FileNotFoundError, read_json, '', opener=open )
        # Improperly formed JSON (ValueError)
        json_text = """error"""
        self.assertRaises( ValueError, read_json,'filename',opener = mock_open(read_data=json_text) )

