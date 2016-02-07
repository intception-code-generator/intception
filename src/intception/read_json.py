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

# Open JSON input file
def read_json(filename,opener=open):
    """
    Opens a JSON file, containing integral information, and uses json module
    to place return a dictionary of json data.
    
    filename:   name of json file
    opener:     name of function used to open file (syntax should be the same as 
                the built-in Python open() function) -- useful for testing.
    """
    import json
    try:
        fd = opener(filename,'r')
    except FileNotFoundError:
        raise FileNotFoundError("Input file not found. Exiting.")
    try:
        data = json.load(fd)
        fd.close()
        return data
    except ValueError:
        raise ValueError("Improperly formed JSON file. See traceback below.")


