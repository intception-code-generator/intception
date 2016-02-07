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
"""
Module providing a simple timer class. This can be used as a context manager using 
Python's "with" statement.
"""
import time

class timer:
    """
    Provides __enter__ and __exit__ methods for timing of events inside a "with" block.
    """
    def __enter__(self):
        """Set start_time on entering block."""
        self.start_time = time.clock()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Set end_time on exiting block."""
        self.end_time   = time.clock()
        self.time_taken = self.end_time - self.start_time
        return None
