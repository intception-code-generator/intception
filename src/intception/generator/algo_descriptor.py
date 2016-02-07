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
Provides the algo_descriptor class, which allows broad algorithm classes to be unambiguously 
identified using a hash (e.g. in a dictionary).
"""

from intception import classtools

class algo_descriptor(classtools.classtools):
    def __init__(self,is_contracted=None,hrr_present=None,aux_index_present=None):
        """
        Key algorithm properties must be set here.

        If additional algorithm properties are added, they must default to None, so that
        algorithm types for which they are irrelevant can still be uniquely described
        by an instance of algo_descriptor. 

        Algorithm properties should be added to tuple_to_hash in __hash__ to ensure that
        algorithms hash to a unique value based on algorithm attributes.
        """
        self.is_contracted     = self.set_attrib(is_contracted,bool)
        self.hrr_present       = self.set_attrib(hrr_present,bool)
        self.aux_index_present = self.set_attrib(aux_index_present,bool)

    def set_attrib(self,attrib,allowed_type):
        """
        Return attrib, if attrib is of allowed type or None.
        """
        if attrib != None:
            assert isinstance(attrib,allowed_type),"attrib must be of type "+str(allowed_type)
        return attrib

    def __hash__(self):
        """
        Returns a hash of object attributes. 
        This must ONLY DEPEND ON OBJECT VALUES, since the hash function
        will be used to match objects with the same attributes but different
        id's.
        """
        tuple_to_hash = (\
                        self.is_contracted,\
                        self.hrr_present,\
                        self.aux_index_present\
                        )
        return hash(tuple_to_hash)

    def __eq__(self,other):
        """ 
        Tests for equality of objects by checking values of attributes are equal.
        """
        return self.__dict__ == other.__dict__
