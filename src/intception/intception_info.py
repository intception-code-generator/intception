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
Provides the prog_unit_info class, instances of which carry information relating to a 
program unit.

Also provides an instance of this class containing information about the
current version of Intception.
"""
from intception import classtools

class prog_unit_info(classtools.classtools):
    """
    Contains information relating to a program unit, e.g. Intception as a whole,
    or an associated code generator.
    """
    def __init__(self, unit_name_str, unit_desc_str, version_str, author_str, license_str, \
                 extra_str_list=None):
        """
        unit_name_str:   name of the program unit, e.g. "Intception".
        unit_desc_str:   descriptive phrase for program unit, e.g. 
                         "Molecular integral code generator".
        version_str:     string describing version of unit (if applicable).
        author_str:      string describing author or authors of unit.
        license_str:     string describing license of unit.
        extra_str:       optional list of strings containing additional
                         information (one item per line output)
        """
        self._unit_name_str = unit_name_str
        self._unit_desc_str = unit_desc_str
        self._version_str   = version_str
        self._author_str    = author_str
        self._license_str   = license_str
        self._extra_str_list= extra_str_list


    def unit_name_str(self):
        return self._unit_name_str

    def unit_desc_str(self):
        return self._unit_desc_str

    def version_str(self):
        return self._version_str

    def author_str(self):
        return self._author_str

    def license_str(self):
        return self._license_str

    def extra_str_list(self):
        return self._extra_str_list

    def formatted_info_str(self,stdout_width):
        """
        Returns formatted program unit information which will fit within
        stdout_width. 
        Useful for providing information to user. 
        
        If any formatted string exceeds provided stdout_width, then an 
        AssertionError exception is raised.
        """
        out = []
        unit_str        = self.unit_name_str() + ": " + \
                             self.unit_desc_str()
        unit_author_str = "Developed by "+ self.author_str()
        unit_version_str= "Version: "    + self.version_str()
        unit_license_str= "License: "    + self.license_str()
        unit_extra_str_list  = self.extra_str_list()
        # Ensure all strings fit within the stdout
        s_n_tuple_list = [ ( unit_str, 'unit_str'),\
                      ( unit_author_str, 'unit_author_str' ),\
                      ( unit_version_str, 'unit_version_str' ),\
                      ( unit_license_str, 'unit_license_str' ),\
                      ]
        if unit_extra_str_list != None: 
            iline = 1
            for line in unit_extra_str_list:
                s_n_tuple_list.append( ( line, 'unit_extra_str line' + str(iline)) )
                iline += 1
        for s, n in s_n_tuple_list:
                assert len( s ) < stdout_width, \
                        n+" must be less than stdout_width ("+str(stdout_width)+")."
        # Build list of centered lines to output
        out.append( unit_str.center(stdout_width) )
        out.append( unit_author_str.center(stdout_width) )
        out.append( unit_version_str.center(stdout_width) )
        out.append( unit_license_str.center(stdout_width) )
        if unit_extra_str_list != None:
            for line in unit_extra_str_list:
                out.append( line.center(stdout_width) )
        return '\n'.join( out )

# Create instance of prog_unit_info class for Intception as a whole
intception_info = prog_unit_info(\
                    unit_name_str = "Intception",\
                    unit_desc_str = "Molecular integral code generator",\
                    version_str   = "1.0",\
                    author_str    = "James C. Womack",\
                    license_str   = "GNU General Public License version 3", \
                    extra_str_list= [  "", # blank line
                            "This program comes with ABSOLUTELY NO WARRANTY.", 
                            "This is free software, and you are welcome to redistribute it",
                            "under certain conditions. See the LICENSE file for details."]
                )

