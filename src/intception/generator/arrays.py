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
from intception import classtools
from intception.dsl import *

class general_array(classtools.classtools):    
    """Parent class with methods and attributes common to representations of arrays."""
    def __init__(self,pointer,array_index,dimensions,getitem_list = []):
        assert isinstance(pointer,dsl_pointer), "pointer must be dsl_pointer"
        assert isinstance(getitem_list,list), "getitem_list must be Python list"
        self._pointer = pointer
        if array_index != None:
            self.set_array_index(array_index)
        else:
            self._array_index = None
        if dimensions != None: 
            self.set_dimensions(dimensions)
        else:
            self._dimensions = None
        self._getitem_list = getitem_list

    def pointer(self):
        return self._pointer

    def array_index(self):
        return self._array_index

    def dimensions(self):
        return self._dimensions

    def getitem_list(self):
        return self._getitem_list

    def set_array_index(self,array_index):
        assert isinstance(array_index,dsl_scalar), "array_index must be dsl_scalar"
        self._array_index = array_index 

    def set_dimensions(self,dimensions):
        assert isinstance(dimensions,list), "dimensions must be Python list"
        self._dimensions = dimensions

    def __getitem__(self,key):
        getitem_list = self.getitem_list()[:]
        getitem_list.append(key)
        return general_array( self.pointer(), self.array_index(), self.dimensions(), getitem_list )

    def assign(self,expr):
        """
        Returns a dsl_binop object assigning the array to an expression.

        Typically this will be used in combination with __getitem__ so that an array element
        can be assigned.
        """
        return dsl_binop(op_assign,self,expr)

class integral_array(general_array):
    """Wrapper class for dsl_scalar objects which refer to an array used to store
    integrals."""
    def __init__(self,pointer,array_index,getitem_list = []):
        general_array.__init__(self,pointer,array_index,[None],getitem_list)
        # the dimensions list has a single entry of None, because while the number of 
        # dimensions is known (1, len( dimensions )), the length of the single dimension
        # is undefined
        # array_index_dict allows additional dsl_scalar objects to be attached to the
        # integral_array object
        # typically the the keys will be integer auxiliary index change (0,1,2) values
        # and the values will be tuples containing and array_index dsl_scalar, and also
        # an array_index_start dsl_scalar, with the array_index_start being used to store
        # the first element of a "m-layer" and the array_index being incremented withing
        # that layer
        # NB. the array_index_dict is also used for the 2-layer HRR algorithm
        self._array_index_dict = {}
        # special_array_index_dict allows additional dsl_scalar objects to be attached to
        # the integral_array object with special purposes, the keys will typically be
        # a string keyword describing the special array_index and the value will be the
        # dsl_scalar object (not a tuple).
        # e.g. for placing the output of the base_function at a custom location in the
        # work array, key = 'base', value = dsl_scalar('iwrkb',vartype='int')
        self._special_array_index_dict = {}
        array_index_start = dsl_scalar( array_index.name()+'_start', vartype='int' )
        self.add_array_index( 0, array_index, array_index_start )
        self._array_offset = dsl_integer_index('offset',0) # integer offset
        # Set current_index to a dsl_index object
        # that is currently being incremented (e.g. in a VRR or HRR).
        self._current_index = None 
        # length_expr can optionally be set to a src_dsl_* expression which describes the
        # 1D length of the array (e.g. for use in memory allocation during runtime)
        self._length_expr = None
        # windex_list can optionally be set to an ordered list of index_wrapper objects
        # which contain the array indexing information for each integral object
        self._windex_list = None

    def __getitem__(self,key):
        """
        Overloads the __getitem__ method of general_array to be specific to integral_array.
        """
        getitem_list = self.getitem_list()[:]
        getitem_list.append(key)
        return integral_array( self.pointer(), self.array_index(), getitem_list )

    def add_array_index(self,key,array_index,array_index_start):
        """
        Add an additional tuple of array index variables to self._array_index_dict.

        In typical use, key would be a value for a change in an auxiliary index, 
        allowing different array_index dsl_scalar objects to refer to different
        auxiliary index "layers" in an array.

        The values will be tuples containing and array_index dsl_scalar, and also
        an array_index_start dsl_scalar, with the array_index_start being used to store
        the first element of a "m-layer" and the array_index being incremented withing
        that layer
        """
        assert isinstance(array_index,dsl_scalar), "array_index must be dsl_scalar"
        assert isinstance(array_index_start,dsl_scalar), "array_index_start must be dsl_scalar"
        self._array_index_dict[key] = ( array_index, array_index_start )

    def add_special_array_index(self,key,array_index):
        """
        Sets the optional array_index_base attribute.

        array_index_base should be a dsl_scalar, and should be used to refer to the 
        start index at which the output of a base function is placed.
        """
        assert isinstance(array_index,dsl_scalar), "array_index must be dsl_scalar"

        self._special_array_index_dict[key] = array_index


    def special_array_index_dict(self):
        return self._special_array_index_dict

    def array_index_dict(self):
        return self._array_index_dict

    def array_offset(self):
        return self._array_offset

    def current_index(self):
        return self._current_index

    def length_expr(self):
        return self._length_expr

    def windex_list(self):
        return self._windex_list

    def set_offset(self,value):
        self._array_offset.set_value(value)

    def set_current_index(self,index):
        self._current_index = index

    def set_length_expr(self,expr):
        self._length_expr = expr

    def set_windex_list(self,windex_list):
        self._windex_list = windex_list

class data_array(general_array):
    """
    Wrapper class for dsl_scalar objects which refer to an array used to store
    data, such as integer values for Cartesian components.
    """
    def __init__(self,pointer,dimensions,data,getitem_list = [],support_data=None):
        """
        pointer:  dsl_pointer object which represents array name and type
        dimensions: Python list with len() == number of dimensions, and each element
                    spec.fying the length of that dimension.
        data:   Python list with dimensionality and length equal to that given in 
                dimensions. Carries literal data to be output in source code.
                There is a special control string for inserting linebreaks in output
                data arrays. If an array element is "#linebreak#", then when it is 
                output by definition_out(), this is replaced by a line break in the source
                code.
        getitem_list: list containing keys from __getitem__ calls, which will be used
                when src() is called.
        support_data: optional attribute for additional data (may be in any format,
                recommend using a dictionary)
        """
        general_array.__init__(self,pointer,None,dimensions,getitem_list)
        if data != None: 
            self.set_data(data)
        else:
            self._data = None
        if support_data != None: 
            self.set_support_data(support_data)
        else:
            self._support_data = None

    def support_data(self):
        return self._support_data

    def data(self):
        return self._data

    def set_support_data(self,support_data):
        self._support_data = support_data

    def set_data(self,data):
        assert isinstance(data,list), "data must be Python list"
        self._data = data 

    def __getitem__(self,key):
        """
        Overloads the __getitem__ method of general_array to be specific to integral_array.
        """
        getitem_list = self.getitem_list()[:]
        getitem_list.append(key)
        return data_array( self.pointer(), self.dimensions(), self.data(),\
                                     getitem_list, self.support_data() )

    def declaration_out(self,p,extern=True):
        """
        Outputs a declaration for the data array, as would usually be placed in a C header file.
        
        p:      cprinter object for output of source code
        extern: if True, variable is declared with extern keyword, so it can be shared among 
                source files.
        """
        src_self = p.converter.expr_full_convert( self )
        src_pointer = p.converter.expr_full_convert( self.pointer() )
        out = []
        if extern == True:
            out.append('extern')
            out.append(' ')
        if self.pointer().is_constant() == True:
            out.append('const')
            out.append(' ')
        out.append( src_pointer.dobj().type() ) 
                                            # no asterisk, as we dereference using array
                                            # style (square brackets)
        out.append( ' ' )
        out.append( src_self.src( alt_getitem_list = self.dimensions() ) )
        p.out( ''.join(out) )

    def definition_out(self,p):
        """
        Outputs a definition of the data array (i.e. assigning data to the array), as would
        usually be placed in a C source file.

        p:      cprinter object for output of source code
        """

        def nested_list_out(out,l):
            """
            Recursively output a list of strings that can be joined to represent 
            nested Python lists in C source code format.

            out:    Python list into which strings to be joined are placed.
            l:      Python (nested) list which will be parsed.
            """
            out.append('{ ')
            for item in l:
                if isinstance(item,list):
                    nested_list_out( out, item )
                    out.append(', ')
                elif item == '#linebreak#':
                    # special linebreak string should not have a comma following it
                    out.append( str(item) )
                else:
                    out.append( str(item) )
                    out.append(', ')
            if out[-1] == '#linebreak#' and out[-2] == ', ':
                del out[-2:]
            elif out[-1] == ', ':
                del out[-1]
            out.append(' }')

        src_self    = p.converter.expr_full_convert( self )
        src_pointer = p.converter.expr_full_convert( self.pointer() )
        if self.pointer().is_constant() == True:
            const_str = 'const '
        else:
            const_str = ''
        p.out( const_str + src_pointer.dobj().type()+' '+\
                src_self.src( alt_getitem_list = self.dimensions() ) +\
                ' = ', endl='' )
        p.indent()
        out = []
        nested_list_out( out, self.data() )
        line = []
        for item in out:
            if item == '#linebreak#':
                # special linebreak string ends current line 
                p.out( ''.join(line), endl='' )
                line = []
            else:
                line.append(item)
        if len( line ) > 0: p.out(''.join(line))
        p.outdent()
