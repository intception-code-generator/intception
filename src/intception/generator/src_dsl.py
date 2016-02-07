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
Provides a classes which encapsulate dsl_* objects from dsl.py. These src_dsl_* objects
carry with them source-code-specific attributes and methods. Automatic dsl_* to src_dsl_*
conversion of expressions (comprised of binary and unary operations) is performed by the
converter class (converter.py). 

Where possible, dsl_* expressions should be manipulated as dsl_* objects until they are
needed to be output as source code. The cprinter class in src_printer.py performs automatic 
conversion when dsl_* objects are passed to its methods as arguments.

The dsl_binop and dsl_unop objects from the dsl_* representation do not have corresponding
src_dsl_* classes. This is because these objects are sufficiently general, and do not themselves
represent "things" in the output source (i.e. they are not variables, or functions, but are
more abstract).

Since the dsl_* was originally developed with C in mind, the src_dsl_* objects can often use
attributes and methods of the corresponding dsl_* class without modification.
"""


import numbers
from intception.dsl import *
from intception.dsl_extensions import *
from intception import classtools

class src_dsl_object(classtools.classtools):
    """
    Parent class for objects that are used to provide valid src() methods for
    dsl_* objects in DSL expressions.
    """
    def __init__(self,dobj,converter):
        self._dobj      = dobj
        self._converter = converter

    def dobj(self):
        return self._dobj

    def converter(self):
        return self._converter

    def type(self):
        # We can simply use the type string from the dsl_* object, since the dsl_* was
        # defined using C-types. A generator for a different language would need to
        # map dsl_* types to src_dsl_* types.
        return self.dobj().type()

    def name(self):
        # We can simply use the name string from the dsl_* object, since the dsl_* was
        # defined to use C-compatible name strings. 
        # A generator for a different language may need to map dsl_* names to 
        # src_dsl_* names.
        return self.dobj().name()

    # Overload the Python operators (as in dsl_base) to allow src_dsl_* objects to be
    # manipulated in expressions.
    # This avoids having to create new expressions using dsl_* objects and then converting
    # them to src_dsl_* objects.
    # In this generator, we use the original op_dict from dsl.py, since this was setup
    # for C source code generation. 
    # Binary operations built into Python
    def __add__(self,other):
        return dsl_binop(op_dict['+'], self, other ) 
    def __radd__(self,other):
        return dsl_binop(op_dict['+'], other, self ) 
    def __sub__(self,other):
        return dsl_binop(op_dict['-'], self, other ) 
    def __rsub__(self,other):
        return dsl_binop(op_dict['-'], other, self ) 
    def __mul__(self,other):
        return dsl_binop(op_dict['*'], self, other ) 
    def __rmul__(self,other):
        return dsl_binop(op_dict['*'], other, self ) 
    def __truediv__(self,other):
        return dsl_binop(op_dict['/'], self, other ) 
    def __rtruediv__(self,other):
        return dsl_binop(op_dict['/'], other, self ) 
    def __floordiv__(self,other):
        return dsl_binop(op_dict['//'], self, other ) 
    def __rfloordiv__(self,other):
        return dsl_binop(op_dict['//'], other, self ) 
    def __pow__(self,other):
        return dsl_binop(op_dict['**'], self, other ) 
    def __rpow__(self,other):
        return dsl_binop(op_dict['**'], other, self ) 
    # Non-overloaded binary operations
    def assign(self,expr):
        return dsl_binop(op_assign,self,expr)
    # Unary operations built into Python
    def __neg__(self):
        return dsl_unop(op_dict['u-'], self )
    def __pos__(self):
        return dsl_unop(op_dict['u+'], self )
    def __abs__(self):
        return dsl_unop(op_dict['abs'], self )

class src_dsl_index(src_dsl_object):
    """
    Instances of this wrap around individual dsl_index objects in a
    DSL expression and offer a src() method which outputs correctly outputs
    dsl_index objects in source code.
    """
    def __init__(self,dindex,converter):
        src_dsl_object.__init__(self,dindex,converter)
        assert isinstance(dindex,dsl_index),\
                'dindex must be dsl_index object'
        self._cc_array = self.converter().wintegral().cc_array()

    def src(self):
        """Returns dsl_index object for output in source code.
        
        The integral_wrapper object associated with converter is used to determine
        what form the source code output should take.
        
        If the index is the "current_index" associated with the integral_wrapper
        work_array, then the value of that index is always output. If the index
        is not the "current_index", then, for a Cartesian index, the cc_array
        associated with the integral_wrapper object will be referenced. 
        
        Implementation notes:
         * For a non-Cartesian index, no behaviour is defined where this index is not
           the "current_index".
        """
        current_index = self.converter().wintegral().work_array().current_index()
        dindex = self.dobj()
        converter = self.converter()
        wintegral = self.converter().wintegral()
        try:
            windex    = wintegral.windex_index_dict()[ dindex.name() ]
        except KeyError:
            windex    = wintegral.windex_aux_dict()[ dindex.name() ]
        cc_array  = self.cc_array()
        direction = self.converter().direction()
        assert hasattr(dindex,'is_cartesian'),"index/aux must have is_cartesian() method"
        if dindex.is_cartesian() == True:
            # dindex is Cartesian, use direction.
            if dindex is current_index:
                if dindex[direction.index()].value() == 0:
                    return src( dsl_zero() )
                else:
                    # recursively call src() of component dsl_index object
                    return src( dindex[ direction.index() ].value() )
            else:
                src_cc_array_ref = src_general_array(\
                        cc_array[ windex.array_index() ][ direction.index() ], converter )
                return src( src_cc_array_ref )
        elif dindex.is_cartesian() == False:
            # dindex is not Cartesian, direction irrelevant.
            # [ Current implementation restriction ]
            assert dindex is current_index,\
                    "src() output for non-Cartesian indexes is not defined where index "+\
                    "is not the current_index of integral_wrapper.work_array object"
            if dindex is current_index:
                if dindex.value() == 0:
                    return src( dsl_zero() )
                else:
                    return src( dindex.value() )

    def cc_array(self):
        return self._cc_array

class src_dsl_pointer(src_dsl_object):
    """
    Instances of this wrap around individual dsl_pointer objects in a
    DSL expression and offer a src() method which outputs correctly outputs
    dsl_variable objects in source code.
    """
    def __init__(self,dpointer,converter):
        src_dsl_object.__init__(self,dpointer,converter)
        assert isinstance(dpointer,dsl_pointer),\
                'dpointer must be dsl_pointer object'

    def src(self,set_array_index = True):
        """Returns string corresponding to dsl_pointer object for output in source code"""
        converter = self.converter()
        out = []
        out.append( src( self.name() ) )
        if converter.contraction_info() != None and set_array_index == True:
            # Contraction required, and array_index setting requested
            # Determine corresponding array index from contraction_obj
            try:
                array_index = converter.contraction_info().array_index_dict[ self.dobj() ]
                out.append( '[ ' )
                out.append( src( array_index ) )
                out.append( ' ]' )
            except KeyError:
                # If there is no key corresponding to self, then skip adding an array index
                pass
        return src( ''.join(out) )

    def __getitem__(self,key):
        """
        Returns a new src_dsl_pointer with __getitem__ called on the encapsulated
        dsl_pointer object.
        """
        dpointer = self.dobj()
        converter = self.converter()
        return converter.expr_full_convert( dpointer[key] )

    def declaration(self):
        """Returns string corresponding to dsl_pointer for output in source code. 
        This is specific to declaration lines, where direction is not relevant."""
        if self.dobj().is_constant() == True:
            const_str = 'const '
        else:
            const_str = ''
        return const_str + src( self.type() ) + ' ' + src( self.name() )

    #def src(self,alt_getitem_list = None):
    #    """Returns string corresponding to dsl_pointer object for output in source code"""
    #    dpointer  = self.dobj()
    #    converter = self.converter() 
    #    if alt_getitem_list != None:
    #        getitem_list = alt_getitem_list
    #    else:
    #        getitem_list = self._getitem_list[:]
    #    out = []
    #    pointer_name = dpointer.name()
    #    out.append( src( pointer_name ) )
    #    for key in getitem_list:
    #        src_key = converter.expr_to_src_convert( key )
    #        #src_key = key
    #        out.append( '[ ' )
    #        out.append( src( src_key ) )
    #        out.append( ' ]' )
    #    return ''.join(out)

    #def __getitem__(self,key):
    #    """Returns a new src_dsl_pointer with __getitem__ called on the encapsulated
    #    dsl_pointer object."""
    #    # A pointer object can only have a single dimension:
    #    assert len( self._getitem_list ) == 0, "src_dsl_pointer objects can only have "+\
    #                                           "a single array index"
    #    new_getitem_list = self._getitem_list[:]
    #    new_getitem_list.append( key )
    #    return src_dsl_pointer( self.dobj(), self.converter(), new_getitem_list )

    #def declaration(self):
    #    """Returns string corresponding to dsl_pointer for output in source code. 
    #    This is specific to declaration lines, where direction is not relevant."""
    #    dpointer = self.dobj()
    #    return src( self.type() ) + ' ' + src( self )

    def assign(self,src_expr):
        """Returns a dsl_binop object corresponding to the assignment of the
        src_dsl_pointer to src_expr."""
        return dsl_binop( op_assign, self, src_expr )

    def type(self):
        # Modify dsl_* type to add an asterisk.
        return src_dsl_object.type(self)+' *'

class src_dsl_variable(src_dsl_object):
    """
    Instances of this wrap around individual dsl_variable objects in a
    DSL expression and offer a src() method which outputs correctly outputs
    dsl_variable objects in source code.
    """
    def __init__(self,dvariable,converter):
        src_dsl_object.__init__(self,dvariable,converter)
        assert isinstance(dvariable,dsl_variable),\
                'dvariable must be dsl_variable object'

    def src(self,set_direction=True):
        """Returns string corresponding to dsl_variable object for output in source code"""
        dvariable = self.dobj()
        direction = self.converter().direction()
        if dvariable.is_cartesian() == True and set_direction == True:
            return src( self[ direction.index() ] )
        else:
            return src( self.name() )

    def __getitem__(self,key):
        """
        Returns a new src_dsl_variable with __getitem__ called on the encapsulated
        dsl_variable object.
        """
        dvariable = self.dobj()
        converter = self.converter()
        return converter.expr_full_convert( dvariable[key] )

    def declaration(self):
        """Returns string corresponding to dsl_variable for output in source code. 
        This is specific to declaration lines, where direction is not relevant."""
        dvariable = self.dobj()
        direction = self.converter().direction()
        if self.dobj().is_constant() == True:
            const_str = 'const '
        else:
            const_str = ''
        if len( dvariable ) != 1:
            return const_str + src( self.type() ) + ' ' + src( self[ len(dvariable) ] ) 
        else:
            return const_str + src( self.type() ) + ' ' + src( self.name() )

    def assign(self,src_expr):
        """Returns a dsl_binop object corresponding to the assignment of the
        src_dsl_variable to src_expr."""
        return dsl_binop( op_assign, self, src_expr )

class src_dsl_value(src_dsl_object):
    """
    Instances of this wrap around individual dsl_value objects in a
    DSL expression and offer a src() method which outputs correctly outputs
    dsl_value objects in source code.
    """
    def __init__(self,dvalue,converter):
        src_dsl_object.__init__(self,dvalue,converter)
        assert isinstance(dvalue,dsl_value),\
                'dindex must be dsl_value object'
        self._dvalue = dvalue # dsl_value object (from DSL expr)

    def src(self):
        """Returns dsl_value object for output in source code"""
        dvalue = self.dobj()
        return src( dvalue.value() )

    def value(self):
        """Returns the the value of the src_dsl_value object."""
        return self.dobj().value()

class src_dsl_function(src_dsl_object):
    """
    Instances of this wrap around individual dsl_function objects in a
    DSL expression and offer a src() method which outputs correctly outputs
    dsl_function objects in source code. call() and prototype() methods are also
    provided.
    """
    def __init__(self,dfunction,converter): 
        src_dsl_object.__init__(self,dfunction,converter)
        assert isinstance(dfunction,dsl_function),\
                'dindex must be dsl_value object'
        self._dfunction = dfunction # dsl_function object (from DSL expr)

    def src(self):
        """Returns dsl_function object for output in source code.
        
        This is the function call, with default arguments.
        """
        return self.call() 

    def prototype(self, const_dict = {}):
        """Returns dsl_function prototype for output in source code"""
        dfunction = self.dobj()
        converter = self.converter()
        out = []
        out.append( src( self.type() ) )
        out.append( ' ' )
        out.append( src( self.name() ) )
        args = dfunction.args()
        if len( args ) > 0 :
            out.append( '( ' )
            for arg in args :
                src_arg = converter.expr_full_convert( arg )
                length = len( src_arg.dobj() )
                # Check if a constant via const_dict or is_constant attribute
                if arg in const_dict.keys() and const_dict[arg] == True:
                    out.append( 'const' )
                    out.append( ' ' )
                elif arg in const_dict.keys() and const_dict[arg] == False:
                    # const_dict[arg] == False cannot override is_constant == True
                    assert src_arg.dobj().is_constant() == False
                elif not arg in const_dict.keys():
                    if src_arg.dobj().is_constant() == True:
                        # No const_dict entry, but dsl_* object has is_constant == True
                        out.append( 'const' )
                        out.append( ' ' )
                else:
                    raise Exception('Unexpected behaviour of const_dict and is_constant attribute')
                if length == 1: # scalar argument (dsl_scalar)
                    out.append( src_arg.type() )
                    out.append( ' ' )
                    out.append( src_arg.name() )
                    out.append( ', ' )
                elif length > 1: # array argument (dsl_position)
                    out.append( src_arg.type() )
                    out.append( ' ' )
                    out.append( src_arg[length].src( set_direction = False ) )
                    out.append( ', ' )
                else:
                    raise Exception( src_arg.name()+' length not understood' )
            del out[-1]
            out.append( ' )' )
        else:
            out.append( '()' )
        return ''.join(out)

    def call(self,_args=None):
        """Returns dsl_function call for output in source code"""
        dfunction = self.dobj()
        converter = self.converter()
        out = []
        out.append( src( self.name() ) )
        if _args == None: 
            args = dfunction.args()
        else:
            args = _args
        if len(args) > 0:
            out.append( '( ' )
            for arg in args:
                src_arg = converter.expr_full_convert( arg )
                if isinstance(src_arg,src_dsl_variable) and\
                        src_arg.dobj().is_cartesian() and\
                        arg.dobj().component_of() == None:
                    # In a function call, we want to pass the first element of 
                    # a Cartesian src_dsl_variable, when a specific component of that 
                    # variable has not been requested
                    out.append( src_arg.src( set_direction = False ) )
                else:
                    out.append( src( src_arg ) )
                out.append( ', ' )
            del out[-1]
            out.append( ' )' )
        else:
            out.append( '()' )
        return ''.join( out )


class src_dsl_integral(src_dsl_object):
    """
    Instances of this wrap around individual dsl_integral objects in a
    DSL expression and offer a src() method which outputs correctly indexes
    integrals in that object.
    """
    def __init__(self,dintegral,converter):
        src_dsl_object.__init__(self,dintegral,converter)
        self._work_array = self.converter().wintegral().work_array() # integral_array object
        self._jump_array = self.converter().wintegral().jump_array()
        self._cc_array   = self.converter().wintegral().cc_array()

    def src(self):
        """Returns correctly index array element for integral"""
        dintegral  = self.dobj()
        converter  = self.converter()
        direction  = self.converter().direction()
        work_array = self.work_array()
        jump_array = self.jump_array()
        work_array_offset     = work_array.array_offset()
        jump_array_min_change = jump_array.support_data()['min change']
        current_dindex  = work_array.current_index()
        current_windex  = converter.wintegral().windex_index_dict()[ current_dindex.name() ]
        loop_index_list = dintegral.index_list()[:]
        # Remove index currently being incremented
        loop_index_list.remove( current_dindex )
        if len( dintegral.aux_list() ) > 0 : 
            # If auxiliary index present, need to select correct array_index object
            assert len( dintegral.aux_list() ) == 1,\
                    "Only one auxiliary index per integral class currently supported."
            daux = dintegral.aux_list()[0]
            # Use increment from corresponding dsl_binop to determine the correct
            # work_array array_index object to use
            aux_binop = dintegral.aux_binop( daux )
            change = aux_binop.right() 
            if change > 0: assert aux_binop.op() is op_add,\
                            "Only increments allowed in auxiliary index."
            # tuple[0] is array_index, tuple[1] is array_index_start
            work_array_index = work_array.array_index_dict()[change][0]
        else:
            # If no auxiliary index present, only one possible array_index object
            work_array_index = work_array.array_index()

        # Create expr for indexing of work array element
        work_array_index_binop = work_array_index

        # Indexing for index currently being looped over (inside subroutine)
        # can be explicit.
        index_binop = dintegral.index_binop( current_dindex )
        assert index_binop.op() is op_add or index_binop.op() is op_sub,\
                "only addition and subtraction allowed for index change"
        if index_binop.op() is op_add:
            change = index_binop.right()
        elif index_binop.op() is op_sub:
            change = -1 * index_binop.right()
        jump = index_binop.left().jump( change, direction ) +\
                work_array_offset.value() + current_dindex.index() - 1
        if jump > 0:
            work_array_index_binop = dsl_binop(op_add, work_array_index_binop, current_windex.work_array_skip() * jump )
        elif jump < 0:
            work_array_index_binop = dsl_binop(op_sub, work_array_index_binop, current_windex.work_array_skip() * jump )
        elif jump == 0:
            pass # jump from work_array().array_index() is zero, so work_array_index_binop
                 # does not change
        # Indexing for remaining Cartesian indices requires use of jump array
        # as Cartesian component is set at runtime.
        for obj in loop_index_list:
            try:
                obj_binop = dintegral.index_binop( obj )
                windex = converter.wintegral().windex_index_dict()[ obj.name() ] 
            except KeyError:
                obj_binop = dintegral.aux_binop( obj )
                windex = converter.wintegral().windex_aux_dict()[ obj.name() ]
            assert hasattr(obj,'is_cartesian'),"index/aux must have is_cartesian() method"
            if obj.is_cartesian() == True:
                # Use jump_array
                if obj_binop.right() >= 1:
                    if obj_binop.op() == op_add:
                        change_value = abs( jump_array_min_change ) + obj_binop.right() - 1
                    elif obj_binop.op() == op_sub:
                        change_value = abs( jump_array_min_change ) - 1 * obj_binop.right()
                    else:
                        raise Exception('operator '+str(obj_binop.op())+' not allowed.')
                    jump = jump_array[ windex.array_index() ][ change_value ][ direction.index() ]
                    # always add this number because the jump array elements carry the sign with them
                    src_jump_ref = src_general_array( jump, converter )
                    work_array_index_binop = dsl_binop(op_add, work_array_index_binop, windex.work_array_skip() * src_jump_ref ) 
            elif obj.is_cartesian() == False:
                # Change is simply integer change
                if obj_binop.right() >= 1:
                    if obj_binop.op() == op_add:
                        change_value = obj_binop.right()
                    elif obj_binop.op() == op_sub:
                        change_value = -1 * obj_binop.right()
                    else:
                        raise Exception('operator '+str(obj_binop.op())+' not allowed.')
                    jump = change_value 
                    work_array_index_binop = dsl_binop(op_add, work_array_index_binop, windex.work_array_skip() * jump ) 

        src_work_array_ref = src_general_array( work_array[ work_array_index_binop ],\
                                                converter )
        return src( src_work_array_ref )
#        return src( work_array[ work_array_index_binop ] )

    def work_array(self):
        return self._work_array

    def jump_array(self):
        return self._jump_array

    def cc_array(self):
        return self._cc_array

class src_general_array:
    """
    Instances of this wrap around individual general_array objects in a
    DSL expression and offer a src() method which outputs correctly outputs
    general_array (and derived types) in source code.

    Does not inherit from src_dsl_object, since the behaviour of general_array
    is different to dsl_* objects.
    """
    def __init__(self,array_obj,converter):
        self._array_obj = array_obj
        self._converter = converter

    def array_obj(self):
        return self._array_obj

    def converter(self):
        return self._converter

    def src(self,alt_getitem_list = None):
        """
        Returns string corresponding to general_array (or derived type) object in 
        C output format.

        Encapsulation in the src_* representation allows a converter object to be 
        associated with the general_array object, so that conversions of dsl_*
        objects in the getitem_list of the object can be done correctly.
        """
        dpointer  = self.array_obj().pointer()
        converter = self.converter() 
        if alt_getitem_list != None:
            getitem_list = alt_getitem_list
        else:
           getitem_list = self.array_obj().getitem_list()
        if self.array_obj().dimensions() != None:
            assert len( getitem_list ) == len( self.array_obj().dimensions() ),\
                    "array reference must have correct dimensionality"
        out = []
        array_name = dpointer.name()
        out.append( src( array_name ) )
        for key in getitem_list:
            # TODO
            # This is a bottleneck, since everytime an integral is output (many times)
            # the expr_full_converter must be run on the indexing expression.
            # /TODO
            src_key = converter.expr_full_convert( key )
            #src_key = key 
            out.append( '[ ' )
            out.append( src( src_key ) )
            out.append( ' ]' )
        return ''.join(out)

