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
from intception.dsl import *
from intception.generator.arrays import general_array
from intception.dsl_extensions import *
#from intception.printer import printer

class wrapper(classtools.classtools):
    """
    Parent class containing methods and attributes general to all wrappers.
    """
    def __init__(self,gen):
        self._generator = gen # generator object this wrapper is associated with
    
    def generator(self):
        return self._generator

class index_wrapper(wrapper):
    """
    Class that wraps around a dsl_index object (or derived type), to associate
    source code-specific variables with the abstract dsl_index object.

    This allows indexing variables to be attached to a given dsl_index.
    """
    def __init__(self,index,gen):
        wrapper.__init__(self,gen)
        const_args = index.constant_attachments()
        self._index = index
        # output_array_skip and output_array_length are always passed as arguments to main integral function, 
        # so may be considered a const
        self._output_array_skip  = dsl_scalar( 'iskip'+index.name(), vartype= 'int', constant = const_args )
        self._output_array_length= dsl_scalar( 'n'+index.name(), vartype = 'int', constant = const_args )
        #
        self._work_array_skip  = dsl_scalar( 'iwskip'+index.name(), vartype= 'int' )
        self._work_array_offset = dsl_scalar( 'iwoff'+index.name(), vartype= 'int' )
        self._work_array_length = dsl_scalar( 'nw'+index.name(), vartype = 'int' )
        self._loop_length = dsl_scalar( 'nloop'+index.name(), vartype = 'int' )
        self._array_index = dsl_scalar( 'idx'+index.name(), vartype = 'int' )
        if index.is_cartesian() == True:
            # If index is Cartesian, use the usual 'l' and 'lmax' naming convention
            # index_value is always used as an argument to the main function, so may be considered a const
            self._index_value = dsl_scalar( 'l'+index.name(), vartype = 'int', constant = const_args )
            #
            self._index_max   = dsl_scalar( 'lmax'+index.name(), vartype = 'int' )
            self._index_loop_max = dsl_scalar( 'lmax_loop'+index.name(), vartype = 'int' )
        elif index.is_cartesian() == False:
            # ... otherwise use a more generic naming convention
            # index_value is always used as an argument to the main function, so may be considered a const
            self._index_value = dsl_scalar( 'ival'+index.name(), vartype = 'int', constant = const_args )
            #
            self._index_max   = dsl_scalar( 'imax'+index.name(), vartype = 'int' )
            self._index_loop_max = dsl_scalar( 'imax_loop'+index.name(), vartype = 'int' )
        # skip dict allows additional skip variables to be attached (e.g. for HRRs)
        # length dict allows additional counter variables to be attached (e.g. for HRRs)
        self._skip_dict = {}
        self._length_dict = {}
        # array_index_dict allows additional array_index variables to be attached (e.g. for contractions)
        self._array_index_dict = {}
        # the contraction array and spherical transformation array are per-index general_array objects
        # contract_array_X is always passed as an argument to the main integral function, so can use const
        dcontract_array = dsl_pointer( name = 'contract_array_'+index.name(), vartype = 'double', constant = True )
        icontract       = dsl_scalar( name = 'icontract'+index.name(), vartype = 'int' )
        self._contract_array  = general_array( pointer = dcontract_array,\
                                          array_index = icontract, dimensions = None)
        # sph_trans_array_X is always passed as an argument to the main integral function, so can use const
        dsph_trans_array = dsl_pointer( name = 'sph_trans_array_'+index.name(), vartype = 'double', constant = True )
        isph_trans       = dsl_scalar( name = 'isph_trans'+index.name(), vartype = 'int' )
        self._sph_trans_array  = general_array( pointer = dsph_trans_array,\
                                              array_index = isph_trans, dimensions = None)

    def add_skip_length(self,key,skip_var,length_var):
        assert isinstance(skip_var,dsl_variable),"skip must be a dsl_variable"
        assert isinstance(length_var,dsl_variable),"length must be a dsl_variable"
        self._skip_dict[key] = skip_var
        self._length_dict[key] = length_var

    def add_array_index(self,key,array_index_var):
        assert isinstance(array_index_var,dsl_variable),"array_index must be a dsl_variable"
        self._array_index_dict[key] = array_index_var

    def index(self):
        return self._index 

    def output_array_skip(self):
        return self._output_array_skip

    def work_array_skip(self):
        return self._work_array_skip

    def work_array_offset(self):
        return self._work_array_offset

    def work_array_length(self):
        return self._work_array_length

    def array_index(self):
        return self._array_index

    def output_array_length(self):
        return self._output_array_length

    def loop_length(self):
        return self._loop_length
    
    def index_value(self):
        return self._index_value

    def index_max(self):
        return self._index_max

    def index_loop_max(self):
        return self._index_loop_max

    def skip_dict(self):
        return self._skip_dict

    def length_dict(self):
        return self._length_dict

    def array_index_dict(self):
        return self._array_index_dict

    def contract_array(self):
        return self._contract_array

    def sph_trans_array(self):
        return self._sph_trans_array

class rr_wrapper(wrapper):
    """
    Class that wraps around a dsl_rr object for a general RR. Should be inherited by classes
    for specific types of RR.
    """
    def __init__(self,rr,wintegral,gen):
        wrapper.__init__(self,gen)
        self._rr = rr
        self._wintegral = wintegral
        self._name = None
        self._rr_function = None

    def rr(self):
        return self._rr

    def name(self):
        return self._name

    def rr_function(self):
        return self._rr_function

    def set_rr_function(self,rr_function):
        assert isinstance(rr_function,function_wrapper),\
            "f_rr must be a function_wrapper object"
        self._rr_function = rr_function

class hrr_wrapper(rr_wrapper):
    """
    Class that wraps around a dsl_rr object for a HRR, to associate source code-specific
    variables with the abstract dsl_rr object.
    """
    def __init__(self,rr,move_from_index,move_to_index,loop_index_list,wintegral,gen):
        rr_wrapper.__init__(self,rr,wintegral,gen)
        self._move_from_index = move_from_index
        self._move_to_index   = move_to_index
        self._loop_index_list = loop_index_list

        # Variables (dsl_* ) related to calling the HRR subroutine
        # arg_list = [ Rfromto, lmax_move_from, l_move_from, loop_l_move_to, \
        #              iskip_move_from0, iskip_move_to0,\
        #              iskip_move_from1, iskip_move_to1,\
        #              dwork_array, iwrk0, iwrk1 ]
        dposition_list = expr_list_objects( rr.expr(), dsl_position )
        assert len(dposition_list) == 1, "Only a single dsl_position object should be in rr.expr()"
        self._Rfromto = dposition_list[0]
        move_from_windex = self._wintegral.windex_index_dict()[ move_from_index.name() ] 
        move_to_windex = self._wintegral.windex_index_dict()[ move_to_index.name() ] 
        self._iskip_move_from0 = move_from_windex.work_array_skip()
        self._iskip_move_to0   = move_to_windex.work_array_skip()
        self._iskip_move_from1 = move_from_windex.skip_dict()['tmp']
        self._iskip_move_to1   = move_to_windex.skip_dict()['tmp']

        
        self._l_move_from    = move_from_windex.index_value()
        self._loop_l_move_to = move_to_windex.index_loop_max()
        self._lmax_move_from = move_to_windex.index_value() + move_from_windex.index_value() -\
                                move_to_windex.index_loop_max()


        # Create a list of variables related to calling the HRR subroutine so these can
        # be added to the local_variable_list for an integral_wrapper
        self._hrr_specific_arg_list = [ self._Rfromto,\
                    self._lmax_move_from, self._l_move_from, self._loop_l_move_to,\
                    self._iskip_move_from0, self._iskip_move_to0,\
                    self._iskip_move_from1, self._iskip_move_to1 ]
        # NB. dwork_array, iwrk0, iwrk1 are created independently of a hrr_wrapper object
        #     so are not added to this list
        if self.generator().generator_options().contracted() == True:
            # Contracted integral additional arguments
            self._ncont_move_from       = move_from_windex.length_dict()['cont']
            self._ncont_move_to          = move_to_windex.length_dict()['cont']
            self._iskip_cont_move_from0 = move_from_windex.skip_dict()['contwork']
            self._iskip_cont_move_to0   = move_to_windex.skip_dict()['contwork']
            self._iskip_cont_move_from1 = move_from_windex.skip_dict()['conttmp']
            self._iskip_cont_move_to1   = move_to_windex.skip_dict()['conttmp']
            contracted_arg_list = [ self._ncont_move_from, self._ncont_move_to, \
                                    self._iskip_cont_move_from0, self._iskip_cont_move_to0, \
                                    self._iskip_cont_move_from1, self._iskip_cont_move_to1 ]
            for var in contracted_arg_list:
                self._hrr_specific_arg_list.append( var )

        # Generate unique name for RR
        self._name = self.generate_unique_name()
        # Create list of dsl_integral objects in simplified expression
        self._dintegral_list = expr_list_objects( rr.expr(), dsl_integral )

    def generate_unique_name(self):
        """
        Creates a unique name for the RR, which can then be used as a function
        name in source code. 

        The naming convention is

        self.rr().name()+move_from_index.name()+move_to_index.name()+'_'+
            [ i.name() for i in self.loop_index_list ]
        """
        funcname = []
        funcname.append( self.rr().name() )
        funcname.append( self.move_from_index().name() )
        funcname.append( self.move_to_index().name() )
        funcname.append( '_')
        for index in self.loop_index_list():
            funcname.append( index.name() )
        if len( self.loop_index_list() ) == 0:
            del funcname[-1]
        return ''.join( funcname )

    def move_from_index(self):
        return self._move_from_index

    def move_to_index(self):
        return self._move_to_index

    def changing_index(self):
        # Note changing_index is the same as move_to_index
        return self._move_to_index

    def loop_index_list(self):
        return self._loop_index_list

    def dintegral_list(self):
        return self._dintegral_list

    ### Integral class HRR-specific instances of general HRR subroutine arguments ###
    def hrr_specific_arg_list(self):
        return self._hrr_specific_arg_list
    
    def contracted_hrr_specific_arg_list(self):
        return self._contracted_hrr_specific_arg_list

    def Rfromto(self):
        return self._Rfromto

    def lmax_move_from(self):
        return self._lmax_move_from

    def l_move_from(self):
        return self._l_move_from

    def loop_l_move_to(self):
        return self._loop_l_move_to

    def iskip_move_from0(self):
        return self._iskip_move_from0

    def iskip_move_to0(self):
        return self._iskip_move_to0

    def iskip_move_from1(self):
        return self._iskip_move_from0

    def iskip_move_to1(self):
        return self._iskip_move_to1

class trr_wrapper(rr_wrapper):
    """
    Class that wraps around a dsl_rr object for a TRR, to associate source code-specific
    variables with the abstract dsl_rr object.
    """
    pass # Not implemented
    
class vrr_wrapper(rr_wrapper):
    """
    Class that wraps around a dsl_rr object for a VRR, to associate source code-specific
    variables with the abstract dsl_rr object.
    """
    def __init__(self,rr,integral_index_list,loop_index_list,wintegral,gen):
        rr_wrapper.__init__(self,rr,wintegral,gen)
        self._loop_index_list = loop_index_list
        assert self._rr.rrtype() == 'vrr', 'Only VRRs allowed in vrr_wrapper'
        self._unrolled_index = self._rr.changing_index()
        self._changing_index = self._rr.changing_index() 
        # Simplify rr.expr() based on loop_index_list
        expr = self._rr.expr()
        for index in integral_index_list:
            if not index in self._loop_index_list and\
               not index is self._unrolled_index:
                expr = expr_find_and_replace( expr, index, dsl_zero() )
        self._simplified_expr = expr_simplify( expr ) 
        # Generate unique name for RR
        self._name = self.generate_unique_name()
        # Create list of dsl_integral objects in simplified expression
        self._dintegral_list = expr_list_objects( self.simplified_expr(), dsl_integral )
        # Create dictionary of increment/decrement values for each integral index for the
        # whole set of RRs for each auxiliary index increment/decrement value
        # Store the decrement / increment pairs as tuples ( index_decrement, aux_increment )
        # Create dictionary of increment values for each auxiliary index for the
        # whole set of RRs
        index_change_dict = {} 
        for dindex in self._wintegral.integral().index_list():
            index_change_dict[dindex.name()] = []
        aux_change_dict = {}
        for daux in self._wintegral.integral().aux_list():
            aux_change_dict[daux.name()] = []
        for dintegral in self.dintegral_list():
            if len( dintegral.aux_list() ) > 0:
                for daux in dintegral.aux_list():
                    binop = dintegral.aux_binop( daux )
                    ival  = binop.right()
                    if binop.op() is op_add:   x = 1
                    elif binop.op() is op_sub: x = -1
                    else: raise Exception( 'Only op_add and op_sub allowed' )
                    aux_change = x * ival
                    aux_change_dict[ daux.name() ].append( aux_change )
                    for dindex in dintegral.index_list():
                        binop = dintegral.index_binop( dindex )
                        ival  = binop.right()
                        if binop.op() is op_add:   x = 1
                        elif binop.op() is op_sub: x = -1
                        else: raise Exception( 'Only op_add and op_sub allowed' )
                        change = x * ival
                        index_change_dict[ dindex.name() ].append( ( change, aux_change ) )
            else:
                    for dindex in dintegral.index_list():
                        binop = dintegral.index_binop( dindex )
                        ival  = binop.right()
                        if binop.op() is op_add:   x = 1
                        elif binop.op() is op_sub: x = -1
                        else: raise Exception( 'Only op_add and op_sub allowed' )
                        change = x * ival
                        index_change_dict[ dindex.name() ].append( ( change, None ) )

        self._index_change_dict = index_change_dict
        self._aux_change_dict = aux_change_dict
        # unrolled_index_lmax is not set when __init__() is run, but when the 
        # integral evaluation algorithm is setup in integral_wrapper.setup_algorithm
        self._unrolled_index_lmax = None
                
    def generate_unique_name(self):
        """
        Creates a unique name for the RR, which can then be used as a function
        name in source code. 

        The naming convention is

        self.rr().name()+unrolled_index.name()+'_'+
            [ i.name() for i in self.integral_index_list ]
        """
        funcname = []
        funcname.append( self.rr().name() )
        funcname.append( self.changing_index().name() )
        funcname.append( '_')
        for index in self.loop_index_list():
            funcname.append( index.name() )
        if len( self.loop_index_list() ) == 0:
            del funcname[-1]
        return ''.join( funcname )

    def set_unrolled_index_lmax(self,lmax):
        assert isinstance(lmax,int), "unrolled_index_lmax must be an integer"
        self._unrolled_index_lmax = lmax

    def loop_index_list(self):
        return self._loop_index_list

    def simplified_expr(self):
        return self._simplified_expr

    def unrolled_index(self):
        return self._unrolled_index

    def changing_index(self):
        return self._changing_index

    def dintegral_list(self):
        return self._dintegral_list

    def index_change_dict(self):
        return self._index_change_dict

    def aux_change_dict(self):
        return self._aux_change_dict

    def unrolled_index_lmax(self):
        return self._unrolled_index_lmax


class special_function_wrapper(wrapper):
    """
    Class that wraps around a special function, encapsulating all information necessary
    to produce source code for the evaluation of the special function.

    In this context a "special function" is a non-standard function that is evaluated
    in an integral base class. It may depend on an auxiliary index.

    The special function is described by the encapsulated function_wrapper object.
    Additional supporting data may be added in classes that inherit from this.
    """
    def __init__(self,function_wrapper_obj,gen):
        wrapper.__init__(self,gen)
        self._function_wrapper_obj = function_wrapper_obj

    def function_wrapper_obj(self):
        return self._function_wrapper_obj

    def special_function_call_out(self, p, array_pointer, array_skip, aux_min, aux_max, x_arg ):
        """
        Standard interface for output of a special function call.

        Calls the appropriate function_wrapper.call_out() method, with arguments.
        """
        #assert isinstance(p, printer), "p must be a printer object"
        self.function_wrapper_obj().call_out( p, args = [ array_pointer, array_skip,\
                                                          aux_min, aux_max, x_arg ] )

class boys_wrapper(special_function_wrapper):
    """
    Class that wraps around a boys_f object, encapsulating all information necessary
    to produce source code for the evaluation of the Boys function (including 
    arrays of tabulated date).

    This class relies on the interface and implementation of boys_f and boys_function
    in boys.py. An alternative Boys function implementation would require a different
    wrapper class.
    """

    def __init__(self,function_wrapper_obj,oof_data_array,smallx_data_array,\
                      mediumx_data_array,gen):
        special_function_wrapper.__init__(self,function_wrapper_obj,gen)
        self._oof_data_array = oof_data_array
        self._smallx_data_array = smallx_data_array
        self._mediumx_data_array = mediumx_data_array

    def oof_data_array(self):
        return self._oof_data_array

    def smallx_data_array(self):
        return self._smallx_data_array

    def mediumx_data_array(self):
        return self._mediumx_data_array

    def special_function_call_out(self, p, array_pointer, array_skip, aux_min, aux_max, x_arg ):
        """
        Overload special_function_call_out specifically for the built-in Boys function
        implementation.
        """
        #assert isinstance(p, printer), "p must be a printer object"
        self.function_wrapper_obj().call_out( p, args = [ array_pointer, array_skip, aux_min, aux_max, x_arg ] )



class function_wrapper(wrapper):
    """
    Class that wraps around a dsl_function object, allowing the Python code 
    necessary to generate the source code for output to be associated with 
    an abstract dsl_function object.

    The function for generating the code should be able to generate all the code
    associated with the dsl_function object using only a integral_wrapper object.
    """
    def __init__(self,f_dsl,f_call,f_prototype,f_source,gen):
        """
        f_dsl:          dsl_function object which describes abstract details
                        of the function which will eventually be translated to
                        source code.
        f_prototype:    Python function which generates a function prototype
                        for output. It should take two arguments:
                        a integral_wrapper object and a printer object.
        f_source:       Python function which generates source code which 
                        describes the behaviour of the function defined in 
                        dsl_function. It should have the same two arguments
                        as f_prototype.
        """
        wrapper.__init__(self,gen)
        self._f_dsl = f_dsl
        self._f_call = f_call
        self._f_prototype = f_prototype
        self._f_source = f_source

    def f_dsl(self):
        return self._f_dsl

    def f_call(self):
        return self._f_call 

    def f_prototype(self):
        return self._f_prototype

    def f_source(self):
        return self._f_source

    def f_const_dict(self):
        return self._f_const_dict

    def call_out(self,p,args=None):
        """
        Calls the callable object self._f_call, which generates source code based on
        an integral_wrapper object. 
        
        All data required to generate the source code is carried with the callable object.

        p:      printer object for output of source code
        """
        self._f_call(p,args)

    def prototype_out(self,p):
        """
        Calls the callable object self._f_prototype, which generates source code based on
        an integral_wrapper object. 
        
        All data required to generate the source code is carried with the callable object.

        p:      printer object for output of source code
        """
        self._f_prototype(p)

    def source_out(self,p):
        """
        Calls the function self._f_source, which generates source code based on
        an integral_wrapper object.
        
        All data required to generate the source code is carried with the callable object.

        p:      printer object for output of source code
        """
        self._f_source(p)

    class callable_function_call:
        """
        Instances are callable objects carrying with them all the information
        necessary to generate a function call. 
        
        This is a parent class which should be inherited by child classes for 
        each type of function call desired.
        """
        def __init__(self,f_dsl, gen):
            assert isinstance(f_dsl,dsl_function),"f_dsl must be a dsl_function object"
            self._f_dsl = f_dsl
            self._generator = gen
        def generator(self):
            return self._generator
        def __call__(self,p,args=None):
            """
            args:    specify an alternative set of arguments to the normal set contained
                     in self._f_dsl
            """
            #assert isinstance(p,printer),"p must be a printer object"
            f = self._f_dsl
            p.funccall( f, args )

    class callable_prototype:
        """
        Instances are callable objects carrying with them all the information
        necessary to generate a function prototype. 
        
        This is a parent class which should be inherited by child classes for 
        each type of function prototype desired.
        """
        def __init__(self,f_dsl, gen, const_dict = {}):
            assert isinstance(f_dsl,dsl_function),"f_dsl must be a dsl_function object"
            self._f_dsl = f_dsl
            self._generator = gen
            # const_dict is an optional argument which allows arguments attached to the
            # f_dsl to be given the const qualifier in the function prototype.
            # This is useful where, for example, we want to assert that a variable is
            # const inside a function, but the corresponding dsl_* object may not have
            # the is_constant attribute set.
            # If both is_constant is set and a const_dict entry is present, the const_dict
            # entry may override a False is_constant attribute. If is_constant == True, then
            # a corresponding const_dict entry must agree.
            self._const_dict = const_dict
        def generator(self):
            return self._generator
        def const_dict(self):
            return self._const_dict
        def __call__(self,p):
            #assert isinstance(p,printer),"p must be a printer object"
            f = self._f_dsl
            p.funcprototype( f, self.const_dict() )
        
    class callable_source:
        """
        Instances are callable objects carrying with them all the information
        necessary to generate executable source code for a given function..
        
        This is a parent class which should be inherited by child classes for 
        each type of function prototype desired.
        """
        def __init__(self,f_dsl, gen, const_dict = {}):
            assert isinstance(f_dsl,dsl_function),"f_dsl must be a dsl_function object"
            self._f_dsl = f_dsl
            self._generator = gen
            # const_dict is an optional argument which allows arguments attached to the
            # f_dsl to be given the const qualifier in the function prototype.
            # This is useful where, for example, we want to assert that a variable is
            # const inside a function, but the corresponding dsl_* object may not have
            # the is_constant attribute set.
            # If both is_constant is set and a const_dict entry is present, the const_dict
            # entry may override a False is_constant attribute. If is_constant == True, then
            # a corresponding const_dict entry must agree.
            self._const_dict = const_dict

        def generator(self):
            return self._generator
        def const_dict(self):
            return self._const_dict
        def __call__(self,p):
            """
            The function source code output will be specific to each child object,
            so this returns NotImplemented.
            """
            #assert isinstance(p,printer),"p must be a printer object"
            return NotImplemented


