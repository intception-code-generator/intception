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
Provides the integral_wrapper class, which wraps around a single dsl_integral object 
representing an abstract integral type.

For each dsl_integral object there can be a corresponding integral_wrapper object, which has associated with it all the fundamental information needed to output source code for 
evaluation of the integral type described in the dsl_integral object, e.g. printer objects, 
function_wrapper objects for the main integral evaluation function and supporting functions,
header file dependencies etc.

The corresponding generator object (in generator.py) creates integral_wrapper objects using 
dsl_integral objects and a set of customisable options relating to the general form of the 
source code. The integral_wrapper object is then parsed and the appropriate output methods 
called in the generator.integral_class_out method.

classes:
    integral wrapper    --  wraps around a dsl_integral object and contains methods and 
                            data necessary to output corresponding C source code for the 
                            evaluation of the integral type represented by the dsl_integral 
                            object.
"""

# General intception modules
from intception.dsl import *
from intception.dsl_extensions import *
from intception.boys import boys_function, boys_f
#from intception.printer import printer

# Generator-specific modules
from intception.generator.src_dsl import *
from intception.generator.supporting_functions import expr_recursive_search, determine_assignment_order
from intception.generator.wrappers import wrapper, index_wrapper, function_wrapper, vrr_wrapper, \
                                          hrr_wrapper, trr_wrapper, boys_wrapper
from intception.generator.arrays import general_array, integral_array, data_array
from intception.generator.algo_descriptor import algo_descriptor

# function_wrapper and callable_* derived classes
from intception.generator.callables import integral_wrapper_main_function
from intception.generator.callables import integral_wrapper_work_array_size_function
from intception.generator.callables import integral_wrapper_rr_function
from intception.generator.callables import integral_wrapper_base_function
import intception.generator.callables.main_source_algorithms as main_source_algorithms

# testing
from intception.generator.src_printer import cprinter
from intception.generator.converter import converter

class integral_wrapper(wrapper):
    """
    Class that wraps around a dsl_integral object to provide non-abstract
    attributes and methods related directly to code generation.

    This is specific to the generator class in generator.py and carries language and 
    implementation-specific attributes/methods. This integral_wrapper class is specific
    to the generation of C source code to the C99 standard.

    Instances of integral_wrapper carry function_wrapper objects that expose the 
    function_wrapper.prototype_out and function_wrapper.source_out methods for output
    of source code. 
    
    Each integral class has: 
        - a "main" function
            the function that outputs an array of evaluated integrals (double precision 
            numbers), taking information such as exponents, centres and angular momentum 
            as arguments.
        - a "base" function
            a function, called by the main function, which evaluates the zero-angular 
            momentum case of the integral class
        - "supporting" functions
            functions which are called by the main function to evaluate integrals, such 
            as a functions containing unrolled inner loops of VRRs.
        - work array size function
            a function which outputs an integer representing the maximum size of work array 
            required for a set of inputs like the angular momentum of the shell -- the 
            maximum of calls to work array size functions for each integral class to be 
            evaluated can be used to determine the size of work array required for all 
            integral classes.
    
    In addition, the integral_wrapper object has associated with is printer objects (from 
    printer.py) for the output code to a source file (.c) and header file (.h), filenames 
    for these, and a list of header files on which the integral class source code depends 
    (e.g. a Boys function header file).
    """
    def __init__(self,integral,full_name,source_filename,header_filename,gen):
        assert isinstance(integral,dsl_integral),\
                "integral must be a dsl_integral object"
        assert isinstance(source_filename,str), "source_filename must be a string"
        assert isinstance(header_filename,str), "header_filename must be a string"
        wrapper.__init__(self,gen)
        self._integral = integral
        self._source_filename = source_filename
        self._header_filename = header_filename

        # Full integral name used for source files and program units associated with this integral
        assert isinstance(full_name,str), "full_name must be a string"
        self._full_name = full_name

        # List of additional header files that are required
        self._header_dependency_list = []

        # printers should be set immediately before output is required and closed
        # immediately after output is completed
        self._source_file    = None
        self._source_printer = None
        self._header_file    = None
        self._header_printer = None

        # Create dsl_direction variable
        self._direction = dsl_direction('x') 

        # Create contraction_info object with empty attributes
        self._contraction_info_obj = integral_wrapper.contraction_info( {}, {} )
        
        # Create converter object (converts dsl_* to src_dsl*)
        self._converter = converter( self, self.direction(), self.contraction_info_obj(),\
                                     self.generator() )

        # Initialize empty source code-related attributes
        self._main_function = None
        self._base_function = None
        self._algorithm_obj = None
        self._work_array_size_function = None
        self._support_function_list = []
        self._integral_argument_list = [] # list of dsl_variable objects which are arguments
                                          # to the main integral evaluation function
        self._auto_local_variable_list = [] # list of dsl_variable objects in assignment order
                                       # such that interdependencies are satisfied 
                                       # for automatic assignment
        self._auto_prim_loop_variable_list = [] # list of dsl_variable objects which 
                                      # (for contracted integrals) are ordered to satisfy
                                      # interdependencies and are automatically assigned 
                                      # inside a loop over primitive exponents (empty for
                                      # primitive integrals)
        self._auto_prim_loop_depends_dict = None # dictionary where dsl_* objects in
                                      # auto_prim_loop_variable_list are keys, and 
                                      # values are sets containing dsl_* objects
                                      # corresponding to the per-primitive integral
                                      # data they depend on (i.e. exponents)
        self._manual_local_variable_list = [] # list of dsl_variable objects which
                                       # must be declared in the main function, but are
                                       # to be manually assigned
        self._wrr_list = [] # list of vrr_wrapper objects
        self._windex_index_list = [] # list of index_wrapper objects (integral index)
        self._windex_index_dict = {} # dict of index_wrapper objects with index.name() keys
        self._windex_aux_list = [] # list of index_wrapper objects (auxiliary index) 
        self._windex_aux_dict = {} # dict of index_wrapper objects with aux.name() keys
        self._work_array = None
        self._output_array = None
        self._cc_array = None
        self._jump_array = None
        self.setup()

    def add_header_dependency(self,header):
        """
        Adds an additional header dependency.

        header:     string with header filename, e.g. name.h

        Additional header dependencies are output at the top of integral source files
        as 

        #include"name.h"
        """
        assert isinstance(header,str), "header must be a string containing the .h filename"
        self._header_dependency_list.append(header)
       
    def setup(self):
        """
        Calls method to simplify RR expressions for eventual output.
        Calls the setup_* methods in order to set object attributes necessary for 
        source code output.
        """
        # Setup index_wrapper objects
        self.setup_windex_attribs()
        # Setup vrr_wrapper objects
        self.setup_wrr_list()
        self.setup_integral_arrays()
        self.setup_transform_arrays()
        self.setup_algorithm()
        #print( "\nIntegral: ", self.integral().name(), "\n" )
        #print( "\nIntegral: ", self.full_name(), "\n" )
        #print( self.algorithm_obj().algorithm_str() )
        self.setup_integral_array_indexing()
        self.setup_work_array_size_function()
        self.setup_data_arrays()
        self.setup_integral_argument_list()
        self.setup_local_variable_lists()
        self.setup_main_function()
        self.setup_rr_functions()
        self.setup_transform_functions()
        self.setup_base_function()

    def setup_windex_attribs(self):
        """
        Creates index_wrapper objects for each dsl_index object in self.integral().index_list()
        and self.integral().aux_list() which carries information about array indexing for these
        dsl_index objects.

        Adds exponents of dsl_cartesian_gaussian objects to 
        self.contraction_info_obj().contraction_convert_dict 
        so these may be converted to dsl_pointer objects when used in integral function
        arguments and also when used in expressions of derived variabless.

        Adds entries to self.contraction_info_obj().array_index_dict, so that the dsl_pointer
        objects referring to the arrays of exponents can be related to their corresponding
        dsl_scalar array index objects.
        """
        for i in self.integral().index_list():
            windex = index_wrapper( i, gen = self.generator() )
            self._windex_index_list.append( windex )
            self._windex_index_dict[ i.name() ] = windex
        for a in self.integral().aux_list():
            windex = index_wrapper( a, gen = self.generator() )
            self._windex_aux_list.append( windex )
            self._windex_aux_dict[ a.name() ] = windex
        # Detect if HRR operation in self.integral().rr_list() or contracted integrals required
        # and if so...
        # add additional skip and length attributes to index_wrapper to integral index index_wrapper objects
        hrr_present = False
        contract_present = False
        sph_trans_present = False
        for rr in self.integral().rr_list():
            if rr.rrtype() == 'hrr':
                hrr_present = True
                break
        if self.generator().generator_options().contracted() == True:
            contract_present = True
        if self.generator().generator_options().spherical_transformed() == True:
            sph_trans_present = True
        if hrr_present == True or contract_present == True or sph_trans_present == True:
            # Add additional skip and length variables to each integral index
            # Reused in HRR, contraction and spherical transform sequences
            # Add additional array_index variable
            for windex in self.windex_index_list():
                index_name = windex.index().name()
                windex.add_skip_length('tmp', dsl_scalar('iskip'+index_name+'_tmp',vartype='int'),\
                                              dsl_scalar('n'+index_name+'_tmp',vartype='int') )
                windex.add_array_index('tmp', dsl_scalar('i'+index_name+'_tmp',vartype='int') )
        if sph_trans_present == True:
            assert contract_present == True, \
                    "Spherical transformation is only supported for contracted integral classes."
            for windex in self.windex_index_list():
                index_name = windex.index().name()
                # isph_transX is used to locate the starting location of the array of spherical transformation
                # coefficients for index X
                windex.add_array_index('sph_trans', dsl_scalar('isph_trans'+index_name,vartype='int') )
        if contract_present == True:
            # If contracted integrals are required, add entries to i
            # self._contraction_convert_dict
            contraction_convert_dict = {}
            array_index_dict = {}
            for windex in self.windex_index_list():
                assert isinstance( windex.index(), dsl_cartesian_gaussian ),\
                        "If contracted == True, then only dsl_cartesian_gaussian indexes "+\
                        "are currently supported."
                # exponents converted from dsl_scalar to dsl_pointer, to allow for array
                # indexing
                assert isinstance( windex.index().exp(), dsl_scalar ),\
                        "Exponent expected to be instance of dsl_scalar."
                index_name = windex.index().name()
                dscalar_exponent  = windex.index().exp()
                dpointer_exponent = dsl_pointer( name = dscalar_exponent.name(),\
                                                 vartype = dscalar_exponent.type(), \
                                                 constant = dscalar_exponent.is_constant() )
                contraction_convert_dict[ dscalar_exponent ] = dpointer_exponent
                # If contracted integrals are required, add additional skip and length attributes
                # for number or primitives (i.e. number of exponents) and number of contractions
                # Also add additional array_index variable for indexing of exponents
                # Since iskipX_cont, nX_prim, nX_cont will be passed as arguments to the
                # main integral function, they can be set as const
                # iskipX_prim is a local variable, so should not be set as const
                windex.add_skip_length('prim', dsl_scalar('iskip'+index_name+'_prim',vartype='int', ),\
                                       dsl_scalar('n'+index_name+'_prim',vartype='int', constant = True) )
                windex.add_skip_length('cont', dsl_scalar('iskip'+index_name+'_cont',vartype='int', constant = True ),\
                                       dsl_scalar('n'+index_name+'_cont',vartype='int', constant = True) )
                windex.add_skip_length('contwork', dsl_scalar('iwskip'+index_name+'_cont',vartype='int' ),\
                                       dsl_scalar('nw'+index_name+'_cont',vartype='int') )
                windex.add_skip_length('conttmp', dsl_scalar('iskip'+index_name+'_cont_tmp',vartype='int' ),\
                                       dsl_scalar('n'+index_name+'_cont_tmp',vartype='int') )
                # iprimX is used to locate the correct exponent in an array of exponents for index X
                windex.add_array_index('prim', dsl_scalar('iprim'+index_name,vartype='int') )
                windex.add_array_index('cont', dsl_scalar('icont'+index_name,vartype='int') )
                # icontractX is used to locate the starting location of the array of contraction coefficients
                # for index X
                windex.add_array_index('contract', dsl_scalar('icontract'+index_name,vartype='int') )
                # Relate dsl_pointer objects for exponent arrays to iprim array_indexes
                array_index_dict[ dpointer_exponent ] = windex.array_index_dict()[ 'prim' ]
            self.contraction_info_obj().contraction_convert_dict = contraction_convert_dict
            self.contraction_info_obj().array_index_dict = array_index_dict
        if len( self.windex_aux_list() ) > 0:
            assert len( self.windex_aux_list() ) == 1, \
                    "Only 0 or 1 auxiliary indexes allowed per integral class."
            # Add additional length to index_wrapper object for auxiliary index, to represent
            # minimum value of auxiliary index.
            windex = self.windex_aux_list()[0]
            index_name = windex.index().name()
            windex.add_skip_length('loop min', dsl_scalar('dummy'), \
                    dsl_scalar('nloop'+index_name+'_min', vartype = 'int' ) )

    def setup_wrr_list(self):
        """
        Creates vrr_wrapper objects for each dsl_rr object in self.integral().rr_list()
        which carries information about the surrounding loop structure over indices 
        (non-auxiliary) and a simplified RR expression based on this.

        Implementation notes:
        * Allowed RRs are currently restricted to VRR and HGP-type HRR.
        """
        # [ Current implementation restriction ]
        for rr in self.integral().rr_list():
            assert rr.rrtype() in [ 'vrr', 'hrr' ], "Only VRR and HRR supported currently"

        integral_index_list = self.integral().index_list()
        rr_index_order_list = [] # order in which indexes are acted upon
        for rr in self.integral().rr_list():
            rr_index_order_list.append( rr.changing_index() )
        for rr in self.integral().rr_list():
            if rr.rrtype() == 'vrr':
                unrolled_index = rr.changing_index()
                i = list( reversed( rr_index_order_list ) ).index( unrolled_index )
                loop_index_list = list( reversed( rr_index_order_list ) )[i+1:]
                wrr = vrr_wrapper( rr, integral_index_list, loop_index_list,\
                        wintegral = self, gen = self.generator() )
            elif rr.rrtype() == 'hrr':
                move_to_index = rr.changing_index()
                # Identify move_from index using HRR expression
                # (this should be the only incremented index)
                dintegral_list = expr_list_objects( rr.expr(), dsl_integral )
                move_from_index = None
                for dintegral in dintegral_list:
                    for dindex in dintegral.index_list():
                        dindex_binop = dintegral.index_binop( dindex )
                        if dindex_binop.op() is op_add and \
                            dindex_binop.right() > 0:
                                assert dindex_binop.right() == 1,\
                                    "Increment in integral index in dependency for HRR should be 1"
                                assert move_from_index == None,\
                                    "HRR should only have one incremented index in dependency integral"
                                move_from_index = dindex
                
                i = list( reversed( rr_index_order_list ) ).index( move_to_index )
                loop_index_list = list( reversed( rr_index_order_list ) )[i+1:]
                loop_index_list.remove( move_from_index )
                wrr = hrr_wrapper( rr, move_from_index, move_to_index, loop_index_list,\
                        wintegral = self, gen = self.generator() )

            self._wrr_list.append( wrr )
        
    def setup_work_array_size_function(self):
        """
        Create a function_wrapper object corresponding to a function which returns an integer
        for the maximum length of the work_array required for this integral class. 
        When multiple integral classes are required, then the work_array_size functions for
        each integral class can be called and the maximum output can be used to set the size of the
        work array.
        """

        # Short names for callable_* derived classes
        work_array_size_call      = integral_wrapper_work_array_size_function.work_array_size_call
        work_array_size_prototype = integral_wrapper_work_array_size_function.work_array_size_prototype
        # The work_array_size_source class is a derived class in callables.main_source_algorithms
        # which is determined in setup_algorithm
        work_array_size_source    = self.algorithm_obj().algo_module.algo_work_array_size_source
        # Create function_wrapper object from work_array_size_prototype() and
        # work_array_size_source()
       
        dvar_list = []
        darg_list = []
        for windex in self.windex_index_list():
            darg_list.append( windex.index_value() )
        # For contracted integrals, work array size depends on number of primitives
        # and number of contractions
        if self.generator().generator_options().contracted() == True:
            for windex in self.windex_index_list():
                darg_list.append( windex.length_dict()['prim'] )
            for windex in self.windex_index_list():
                darg_list.append( windex.length_dict()['cont'] )

        # Setup a const_dict object, to assert that all variables passed in should be
        # const
        const_dict = {}
        for arg in darg_list:
            const_dict[arg] = True

        #funcname = self.integral().name()+'_'+self.work_array().pointer().name()+'_size' 
        funcname = self.full_name()+'_'+self.work_array().pointer().name()+'_size' 
        f_dsl = dsl_function(funcname,args = darg_list, vartype="int")
        f_call = work_array_size_call( f_dsl, gen = self.generator() )
        f_prototype = work_array_size_prototype( f_dsl, gen = self.generator(), const_dict = const_dict )
        f_source = work_array_size_source( f_dsl, wintegral = self, \
                local_variable_list = dvar_list, gen = self.generator(), const_dict = const_dict )
        wf = function_wrapper(f_dsl,f_call,f_prototype,f_source, gen = self.generator() )
        self._work_array_size_function = wf
        self._support_function_list.append( wf )

    class contraction_info(classtools.classtools):
        """
        Class encapsulating information about contractions.
        """
        def __init__(self, contraction_convert_dict = None, array_index_dict = None ):
            # contraction_convert_dict contains entries which map instances of dsl_* objects 
            # to instances of alternative classes. The dsl_* objects are the keys (these must 
            # be hashable) and the new objects are values. 
            # If self.generator().generator_options().contracted() == True, then this dictionary
            # is used to convert dsl_* objects to appropriate objects of other classes 
            # (e.g. dsl_scalar Gaussian exponents to dsl_pointer, for array access)
            self.contraction_convert_dict    = contraction_convert_dict
            # array_index_dict contains entries which map instances of objects to array indexes 
            # e.g. if a dsl_pointer is used to refer to an array of exponents xa[:] for index a then
            # the dictionary would map the dsl_pointer object (key) to a dsl_scalar object (value) for
            # array index iprima.
            self.array_index_dict = array_index_dict

        def contraction_convert_list(self,var_list):
            """
            Returns a new list, where objects matching keys in self.contraction_convert_dict
            are replaced by the corresponding values in self.contraction_convert_dict and
            the original ordering of the list is preserved.

            var_list:   expected to be an ordered iterable (generally a Python list)
            """
            new_list = []
            for var in var_list:
                try:
                    new_var = self.contraction_convert_dict[ var ]
                except KeyError:
                    # If arg not a key in contraction_convert_dict, do no conversion!
                    new_var = var
                new_list.append( new_var )
            return new_list


    class algorithm(classtools.classtools):
        """
        Stores data relevant to the integral evaluation algorithm.

        Note on ordering:
            * Lists of classes descring operation types, e.g.
                self.wvrr_group_list
                self.early_transform_list
                self.whrr_group_list
              ... are all ordered in order of execution in the generated code, i.e.
              the first list entry is executed first.
            * Ordered lists contained within wvrr_group, whrr_group, contraction objects 
              which contain indexes are ordered in "loop output" order, so the last item 
              would be relevant to the innermost loop, while the first item is relevant
              to the outermost loop.
        """
        def __init__(self, wvrr_group_list, early_transform_list, whrr_group_list, late_transform_list, \
                     m_min,m_max,m_layers,algo_module):
            self.wvrr_group_list  = wvrr_group_list
            self.early_transform_list = early_transform_list
            self.whrr_group_list  = whrr_group_list
            self.late_transform_list = late_transform_list
            self.m_min = m_min       # minimum value of auxiliary index required for zero-angular momentum case
            self.m_max = m_max       # maximum value of auxiliary index
            self.m_layers = m_layers # maximum number of layers required for auxiliary index
            # algo_module is one of the modules in 
            # intception.generator.callables.main_function_algorithms
            self.algo_module = algo_module


        class operation(classtools.classtools):
            def __init__(self, optype,\
                        start_order_list, end_order_list,\
                        start_index_list, end_index_list,\
                        start_cont_dict,  end_cont_dict,\
                        start_min_l_expr_dict, end_min_l_expr_dict,\
                        start_max_l_expr_dict, end_max_l_expr_dict,\
                        start_ncc_type_dict, end_ncc_type_dict):
                # optype is a string description of the type of operation, e.g.
                # 'vrr', 'hrr', 'contraction'
                self.optype = optype
                # start_order_list and end_order_list are lists of indexes (including 
                # auxiliary indexes) in the order they appear in the work_array before
                # an operation and after completion of the operation.
                # Even if indexes have not been incremented (0 angular momentum)
                # they are still listed, since primitive exponent/contraction indexes may
                # be present
                # The last entry is the "fastest changing" index, and the first is the
                # "slowest changing" index.
                self.start_order_list = start_order_list
                self.end_order_list   = end_order_list
                # start_index_list and end_index_list are the list of indexes that have been
                # incremented before the operation, and the list of indexes that
                # have been incremented after completion of the operation (typically
                # end_index_list will have 1 additional index present that was incremented for a
                # RR operation).
                self.start_index_list  = start_index_list
                self.end_index_list    = end_index_list
                # start_cont_dict and end_cont_dict contain strings which describe whether an 
                # integral index is contracted or primitive before and after the completion 
                # of the operation.
                self.start_cont_dict   = start_cont_dict
                self.end_cont_dict     = end_cont_dict
                # start_{min,max}_l_expr_dict and end_{min,max}_l_expr_dict contain expressions, 
                # in terms of user input angular momentum values, for the angular momentum in 
                # the block of useful integral intermediates before the operation, and 
                # after completion of the operation, the keys are the same as the 
                # {start,end}_index_list elements, and the values are determined by the algorithm
                # i.e. is more angular momentum required in an index at start to successfully
                # obtain the required intermediates at the end?
                # These effectively characterise the range of angular momentum in integral intermediates
                # that are necessary for the RR operation to produce the desired outcome.
                self.start_min_l_expr_dict = start_min_l_expr_dict
                self.end_min_l_expr_dict   = end_min_l_expr_dict
                self.start_max_l_expr_dict = start_max_l_expr_dict
                self.end_max_l_expr_dict   = end_max_l_expr_dict
                # {start,end}_ncc_type_dict contain strings that describe the type of expression that
                # is required to express the number of (Cartesian) components for each integral
                # index in the intermediate (work) array before (start) the RR operation and after
                # (end) the RR operation. Current strings:
                # 'int'  : integer index, so number of components is equal to range of values
                # 'ncc'  : (l+1)*(l+2)/2 (only a single value of l)
                # 'ncc0' : (l+1)*(l+2)*(l+3)/6 (range of l values from 0..l)
                # 'ncc0-ncc0' : 
                #          (l_max+1)*(l_max+2)*(l_max+3)/6 - l_min*(l_min+1)*(l_min+2)/6
                # 'nsph' : 2l + 1
                # Keys are index.name() for integral indexes.
                self.start_ncc_type_dict = start_ncc_type_dict
                self.end_ncc_type_dict   = end_ncc_type_dict

        class contraction(operation):
            def __init__(self,index_to_contract,loop_index_list,\
                         start_order_list, end_order_list,\
                         start_index_list, end_index_list,\
                         start_cont_dict,  end_cont_dict,\
                         start_min_l_expr_dict, end_min_l_expr_dict,\
                         start_max_l_expr_dict, end_max_l_expr_dict,\
                         start_ncc_type_dict, end_ncc_type_dict):
                # Use __init__ from parent class
                integral_wrapper.algorithm.operation.__init__(self,'contraction',\
                         start_order_list, end_order_list,\
                         start_index_list, end_index_list,\
                         start_cont_dict,  end_cont_dict,\
                         start_min_l_expr_dict, end_min_l_expr_dict,\
                         start_max_l_expr_dict, end_max_l_expr_dict,\
                         start_ncc_type_dict, end_ncc_type_dict)
                # Contraction-specific attributes
                # Index which will be contracted
                self.index_to_contract = index_to_contract
                # Indexes which must be looped over outside of subroutine, i.e.
                # indexes that are not yet contracted
                self.loop_index_list   = loop_index_list
                # To be set later...
                self.wfunction = None # function_wrapper object for calling contract

        class spherical_transform(operation):
            def __init__(self,index_to_transform,loop_index_list,\
                         start_order_list, end_order_list,\
                         start_index_list, end_index_list,\
                         start_cont_dict,  end_cont_dict,\
                         start_min_l_expr_dict, end_min_l_expr_dict,\
                         start_max_l_expr_dict, end_max_l_expr_dict,\
                         start_ncc_type_dict, end_ncc_type_dict):
                # Use __init__ from parent class
                integral_wrapper.algorithm.operation.__init__(self,'spherical transform',\
                         start_order_list, end_order_list,\
                         start_index_list, end_index_list,\
                         start_cont_dict,  end_cont_dict,\
                         start_min_l_expr_dict, end_min_l_expr_dict,\
                         start_max_l_expr_dict, end_max_l_expr_dict,\
                         start_ncc_type_dict, end_ncc_type_dict)
                # Spherical transform-specific attributes
                # Index which will be contracted
                self.index_to_transform = index_to_transform
                # Indexes which must be looped over outside of subroutine, i.e.
                # indexes that are not yet contracted
                self.loop_index_list   = loop_index_list
                # To be set later...
                self.wfunction = None # function_wrapper object for calling sph_trans

        class wrr_group(operation):
            def __init__(self, rrtype, wrr_list,\
                        start_order_list, end_order_list,\
                        start_index_list, end_index_list,\
                         start_cont_dict,  end_cont_dict,\
                        start_min_l_expr_dict, end_min_l_expr_dict,\
                        start_max_l_expr_dict, end_max_l_expr_dict,\
                        start_ncc_type_dict, end_ncc_type_dict):
                # Use __init__ from parent class
                integral_wrapper.algorithm.operation.__init__(self,rrtype,\
                         start_order_list, end_order_list,\
                         start_index_list, end_index_list,\
                         start_cont_dict,  end_cont_dict,\
                         start_min_l_expr_dict, end_min_l_expr_dict,\
                         start_max_l_expr_dict, end_max_l_expr_dict,\
                         start_ncc_type_dict, end_ncc_type_dict)
                # RR-specific attributes
                # List of RR operations
                self.wrr_list = wrr_list
                
        class whrr_group(wrr_group):
            def __init__(self, wrr,\
                         start_order_list, end_order_list,\
                         start_index_list, end_index_list,\
                         start_cont_dict,  end_cont_dict,\
                         start_min_l_expr_dict, end_min_l_expr_dict,\
                         start_max_l_expr_dict, end_max_l_expr_dict,\
                         start_ncc_type_dict, end_ncc_type_dict):
                integral_wrapper.algorithm.wrr_group.__init__(self, 'hrr', [ wrr ],
                         start_order_list, end_order_list,\
                         start_index_list, end_index_list,\
                         start_cont_dict,  end_cont_dict,\
                         start_min_l_expr_dict, end_min_l_expr_dict,\
                         start_max_l_expr_dict, end_max_l_expr_dict,\
                         start_ncc_type_dict, end_ncc_type_dict)
                ### Additional HRR-specific attributes ###
                self.wrr = wrr

        class wvrr_group(wrr_group):
            def __init__(self, wrr_list, index_loop_max_expr_dict, loop_length_expr_dict, \
                         start_order_list, end_order_list,\
                         start_index_list, end_index_list, \
                         start_l_expr_dict, end_l_expr_dict, \
                         auxiliary_index, start_m, end_m, m_layers ):
                # {start,end}_min_l_expr_dict for VRR operations always start from l = 0.
                start_min_l_expr_dict = {}
                end_min_l_expr_dict   = {}
                for k in start_l_expr_dict.keys():
                    start_min_l_expr_dict[ k ] = 0
                for k in end_l_expr_dict.keys():
                    end_min_l_expr_dict[ k ] = 0
                # {start,end}_ncc_type_dict for VRR operations will always be set to 'ncc0'
                start_ncc_type_dict = {}
                end_ncc_type_dict = {}
                for index in start_index_list:
                    assert index.is_cartesian() == True,\
                            "only Cartesian integral indexes supported currently"
                    start_ncc_type_dict[ index.name() ] = 'ncc0'
                for index in end_index_list:
                    assert index.is_cartesian() == True,\
                            "only Cartesian integral indexes supported currently"
                    end_ncc_type_dict[ index.name() ] = 'ncc0'
                # {start,end}_cont_dict for VRR operations are always set to 'prim'
                start_cont_dict = {}
                end_cont_dict   = {}
                for index in end_order_list:
                    start_cont_dict[ index.name() ] = 'prim'
                    end_cont_dict[ index.name() ]   = 'prim'

                # Use __init__ from parent class
                integral_wrapper.algorithm.wrr_group.__init__(self, 'vrr', wrr_list,
                         start_order_list, end_order_list,\
                         start_index_list, end_index_list,\
                         start_cont_dict,  end_cont_dict,\
                         start_min_l_expr_dict, end_min_l_expr_dict,\
                         start_l_expr_dict, end_l_expr_dict,\
                         start_ncc_type_dict, end_ncc_type_dict)
                ### Additional VRR-specific attributes ###
                # index_loop_max_expr refers to the expression for the maximum index value
                # which may depend on the iteration of the loop over auxiliary index for VRR 
                # type groups
                self.index_loop_max_expr_dict = index_loop_max_expr_dict
                # loop_length_expr refers to the length of the loop over a particular index
                # while the RR operation is performed in the innermost loop
                self.loop_length_expr_dict = loop_length_expr_dict
                # auxiliary_index is a dsl_integer_index object referring to an auxiliary index
                # if this must also be iterated over in the RR operation group -- if not, set to
                # None. Only a single auxiliary index is supported.
                self.auxiliary_index   = auxiliary_index
                # start_m and end_m are dsl expressions for the start and end values of the 
                # auxiliary index in terms of user set angular momentum values (determined at
                # generated-code runtime), set to None is auxiliary index does not need
                # incrementing in this RR group).
                self.start_m = start_m
                self.end_m   = end_m
                # m_layers is an integer which specifies the number of "layers" required in the
                # VRR evaluation algorithm, set to None if auxiliary index does not need 
                # incrementing in this RR group).
                # m_layers == 1 allows the use of an implied auxiliary index in some circumstances
                #               i.e. where wvrr_group is the last operation object in wvrr_group_list
                #               and where the wvrr in wrr_list is the last VRR operation to be 
                #               executed.
                # m_layers == 2 enfores the use of a two layer algorithm for all VRRs in the group
                self.m_layers = m_layers

        def algorithm_str(self):
            """
            Returns a string containing a human-readable summary of the algorithm for debugging.
            """
            def operation_out(operation_obj):
                nonlocal out
                out.append("start_order_list = "+\
                        str([ i.name() for i in operation_obj.start_order_list ] ))
                out.append("end_order_list   = "+\
                        str([ i.name() for i in  operation_obj.end_order_list ] ))
                out.append("start_index_list = "+\
                        str([ i.name() for i in operation_obj.start_index_list ] ))
                out.append("end_index_list   = "+\
                        str([ i.name() for i in  operation_obj.end_index_list ] ))
                out.append("start_cont_dict = "+str(\
                [ (k, src(v)) for k, v in operation_obj.start_cont_dict.items() ] ) )
                out.append("end_cont_dict = "+str(\
                [ (k, src(v)) for k, v in operation_obj.end_cont_dict.items() ] ) )
                out.append("start_min_l_expr_dict = "+str(\
                [ (k, src(v)) for k, v in operation_obj.start_min_l_expr_dict.items() ] ) )
                out.append("start_max_l_expr_dict = "+str(\
                [ (k, src(v)) for k, v in operation_obj.start_max_l_expr_dict.items() ] ) )
                out.append("start_ncc_type_dict = "+str(\
                [ (k, src(v)) for k, v in operation_obj.start_ncc_type_dict.items() ] ) )
                out.append("end_min_l_expr_dict = "+str(\
                [ (k, src(v)) for k, v in operation_obj.end_min_l_expr_dict.items() ] ) )
                out.append("end_max_l_expr_dict = "+str(\
                [ (k, src(v)) for k, v in operation_obj.end_max_l_expr_dict.items() ] ) )
                out.append("end_ncc_type_dict = "+str(\
                [ (k, src(v)) for k, v in operation_obj.end_ncc_type_dict.items() ] ) )

            out = []
            out.append("*** Algorithm summary ***")
            out.append("\nAuxiliary index algorithm information:\n")
            out.append("m_min    = "+str(self.m_min ))
            out.append("m_max    = "+str(self.m_max ))
            out.append("m_layers = "+str(self.m_layers))
            out.append("\nRR groups:\n")
            i = 1
            for wrr_group in self.wvrr_group_list:
                out.append("VRR group"+str(i))
                out.append("auxiliary_index = "+str(wrr_group.auxiliary_index))
                out.append("start_m  = "+str( wrr_group.start_m))
                out.append("end_m    = "+str( wrr_group.end_m))
                out.append("m_layers = "+str( wrr_group.m_layers))
                out.append("index_loop_max_expr_dict = "+str(
                    [ (k, src(v)) for k, v in wrr_group.index_loop_max_expr_dict.items() ] ) )
                out.append("loop_length_expr_dict = "+str(
                    [ (k, src(v)) for k, v in wrr_group.loop_length_expr_dict.items() ] ) )
                operation_out( wrr_group )
                i += 1
            i = 1
            for operation in self.early_transform_list:
                out.append( "Early transformation "+str(i) )
                if isinstance( operation, integral_wrapper.algorithm.contraction ):
                    out.append( "Contraction" )
                    out.append( "index_to_contract = "+str( operation.index_to_contract.name() ) )
                elif isinstance( operation, integral_wrapper.algorithm.spherical_transform ):
                    out.append( "Spherical transform" )
                    out.append( "index_to_transform = "+str( operation.index_to_transform.name() ) )
                out.append( "loop_index_list   = "+str([ i.name() for i in operation.loop_index_list ] ))
                operation_out( operation )
                i +=1
            i = 1
            for wrr_group in self.whrr_group_list:
                out.append("HRR group"+str(i))
                operation_out( wrr_group )
                i += 1
            i = 1
            for operation in self.late_transform_list:
                out.append( "Late transformation "+str(i) )
                if isinstance( operation, integral_wrapper.algorithm.contraction ):
                    out.append( "Contraction" )
                    out.append( "index_to_contract = "+str( operation.index_to_contract.name() ) )
                elif isinstance( operation, integral_wrapper.algorithm.spherical_transform ):
                    out.append( "Spherical transform" )
                    out.append( "index_to_transform = "+str( operation.index_to_transform.name() ) )
                out.append( "loop_index_list   = "+str([ i.name() for i in operation.loop_index_list ] ))
                operation_out( operation )
                i +=1

            return '\n'.join(out)

    def setup_algorithm(self):
        """
        Create an algorithm object which contains information relating
        to the algorithm used.

        Implementation notes:
        * Currently only a single HRR operation per integral class is supported.
        """
        # Attributes to be set for entire algorithm object:
        wvrr_group_list = []
        whrr_group_list = []
        m_min = None
        m_max = None
        m_layers = None
        # Dictionary of maximum angular momentum required in unrolled indexes, set for each vrr_wrapper
        unrolled_index_lmax_dict = {}
        # Attributes for final wvrr_group (the last to be executed at runtime):
        # (if whrr_group objects are setup before the wvrr_group, then these will depend on the HRR
        # operations)
        index_loop_max_expr_dict = {}
        loop_length_expr_dict = {}
        start_l_expr_dict = {}
        end_l_expr_dict = {}
        # Dictionary for relating a VRR operation to the index it increments
        wvrr_dict = {}
        # Lists of rr_wrapper objects
        wvrr_list = [] 
        whrr_list = []
        # List for keeping track of ordering of indexes in work_array
        index_order_list = []
        # Some lists for keeping track of indexes incremented by VRR and HRR operations and overall
        vrr_group_index_list = []
        hrr_group_index_list = []
        incremented_index_list = []
        # Temporary lists of w*rr_group objects, which will be reversed when set for algorithm object
        tmp_whrr_group_list = []
        tmp_wvrr_group_list = []
        # Temporary list of contraction and spherical_transform objects
        tmp_early_transform_list = []
        tmp_late_transform_list = []

        # Sort VRR and HRR operations into separate lists
        hrr_found = False
        for wrr in self.wrr_list():
            if wrr.rr().rrtype() == 'vrr':
                assert hrr_found == False, "HRRs cannot precede VRRs"
                wvrr_dict[ wrr.changing_index().name() ] = wrr # maps changing_index to wrr
                wvrr_list.append( wrr )
                vrr_group_index_list.append( wrr.changing_index() )
                incremented_index_list.append( wrr.changing_index() )
            elif wrr.rr().rrtype() == 'hrr':
                assert hrr_found == False, "Currently only a single HRR operation per "+\
                                            "integral class is supported"
                assert wrr.move_from_index() in vrr_group_index_list, \
                        "Angular momentum can only be shifted from indexes incremented using a VRR, "+\
                        "not a HRR."
                hrr_found = True
                whrr_list.append( wrr )
                hrr_group_index_list.append( wrr.changing_index() )
                incremented_index_list.append( wrr.move_to_index() )
        # Initialize index_loop_max_expr_dict with the set of index_value() variables
        # corresponding to the final output integral class and initialize
        # loop_length_expr_dict with corresponding nloop expressions
        # This is the correct state for a simple VRR-only algorithm, without auxiliary indexes
        for windex in self.windex_index_list():
            assert windex.index() in incremented_index_list,\
                    "All integral indexes must be incremented by some type of RR"
            if windex.index().is_cartesian() == True:
                # Cartesian index
                loop_lmax = windex.index_value()
                unrolled_index_lmax = self.generator().generator_options().global_user_lmax()
                nloop = (loop_lmax+1)*(loop_lmax+2)*(loop_lmax+3) / 6
                index_loop_max_expr_dict[ windex.index().name() ] = loop_lmax
                unrolled_index_lmax_dict[ windex.index().name() ] = unrolled_index_lmax
                loop_length_expr_dict[ windex.index().name() ] = nloop
                start_l_expr_dict[ windex.index().name() ] = loop_lmax
                end_l_expr_dict[ windex.index().name() ]   = loop_lmax
            elif windex.index().is_cartesian() == False:
                # non-Cartesian index
                assert windex.index().is_cartesian() == True,\
                        "Only Cartesian integral indexes are currently supported"

        # Determine ordering of primitive exponent and contraction indexes.
        # * For VRR only or VRR + single HRR algorithm, the ordering of these always stays
        #   the same.
        # * Each primitive exponent/contraction index is attached to the corresponding
        #   Cartesian component index, so it forms a compound index ncc * nprim or 
        #   ncc * ncont.
        # * The ordering of indexes throughout is determined by the required ordering for
        #   the last RR operation (for multiple HRRs, using current HRR routine, this would
        #   be more complex, since rearrangement of indexes would be necessary).
        # This data will be stored in each wrr_group and contraction object attached to
        # the algorithm object.
        if hrr_found == True:
            assert len( whrr_list ) == 1, "Only a single HRR per integral class currently "+\
                                          "supported."
            # Ordering determined by single HRR operation
            whrr = whrr_list[0]
            index_order_list = incremented_index_list[:]
            # HRR requires that move_from_index is the fastest changing index and
            # move_to_index is the second fastest changeing index:
            # Rearrange index_order_list to reflect this
            index_order_list.remove( whrr.move_from_index() )
            index_order_list.remove( whrr.move_to_index() )
            index_order_list.append( whrr.move_to_index() )
            index_order_list.append( whrr.move_from_index() )
            # Reorder incremented_index_list to reflect order in index_order_list
            incremented_index_list = index_order_list[:]

        elif hrr_found == False:
            assert len( whrr_list ) == 0, "hrr_found variable disagrees with length of "+\
                                          "whrr_list."
            # Ordering determined by order of VRR operations set by user
            index_order_list = incremented_index_list[:]

        ### HRRs ###
        # Process HRRs before VRRs: the processing of these is independent of whether an auxiliary 
        # index is involved VRRs.
        # For each HRR, add a whrr_group, remove the move_to_index from incremented_index_list.
        # Modify the lists and dictionaries for the final wvrr_group object
        for wrr in reversed( whrr_list ): 
            # Temporary dictionaries created so each whrr_group has a separate instance
            tmp_start_cont_dict = {}
            tmp_end_cont_dict = {}
            tmp_start_min_l_expr_dict = {}
            tmp_end_min_l_expr_dict = {}
            tmp_start_max_l_expr_dict = {}
            tmp_end_max_l_expr_dict = {}
            tmp_start_ncc_type_dict = {}
            tmp_end_ncc_type_dict = {}
            # {start,end}_order_list do not change (single HRR algorithm)
            assert len( whrr_list ) <= 1,\
                    "Currently only 1 HRR operation per integral class is supported."
            tmp_start_order_list = index_order_list[:]
            tmp_end_order_list = index_order_list[:]
            # tmp_{start,end}_cont_dict are initially all set to 'prim', but
            # are modified later if contractions are necessary
            for index in tmp_end_order_list:
                tmp_start_cont_dict[ index.name() ] = 'prim'
                tmp_end_cont_dict[ index.name() ]   = 'prim'
            # Only setup the "max" l_expr dictionaries in this reversed loop. 
            # The "min" l_expr dictionaries can be setup in a subsequent forward loop.
            # The ncc_type dictionaries will also be setup in a subsequent forward loop.
            end_index_list = incremented_index_list[:] 
            # permute indexes in list so that wrr.move_from_index() is the final element in
            # end_index_list
            end_index_list.remove( wrr.move_from_index() )
            end_index_list.append( wrr.move_from_index() )
            # move_to_index has been incremented in HRR
            incremented_index_list.remove( wrr.move_to_index() )
            start_index_list = end_index_list[:]
            start_index_list.remove( wrr.move_to_index() )
            # remove incremented index for set of indexes of start intermediate array
            del start_l_expr_dict[ wrr.move_to_index().name() ]
            max_l_mf_expr = end_l_expr_dict[ wrr.move_from_index().name() ]
            max_l_mt_expr = end_l_expr_dict[ wrr.move_to_index().name() ]
            # modify {start,end}_l_expr_dict objects (will be used by final wvrr_group)
            start_l_expr_dict[ wrr.move_from_index().name() ] = max_l_mf_expr + max_l_mt_expr
            end_l_expr_dict[ wrr.move_from_index().name() ]   = max_l_mf_expr
            end_l_expr_dict[ wrr.move_to_index().name() ]     = max_l_mt_expr
            for index in start_index_list: 
                tmp_start_max_l_expr_dict[ index.name() ] =\
                        expr_deepish_copy( start_l_expr_dict[ index.name() ] )
            for index in end_index_list:
                tmp_end_max_l_expr_dict[ index.name() ] =\
                        expr_deepish_copy( end_l_expr_dict[ index.name() ] )
            # modify index_loop_lmax_expr_dict (will be used by final wvrr_group)
            loop_lmax_mf_expr = index_loop_max_expr_dict[ wrr.move_from_index().name() ] 
            loop_lmax_mt_expr = index_loop_max_expr_dict[ wrr.move_to_index().name() ]
            loop_lmax = loop_lmax_mf_expr + loop_lmax_mt_expr
            index_loop_max_expr_dict[ wrr.move_from_index().name() ] = loop_lmax
            # modify unrolled_index_lmax_dict (unrolled_index_lmax will be set for each vrr_wrapper)
            unrolled_index_lmax_mf = unrolled_index_lmax_dict[ wrr.move_from_index().name() ]
            unrolled_index_lmax_mt = unrolled_index_lmax_dict[ wrr.move_to_index().name() ]
            unrolled_index_lmax = unrolled_index_lmax_mf + unrolled_index_lmax_mt
            unrolled_index_lmax_dict[ wrr.move_from_index().name() ] = unrolled_index_lmax
            # modify the loop_length_expr_dict (will be used by final wvrr_group)
            nloop = (loop_lmax+1)*(loop_lmax+2)*(loop_lmax+3) / 6
            loop_length_expr_dict[ wrr.move_from_index().name() ] = nloop

            # create a whrr_group and add to tmp_whrr_group_list (reverse order list)
            tmp_whrr_group_list.append( integral_wrapper.algorithm.whrr_group( wrr, \
                    tmp_start_order_list, tmp_end_order_list,\
                    start_index_list, end_index_list,\
                    tmp_start_cont_dict, tmp_end_cont_dict,\
                    tmp_start_min_l_expr_dict, tmp_end_min_l_expr_dict,\
                    tmp_start_max_l_expr_dict, tmp_end_max_l_expr_dict,\
                    tmp_start_ncc_type_dict,   tmp_end_ncc_type_dict ) )

            # end_max_l_expr_dict becomes previous start_max_l_expr_dict
            end_l_expr_dict = {}
            for k, v in start_l_expr_dict.items():
                end_l_expr_dict[ k ] = expr_deepish_copy( v )

        # Now run forward through list of whrr_group objects and set {start,end}_min_l_expr_dict
        # and {start,end}_ncc_type_dict values
        tmp_start_min_l_expr_dict = {}
        for k, v in start_l_expr_dict.items():
            tmp_start_min_l_expr_dict[ k ] = 0
        for whrr_group in reversed( tmp_whrr_group_list ):
            wrr = whrr_group.wrr
            # Copy tmp_start_min_l_expr_dict to whrr_group local start_min_l_expr_dict dictionary
            for k, v in tmp_start_min_l_expr_dict.items():
                whrr_group.start_min_l_expr_dict[ k ] = v
            # At the end of a HRR operation, we have reduced the angular momentum range for
            # move_from and move_to index to only the particular angular momentum value of interest
            whrr_group.end_min_l_expr_dict[ wrr.move_from_index().name() ] = \
                    whrr_group.end_max_l_expr_dict[ wrr.move_from_index().name() ]
            whrr_group.end_min_l_expr_dict[ wrr.move_to_index().name() ] = \
                    whrr_group.end_max_l_expr_dict[ wrr.move_to_index().name() ]
            for index in wrr.loop_index_list():
                assert len( tmp_whrr_group_list ) == 1,\
                    "Only one HRR operation per integral class is currently supported."
                # Since only one HRR operation is supported, we can truncate all the
                # remaining integral indexes (already incremented by VRRs) to so the 
                # only Cartesian components that are kept are those of the desired 
                # value of l. 
                # TODO for multiple HRRs, this would need modification, since later 
                # operations would require larger angular momentum ranges in some
                # already-incremented indexes.
                whrr_group.end_min_l_expr_dict[ index.name() ] = \
                        whrr_group.end_max_l_expr_dict[ index.name() ]
            # Determine ncc_type from values in start_min_l_expr_dict and set for whrr_group dictionaries
            for min_l_dict, max_l_dict, ncc_dict in \
                    [ ( whrr_group.start_min_l_expr_dict, whrr_group.start_max_l_expr_dict,\
                        whrr_group.start_ncc_type_dict ), \
                      ( whrr_group.end_min_l_expr_dict, whrr_group.end_max_l_expr_dict,\
                        whrr_group.end_ncc_type_dict ) ]:
                for k in min_l_dict.keys():
                    min_l = min_l_dict[k]
                    max_l = max_l_dict[k]
                    if min_l is max_l:
                        # min and max l values the same, use 'ncc'
                        ncc_dict[ k ] = 'ncc'
                    elif min_l == 0 and not min_l is max_l:
                        # min l = 0, max l != 0, use 'ncc0'
                        ncc_dict[ k ] = 'ncc0'
                    else:
                        raise Exception( "combination of start_l and end_l expressions not understood" )
            # Set tmp_start_min_l_expr_dict equal to end_min_l_expr_dict for whrr_group in this iteration
            # so that the end point of the previous whrr_group is set as the starting point for the next
            tmp_start_min_l_expr_dict = {}
            for k, v in whrr_group.end_min_l_expr_dict.items():
                tmp_start_min_l_expr_dict[k] = v
                        

        # Set the unrolled_index_lmax variables for each vrr_wrapper object now that
        # they have been determined by the HRRs that follow the VRR group.
        for wvrr in wvrr_list:
            wvrr.set_unrolled_index_lmax(\
                    unrolled_index_lmax_dict[ wvrr.unrolled_index().name() ] )

        ### VRRs ###
        if len( self._windex_aux_list ) < 1:
            # If auxiliary index not present in VRRs, use simple algorithm
            # index order is unchanged by a VRR operation
            tmp_start_order_list = index_order_list[:]
            tmp_end_order_list = index_order_list[:]
            # For non-auxiliary index list, we only have one wrr_group,
            # and thus all integral indexes must be incremented in that
            # one group. The remaining indexes in incremented_index_list 
            # are incremented in this group, and the index_loop_max_expr_dict
            # and loop_length_expr_dict initialized above (and possibly 
            # modified for each HRR are used in this wrr_group).
            tmp_wvrr_group_list.append( integral_wrapper.algorithm.wvrr_group( \
                    wvrr_list, index_loop_max_expr_dict, loop_length_expr_dict,\
                    tmp_start_order_list, tmp_end_order_list, \
                    [], incremented_index_list, start_l_expr_dict, end_l_expr_dict, \
                    None, None, None, None ) )
        elif len( self.windex_aux_list() ) > 0:
            # index order is unchanged by a VRR operation
            tmp_start_order_list = index_order_list[:]
            tmp_end_order_list = index_order_list[:]
            group_wvrr_list = wvrr_list
            group_whrr_list = whrr_list
            # If auxiliary index present in VRRs, determine algorithm based on set of 
            # VRRs provided. Gather required information from wrr_list.
            #single_layer = False
            wvrr_aux_changes_set = set()
            # Implementation restriction:
            # * Allow only one m_layers == 1 VRR operation, which must be the final VRR operation
            # ...this is because an implied auxiliary index VRR only gives the correct m = 0 integrals
            # for the desired angular momentum value l, for integrals with angular momentum < l, the
            # auxiliary index is incremented. Currently detection of whether this would be problematic
            # for subsequent VRR operations is not implemented.
            waux = self.windex_aux_list()[0]
            daux = waux.index()
            wvrr = wvrr_list[-1]
            aux_change_dict = {}   # a set of changes of auxiliary indexes for each index
            aux_change_set = set() # a set of changes of auxiliary indexes for ALL indexes
            for index in vrr_group_index_list:
                # Create list of change values for auxiliary index for indexes that are
                # incremented in VRR sequence.
                aux_change_dict[ index.name() ] = set()
                for tup in wvrr.index_change_dict()[index.name()]:
                    aux_change_dict[ index.name() ].add( tup[1] )
                    aux_change_set.add( tup[1] )
            if not 0 in aux_change_set and 1 in aux_change_set and\
                    self.generator().generator_options().contracted() == False:
                # Create info message for user to inform that the implied auxiliary index algorithm is
                # not implemented for primitive integrals
                message_str = "Implied auxiliary index algorithm not implemented for primitive integrals "+\
                              "(see e.g. Ahlrichs, R. Phys. Chem. Chem. Phys. 6, 5119–5121 (2004) ). "+\
                              "You will not be able to take full advantage of implied auxiliary indexes for "+\
                              "integral "+self.full_name()+"."
                              #"integral "+self.integral().name()+"."
                info_str = str(__name__)+", integral_wrapper.setup_algorithm"
                self.generator().add_message( message_str, info_str )

            if not 0 in aux_change_set and 1 in aux_change_set and\
                    self.generator().generator_options().contracted() == True:
                # Note: Currently only available for contracted integral algorithms
                # No integral dependencies have an unincremented auxiliary index.
                # --> Can use a single-layer (implied auxiliary index) algorithm for final wvrr_group
                for c in aux_change_set: assert c < 2, "Auxiliary index max increment is +1."
                m_layers = 1
                changing_index = wvrr.changing_index()
                changing_windex = self.windex_index_dict()[ changing_index.name() ]
                start_m = index_loop_max_expr_dict[ changing_index.name() ]
                end_m = 0
                tmp_group_wvrr_list = group_wvrr_list[:]
                tmp_vrr_group_index_list = vrr_group_index_list[:]
                tmp_vrr_group_index_list.remove( changing_index )
                tmp_group_wvrr_list.remove( wvrr_dict[ changing_index.name() ] )
                # The increment in the auxiliary index must be +1 for all indexes, regardless
                # of how much they are changed in each step
                for index in wvrr.loop_index_list() + [ wvrr.changing_index() ]:
                    for index_decrement, aux_increment in wvrr.index_change_dict()[ index.name() ]:
                        assert aux_increment == 1, "increments in auxiliary index should be +1 for "+\
                                                   "single layer algorithm"
                # The size of the "layer" for each increment of auxiliary index can only be 
                # reduced in a single direction (integral index).
                tmp_index_loop_max_expr_dict = {}
                tmp_loop_length_expr_dict = {}
                loop_max_var = changing_windex.index_loop_max()
                loop_max_expr = start_m - waux.array_index()
                tmp_index_loop_max_expr_dict[ changing_index.name() ] = loop_max_expr
                tmp_loop_length_expr_dict[ changing_index.name() ] =\
                        (loop_max_var+1)*(loop_max_var+2)*(loop_max_var+3)/6
                for index in wrr.loop_index_list():
                    windex = self.windex_index_dict()[ index.name() ] # decremented index
                    loop_max_var  = windex.index_loop_max()
                    loop_max_expr = index_loop_max_expr_dict[ index.name() ] 
                    tmp_index_loop_max_expr_dict[ index.name() ] = loop_max_expr
                    tmp_loop_length_expr_dict[ index.name() ] =\
                            (loop_max_expr+1)*(loop_max_expr+2)*(loop_max_expr+3)/6
                # Set tmp_start_l_expr_dict and tmp_end_l_expr_dict dictionaries
                # (expr values do not change, but number of elements does)
                tmp_start_l_expr_dict = {}
                tmp_end_l_expr_dict = {}
                for start_index in tmp_vrr_group_index_list:
                   tmp_start_l_expr_dict[ start_index.name() ] = \
                             expr_deepish_copy( start_l_expr_dict[ start_index.name() ] )
                for end_index in vrr_group_index_list:
                    tmp_end_l_expr_dict[ end_index.name() ] = \
                           expr_deepish_copy( end_l_expr_dict[ end_index.name() ] )

                tmp_wvrr_group_list.append( integral_wrapper.algorithm.wvrr_group( \
                    [ wvrr_dict[ changing_index.name() ] ], tmp_index_loop_max_expr_dict, tmp_loop_length_expr_dict,\
                    tmp_start_order_list, tmp_end_order_list, \
                    tmp_vrr_group_index_list, vrr_group_index_list, \
                    tmp_start_l_expr_dict, tmp_end_l_expr_dict, \
                    waux.index() , start_m, end_m, m_layers ) )

                vrr_group_index_list = tmp_vrr_group_index_list
                group_wvrr_list = tmp_group_wvrr_list
                # Set m_min for 1-layer implied auxiliary index algorithm
                m_min = start_m
            else:
                # No single layer algorithm possible
                start_m = None
                end_m = None
                # Set m_min value for 2-layer algorithm
                m_min = 0

            # 2-layer algorithm
            m_layers = 2

            # Set up VRR wrr_group objects
            daux = self.windex_aux_list()[0].index()
            while len( vrr_group_index_list ) > 0:
                #print( "group_index_list", [i.name() for i in group_index_list ] )
                #print( "group_wrr_list", [ w.name() for w in group_wrr_list ] )
                tmp_vrr_group_index_list = vrr_group_index_list[:]
                tmp_group_wvrr_list = group_wvrr_list[:]
                smallest_decrement_dict = {}
                non_zero_decrement_dict = {}
                for index in vrr_group_index_list:
                    smallest_decrement_dict[ index.name() ] = None
                for wrr in group_wvrr_list:
                    # Smallest decrement in each index for incremented auxiliary index
                    for index in wrr.loop_index_list() + [ wrr.changing_index() ]:
                        for index_decrement, aux_increment in wrr.index_change_dict()[ index.name() ]:
                            if aux_increment == None: # no auxiliary index present
                                auxin = 0
                            else: 
                                auxin = aux_increment
                            assert auxin <= 1, "increments in auxiliary index greater than 1 "\
                                                      +"not supported"
                            if auxin > 0:
                                if smallest_decrement_dict[ index.name() ] == None:
                                    smallest_decrement_dict[ index.name() ] = index_decrement
                                elif smallest_decrement_dict[ index.name() ] < index_decrement:
                                    smallest_decrement_dict[ index.name() ] = index_decrement
                for index_name, smallest_decrement in list( smallest_decrement_dict.items() ):
                    if smallest_decrement != 0:
                        index = self.windex_index_dict()[ index_name ].index()
                        non_zero_decrement_dict[ index_name ] = smallest_decrement
                        # as the auxiliary index is incremented, this index will go to zero
                        tmp_vrr_group_index_list.remove( index )
                        tmp_group_wvrr_list.remove( wvrr_dict[ index.name() ] )
                assert len(non_zero_decrement_dict.values() ) > 0,\
                        "There must be at least one non-zero decrement in an index where the "+\
                        "auxiliary index is incremented."
                #print( "smallest", smallest_decrement_dict )
                #print( "non-zero", non_zero_decrement_dict )
                if len( non_zero_decrement_dict.values() ) > 1:
                    # The size of the "layer" for each increment of auxiliary index can be reduced 
                    # in multiple directions (integral indexes).
                    raise Exception("The simultaneous decrement of multiple integral indexes for "+\
                                    "incremented auxiliary index is not supported.")
                elif len(non_zero_decrement_dict.values() ) == 1:
                    # The size of the "layer" for each increment of auxiliary index can only be 
                    # reduced in a single direction (integral index).
                    tmp_index_loop_max_expr_dict = {}
                    tmp_loop_length_expr_dict = {}
                    index_name = list( non_zero_decrement_dict.keys() )[0]
                    windex = self.windex_index_dict()[ index_name ] # decremented index
                    if end_m == None:
                        end_m = 0
                    else:
                        end_m = start_m
                    if start_m == None:
                        #start_m = windex.index_value()
                        start_m = index_loop_max_expr_dict[ index_name ]
                    else:
                        #start_m = end_m + windex.index_value() 
                        start_m = end_m + index_loop_max_expr_dict[ index_name ]
                    loop_max_var  = windex.index_loop_max()
                    #loop_max_expr = windex.index_value() - waux.array_index()
                    loop_max_expr = start_m - waux.array_index()
                    tmp_index_loop_max_expr_dict[ index_name ] = loop_max_expr
                    tmp_loop_length_expr_dict[ index_name ] =\
                            (loop_max_var+1)*(loop_max_var+2)*(loop_max_var+3)/6
                    #print( index_loop_max_expr_dict )
                    #print( loop_length_expr_dict )
                    for index in wrr.loop_index_list():
                        windex = self.windex_index_dict()[ index.name() ] # decremented index
                        loop_max_var  = windex.index_loop_max()
                        #loop_max_expr = windex.index_value()
                        loop_max_expr = index_loop_max_expr_dict[ index.name() ] 
                        tmp_index_loop_max_expr_dict[ index.name() ] = loop_max_expr
                        #tmp_loop_length_expr_dict[ index.name() ] =\
                        #        (loop_max_var+1)*(loop_max_var+2)*(loop_max_var+3)/6
                        tmp_loop_length_expr_dict[ index.name() ] =\
                                (loop_max_expr+1)*(loop_max_expr+2)*(loop_max_expr+3)/6
                    # Set tmp_start_l_expr_dict and tmp_end_l_expr_dict dictionaries
                    # (expr values do not change, but number of elements does)
                    tmp_start_l_expr_dict = {}
                    tmp_end_l_expr_dict = {}
                    for start_index in tmp_vrr_group_index_list:
                       tmp_start_l_expr_dict[ start_index.name() ] = \
                                 expr_deepish_copy( start_l_expr_dict[ start_index.name() ] )
                    for end_index in vrr_group_index_list:
                        tmp_end_l_expr_dict[ end_index.name() ] = \
                               expr_deepish_copy( end_l_expr_dict[ end_index.name() ] )
                    tmp_wvrr_group_list.append( integral_wrapper.algorithm.wvrr_group( \
                        group_wvrr_list, tmp_index_loop_max_expr_dict, tmp_loop_length_expr_dict,\
                        tmp_start_order_list, tmp_end_order_list, \
                        tmp_vrr_group_index_list, vrr_group_index_list, \
                        tmp_start_l_expr_dict, tmp_end_l_expr_dict, \
                        waux.index() , start_m, end_m, m_layers ) )
                    vrr_group_index_list = tmp_vrr_group_index_list
                    group_wvrr_list   = tmp_group_wvrr_list
            
            # Set m_max value
            m_max = start_m

        ### Early transform (pre-HRR) and late transform (post-HRR) ###
        if self.generator().generator_options().contracted() == True or\
                self.generator().generator_options().spherical_transformed() == True:
            # Implementation restriction: 
            # * Currently, spherical transformations are only supported for 
            #   contracted integrals
            # * Currently, spherical transformations always follow contractions, 
            #   so we can assume that angular momentum truncation is always done
            #   by a contraction, rather than a spherical transformation.
            assert self.generator().generator_options().contracted() == True, \
                    "Spherical transformation is only supported for contracted integral classes."
            final_wvrr_group = tmp_wvrr_group_list[0]
            # Order of integral indexes is unchanged during transforms
            tmp_start_order_list = index_order_list[:]
            tmp_end_order_list   = index_order_list[:]
            # Integral indexes not incremented during transforms
            tmp_start_index_list = final_wvrr_group.end_index_list[:]
            tmp_end_index_list   = final_wvrr_group.end_index_list[:] 
            # At the start of the early transform sequence, all indexes are
            # uncontracted (primitive)
            tmp_start_cont_dict  = {}
            tmp_end_cont_dict    = {}
            for index in tmp_end_order_list:
                tmp_start_cont_dict[ index.name() ] = 'prim'
                tmp_end_cont_dict[ index.name() ]   = 'prim'
            # In the first contraction, we can convert all ncc0 indexes that
            # will not be move_from, or move_to indexes in subsequent HRR operations
            # to ncc (i.e. remove unecessary angular momentum components).
            tmp_start_min_l_expr_dict = {}
            tmp_end_min_l_expr_dict   = {}
            tmp_start_max_l_expr_dict = {}
            tmp_end_max_l_expr_dict   = {}
            tmp_start_ncc_type_dict   = {}
            tmp_end_ncc_type_dict     = {}
            # We assume all indexes will undergo contraction
            index_to_contract_list = list( reversed( index_order_list[:] ) )
            # All indexes which have been incremented by VRRs and do not feature in any
            # subsequent HRR operations can be spherically transformed in the
            # early_transform sequence
            if self.generator().generator_options().spherical_transformed() == True:
                # We assume that all indexes will (eventually) undergo spherical
                # transformation
                index_to_sph_trans_list = list( reversed( index_order_list[:] ) )
            else:
                index_to_sph_trans_list = []
            early_sph_trans_list = index_to_sph_trans_list[:]
            late_sph_trans_list  = index_to_sph_trans_list[:]
            for whrr_group in tmp_whrr_group_list:
                move_from_index = whrr_group.wrr.move_from_index()
                move_to_index = whrr_group.wrr.move_to_index()
                if move_from_index in early_sph_trans_list:
                    early_sph_trans_list.remove( move_from_index )
                if move_to_index in early_sph_trans_list:
                    early_sph_trans_list.remove( move_to_index )
            for index in early_sph_trans_list:
                if index in late_sph_trans_list:
                    late_sph_trans_list.remove( index )
            assert len( index_to_contract_list ) > 0,\
                    "For contracted integrals, at least one index should be contracted."
            # First contraction is special -- discards unneeded integrals
            # First contraction starts with state of final VRR
            for con_d, vrr_d in \
                    [ (tmp_start_min_l_expr_dict, final_wvrr_group.end_min_l_expr_dict ),\
                    (tmp_start_max_l_expr_dict, final_wvrr_group.end_max_l_expr_dict ) ]:
                for index in tmp_start_index_list:
                    con_d[ index.name() ] = expr_deepish_copy( vrr_d[index.name()] )
            for index in tmp_start_index_list:
                tmp_start_ncc_type_dict[ index.name() ] =\
                        final_wvrr_group.end_ncc_type_dict[ index.name() ]

            # First contraction truncates the cartesian component range of all indexes
            # that are not going to be move_from or move_to indexes in subsequent HRRs
            # It also truncates the cartesian component range of indexes from which
            # angular momentum is moved from in HRRs, e.g. is moving from index b to a
            # end_max_l_expr for b should be la + lb, and end_min_l_expr is 0 for 
            # the preceding VRR phase. The first ncc0(lb-1) integrals can be discarded,
            # since we only care about integrals with angular momentum lb in the end
            # result.
            move_from_set = set()
            move_to_set = set()
            for whrr_group in tmp_whrr_group_list:
                move_from_set.add( whrr_group.wrr.move_from_index() )
                move_to_set.add( whrr_group.wrr.move_to_index() )
            # Indexes to truncate so end_min_l_expr == end_max_l_expr
            index_to_truncate_list = [ i for i in index_to_contract_list \
                                       if not i in move_to_set and not i in move_from_set ]
            # Indexes to truncate where the first ncc0(lX-1) Cartesian components can
            # be discarded (i.e. indexes that will be move_from_index in subsequent HRRs).
            mf_truncate_list  = [ i for i in move_from_set ]
            for index in tmp_start_index_list:
                if index in index_to_truncate_list:
                   # Indexes to be truncated have end_min_l_expr == end_max_l_expr, ncc_type == 'ncc'
                    tmp_end_min_l_expr_dict[ index.name() ] =\
                            expr_deepish_copy( final_wvrr_group.end_max_l_expr_dict[ index.name() ] )
                    tmp_end_max_l_expr_dict[ index.name() ] =\
                            expr_deepish_copy( final_wvrr_group.end_max_l_expr_dict[ index.name() ] )
                    tmp_end_ncc_type_dict[ index.name() ] = 'ncc' 

                elif index in mf_truncate_list:
                    # Indexes to be truncated by removing first ncc0(lX-1) CCs.
                    windex = self.windex_index_dict()[ index.name() ]
                    tmp_end_min_l_expr_dict[ index.name() ] = windex.index_value()
                    tmp_end_max_l_expr_dict[ index.name() ] =\
                            expr_deepish_copy( final_wvrr_group.end_max_l_expr_dict[ index.name() ] )
                    tmp_end_ncc_type_dict[ index.name() ] = 'ncc0-ncc0'

                else:
                    # Indexes that cannot be truncated, remain unchanged
                    tmp_end_min_l_expr_dict[ index.name() ] =\
                            expr_deepish_copy( final_wvrr_group.end_min_l_expr_dict[ index.name() ] )
                    tmp_end_max_l_expr_dict[ index.name() ] =\
                            expr_deepish_copy( final_wvrr_group.end_max_l_expr_dict[ index.name() ] )
                    tmp_end_ncc_type_dict[ index.name() ] =\
                            final_wvrr_group.end_ncc_type_dict[ index.name() ]

            # First transform is a contraction (all contractions always precede
            # spherical transformations) which truncates the angular momentum 
            # range before HRR sequence
            index = index_to_contract_list[0]
            # Modify tmp_end_index_dict so that index is shown as contracted
            tmp_end_cont_dict[ index.name() ] = 'cont'
            # Setup loop_index_list based on order of indexes in work array
            pos = tmp_end_order_list.index( index )
            loop_index_list = tmp_end_order_list[ :pos ]
            # Create contraction object and append to tmp_early_transform_list for first index to contract
            tmp_early_transform_list.append(\
                    integral_wrapper.algorithm.contraction(index,loop_index_list,\
                                    tmp_start_order_list, tmp_end_order_list,\
                                    tmp_start_index_list, tmp_end_index_list,\
                                    tmp_start_cont_dict,  tmp_end_cont_dict,\
                                    tmp_start_min_l_expr_dict, tmp_end_min_l_expr_dict,\
                                    tmp_start_max_l_expr_dict, tmp_end_max_l_expr_dict,\
                                    tmp_start_ncc_type_dict, tmp_end_ncc_type_dict) )

            # Propagate contraction status through subsequent whrr_group operations
            for whrr_group in tmp_whrr_group_list:
                whrr_group.start_cont_dict[ index.name() ] = 'cont' 
                whrr_group.end_cont_dict[ index.name() ]   = 'cont' 

            # Modify whrr_group ncc_type, min_l_expr and max_l_expr values for 
            # indexes that are truncated in the first contraction.
            # This is necessary because the setup procedure for whrr_groups assumes 
            # no contractions, so the first whrr_group truncates from ncc0 --> ncc.
            for whrr_group in tmp_whrr_group_list:
                for index in index_to_truncate_list:
                    whrr_group.start_min_l_expr_dict[ index.name() ] = \
                            expr_deepish_copy(tmp_end_min_l_expr_dict[ index.name() ] )
                    whrr_group.start_max_l_expr_dict[ index.name() ] = \
                            expr_deepish_copy(tmp_end_max_l_expr_dict[ index.name() ] )
                    whrr_group.end_min_l_expr_dict[ index.name() ] = \
                            expr_deepish_copy(tmp_end_min_l_expr_dict[ index.name() ] )
                    whrr_group.end_max_l_expr_dict[ index.name() ] = \
                            expr_deepish_copy(tmp_end_max_l_expr_dict[ index.name() ] )
                    whrr_group.start_ncc_type_dict[ index.name() ] = \
                                    tmp_end_ncc_type_dict[ index.name() ]
                    whrr_group.end_ncc_type_dict[ index.name() ] = \
                                    tmp_end_ncc_type_dict[ index.name() ]

            # For the whrr_groups preceding the whrr_group where an index in mf_truncate_list
            # is the move_from_index, modify start_ and end_ values of l_expr and ncc_type
            # For the whrr_group where index is move_from_index, only modify start_ values.
            for index in mf_truncate_list:
                for whrr_group in reversed( tmp_whrr_group_list ):
                    whrr_group.start_min_l_expr_dict[ index.name() ] = \
                            expr_deepish_copy( tmp_end_min_l_expr_dict[ index.name() ] )
                    whrr_group.start_max_l_expr_dict[ index.name() ] = \
                            expr_deepish_copy( tmp_end_max_l_expr_dict[ index.name() ] )
                    whrr_group.start_ncc_type_dict[index.name() ] = \
                            tmp_end_ncc_type_dict[ index.name() ]
                    if index is whrr_group.wrr.move_from_index():
                        break
                    else:
                        whrr_group.end_min_l_expr_dict[ index.name() ] = \
                                expr_deepish_copy( tmp_end_min_l_expr_dict[ index.name() ] )
                        whrr_group.end_max_l_expr_dict[ index.name() ] = \
                                expr_deepish_copy( tmp_end_max_l_expr_dict[ index.name() ] )
                        whrr_group.end_ncc_type_dict[index.name() ] = \
                                tmp_end_ncc_type_dict[ index.name() ]


            # Remaining contractions (preceding any spherical transformations)
            # Create per-loop cont_dict
            loop_cont_dict = {}
            for index in tmp_end_order_list:
                loop_cont_dict[ index.name() ] = tmp_end_cont_dict[ index.name() ]
            i = 1
            if len( index_to_contract_list ) > 1:
                for index in index_to_contract_list[1:]:
                    i += 1
                    # All subsequent contractions have the same attributes as the first contraction
                    # Make copies of attributes for first contraction (this avoids potential
                    # issues with the same object being an attribute to multiple operation
                    # objects).
                    tmp1_start_order_list = tmp_start_order_list[:]
                    tmp1_end_order_list   = tmp_end_order_list[:]
                    tmp1_start_index_list = tmp_start_index_list[:]
                    tmp1_end_index_list   = tmp_end_index_list[:]
                    tmp1_start_cont_dict  = {}
                    tmp1_end_cont_dict    = {}
                    tmp1_start_min_l_expr_dict = {}
                    tmp1_end_min_l_expr_dict   = {}
                    tmp1_start_max_l_expr_dict = {}
                    tmp1_end_max_l_expr_dict   = {}
                    tmp1_end_max_l_expr_dict   = {}
                    tmp1_start_ncc_type_dict   = {}
                    tmp1_end_ncc_type_dict     = {}
                    # The l_expr and ncc_type state at the end of the first contraction is
                    # unchanged for subsequent contractions
                    for tmp1_d, tmp_d in \
                                [ ( tmp1_start_cont_dict, loop_cont_dict ),\
                                  ( tmp1_end_cont_dict, loop_cont_dict ),\
                                  ( tmp1_start_min_l_expr_dict, tmp_end_min_l_expr_dict ),\
                                  ( tmp1_end_min_l_expr_dict, tmp_end_min_l_expr_dict ),\
                                  ( tmp1_start_max_l_expr_dict, tmp_end_max_l_expr_dict ),\
                                  ( tmp1_end_max_l_expr_dict, tmp_end_max_l_expr_dict ),\
                                  ( tmp1_start_ncc_type_dict, tmp_end_ncc_type_dict ),\
                                  ( tmp1_end_ncc_type_dict, tmp_end_ncc_type_dict ) ] :
                        for k, v in tmp_d.items():
                            tmp1_d[ k ] = v
                    # Update loop_cont_dict and end_cont_dict to reflect that index is
                    # contracted in this contraction operation
                    loop_cont_dict[ index.name() ] = 'cont'
                    tmp1_end_cont_dict[ index.name() ] = 'cont'
                    
                    # Propagate contraction status through subsequent whrr_group operations
                    for whrr_group in tmp_whrr_group_list:
                        whrr_group.start_cont_dict[ index.name() ] = 'cont' 
                        whrr_group.end_cont_dict[ index.name() ]   = 'cont' 
                    # Create loop_index_list from indexes which come "after" index being contracted
                    pos = tmp_end_order_list.index( index )
                    loop_index_list = tmp_end_order_list[ :pos ]
                    # Add contraction object to list of contractions
                    tmp_early_transform_list.append(\
                            integral_wrapper.algorithm.contraction(index,loop_index_list,\
                                        tmp1_start_order_list, tmp1_end_order_list,\
                                        tmp1_start_index_list, tmp1_end_index_list,\
                                        tmp1_start_cont_dict,  tmp1_end_cont_dict,\
                                        tmp1_start_min_l_expr_dict, tmp1_end_min_l_expr_dict,\
                                        tmp1_start_max_l_expr_dict, tmp1_end_max_l_expr_dict,\
                                        tmp1_start_ncc_type_dict, tmp1_end_ncc_type_dict) )
            else:
                # Only one index to contract -- need to initialize
                # tmp1_end_ncc_type_dict for any subsequent spherical transformations
                tmp1_end_ncc_type_dict = {}
                for k, v in tmp_end_ncc_type_dict.items():
                    tmp1_end_ncc_type_dict[ k ] = v
            
            # Constrain spherical_transformations to always follow contractions 
            # (for simplicity)
            loop_ncc_type_dict = {}
            for k, v in tmp1_end_ncc_type_dict.items():
                loop_ncc_type_dict[ k ] = v

            # At end of sequence of transformations, loop_cont_dict contains the 
            # state of the final transformation (should be all contracted)
            for k, v in loop_cont_dict.items():
                assert v == 'cont', "After final contraction, all indexes should be contracted."
            if self.generator().generator_options().spherical_transformed() == True:
                i = 0
                for index in early_sph_trans_list:
                    i += 1
                    # Make copies of attributes for contractions (this avoids potential
                    # issues with the same object being an attribute to multiple operation
                    # objects).
                    tmp1_start_order_list = tmp_start_order_list[:]
                    tmp1_end_order_list   = tmp_end_order_list[:]
                    tmp1_start_index_list = tmp_start_index_list[:]
                    tmp1_end_index_list   = tmp_end_index_list[:]
                    tmp1_start_cont_dict = {}
                    tmp1_end_cont_dict = {}
                    tmp1_start_min_l_expr_dict = {}
                    tmp1_end_min_l_expr_dict   = {}
                    tmp1_start_max_l_expr_dict = {}
                    tmp1_end_max_l_expr_dict   = {}
                    tmp1_start_ncc_type_dict   = {}
                    tmp1_end_ncc_type_dict     = {}
                    # The l_expr at the end of the spherical transformation sequence is
                    # unchanged for subsequent contractions
                    for tmp1_d, tmp_d in \
                                [ ( tmp1_start_cont_dict, loop_cont_dict ),\
                                  ( tmp1_end_cont_dict, loop_cont_dict ),\
                                  ( tmp1_start_min_l_expr_dict, tmp_end_min_l_expr_dict ),\
                                  ( tmp1_end_min_l_expr_dict, tmp_end_min_l_expr_dict ),\
                                  ( tmp1_start_max_l_expr_dict, tmp_end_max_l_expr_dict ),\
                                  ( tmp1_end_max_l_expr_dict, tmp_end_max_l_expr_dict ), \
                                  ( tmp1_start_ncc_type_dict, loop_ncc_type_dict ) ] :
                        for k, v in tmp_d.items():
                            tmp1_d[ k ] = v
                    # sph_trans operations modify the ncc type
                    loop_ncc_type_dict[ index.name() ] = 'nsph'
                    for k, v in loop_ncc_type_dict.items():
                        tmp1_end_ncc_type_dict[ k ] = v
                    # Create loop_index_list from indexes which come "after" index being contracted
                    pos = tmp1_end_order_list.index( index )
                    loop_index_list = tmp1_end_order_list[ :pos ]
                    # Add contraction object to list of transforms
                    tmp_early_transform_list.append(\
                            integral_wrapper.algorithm.spherical_transform(index,loop_index_list,\
                                        tmp1_start_order_list, tmp1_end_order_list,\
                                        tmp1_start_index_list, tmp1_end_index_list,\
                                        tmp1_start_cont_dict,  tmp1_end_cont_dict,\
                                        tmp1_start_min_l_expr_dict, tmp1_end_min_l_expr_dict,\
                                        tmp1_start_max_l_expr_dict, tmp1_end_max_l_expr_dict,\
                                        tmp1_start_ncc_type_dict, tmp1_end_ncc_type_dict) )
                    # Propagate ncc_type change through subsequent HRR operations
                    for whrr_group in tmp_whrr_group_list:
                        whrr_group.start_ncc_type_dict[ index.name() ] = 'nsph'
                        whrr_group.end_ncc_type_dict[ index.name() ] = 'nsph'
    
                if len( late_sph_trans_list ) > 0:
                    assert len( tmp_whrr_group_list ) > 0, \
                            "late_sph_trans_list should be empty if no HRR operations are required."
                    # Late transformations are those that follow one or more HRR operations.
                    # Currently, these will always be spherical transformations, as
                    # contractions can always precede a HRR operation.
                    # Once a HRR operation has moved angular momentum from the move_from_index
                    # to the move_to_index, then both can be spherical transfomed, as they are no
                    # longer involved in any further Cartesian RRs.
                    ihrr = 0
                    for whrr_group in reversed( tmp_whrr_group_list ): 
                        move_from_index = whrr_group.wrr.move_from_index()
                        move_to_index   = whrr_group.wrr.move_to_index()
                        assert move_from_index in late_sph_trans_list, \
                                "move_from_index for all HRR operations must be in late_sph_trans_list"
                        assert move_to_index in late_sph_trans_list, \
                                "move_to_index for all HRR operations must be in late_sph_trans_list"
                        # Make copies of attributes for HRR operation (this avoids potential
                        # issues with the same object being an attribute to multiple operation
                        # objects). These attributes are not modified by spherical transformation.
                        tmp1_start_order_list = whrr_group.start_order_list[:]
                        tmp1_end_order_list   = whrr_group.end_order_list[:]
                        tmp1_start_index_list = whrr_group.start_index_list[:]
                        tmp1_end_index_list   = whrr_group.end_index_list[:]
                        if late_sph_trans_list.index( move_from_index ) < late_sph_trans_list.index( move_to_index ):
                            per_hrr_group_index_list = [ move_from_index, move_to_index ]
                        else:
                            per_hrr_group_index_list = [ move_to_index, move_from_index ]
                        # ncc_type changes after spherical transformation (ncc --> nsph), so copy these
                        # to temporary dictionary objects that are modified inside each loop over
                        # index to transform
                        tmp_end_ncc_type_dict = {}
                        for k, v in whrr_group.end_ncc_type_dict.items():
                            tmp_end_ncc_type_dict[ k ] = v
                        for index in per_hrr_group_index_list:
                            # The l_expr at the end of the spherical transformation sequence is
                            # unchanged for subsequent contractions
                            tmp1_start_min_l_expr_dict = {}
                            tmp1_end_min_l_expr_dict   = {}
                            tmp1_start_max_l_expr_dict = {}
                            tmp1_end_max_l_expr_dict   = {}
                            tmp1_start_ncc_type_dict   = {}
                            tmp1_end_ncc_type_dict     = {}
                            tmp1_start_cont_dict       = {}
                            tmp1_end_cont_dict         = {}
                            for tmp1_d, tmp_d in \
                                        [ ( tmp1_start_cont_dict, whrr_group.end_cont_dict ),\
                                          ( tmp1_end_cont_dict, whrr_group.end_cont_dict ),\
                                          ( tmp1_start_min_l_expr_dict, whrr_group.end_min_l_expr_dict ),\
                                          ( tmp1_end_min_l_expr_dict, whrr_group.end_min_l_expr_dict ),\
                                          ( tmp1_start_max_l_expr_dict, whrr_group.end_max_l_expr_dict ),\
                                          ( tmp1_end_max_l_expr_dict, whrr_group.end_max_l_expr_dict ), \
                                          ( tmp1_start_ncc_type_dict, tmp_end_ncc_type_dict ) ] :
                                for k, v in tmp_d.items():
                                    tmp1_d[ k ] = v
                            # sph_trans operations modify the ncc type
                            tmp_end_ncc_type_dict[ index.name() ] = 'nsph'
                            for k, v in tmp_end_ncc_type_dict.items():
                                tmp1_end_ncc_type_dict[ k ] = v
                            # Create loop_index_list from indexes yet to be spherical transformed
                            pos = tmp1_end_order_list.index( index )
                            loop_index_list = tmp1_end_order_list[:pos]
                            # Add operation object to list of contractions
                            tmp_late_transform_list.append(\
                                    integral_wrapper.algorithm.spherical_transform(index,loop_index_list,\
                                                tmp1_start_order_list, tmp1_end_order_list,\
                                                tmp1_start_index_list, tmp1_end_index_list,\
                                                tmp1_start_cont_dict,  tmp1_end_cont_dict,\
                                                tmp1_start_min_l_expr_dict, tmp1_end_min_l_expr_dict,\
                                                tmp1_start_max_l_expr_dict, tmp1_end_max_l_expr_dict,\
                                                tmp1_start_ncc_type_dict, tmp1_end_ncc_type_dict) )
                            # Propagate ncc_type change through subsequent whrr_group objects
                            for whrr_group in list( reversed( tmp_whrr_group_list ) )[ihrr+1:]:
                                whrr_group.start_ncc_type_dict[ index.name() ] = 'nsph'
                                whrr_group.end_ncc_type_dict[ index.name() ] = 'nsph'
                        ihrr += 1
    
        ### algo_module selector ###
        # Determines the algorithm used in outputting the main function (called when evaluating
        # a particular integral class
        # Set boolean variables
        is_contracted     = False
        hrr_present       = False
        aux_index_present = False
        if self.generator().generator_options().contracted() == True:
            is_contracted = True
        if len( tmp_whrr_group_list ) > 0:
            hrr_present = True
        if len( self.windex_aux_list() ) >0:
            aux_index_present = True
        # Create algorithm_descriptor object specific to this integral class
        integral_algo_desc = algo_descriptor( is_contracted = is_contracted, \
                                          hrr_present = hrr_present, \
                                          aux_index_present = aux_index_present )
        # Import dictionary of known algorithm_descriptor objects
        algo_desc_dict = main_source_algorithms.algo_desc_dict
        # Use hash of algo_descriptor to identify correct algorithm module
        try:
            algo_module    = algo_desc_dict[ integral_algo_desc ] 
        except KeyError:
            # If the particular algo_descriptor is not matched, then throw an exception
            print( "" )
            print( "ERROR:" )
            print( "The algorithm with the following algo_descriptor unsupported:")
            print( integral_algo_desc.debug_str() )
            print( "" )
            raise

        ### Setup algorithm object ###
        wvrr_group_list  = list( reversed( tmp_wvrr_group_list ) )
        early_transform_list = tmp_early_transform_list[:] 
        whrr_group_list  = list( reversed( tmp_whrr_group_list ) )
        late_transform_list  = tmp_late_transform_list[:]
        self._algorithm_obj = integral_wrapper.algorithm( wvrr_group_list, early_transform_list, \
                whrr_group_list, late_transform_list, m_min, m_max, m_layers, algo_module )
        return

    def setup_integral_arrays(self):
        """
        Create an integral_array object representing the work (intermediate) and
        output arrays which store integral values.
        """
        ### Work array ###
        # Create associated DSL objects
        work_array = dsl_pointer(name = 'work_array', vartype = 'double' )
        work_array_index = dsl_scalar( name = 'iwrk0', vartype = 'int' )

        # Generate DSL expression for array size
        pass

        # Create integral_array object
        self._work_array = integral_array( work_array, work_array_index )

        ### Output array ### 
        # Create associated DSL objects
        output_array = dsl_pointer(name = 'output_array', vartype = 'double' )
        output_array_index = dsl_scalar( name = 'iopt', vartype = 'int' )

        # Generate DSL expression for array size
        pass 

        # Create integral_array object
        self._output_array = integral_array( output_array, output_array_index )

    def setup_transform_arrays(self):
        """
        Create a general_array object representing the array of contraction coefficients
        (if necessary) and set self._contract_array accordingly.
        """
        ### DEPRECATED ###
        # contract_array and sph_trans_array objects are associated with individual index_wrapper
        # objects, rather than for the entire integral_wrapper object.
        ### 
        pass

    def setup_integral_array_indexing(self):
        """
        Initializes variables related to indexing based on the integral indexes and auxiliary
        indexes by setting expressions for these.
        """
        algorithm_obj = self.algorithm_obj()

        # Create instance of integral_array_indexing_setup class for selected algorithm module
        array_indexing_setup_obj = algorithm_obj.algo_module.integral_array_indexing_setup(\
                                    self, self.generator() )
        # Call algorithm-specific routine for array indexing setup (sets up work_array integral_array
        # object).
        array_indexing_setup_obj()
        
        # Set work_array().array_index() initial value (zero)
        #self.work_array().array_index().set_expr( 0 )

        ## Go through RR list to determine order of indexes in work array (for each auxiliary index
        ## combination)
        #work_array_windex_index_list = []
        #vrr_windex_index_list = []
        #hrr_windex_index_list = []
        ## (see definition in vrr_wrapper class)
        #work_array_windex_aux_list = []
        #wvrr_list = []
        #whrr_list = []
        ## Separate vrr_wrapper and hrr_wrapper objects
        #hrr_found = False
        #for wrr in self.wrr_list():
        #    assert wrr.rr().rrtype() in [ 'vrr', 'hrr' ], \
        #        "Only VRR and HRR-type RRs supported currently"
        #    if wrr.rr().rrtype() == 'vrr':
        #        assert hrr_found == False, "VRRs may not follow HRRs in algorithm"
        #        wvrr_list.append(wrr)
        #        windex = self.windex_index_dict()[ wrr.changing_index().name() ]
        #        work_array_windex_index_list.append( windex )
        #        vrr_windex_index_list.append( windex )
        #    elif wrr.rr().rrtype() == 'hrr':
        #        hrr_found = True
        #        whrr_list.append(wrr)
        #        windex = self.windex_index_dict()[ wrr.move_to_index().name() ]
        #        work_array_windex_index_list.append( windex )
        #        hrr_windex_index_list.append( windex )

        ## Any additional integral indexes not captured from wrr_list
        #for windex in self.windex_index_list():
        #    assert windex in work_array_windex_index_list,\
        #            "Currently, all integral indexes must be incremented in an RR"
        #    #if not windex in self.windex_index_list():
        #    #    work_array_windex_index_list.append( windex )
        ## * Get (unordered) list of auxiliary indexes and add to list
        ## * Create new array_index_dict entries for each change value of 
        ##   auxiliary index
        #all_wvrr_aux_change_dict = {}
        #for windex in self.windex_aux_list():
        #    work_array_windex_aux_list.append( windex )
        #    all_wvrr_aux_change_dict[ windex.index().name() ] = set()
        #for wvrr in wvrr_list:
        #    for windex in work_array_windex_aux_list:
        #        for change in wvrr.aux_change_dict()[ windex.index().name() ]:
        #            all_wvrr_aux_change_dict[ windex.index().name() ].add( change )
        ## Add array_index for placement of base function output in work_array
        ## The following code assumes only a single auxiliary index
        #work_array_index_base = \
        #                dsl_scalar( name = 'iwrkb', vartype = 'int' )
        #self.work_array().add_special_array_index('base',work_array_index_base )
        #assert len( work_array_windex_aux_list ) <= 1,\
        #        "Only 0 or 1 auxiliary indexes per integral class supported"
        #if len(work_array_windex_aux_list) > 0:
        #    aux_windex = work_array_windex_aux_list[0]
        #    aux_change_list = list( all_wvrr_aux_change_dict[ aux_windex.index().name() ] )
        #    # change = 0 automatically set in array_index_dict for integral_array objects
        #    if 0 in aux_change_list: aux_change_list.remove(0)
        #    for change in aux_change_list:
        #        work_array_index = \
        #                dsl_scalar( name = 'iwrk'+str(change), vartype = 'int' )
        #        work_array_index_start = \
        #                dsl_scalar( name = 'iwrk'+str(change)+'_start', vartype = 'int' )
        #        self.work_array().add_array_index( change, work_array_index,\
        #                work_array_index_start )
        ## Add special_array_index for HRR "layer" or contraction "layer", if required
        #hrr_present = False
        #contract_present = False
        #if len( self.algorithm_obj().whrr_group_list ) > 0: 
        #    hrr_present = True
        #if self.generator().generator_options().contracted() == True:
        #    contract_present = True
        #if hrr_present == True or contract_present == True:
        #    # if a iwrk1 is already present due to auxiliary index algorithm, use that,
        #    # else create a new iwrk1 array_index tuple
        #    if not 1 in list( self.work_array().array_index_dict().keys() ):
        #        # if iwrk1 not present, then create iwrk1
        #        work_array_index = \
        #                    dsl_scalar( name = 'iwrk1', vartype = 'int' )
        #        work_array_index_start = \
        #                dsl_scalar( name = 'iwrk1_start', vartype = 'int' )
        #        self.work_array().add_array_index( 1, work_array_index,\
        #                work_array_index_start )

        ##print( self.work_array().array_index_dict() )

        ## Set work_array windex_list to integral indexes (ordered) +
        ## auxiliary indexes (unordered)
        #self.work_array().set_windex_list( work_array_windex_index_list +\
        #                                    work_array_windex_aux_list )

    def setup_integral_argument_list(self):
        """
        Create a list of dsl_variable objects which represent the arguments
        for the main integral evaluation function. These will be based on the indexes
        of the integral class and the objects used to pass data in and out of the
        function.
        """
        algorithm_obj = self.algorithm_obj()

        arg_list = []
        
        # Create instance of integral_argument_setup class for selected algorithm module
        arg_setup_obj = algorithm_obj.algo_module.integral_argument_setup( self, self.generator() )

        # Setup integral argument list
        arg_setup_obj( arg_list )
        
        # Set integral_wrapper attribute
        self._integral_argument_list = arg_list

    def setup_local_variable_lists(self):
        """
        Create a list of src_dsl_variable objects ordered based on dependency.
        If the va.expr() contains on vb, then va will be after vb in the list.
        The list is set up so that all such interdependencies are satisfied.
        """

        algorithm_obj = self.algorithm_obj()

        ### VARIABLES FOR MANUAL ASSIGNMENT ###
        manual_list = []
        
        # Create instance of manual_local_variable_setup class for selected algorithm module
        manual_setup_obj = \
                algorithm_obj.algo_module.manual_local_variable_setup( self, self.generator() )

        # Setup manual_set
        manual_setup_obj( manual_list )

        # Set integral_wrapper attribute
        self._manual_local_variable_list = manual_list

        ### VARIABLES FOR AUTOMATIC ASSIGMENT ###
        ordered_list_list = []
        ordered_dict_list = []

        # Create instance of manual_local_variable_setup class for selected algorithm module
        automatic_setup_obj = \
                algorithm_obj.algo_module.automatic_local_variable_setup( self, self.generator() )

        # Setup ordered_list 
        automatic_setup_obj( ordered_list_list, ordered_dict_list )

        # Setup integral_wrapper attribute
        self._auto_local_variable_list     = ordered_list_list[0]
        self._auto_prim_loop_variable_list = ordered_list_list[1]
        # ordered_dict_list[1] contains a dictionary with dsl_variable 
        # objects from ordered_list_list[1] as keys and sets of 
        # per-primitive dsl_variable objects as values, to indicate that
        # the dsl_variable key depends on the set of dsl_variable objects as
        # keys -- this allows automatic assignment of variables inside of
        # loops over primitive indexes
        # ordered_dict_list[0] is empty, as we do not need to relate these
        # dsl_variables to any other object for code output.
        self._auto_prim_loop_depends_dict  = ordered_dict_list[1]

    def setup_main_function(self):
        """
        Create a function_wrapper object corresponding to the main function for 
        generating the integal class in self.integral() (the main function calls
        the supporting functions (function_wrappers for these are in 
        self._support_function_list) to generate the integral class.
        """
        # Short names for callable_* derived classes
        main_call      = integral_wrapper_main_function.main_call
        main_prototype = integral_wrapper_main_function.main_prototype
        # The main_source class is a derived class in callables.main_source_algorithms
        # which is determined in setup_algorithm
        main_source    = self.algorithm_obj().algo_module.algo_main_source

        darg_list = self.integral_argument_list()

        # Setup a const_dict object, to assert that all variables passed in should be
        # const, except self.work_array().pointer() and self.output_array.pointer() 
        # -- this will be expressed in the function prototyp
        const_dict = {}
        non_const_list = [ self.work_array().pointer(), self.output_array().pointer() ]
        for arg in darg_list:
            if not arg in non_const_list :
                const_dict[arg] = True
            else:
                const_dict[arg] = False

        #funcname = self.integral().name()+'_main'
        funcname = self.full_name()+'_main'
        f_dsl = dsl_function(funcname,darg_list,vartype="void")
        f_call = main_call( f_dsl, gen = self.generator() )
        f_prototype = main_prototype( f_dsl, gen = self.generator(), const_dict = const_dict )
        f_source = main_source( f_dsl, wintegral = self, gen = self.generator(), const_dict = const_dict )
        wf = function_wrapper(f_dsl,f_call,f_prototype,f_source, gen = self.generator() )
        self._main_function = wf

    def setup_rr_functions(self):
        """
        For each dsl_rr object in self.integral(), create a function_wrapper
        object that provides methods for generating function prototypes and
        executable source code related to that RR.
        """
        # Short names for callable_* derived classes
        rr_call      = integral_wrapper_rr_function.rr_call
        rr_prototype = integral_wrapper_rr_function.rr_prototype
        vrr_source   = integral_wrapper_rr_function.vrr_source

        # Callable classes containing code necessary to generate prototype and source
        # specific to RR functions

        for wrr in self.wrr_list():
            assert wrr.rr().rrtype() in ['vrr','hrr'], "Only VRR and HRR type RRs supported currently"
            if wrr.rr().rrtype() == 'vrr':
                # Check for single layer algorithm
                # TODO improve inefficient check of algorithm_obj for information on wrr
                single_layer = False
                for wvrr_group in self.algorithm_obj().wvrr_group_list:
                    if wvrr_group.m_layers == 1:
                        assert len( wvrr_group.wrr_list ) == 1, \
                                "Only a single VRR per wvrr_group with m_layers == 1 supported."
                        assert wvrr_group is self.algorithm_obj().wvrr_group_list[-1],\
                                "Only a single wvrr_group with m_layers == 1 is supported and "+\
                                "this must be the final wvrr_group operation."
                        if wvrr_group.wrr_list[0] is wrr:
                            # Current wrr is in m_layers == 1 wvrr_group
                            single_layer = True
                            break
                # If a single layer algorithm modify wrr expression so that integrals
                # appear to have an unchanged auxiliary index.
                simplified_expr = expr_deepish_copy( wrr.simplified_expr() )
                # Get list of dsl_integral objects in simplified_expr
                dsl_integral_list = expr_list_objects( simplified_expr, dsl_integral )
                if single_layer == True:
                    # Modify dsl_integral objects to remove increment in auxiliary index
                    # but only change the copied expression, rather than the original
                    # wrr.simplified_expr().
                    for dintegral in dsl_integral_list:
                        assert len( dintegral.aux_list() ) == 1,\
                                "For single layer auxiliary index algorithm, only a single "+\
                                "auxiliary index per integral class is supported."
                        daux = dintegral.aux_list()[0]
                        new_dintegral = dintegral.int()
                        new_dintegral.add_auxiliary( daux, change=0, op=op_add, update = True )
                        simplified_expr = expr_find_and_replace( simplified_expr, dintegral, new_dintegral )
                else:
                    # Do not modify expression if single layer auxiliary index algorithm not present
                    pass
                # Get list of dsl_variable objects in simplified_expr
                dsl_variable_set = set( expr_list_objects( simplified_expr, dsl_variable ) )
                # Get list of integral index objects (allowed types) other than unrolled_index
                dsl_index_set = set()
                for index_type in self.integral().indextypes():
                    dsl_index_set = dsl_index_set.union( \
                            set( expr_list_objects( simplified_expr, index_type ) ) )
                    dsl_index_set.discard( wrr.unrolled_index() )
                # Get list of auxiliary index objects (allowed types)
                dsl_aux_set = set()
                for aux_type in self.integral().auxtypes():
                    dsl_aux_set = dsl_aux_set.union( \
                            set( expr_list_objects( simplified_expr, aux_type ) ) )
                for integral in dsl_integral_list:
                    for index in integral.index_list():
                        if integral.index_binop(index).right() != 0:
                            # index is incremented/decemented, so we need the 
                            # corresponding loop index and skip variables to be passed.
                            if index != wrr.unrolled_index(): 
                                dsl_index_set.add(index)
                    for aux in integral.aux_list():
                        if integral.aux_binop(aux).right() != 0:
                            # aux is incremented/decemented, so we need the 
                            # corresponding loop index and skip variables to be passed.
                            dsl_aux_set.add(aux)
                    for arg in integral.arg_list():
                        dsl_variable_set.add( arg )
                # Collect arguments for RR function and add indexing arguments
                arg_list = list( dsl_variable_set )
                for x in dsl_index_set:
                    windex = self.windex_index_dict()[ x.name() ]
                    arg_list.append( windex.array_index() )
                    arg_list.append( windex.work_array_skip() )
                unrolled_windex = self.windex_index_dict()[ wrr.unrolled_index().name() ]
                arg_list.append( unrolled_windex.work_array_skip() )
                arg_list.append( unrolled_windex.index_loop_max() )
    
                # Add integral (work) array to list of arguments
                arg_list.append( self.work_array().pointer() )     #work_array
                arg_list.append( self.work_array().array_index() ) #iwrk0
                # Add work_array array_index variables for each aux_index change value
                # If m_layers == 1, then do not add additional array_index variables,
                # since this is unnecessary.
                #
                change_set = set()
                for l in wrr.aux_change_dict().values():
                    for change in l:
                        assert isinstance(change,int) and change >= 0,\
                                "auxiliary index change must be int >= 0"
                        change_set.add(change)
                ordered_change_list = sorted( list( change_set ) )
                if len( ordered_change_list ) > 0 and single_layer == False:
                    istart = 0
                    if ordered_change_list[0] == 0: istart = 1
                    for change in ordered_change_list[istart:]:
                        arg_list.append( self.work_array().array_index_dict()[ change ][0] )

                # Work out the maximum (l) value of unrolled_index required in this RR
                if isinstance( wrr.unrolled_index_lmax(), int ):
                    src_max = wrr.unrolled_index_lmax()
                else: # TEMPORARY code, since unrolled_index_lmax not implemented for 
                      # code with auxiliary indexes yet.
                    src_max = self.generator().generator_options().global_user_lmax()

                # Setup a const_dict object, to assert that all variables passed in should be
                # const, except self.work_array().pointer() -- this will be expressed in the
                # function prototyp
                const_dict = {}
                for arg in arg_list:
                    if not arg is self.work_array().pointer():
                        const_dict[arg] = True
                    else:
                        const_dict[arg] = False

                # Create dsl_function object to represent RR function
                if self.generator().generator_options().parallelism().vectorized() == True:
                    raise Exception("vector output not implemented yet!")
                    # modify arg_list to deal with vectorization
                    pass
                elif self.generator().generator_options().parallelism().vectorized() == False:
                    # Serial output
                    #funcname = self.integral().name()+'_'+wrr.name()
                    funcname = self.full_name()+'_'+wrr.name()
                    f_dsl = dsl_function( funcname, arg_list, vartype="void")
                # Create function_wrapper object for RR function
                f_call = rr_call( f_dsl, gen = self.generator() )
            # /end temporary solution
                f_prototype = rr_prototype( f_dsl, gen= self.generator(), const_dict = const_dict )
                f_source = vrr_source( f_dsl, expr = simplified_expr,\
                                             windex = unrolled_windex,\
                                             windex_src_max = src_max,\
                                             wintegral = self,\
                                             gen = self.generator(),\
                                             const_dict = const_dict )
                wf = function_wrapper(f_dsl,f_call,f_prototype,f_source, gen = self.generator() )
                wrr.set_rr_function( wf )
                # Add to integral_wrapper supporting_function_list, so that code is output
                # with integral class main function
                self._support_function_list.append(wf)
            elif wrr.rr().rrtype() == 'hrr':
                # Create dsl_function object to represent RR function
                # The HRR function is generic, so only one instance of the function needs
                # to be generated, however, since each integral class has a unique set of 
                # objects corresponding to the arguments that the function should use, we
                # need to create a new function_wrapper object with an integral-class specific
                # callable_source object.
                generic_hrr_function_wrapper = self.generator().hrr_function_wrapper_obj()
                # Need to make a integral-class specific arg_list based on arg_list used in generator.py:

                # arg_list = [ Rfromto, lmax_move_from, l_move_from, l_move_to, \
                #              iskip_move_from0, iskip_move_to0,\
                #              iskip_move_from1, iskip_move_to1,\
                #              dwork_array, iwrk0, iwrk1 ]
                move_from_windex = self.windex_index_dict()[ wrr.move_from_index().name() ]
                move_to_windex = self.windex_index_dict()[ wrr.move_to_index().name() ]
                arg_list = wrr.hrr_specific_arg_list()[:]
                # dwork_array     
                arg_list.append( self.work_array().pointer() )
                # iwrk0
                arg_list.append( self.work_array().array_index_dict()[0][0] )
                # iwrk1
                arg_list.append( self.work_array().array_index_dict()[1][0] )

                generic_f_dsl = generic_hrr_function_wrapper.f_dsl()
                specific_f_dsl = dsl_function( generic_f_dsl.name(), arg_list, \
                                                vartype=generic_f_dsl.type() )
                specific_f_call = rr_call( specific_f_dsl, gen = self.generator() )
                generic_f_prototype = generic_hrr_function_wrapper.f_prototype()
                generic_f_source    = generic_hrr_function_wrapper.f_source()
                specific_hrr_function_wrapper = function_wrapper( specific_f_dsl, specific_f_call, \
                                                            generic_f_prototype, generic_f_source, \
                                                            gen = self.generator() )
                                                     
                wf = specific_hrr_function_wrapper
                wrr.set_rr_function( wf )

    def setup_transform_functions(self):
        """
        For each transformation object in self.algorithm_obj() (contraction or 
        spherical_transform) create a function_wrapper object that provides methods 
        for generating function prototypes and executable source code related 
        to that contraction.
        """
        class transform_call(function_wrapper.callable_function_call):
            """
            Instances are callable objects carrying with them all the information
            necessary to generate function calls specific to a contraction.
            """
            pass

        if self.generator().generator_options().contracted() == True or\
           self.generator().generator_options().spherical_transformed() == True:
            # Create a new function_wrapper object by modifying the arg_list from the
            # generic function_wrapper object created in generator.initialize_global_functions()
            # Get list of indexes that have been incremented at end of VRR phase
            incremented_by_vrr_list = self.algorithm_obj().wvrr_group_list[-1].end_index_list
            if len( self.algorithm_obj().whrr_group_list ) > 0:
                incremented_by_vrr_and_hrr_list = self.algorithm_obj().whrr_group_list[-1].end_index_list
            else:
                incremented_by_vrr_and_hrr_list = []
            i = 0
            for t_list, incr_list in [\
                ( self.algorithm_obj().early_transform_list, incremented_by_vrr_list ), \
                ( self.algorithm_obj().late_transform_list, incremented_by_vrr_and_hrr_list ) \
                ]:
                for t in t_list :
                    n_pre_expr  = None
                    arg_list = []
                    # dwork_array     
                    arg_list.append( self.work_array().pointer() )
                    # iwrk0
                    arg_list.append( self.work_array().array_index_dict()[0][0] )
                    # iwrk1
                    arg_list.append( self.work_array().array_index_dict()[1][0] )
                    if isinstance(t,integral_wrapper.algorithm.contraction):
                        assert not t in self.algorithm_obj().late_transform_list, \
                                "contractions must be in early transform sequence"
                        generic_function_wrapper = self.generator().contract_function_wrapper_obj()
                        generic_f_dsl = generic_function_wrapper.f_dsl()
                        generic_f_prototype = generic_function_wrapper.f_prototype()
                        generic_f_source    = generic_function_wrapper.f_source()
                        windex_to_transform = self.windex_index_dict()[ t.index_to_contract.name() ]
                        # Setup arguments for transform function
                        # Generic arg_list:
                        pos = list( reversed( t.end_order_list ) ).index( windex_to_transform.index() )
                        n_pre_expr = None
                        for index in list( reversed( t.end_order_list ) )[ :pos ]:
                            windex = self.windex_index_dict()[ index.name() ]
                            cont_str = t.end_cont_dict[ index.name() ]
                            if n_pre_expr == None:
                                n_pre_expr = windex.work_array_length() * windex.length_dict()[cont_str]
                            else:
                                n_pre_expr = n_pre_expr * windex.work_array_length() * windex.length_dict()[cont_str]
                        if not windex_to_transform.index() in incr_list:
                            tmp_work_array_length = 1
                        else:
                            tmp_work_array_length = windex_to_transform.work_array_length()
                        if n_pre_expr == None:
                            n_pre_expr = tmp_work_array_length
                        else:
                            n_pre_expr = n_pre_expr * tmp_work_array_length
                        #arg_list = [ dwork_array, iwrk0, iwrk1, \
                        #             dcontract_array, \
                        #             n_pre, n_prim, n_cont, iskip_prim ]
                        dtransform_array = windex_to_transform.contract_array().pointer()
                        # Since we only ever call the specific contract functions, we can use
                        # dsl_binop expressions as arguments
                        dnm              = n_pre_expr 
                        dnn              = windex_to_transform.length_dict()['cont']
                        dnk              = windex_to_transform.length_dict()['prim']
                        diskip           = windex_to_transform.skip_dict()['prim']
                        if i == 0: # first transformation
                            diskip = windex_to_transform.skip_dict()['prim']
                        else:
                            diskip = expr_deepish_copy( n_pre_expr )
                    elif isinstance(t,integral_wrapper.algorithm.spherical_transform):
                        generic_function_wrapper = self.generator().sph_trans_function_wrapper_obj()
                        generic_f_dsl = generic_function_wrapper.f_dsl()
                        generic_f_prototype = generic_function_wrapper.f_prototype()
                        generic_f_source    = generic_function_wrapper.f_source()
                        windex_to_transform = self.windex_index_dict()[ t.index_to_transform.name() ]
                        assert windex_to_transform.index() in incr_list, \
                                "index must be incremented before spherical transform"
                        # Setup arguments for transform function
                        # Generic arg_list:
                        pos = list( reversed( t.end_order_list ) ).index( windex_to_transform.index() )
                        n_pre_expr = 1
                        for index in list( reversed( t.end_order_list ) )[ :pos ]:
                            windex = self.windex_index_dict()[ index.name() ]
                            cont_str = t.end_cont_dict[ index.name() ]
                            if n_pre_expr == 1:
                                n_pre_expr = windex.work_array_length()                                 
                                # Multiply by nX_prim/nX_cont, for contracted integral class
                                n_pre_expr = n_pre_expr * windex.length_dict()[cont_str]
                            else:
                                n_pre_expr = n_pre_expr * windex.work_array_length()
                                # Multiply by nX_prim/nX_cont, for contracted integral class
                                n_pre_expr = n_pre_expr * windex.length_dict()[cont_str]
                        # Generic arg_list:
                        #arg_list = [ dwork_array, iwrk0, iwrk1, \
                        #             dsph_trans_array, \
                        #             n_pre, n_cart, n_sph, iskip_cart ]
                        dtransform_array = windex_to_transform.sph_trans_array().pointer()
                        # Since we only ever call the specific contract functions, we can use
                        # dsl_binop expressions as arguments
                        dnm              = n_pre_expr 
                        dnn              = windex_to_transform.work_array_length()
                        dnk              = windex_to_transform.length_dict()['tmp']
                        if i == 0: # first transformation
                            diskip = windex_to_transform.work_array_skip()
                        else:
                            diskip = expr_deepish_copy( n_pre_expr )
                        pos = list( reversed( t.end_index_list ) ).index( windex_to_transform.index() )
                        n_pre_expr = 1
                        for index in list( reversed( t.end_order_list ) )[ :pos ]:
                            windex = self.windex_index_dict()[ index.name() ]
                            cont_str = t.end_cont_dict[ index.name() ]
                            if n_pre_expr == 1:
                                n_pre_expr = windex.work_array_length()
                                # Multiply by nX_prim/nX_cont, for contracted integral class
                                n_pre_expr = n_pre_expr * windex.length_dict()[cont_str]
                            else:
                                n_pre_expr = n_pre_expr * windex.work_array_length()
                                # Multiply by nX_prim/nX_cont, for contracted integral class
                                n_pre_expr = n_pre_expr * windex.length_dict()[cont_str]
                    ### Setup specific arg_list ###
                    # dcontract_array / dsph_trans_array
                    arg_list.append( dtransform_array )
                    # Since we only ever call the specific contract functions, we can use
                    # dsl_binop expressions as arguments
                    # n_pre 
                    arg_list.append( dnm )
                    # nprim / n_cart
                    arg_list.append( dnk )
                    # ncont / n_sph
                    arg_list.append( dnn )
                    # iskip_prim / iskip_cart
                    arg_list.append( diskip )
                    # Create a function_wrapper object with generic callable_prototype and 
                    # callable_source objects, but a callable_function_call object specific
                    # to the integral class and index being contracted.
                    # Attach this function_wrapper to the correponding contraction object in
                    # self.algorithm_obj().early_transform_list
                    ### Setup specific dsl_function and contract_call objects ###
                    specific_f_dsl = dsl_function( generic_f_dsl.name(), arg_list, \
                                                   vartype=generic_f_dsl.type() )
                    specific_f_call = transform_call( specific_f_dsl, gen = self.generator() )
                    ### Create a specific function_wrapper object ###
                    specific_function_wrapper = function_wrapper( specific_f_dsl, specific_f_call, \
                                                                generic_f_prototype, generic_f_source, \
                                                                gen = self.generator() )
                    ### Set contraction.wfunction attribute ###
                    t.wfunction = specific_function_wrapper
                    i += 1

    def setup_base_function(self):
        """
        For the base() expression in self.integral(), create a function_wrapper
        object that provides methods for generating function prototypes and
        executable source code related to the base expression.

        Also detects the presence of special functions which require calling 
        non-standard functions to evaluate (e.g. Boys function) and adds any
        additional necessary header files to self._header_dependency_list.
        """

        # Short names for callable_* derived classes
        base_call      = integral_wrapper_base_function.base_call
        base_prototype = integral_wrapper_base_function.base_prototype
        base_source    = integral_wrapper_base_function.base_source

        arg_set = set()
        base_expr = self.integral().base()
        # Collect dsl_variable arguments for base function
        arg_set = arg_set.union( set( expr_list_objects( base_expr, dsl_variable ) ) )

        # Create list of local variables for base function
        local_set = set()

        # Detect presence of special functions (current implementation only allows Boys function)
        # Detect special function objects an append to list of function_wrapper objects
        special_function_tuple_list = []
        # ... consists of tuples ( special function object, special function wrapper obj )

        ### Default Boys function (built-in) ###
        # These will only be detected if they are of type
        # self.generator().generator_options().default_boys_f()
        default_boys_f = self.generator().generator_options().default_boys_f()
        default_boys_f_list = expr_list_objects( base_expr, default_boys_f )
        if len( default_boys_f_list ) > 0:
            # boys_f object present in base expression, will need to incorporate a call 
            # to the Boys function
            # [ Current implementation restriction ]
            # * Only one instance of default Boys function supported in base expression
            assert len( default_boys_f_list ) <= 1, \
                "Only one instance of the Boys function supported in base expression currently"
            for bf in default_boys_f_list:
                special_function_tuple_list.append( ( bf, self.generator().boys_wrapper_obj() ) )
                # Add additional arguments to base function from boys_f object
                assert bf.m().name() in list( self.windex_aux_dict().keys() ),\
                        "boys function auxiliary index must be part of integral definition"+\
                        " i.e. it must be in windex_aux_dict().keys()"
                bf_aux_windex = self.windex_aux_dict()[ bf.m().name() ]
                bf_arg_dvars = expr_list_objects( bf.x(), dsl_variable )
                bf_local_dvars = []
                bf_local_dvars.append( bf_aux_windex.array_index() )
                dnloop_min = dsl_scalar( bf_aux_windex.loop_length().name()+'_min', vartype = 'int' )
                bf_arg_dvars.append( dnloop_min ) # allows minimum value to be set for Boys function
                                                  # evaluation -- in some cases, may reduce number of
                                                  # operations required (not always)
                bf_arg_dvars.append( bf_aux_windex.loop_length() )
                # work_array_skip unnecessary for 2-layer memory-reduced algorithm
                for dvar in bf_arg_dvars:
                    arg_set.add( dvar )
                for dvar in bf_local_dvars:
                    local_set.add( dvar )
            # Add associated header file to self._header_dependency_list
            self.add_header_dependency( \
                    self.generator().generator_options().boys_function_filename()+'.h')
        
        ### Other special functions ###
        pass # Not currently supported

        # Convert arg_set and local_set to lists
        arg_list = list( arg_set )
        local_list = list( local_set )

        
        # Add integral (work) array to list of arguments (at end)
        arg_list.append( self.work_array().pointer() )
        # Add special base function array_index object (at end)
        work_array_index_base = self.work_array().special_array_index_dict()['base']
        arg_list.append( work_array_index_base )

        # Setup a const_dict object, to assert that all variables passed in should be
        # const, except self.work_array().pointer() -- this will be expressed in the
        # function prototyp
        const_dict = {}
        for arg in arg_list:
            if not arg is self.work_array().pointer():
                const_dict[arg] = True
            else:
                const_dict[arg] = False

        #funcname = self.integral().name()+'_base'
        funcname = self.full_name()+'_base'
        f_dsl = dsl_function( funcname, arg_list, local_list, vartype="void")
        f_call = base_call( f_dsl, gen= self.generator() )
        f_prototype = base_prototype( f_dsl, gen= self.generator(), const_dict = const_dict )
        f_source = base_source( f_dsl, expr = base_expr,\
                                       wintegral = self,\
                                       gen = self.generator(),\
                                       const_dict = const_dict, \
                                       special_function_tuple_list = special_function_tuple_list )
        wf = function_wrapper(f_dsl,f_call,f_prototype,f_source, gen = self.generator() )
        self._support_function_list.append(wf)
        self._base_function = wf

    def setup_data_arrays(self):
        """
        Create data_array objects representing the jump and ccs arrays.
        """
        self._cc_array = self.generator().global_cc_array()
        self._jump_array = self.generator().global_jump_array()
        
    def setup_boys_function(self):
        pass
    
    def algorithm_obj(self):
        return self._algorithm_obj

    def header_dependency_list(self):
        return self._header_dependency_list

    def integral(self):
        return self._integral

    def source_filename(self):
        return self._source_filename

    def header_filename(self):
        return self._header_filename

    def full_name(self):
        return self._full_name

    def source_printer(self):
        return self._source_printer

    def header_printer(self):
        return self._header_printer

    def direction(self):
        return self._direction

    def converter(self):
        return self._converter 

    def base_function(self):
        return self._base_function

    def main_function(self):
        return self._main_function

    def work_array_size_function(self):
        return self._work_array_size_function

    def support_function_list(self):
        return self._support_function_list

    def contraction_info_obj(self):
        return self._contraction_info_obj

    def wrr_list(self):
        return self._wrr_list

    def windex_index_list(self):
        return self._windex_index_list

    def windex_aux_list(self):
        return self._windex_aux_list

    def windex_index_dict(self):
        return self._windex_index_dict

    def windex_aux_dict(self):
        return self._windex_aux_dict

    def work_array(self):
        return self._work_array

    def output_array(self):
        return self._output_array

    def jump_array(self):
        return self._jump_array

    def cc_array(self):
        return self._cc_array

    def auto_local_variable_list(self):
        return self._auto_local_variable_list

    def auto_prim_loop_variable_list(self):
        return self._auto_prim_loop_variable_list

    def auto_prim_loop_depends_dict(self):
        return self._auto_prim_loop_depends_dict

    def manual_local_variable_list(self):
        return self._manual_local_variable_list

    def integral_argument_list(self):
        return self._integral_argument_list

    def set_direction(self,d):
        """
        Updates the Cartesian direction used for RR expression output.
        """
        self._direction.set(d)

    def create_printer(self,file_obj,tab=2,preamble_comment_str=None):
        """
        Set a printer for a file associated with this integral class.

        file_obj:    file, or file-like object for output
        preamble_comment_str:
                     a str which will be automatically output as a comment at the
                     start of the source file
        """
        assert hasattr(file_obj,'write'),\
                "file_obj must be a writable file or file-like object"
        #return printer(endl=';',output_file=file_obj,tab=tab)
        return cprinter(self.converter(),output_file=file_obj,\
                preamble_comment_str=preamble_comment_str)

    def set_source_printer(self,file_obj,preamble_comment_str=None):
        """
        Creates a printer object for output to source file (.c) using the passed in
        file or file-like object. Also sets internal attribute for source file.

        file_obj:    file, or file-like object for output
        preamble_comment_str:
                     a str which will be automatically output as a comment at the
                     start of the source file
        """
        self._source_file = file_obj
        self._source_printer = self.create_printer(self._source_file,\
                preamble_comment_str=preamble_comment_str)

    def set_header_printer(self,file_obj,preamble_comment_str=None):
        """
        Creates a printer object for output to header file (.h) using the passed in
        file or file-like object. Also sets internal attribute for source file.

        file_obj:    file, or file-like object for output
        preamble_comment_str:
                     a str which will be automatically output as a comment at the
                     start of the source file
        """
        self._header_file = file_obj
        self._header_printer = self.create_printer(self._header_file,\
                preamble_comment_str=preamble_comment_str)

