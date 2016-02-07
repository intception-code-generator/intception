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
Provides a generator class and supporting classes for converting Intception dsl_integral
objects into C source code, which adheres to the C99 standard.

classes:
    parallelism         -- defines parallelism / vectorization of output source code
    generator_options   -- encapsulates user-configurable customizations for code generation
    generator           -- methods for generating C99 source code from dsl_integral objects
"""

# General Python modules
import os
import sys
import math
import datetime

# General intception modules
from intception import classtools
from intception.intception_info import prog_unit_info, intception_info
from intception.dsl import *
from intception.dsl_extensions import *
from intception.boys import boys_function, boys_f
from intception.timer import timer

# Generator-specific modules
from intception.generator.src_dsl import *
from intception.generator.supporting_functions import expr_recursive_search, determine_assignment_order
from intception.generator.wrappers import index_wrapper, function_wrapper, vrr_wrapper, \
                                          hrr_wrapper, boys_wrapper
from intception.generator.integral_wrapper import integral_wrapper
from intception.generator.arrays import general_array, integral_array, data_array
from intception.generator.converter import converter
from intception.generator.src_printer import cprinter

# function_wrapper and callable_* derived classes
from intception.generator.callables import generator_sph_trans_function
from intception.generator.callables import generator_contract_function
from intception.generator.callables import generator_hrr_function
from intception.generator.callables import generator_boys_function

# Create instance of prog_unit_info class for C99 code generator
generator_info = prog_unit_info(\
                    unit_name_str = "C code generator",\
                    unit_desc_str = "C99 standard",\
                    version_str   = intception_info.version_str(),\
                    author_str    = intception_info.author_str(),\
                    license_str   = intception_info.license_str() \
                )

# Classes supporting code generation
class parallelism(classtools.classtools):
    """
    Instances of parallelism class carry all the data necessary to define
    the types of parallelism that should be employed by a generator object.

    A parallelism object is associated with generator_options objects, which
    are passed to / created by a generator object. To customise the
    behaviour / output of generator with respect to the parallelism of the output
    source code, users should create a custom parallelism object and pass this
    to a custom generator_options object. This generator_options object should 
    then be passed to generator.__init__() via the opt keyword argument.
    """

    def __init__(self,vectorized = False, vector_batch_size = 4 ):
        """
        Initialize parallelism object.

        vectorized:         Boolean, if True, generate vectorized code.
        vector_batch_size:  length of vectorization batches (innermost loop).
        """
        # vectorized determines whether inner loops in RR functions are
        # vectorized over shells
        self._vectorized = vectorized
        # vector_batch_size is the length of the inner loop where that
        # loop is vectorized
        self._vector_batch_size = vector_batch_size

    def vectorized(self):
        return self._vectorized

    def vector_batch_size(self):
        return self._vector_batch_size

class generator_options(classtools.classtools):
    """
    Instances of generator_options carry data relating to the functioning
    of instances of generator. 
    
    If __init__() is called without arguments, then a set of sane defaults are used, 
    otherwise these can be set to custom values using keyword arguments for __init__().

    When the user creates an instance of generator, an instance of generator_options
    can be provided to customise the behavior of generator. Using a separate class for this
    avoids complexity in the interface to the generator class itself.
    """

    def __init__(self, contracted = None, spherical_transformed = None, \
                    spherical_transform_conditionals = None, \
                    global_user_lmax = None, pllm = None, default_boys_f = None, \
                    global_data_filename = None, global_functions_filename = None, \
                    boys_function_filename = None, \
                    intception_main_filename = None, \
                    use_cblas = None, \
                    cblas_header_filename = None, \
                    cblas_dgemm_min_lengths = None,
                    output_directory = None, \
                    auto_suffix = None, \
                    prim_suffix = None, \
                    cont_suffix = None, \
                    cart_suffix = None, \
                    sph_suffix  = None, \
                    add_suffix_to_global = None, \
                    silent = None, stdout_width = None ):
        """
        Initialize generator_options object.

        If arguments are not present, default values are set.

        contracted:                 whether integrals over primitive indexes or contracted
                                    indexes should be generated (default: True)
        spherical_transformed:      whether Cartesian indexes should be spherically 
                                    transformed (default: True)
        spherical_transform_conditionals:
                                    by default, generated spherical transformed integral code always
                                    does the transformation, regardless of whether l < 2, if True
                                    then conditional blocks are added such that spherical transformation
                                    is done only if l>=2 (default: False)
        global_user_lmax:           maximum angular momentum expected in integral arguments (default: 6)
        pllm:                       parallelism object (default: create instance with default args)
        default_boys_f:             type of objects used in DSL expressions to represent built-in 
                                    Boys function implementation (default: boys_f from boys.py).
        global_data_filename:       base of filename global data arrays are output to (default: 
                                    'global_data')
        global_functions_filename:  base of filename global functions are output to (default:
                                    'global_functions')
        boys_function_filename:     base of filename built-in Boys function code is output to (default:
                                    'boys_function')
        intception_main_filename:   base of filename for the output of the main Intception code, i.e.
                                    where functions acting on the entire library will be output
                                    (default: 'intception_main')
        use_cblas:                  if True, then assume CBLAS interface can be used, and generate 
                                    code which utilizes CBLAS for linear-algebra operations, if False
                                    then generate code which does not rely on CBLAS calls (default: False).
        cblas_header_filename:      header file used in source files where CBLAS calls are made (the directory
                                    containing this header file is assumed to be provided to the compiler,
                                    so only the filename is used in the #include directive, default: 'cblas.h').
        cblas_dgemm_min_lengths:    a Python list containing integers which define the minimum dimensions used in
                                    calling a the cblas_dgemm function, of the form [ m, n, k ], where m, n, k
                                    correspond to the integer arguments for cblas_dgemm (default: [ 3, 3, 10 ] ).
        output_directory:           directory into which Intception source code (which may involve the
                                    creation of subdirectories) is output (default: os.getcwd() )
        auto_suffix:                if True, then all generated filenames and subroutine names have a string 
                                    suffix appended to them which indicates the level of contraction and whether
                                    the indexing of evaluated integrals (default: True )
        prim_suffix:                the suffix appended for primitive (uncontracted) integrals
                                    (default: 'prim')
        cont_suffix:                the suffix appended for contract integrals
                                    (default: 'cont')
        cart_suffix:                the suffix appended for Cartesian integrals
                                    (default: 'cart')
        sph_suffix:                 the suffix appended for spherical integrals
                                    (default: 'sph')
        add_suffix_to_global:       if True, and auto_suffix == True, then the suffixes are added to 
                                    any global functions and source files that are generated
                                    (default: 'False')
        silent:                     if True, intception does not output program status to
                                    sys.stdout (e.g. stage of output processing, timings),
                                    (default: False )
        stdout_width:               integer value for width of status messages to stdout 
                                    (used for nice formatting, default: 80)

        Implementation notes:
            * Options for contraction and spherical transformation are global -- apply to all
              integral classes generated. In theory, each class could carry this option,
              but for the current implementation this applies to all classes.
            * __init__ calls a method of generator_options which enforces 
              limitations based on current implementation.
        """
        attrib_list = []
        # attrib_list items are tuples:
        # ( [object attribute], [__init__ argument], [required object type], [ default value ] )
        self._contracted = self.set_attrib( contracted, True, bool)
        if self.contracted() == True:
            self._spherical_transformed = self.set_attrib( spherical_transformed, True, bool )
        else:
            self._spherical_transformed = self.set_attrib( spherical_transformed, False, bool )
        self._spherical_transform_conditionals = self.set_attrib( spherical_transform_conditionals, False, bool )
        self._global_user_lmax = self.set_attrib( global_user_lmax, 6, int )
        self._parallelism = self.set_attrib( pllm, parallelism(), parallelism )
        # If the user sets a new default_boys_f class, then this class must have an interface
        # which is compatible with the standard boys_f class in boys.py, i.e. with
        # def __init__(self, m, x, [ optional arguments ] )
        # where m is a dsl_integer_index and x is a DSL expression of dsl_variable objects.
        # Implementation notes: 
        # * Currently user-set default_boys_f classes are not supported, as the implementation of
        #   the Boys function relies upon the class "boys_function" in boys.py.
        # * An alternative implementation would require either the same interface 
        #   (via default_boys_f.boys_function_obj() ) or some rewriting of code in setup_boys_function()
        #   in generator to accept an alternative interface.
        self._default_boys_f = self.set_attrib( default_boys_f, boys_f, type ) 
        self._global_data_filename = self.set_attrib( global_data_filename, 'global_data', str )
        self._global_functions_filename = self.set_attrib( global_functions_filename, 'global_functions', str )
        self._boys_function_filename = self.set_attrib( boys_function_filename, 'boys_function', str )
        self._intception_main_filename = self.set_attrib( intception_main_filename, 'intception_main', str )
        self._use_cblas = self.set_attrib( use_cblas, False, bool )
        self._cblas_header_filename = self.set_attrib( cblas_header_filename, 'cblas.h', str )
        self._cblas_dgemm_min_lengths = self.set_attrib( cblas_dgemm_min_lengths, [ 3, 3, 10 ], list )
        assert len( self._cblas_dgemm_min_lengths ) == 3, "cblas_dgemm_min_lengths must be a Python list"
        for i in self._cblas_dgemm_min_lengths: assert isinstance( i, int ), \
                "each element of cblas_dgemm_min_lengths must be an integer"
        self._auto_suffix = self.set_attrib( auto_suffix, True, bool )
        self._prim_suffix = self.set_attrib( prim_suffix, 'prim', str )
        self._cont_suffix = self.set_attrib( cont_suffix, 'cont', str )
        self._cart_suffix = self.set_attrib( cart_suffix, 'cart', str )
        self._sph_suffix  = self.set_attrib( sph_suffix, 'sph', str )
        self._add_suffix_to_global = self.set_attrib( add_suffix_to_global, False, bool )
        self._output_directory = self.set_attrib( output_directory, os.getcwd(), str )
        self._silent = self.set_attrib( silent, False, bool )
        self._stdout_width = self.set_attrib( stdout_width, 80, int )

        # Check implementation restrictions are satisfied
        self.implementation_restriction_check()

    def set_attrib(self,arg,default,allowed_type):
        """
        Return an eiether default or a custom value, checking the type of the attribute.
        
        For arg != None, check whether arg is of allowed_type, and if so, return arg.
        For arg == None, check whether the supplied default is of allowed type, and if so return default.
        """
        if arg != None:
            assert isinstance( arg, allowed_type )
            return arg
        else:
            assert isinstance( default, allowed_type )
            return default

    def implementation_restriction_check(self):
        """
        Run through a set of assertions that attempt to prevent Intception being
        run with option combinations that are not yet supported.
        """
        assert self.default_boys_f() == boys_f,\
                "user-set boys_f classes are not currently supported"
        if self.use_cblas() == True:
            assert self.contracted() == True,\
                "CBLAS use only supported in combination with contracted integrals"
        if self.spherical_transformed() == True:
            assert self.contracted() == True,\
                "spherical transformation only supported for contracted integrals"

    def contracted(self):
        return self._contracted

    def spherical_transformed(self):
        return self._spherical_transformed

    def spherical_transform_conditionals(self):
        return self._spherical_transform_conditionals

    def global_user_lmax(self):
        return self._global_user_lmax

    def parallelism(self):
        return self._parallelism

    def default_boys_f(self):
        return self._default_boys_f

    def intception_main_filename(self):
        return self._intception_main_filename

    def global_data_filename(self):
        return self._global_data_filename
    
    def global_functions_filename(self):
        return self._global_functions_filename

    def boys_function_filename(self):
        return self._boys_function_filename

    def use_cblas(self):
        return self._use_cblas

    def cblas_header_filename(self):
        return self._cblas_header_filename

    def cblas_dgemm_min_lengths(self):
        return self._cblas_dgemm_min_lengths

    def output_directory(self):
        return self._output_directory

    def auto_suffix(self):
        return self._auto_suffix

    def prim_suffix(self):
        return self._prim_suffix

    def cont_suffix(self):
        return self._cont_suffix

    def cart_suffix(self):
        return self._cart_suffix

    def sph_suffix(self):
        return self._sph_suffix

    def add_suffix_to_global(self):
        return self._add_suffix_to_global

    def silent(self):
        return self._silent

    def stdout_width(self):
        return self._stdout_width

class context_timer(timer):
    """
    Class which times "with" block and can output messages formatted nicely to stdout with
    names of operations and completion times.

    NB. Messages will only be nicely formatted if nothing is output to stdout inside
    the "with" block!
    """
    def __init__(self,block_name = "Block",message_width = 80,silent = False):
        """
        block_name:     name of block, output at start of operation
        message_width:  desired width of output (if message has length greater than
                        this, then the output will be split over 2 lines)
        silent:         if True, output to stdout will be suppressed. Object attributes
                        can still be accessed if the "with" statement is used with the
                        "as" keyword.
        """
        self.block_name = block_name
        self.message_width = message_width
        self.silent = silent

    def __enter__(self):
        timer.__enter__(self)
        if not self.silent:
            print( self.block_name,end="",flush=True )
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        timer.__exit__(self, exc_type, exc_value, traceback )
        if not self.silent:
            close_str = "completed in {:.3f} s.".format(self.time_taken)
            if len(self.block_name+close_str) > self.message_width:
                rjust_width = self.message_width
                print( '\n' )
            else:
                rjust_width = self.message_width-len(self.block_name)
            print( close_str.rjust( rjust_width ) )
        return None

# Class for code generator
class generator(classtools.classtools):
    """
    Class providing methods to convert abstract dsl_integral objects into 
    compilable source code.

    Users should create an instance of generator, passing previously created dsl_integral
    objects as arguments. A generator_options object may also be optionally 
    supplied by the user to customise the behaviour of generator (otherwise these are
    automatically created with sane default settings).

    The dsl_integral arguments will be processed by the generator and various wrapper 
    classes representing dsl_* objects in source code form will be created.
    When the user calls the generator.out() method, C source code (C99 standard) for 
    the evaluation of these integral classes will be output according to the options 
    set in the generator_options objects attached to generator.
    """
    def __init__(self,*args,opt = None):
        """
        Create an instance of generator.

        *args:  all non-keyword arguments collected into a tuple to be processed by
                generator.parse_args() method, generally these should be dsl_integral objects
        opt:    optional generator_options object, to allow user customisation of code output --
                if not provided, a default instance of generator_options is created.
        """

        ### Create some object attributes ###
        self._dintegral_list = [] # list of dsl_integral objects
        self._wi_list = [] # list of wrapped integrals
        self._wi_dict = {} # dictionary of wrapper integrals with integral.name() 
                           # as the keys
        self._default_boys_present = False          # set by detect_special_functions
        self._global_functions_present = False      # set during creation of integral_wrappers
        self._hrr_present = False                   # set during creation of integral_wrappers
        self._contract_present = False              # set in initialize_global_functions()
        self._sph_trans_present = False             # set in initialize_global_functions()
        self._boys_wrapper_obj = None               # set in setup_boys_function()
        self._hrr_function_wrapper_obj = None       # set in initialize_global_functions()
        self._contract_function_wrapper_obj = None  # set in initialize_global_functions()
        self._sph_trans_function_wrapper_obj = None # set in initialize_global_functions()
        self._global_functions_header_dependency_list = []
                                                    # set in initialize_global_functions()

        
        # Additional information output after completion of execution if 
        # generator_options.silent() == False, to inform / warn the user.
        self._message_list = []                     # a list of tuples, with the first element
                                                    # containing the message and second element
                                                    # containing the program unit which issued
                                                    # the message

        if opt == None:
            self._generator_options = generator_options() # no opt set by user, use defaults
        else:
            assert isinstance(opt,generator_options), "opt must be of type generator_options"
            self._generator_options = opt
        # Generate suffix for naming routines and files based on contraction level and coordinate system
        # if auto_suffix is True (this is the default)
        if self.generator_options().auto_suffix() == True:
            full_name_suffix_list = []
            # Add suffix for contraction level
            if self.generator_options().contracted() == False:
                # Primitive (uncontracted)
                full_name_suffix_list.append( self.generator_options().prim_suffix() )
            elif self.generator_options().contracted() == True:
                # Contracted
                full_name_suffix_list.append( self.generator_options().cont_suffix() )
            # Add suffix for coordinate system (Cartesian or spherical)
            if self.generator_options().spherical_transformed() == False:
                # Cartesian
                full_name_suffix_list.append( self.generator_options().cart_suffix() )
            elif self.generator_options().spherical_transformed() == True:
                # Spherical
                full_name_suffix_list.append( self.generator_options().sph_suffix() )
            out = []
            for suffix_str in full_name_suffix_list:
                if len( suffix_str ) > 0:
                    out.append( '_' )
                    out.append( suffix_str )
            self._full_suffix = ''.join( out )
        else:
            self._full_suffix = ""
        # global_internal_lmax is updated depending on the required lmax values for integral
        # indexes for processed integral classes
        self._global_internal_lmax = self.generator_options().global_user_lmax()
        # Set names for global data, global functions and Boys function, with optional added
        # suffixes
        self._global_data_name = self.generator_options().global_data_filename()
        self._global_functions_name = self.generator_options().global_functions_filename()
        self._boys_function_name = self.generator_options().boys_function_filename()

        # [ Implementation restriction: ]
        # Adding suffixes to the global data and functions is not supported currently.
        # It is not clear what the best naming convention would be in this case.
        # Recommend that source files with different contraction levels and coordinate systems
        # are output to different directories, to avoid overwriting of global source files (which
        # may result in unexpected behaviour).
        assert self.generator_options().add_suffix_to_global() == False,\
                "Adding a suffix to global source files, data and functions is not supported at present."
        # If requested, add full_suffix to end of name of global functions
        #if self.generator_options().add_suffix_to_global() == True:
        #    self._global_data_name += self.full_suffix()
        #    self._global_functions_name += self.full_suffix()
        #    self._boys_function_name += self.full_suffix()

        # Output some information to stdout
        if not self.generator_options().silent() : 
            stdout_width = self.generator_options().stdout_width()
            print( str( ''.join( [ '-' for i in range(stdout_width) ] ) ) )
            print( intception_info.formatted_info_str( stdout_width ) )
            print( str( ''.join( [ '-' for i in range(stdout_width) ] ) ) )
            print( generator_info.formatted_info_str( stdout_width ) )
            print( str( ''.join( [ '-' for i in range(stdout_width) ] ) ) )
        init_str = "Initializing generator..."
        with context_timer(init_str, self.generator_options().stdout_width(),\
                                     self.generator_options().silent() ) as t: 
            ### Initialize object attributes using args ###
            self.parse_args(args)
            # Initialize global objects (without data)
            self.initialize_global_data()
            self.initialize_global_functions()
            # Detect special functions in dsl_integral base expressions
            # (currently only Boys function supported)
            self.detect_special_functions()
            # If necessary, setup a boys_wrapper object
            self.setup_boys_wrapper()
            # Create integral_wrapper objects for each dsl_integral
            self.setup_integral_wrappers()
            # Set global_internal_lmax based on lmax values for all integral classes
            self.set_global_internal_lmax()
            # Set data for global objects
            self.set_global_data() 
        if not self.generator_options().silent() : 
            print( str( ''.join( [ '-' for i in range(stdout_width) ] ) ) )

        # Assemble preamble comment to be included in all output files
        preamble_list = []
        preamble_list.append( "Source code generated by " )
        preamble_list.append( intception_info.unit_name_str() )
        preamble_list.append( " version " )
        preamble_list.append( intception_info.version_str() )
        preamble_list.append( " using " )
        preamble_list.append( generator_info.unit_name_str() )
        preamble_list.append( " version " )
        preamble_list.append( generator_info.version_str() )
        preamble_list.append( ". " )
        preamble_list.append( datetime.datetime.now().strftime("%c") )
        preamble_list.append( "." )
        self._preamble_comment_str = ''.join( preamble_list )

    def parse_args(self,args):
        """
        Iterate through a tuple of arguments passed to the __init__() method and update corresponding
        object attributes.

        If an arg is dsl_integral, call generator.add_integral for that dsl_integral object.

        args:   tuple of arguments passed to generator.__init__().
        """
        i = 1
        for arg in args:
            if isinstance(arg,dsl_integral):
                self.add_integral(arg)
                
            elif isinstance(arg,parallelism):
                raise Exception("parallelism must be defined in a generator_options object")
            else:
                raise Exception('Argument number '+str(i)+': '+str(arg)+\
                                ' not understood.')
            i +=1

    def decide(self,a):
        """
        Return the index of the first element in an array of integers that is
        greater than 0. 
        
        Used for selecting the Cartesian direction of an RR expression when it
        is output. This is a dumb selector that does not optimize the route 
        through the RR.
        """
        for i in range(len(a)):
            if a[i] > 0: return i
        return 0

    def dintegral_list(self):
        return self._dintegral_list 

    def wi_list(self):
        return self._wi_list

    def generator_options(self):
        return self._generator_options

    def global_cc_array(self):
        return self._global_cc_array

    def global_jump_array(self):
        return self._global_jump_array

    def global_hrr_layer_array(self):
        return self._global_hrr_layer_array

    def global_internal_lmax(self):
        return self._global_internal_lmax
    
    def hrr_present(self):
        return self._hrr_present

    def hrr_function_wrapper_obj(self):
        return self._hrr_function_wrapper_obj

    def set_hrr_present(self,hrr_present):
        assert isinstance( hrr_present, bool ), "hrr_present must be bool"
        self._hrr_present = hrr_present

    def contract_present(self):
        return self._contract_present

    def contract_function_wrapper_obj(self):
        return self._contract_function_wrapper_obj

    def set_contract_present(self,contract_present):
        assert isinstance( contract_present, bool ), "contract_present must be bool"
        self._contract_present = contract_present

    def sph_trans_present(self):
        return self._sph_trans_present
    
    def sph_trans_function_wrapper_obj(self):
        return self._sph_trans_function_wrapper_obj

    def set_sph_trans_present(self,sph_trans_present):
        assert isinstance( sph_trans_present, bool ), "sph_trans_present must be bool"
        self._sph_trans_present = sph_trans_present

    def global_functions_present(self):
        return self._global_functions_present

    def set_global_functions_present(self,global_functions_present):
        assert isinstance( global_functions_present, bool ), "global_functions_present must be bool"
        self._global_functions_present = global_functions_present

    def global_functions_header_dependency_list(self):
        return self._global_functions_header_dependency_list

    def default_boys_present(self):
        return self._default_boys_present

    def boys_wrapper_obj(self):
        return self._boys_wrapper_obj

    def set_default_boys_present(self,default_boys_present):
        assert isinstance( default_boys_present, bool ),"default_boys_present must be bool"
        self._default_boys_present = default_boys_present

    def message_list(self):
        return self._message_list

    def full_suffix(self):
        return self._full_suffix

    def global_data_name(self):
        return self._global_data_name

    def global_functions_name(self):
        return self._global_functions_name

    def boys_function_name(self):
        return self._boys_function_name

    def add_message(self,message_str,info_str):
        """
        Add a message to a list of messages to be printed to stdout after execution.
        The messages are output after the generator.out method is called (if
        generator_options.silent() == False).

        message_str:    informational message for user
        info_str:       context for message (usually program unit providing message)
        """
        assert isinstance( message_str, str ),\
                "message_str must be an instance of str"
        assert isinstance( info_str, str ),\
                "info_str must be an instance of str"
        self._message_list.append( ( message_str, info_str ) )

    def add_integral(self,integral):
        """
        Add a dsl_integral object to self._dintegral_list.
        """
        assert isinstance(integral,dsl_integral), \
                "added integral must be a dsl_integral object"
        assert not integral in self._dintegral_list, \
                "duplicate integral classes not allowed"
        self._dintegral_list.append( integral )

    def preamble_comment_str(self):
        """
        Returns the text of the comment that should be placed at the top of
        every output source file.
        """
        return self._preamble_comment_str

    def initialize_global_data(self):
        """
        Initializes objects representing static global data that is to be shared with 
        all source files, but does not set data.

        cc_array is a data_array object which represents an array of Cartesian 
        component vectors, mapping a single array index to a set of Cartesian 
        components (x, y, z). It has the form ccs[ ncc ][ 3 ], where ncc is the 
        number of Cartesian components, e.g. ccs[ i ][ j ]
        
        int ccs[ 84 ][ 3 ] = 
              { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 2, 0, 0 }, ... }

        jump_array is a data array object which represents an array of 
        'jump' values, which map an increment / decrement in a Cartesian component 
        (x, y, z) to a change in an array index.
        It has the form jump_array[ ncc ][ change ][ 3 ], with change being an 
        integer >= 0. change = 0 is the lowest change value (can be negative), so
        for a situation where we need in/decrements of -2, -1, +1, we would have
        change = 0 (-2), change = 1 (-1), change = 2 (+1). Note that in/decrement = 0 
        is skipped as all 'jump' values woud be zero. E.g. jump[ i ][ j ][ k ], e.g.

        int jump[ 84 ][ 2 ][ 3 ] = 
              { { { 0, 0, 0 }, { 0, 0, 0 } }, 
              { { 0, 0, 0 }, { -1, 0, 0 } }, 
              { { 0, 0, 0 }, { 0, -2, 0 } }, 
              { { 0, 0, 0 }, { 0, 0, -3 } }, 
                ... }

        Note that values of 0 in jump mean that a decrement would result in an 
        unphysical negative angular momentum.

        hrr_layer is a data_array object which represents an array of maximum
        "layer" sizes for a combination of angular momentum values used in the 
        HGP-type HRR, i.e.

            ( a | b ) = ( a - 1i | b + 1i ) - ABi ( a - 1i | b )

        Each "layer" is the number of array elements needed to store the integrals
        with lb..la+lb-lq units of angular momentum on centre B and lq units of
        angular momentum on A. The array contains the maximum "layer" size for all
        "layers" with lq = 0..lb for given combinations of la and lb.

        The hrr_layer array has the form hrr_layer[ la ][ lb ], where la and lb are the
        angular momentum for the final result of the HRR, with angular momentum having
        been shifted from centre B to centre A, e.g. for max la, lb = 6

        int hrr_layer[ 6 ][ 6 ] = 
                { {1, 4, 12, 30, 60, 120, 210 },
                { 4, 12, 30, 60, 120, 210, 350 },
                { 10, 30, 60, 120, 210, 350, 560 },
                { 20, 60, 120, 210, 350, 560, 840 },
                { 35, 105, 210, 350, 560, 840, 1260 },
                { 56, 168, 336, 560, 840, 1260, 1800 },
                { 84, 252, 504, 840, 1260, 1800, 2520 } }

        These are output in source and header files with the base name
        generator.generator_options().global_data_filename() (with optional suffixes).
        """
        ### cc_array ###
        # Create associated DSL objects
        cc_array = dsl_pointer(name = 'ccs', vartype = 'int', constant = True )
        # Create data_array object
        self._global_cc_array = data_array(cc_array,None,None)

        ### jump_array ###
        # Create associated DSL objects
        jump_array = dsl_pointer(name = 'jump', vartype = 'int', constant = True )
        # Create data_array object
        self._global_jump_array = data_array(jump_array,None,None)

        ### hrr_layer_array ###
        # Create associated DSL objects
        hrr_layer_array = dsl_pointer( name = 'hrr_layer', vartype = 'int', constant = True )
        # Create data_array object
        self._global_hrr_layer_array = data_array(hrr_layer_array,None,None)

    def initialize_global_functions(self):
        """
        Initializes objects representing global functions that are to be shared with 
        all source files, if necessary.

        These are output in source and header files with the base name
        generator.generator_options().global_functions_filename() (with optional suffixes).
        """
        ### hrr_function_wrapper_obj ###
        # Short name for callable_* objects
        hrr_call       = generator_hrr_function.hrr_call
        hrr_prototype  = generator_hrr_function.hrr_prototype
        # Select whether to generate HRR function with or without support for 
        # contarcted integrals
        assert isinstance(self.generator_options().contracted(), bool),\
                "generator_options.contracted() must be bool"
        if self.generator_options().contracted() == True:
            hrr_source     = generator_hrr_function.contracted_hrr_source
        elif self.generator_options().contracted() == False:
            hrr_source     = generator_hrr_function.primitive_hrr_source

        ### hrr_function_wrapper_obj ###    
        # Setup function arg_list
        # Create generic argument list for 
        #  (move_from, move_to ) = ( move_from+1,move_to-1 ) + Rfromto * (move_from,move_to-1 ) 
        # This needs to overidden for local calls inside integral class main functions.
        Rfromto         = dsl_position( name = 'Rfromto', vartype = 'double' )
        lmax_move_from  = dsl_scalar( name = 'lmax_move_from', vartype = 'int' )
        l_move_from     = dsl_scalar( name = 'l_move_from', vartype = 'int' )
        loop_l_move_to       = dsl_scalar( name = 'loop_l_move_to', vartype = 'int' )
        iskip_move_from0 = dsl_scalar( name = 'iskip_move_from0', vartype = 'int' )
        iskip_move_to0   = dsl_scalar( name = 'iskip_move_to0', vartype = 'int' )
        iskip_move_from1 = dsl_scalar( name = 'iskip_move_from1', vartype = 'int' )
        iskip_move_to1   = dsl_scalar( name = 'iskip_move_to1', vartype = 'int' )
        dwork_array     = dsl_pointer( name = 'work_array', vartype = 'double' )
        iwrk0           = dsl_scalar( name = 'iwrk0', vartype = 'int' )
        iwrk1           = dsl_scalar( name = 'iwrk1', vartype = 'int' )
        work_array      = integral_array( pointer = dwork_array, array_index = iwrk0 )
        work_array.add_special_array_index( 'hrr', iwrk1 )
        # Arguments for both primitive and contracted version of hrr function
        arg_list = [ Rfromto, lmax_move_from, l_move_from, loop_l_move_to, \
                     iskip_move_from0, iskip_move_to0,\
                     iskip_move_from1, iskip_move_to1 ]
        # Additional arguments for contracted version of hrr function
        if self.generator_options().contracted() == True:
            ncont_move_from = dsl_scalar( name = 'ncont_move_from', vartype = 'int' )
            ncont_move_to = dsl_scalar( name = 'ncont_move_to', vartype = 'int' )
            iskip_cont_move_from0 = dsl_scalar( name = 'iskip_cont_move_from0', vartype = 'int' )
            iskip_cont_move_to0 = dsl_scalar( name = 'iskip_cont_move_to0', vartype = 'int' )
            iskip_cont_move_from1 = dsl_scalar( name = 'iskip_cont_move_from1', vartype = 'int' )
            iskip_cont_move_to1 = dsl_scalar( name = 'iskip_cont_move_to1', vartype = 'int' )
            for var in [ ncont_move_from, ncont_move_to, iskip_cont_move_from0, iskip_cont_move_to0, \
                         iskip_cont_move_from1, iskip_cont_move_to1 ]:
                arg_list.append( var )
        # work_array and array indexes
        for var in [ dwork_array, iwrk0, iwrk1 ]:
            arg_list.append( var )

        if self.generator_options().parallelism().vectorized() == True:
             raise Exception("vector output not implemented yet!")
             # modify arg_list to deal with vectorization
             pass
        elif self.generator_options().parallelism().vectorized() == False:
             # Serial output
             funcname = 'hrr'
             f_dsl = dsl_function( funcname, arg_list, vartype="void")

        # Setup a const_dict object, to assert that all variables passed in should be
        # const, except work_array.pointer() -- this will be expressed in the
        # function prototyp
        const_dict = {}
        for arg in arg_list:
            if not arg is work_array.pointer():
                const_dict[arg] = True
            else:
                const_dict[arg] = False

        # Setup callable_* objects
        f_call = hrr_call( f_dsl, gen = self )
        f_prototype = hrr_prototype( f_dsl, gen = self, const_dict = const_dict )
        if self.generator_options().contracted() == True:
            # Contracted version of hrr function
            f_source = hrr_source( f_dsl, Rfromto, lmax_move_from, \
                                            l_move_from, loop_l_move_to,\
                                            iskip_move_from0, iskip_move_to0, \
                                            iskip_move_from1, iskip_move_to1, \
                                            ncont_move_from,  ncont_move_to, \
                                            iskip_cont_move_from0, iskip_cont_move_to0, \
                                            iskip_cont_move_from1, iskip_cont_move_to1, \
                                            work_array, gen = self, const_dict = const_dict )
        elif self.generator_options().contracted() == False:
            # Primitive version of hrr function
            f_source = hrr_source( f_dsl, Rfromto, lmax_move_from, \
                                            l_move_from, loop_l_move_to,\
                                            iskip_move_from0, iskip_move_to0, \
                                            iskip_move_from1, iskip_move_to1, \
                                            work_array, gen = self, const_dict = const_dict )

        self._hrr_function_wrapper_obj = \
                function_wrapper(f_dsl,f_call,f_prototype,f_source, gen = self )

        ### contract_function_wrapper_obj ###
        # Short name for callable_* objects
        if self.generator_options().contracted() == True:
            contract_call       = generator_contract_function.contract_call
            contract_prototype  = generator_contract_function.contract_prototype
            # Select whether to generate contract function with or without CBLAS
            assert isinstance(self.generator_options().use_cblas(), bool),\
                    "generator_options.use_cblas() must be bool"
            if self.generator_options().use_cblas() == True:
                contract_source = generator_contract_function.cblas_contract_source
                cblas_header_filename = self.generator_options().cblas_header_filename()
                # Add user-set CBLAS header filename to global_functions_header_dependency_list 
                # if not already present
                if not cblas_header_filename in self.global_functions_header_dependency_list():
                    self._global_functions_header_dependency_list.append( cblas_header_filename )
            elif self.generator_options().use_cblas() == False:
                contract_source = generator_contract_function.noblas_contract_source   

            iwrk0           = dsl_scalar( name = 'iwrk0', vartype = 'int', constant = True )
            iwrk1           = dsl_scalar( name = 'iwrk1', vartype = 'int', constant = True )
            dwork_array     = dsl_pointer( name = 'work_array', vartype = 'double', constant = False )
            work_array      = integral_array( pointer = dwork_array, array_index = iwrk0 )
            work_array.add_special_array_index( 'contract', iwrk1 )
            dcontract_array = dsl_pointer( name = 'contract_array', vartype = 'double', constant = True )
            icontract       = dsl_scalar( name = 'icontract', vartype = 'int', constant = True )
            contract_array  = general_array( pointer = dcontract_array,\
                                              array_index = icontract, dimensions = [ None ] )
            n_pre           = dsl_scalar( name = 'n_pre', vartype = 'int', constant = True )
            n_prim          = dsl_scalar( name = 'n_prim', vartype = 'int', constant = True )
            n_cont          = dsl_scalar( name = 'n_cont', vartype = 'int', constant = True )
            iskip_prim      = dsl_scalar( name = 'iskip_prim', vartype = 'int', constant = True )
            arg_list = [ dwork_array, iwrk0, iwrk1, \
                         dcontract_array, \
                         n_pre, n_prim, n_cont, iskip_prim ]
            if self.generator_options().parallelism().vectorized() == True:
                 raise Exception("vector output not implemented yet!")
                 # modify arg_list to deal with vectorization
                 pass
            elif self.generator_options().parallelism().vectorized() == False:
                 # Serial output
                 funcname = 'contract'
                 f_dsl = dsl_function( funcname, arg_list, vartype="void")

            # Setup a const_dict object, to assert that all variables passed in should be
            # const, except work_array.pointer() -- this will be expressed in the
            # function prototype and definition
            const_dict = {}
            for arg in arg_list:
                if not arg is work_array.pointer():
                    const_dict[arg] = True
                else:
                    const_dict[arg] = False

            # Setup callable_* objects
            f_call = contract_call( f_dsl, gen = self )
            f_prototype = contract_prototype( f_dsl, gen = self, const_dict = const_dict )
            f_source = contract_source( f_dsl, \
                                        work_array, contract_array,\
                                        n_pre, n_prim, n_cont, iskip_prim, \
                                        gen = self, const_dict = const_dict )
            self._contract_function_wrapper_obj = \
                    function_wrapper(f_dsl,f_call,f_prototype,f_source, gen = self )
            # Set generator.contract_present() to indicate need for contraction function
            self.set_contract_present(True)
        else:
            self._contract_function_wrapper_obj = None
            # Set generator.contract_present() to indicate no need for contraction function
            self.set_contract_present(False)

        ### sph_trans_function_wrapper_obj ###
        # Short name for callable_* objects
        if self.generator_options().spherical_transformed() == True:
            sph_trans_call       = generator_sph_trans_function.sph_trans_call
            sph_trans_prototype  = generator_sph_trans_function.sph_trans_prototype
            # Select whether to generate sph_trans function with or without CBLAS
            assert isinstance(self.generator_options().use_cblas(), bool),\
                    "generator_options.use_cblas() must be bool"
            if self.generator_options().use_cblas() == True:
                sph_trans_source = generator_sph_trans_function.cblas_sph_trans_source
                cblas_header_filename = self.generator_options().cblas_header_filename()
                # Add user-set CBLAS header filename to global_functions_header_dependency_list
                # if not already present
                if not cblas_header_filename in self.global_functions_header_dependency_list():
                    self._global_functions_header_dependency_list.append( cblas_header_filename )
            elif self.generator_options().use_cblas() == False:
                sph_trans_source = generator_sph_trans_function.noblas_sph_trans_source   

            iwrk0           = dsl_scalar( name = 'iwrk0', vartype = 'int', constant = True )
            iwrk1           = dsl_scalar( name = 'iwrk1', vartype = 'int', constant = True )
            dwork_array     = dsl_pointer( name = 'work_array', vartype = 'double', constant = False )
            work_array      = integral_array( pointer = dwork_array, array_index = iwrk0 )
            work_array.add_special_array_index( 'sph_trans', iwrk1 )
            dsph_trans_array = dsl_pointer( name = 'sph_trans_array', vartype = 'double', constant = True )
            isph_trans       = dsl_scalar( name = 'isph_trans', vartype = 'int', constant = True )
            sph_trans_array  = general_array( pointer = dsph_trans_array,\
                                              array_index = isph_trans, dimensions = [ None ] )
            n_pre           = dsl_scalar( name = 'n_pre', vartype = 'int', constant = True )
            n_cart          = dsl_scalar( name = 'n_cart', vartype = 'int', constant = True )
            n_sph           = dsl_scalar( name = 'n_sph', vartype = 'int', constant = True )
            iskip_cart      = dsl_scalar( name = 'iskip_cart', vartype = 'int', constant = True )
            arg_list = [ dwork_array, iwrk0, iwrk1, \
                         dsph_trans_array, \
                         n_pre, n_cart, n_sph, iskip_cart ]
            if self.generator_options().parallelism().vectorized() == True:
                 raise Exception("vector output not implemented yet!")
                 # modify arg_list to deal with vectorization
                 pass
            elif self.generator_options().parallelism().vectorized() == False:
                 # Serial output
                 funcname = 'sph_trans'
                 f_dsl = dsl_function( funcname, arg_list, vartype="void")

            # Setup a const_dict object, to assert that all variables passed in should be
            # const, except work_array.pointer() -- this will be expressed in the
            # function prototype and definition
            const_dict = {}
            for arg in arg_list:
                if not arg is work_array.pointer():
                    const_dict[arg] = True
                else:
                    const_dict[arg] = False

            # Setup callable_* objects
            f_call = sph_trans_call( f_dsl, gen = self )
            f_prototype = sph_trans_prototype( f_dsl, gen = self, const_dict = const_dict )
            f_source = sph_trans_source( f_dsl, \
                                        work_array, sph_trans_array,\
                                        n_pre, n_cart, n_sph, iskip_cart, \
                                        gen = self, const_dict = const_dict )
            self._sph_trans_function_wrapper_obj = \
                    function_wrapper(f_dsl,f_call,f_prototype,f_source, gen = self )
            # Set generator.sph_trans_present() to indicate need for sph_transion function
            self.set_sph_trans_present(True)
        else:
            self._sph_trans_function_wrapper_obj = None
            # Set generator.sph_trans_present() to indicate no need for sph_transion function
            self.set_sph_trans_present(False)


    def set_global_internal_lmax(self):
        """
        Set self._global_internal_lmax value based on the maximum angular momentum required
        in an integral index for all integral classes.

        integral_wrapper objects for each class of dsl_integral provided by the user should 
        have been created by setup_integral_wrappers prior to the execution of this method.

        This is used in determining the required extents of the precomputed data_array objects
        in set_global_data.
        """
        lmax = self.generator_options().global_user_lmax()
        for wi in self.wi_list():
            for wrr in wi.wrr_list():
                if wrr.rr().rrtype() == 'vrr':
                    # Since VRR operations always come before HRR operations, and the HRR
                    # operations determine the maximum angular momentum required in the VRR
                    # we can use the wrr.unrolled_index_lmax values of VRRs to determine the
                    # highest possible angular momentum value required for all integral
                    # classes.
                    if wrr.unrolled_index_lmax() > lmax:
                        lmax = wrr.unrolled_index_lmax()
        self._global_internal_lmax = lmax

    def set_global_data(self):
        """
        Set static global data that is to be shared with all source files for already
        initialized global data objects.

        cc_array, jump_array, and hrr_layer data are precomputed, and set for the relevant 
        data_array objects. See the docstring for initialize_global_data for details
        of the contents of these arrays.

        These are output in source and header files with the base name
        generator.generator_options().global_data_filename() (with optional suffixes).
        """
        ### cc_array ###
        # Set data
        ccs_per_line = 5 # determines when special linebreak string appended
        cc_array_data = []
        lmax = self.global_internal_lmax()
        tmp_index = dsl_cartesian_gaussian('tmp')
        ncc = 0
        tmp_index.set_angmom([0,0,0])
        while tmp_index.angmom()[0] != lmax + 1:
            ncc += 1 
            cc_array_data.append( tmp_index.angmom()[:] )
            if ncc % ccs_per_line == 0:
                # Insert special linebreak character recognized by data_array.definition_out()
                # but invalid in C code (will be replaced by a line break on source code).
                cc_array_data.append( '#linebreak#' ) 
            tmp_index.increment()
        cc_array_dimensions = [ ncc, 3 ]
        self.global_cc_array().set_dimensions(cc_array_dimensions)
        self.global_cc_array().set_data(cc_array_data)

        ### jump_array ###
        # Determine maximum and minimum changes in Cartesian indexes for all integrals
        jumps_per_line = 1 # determines when special linebreak string appended
        changes_set = set()
        for wi in self.wi_list():
            dintegral = wi.integral()
            for rr in dintegral.rr_list():
                dintegral_list = expr_list_objects( rr.expr(), dsl_integral )
                for dintegral in dintegral_list:
                    for index in dintegral.index_list():
                        if index.is_cartesian() == True:
                            index_binop = dintegral.index_binop( index )
                            if index_binop.op() == op_add:
                                change = index_binop.right()
                            elif index_binop.op() == op_sub:
                                change = -1 * index_binop.right() 
                            else:
                                raise Exception( 'operator '+index_binop.op()+\
                                                 ' not understood.')
                            if change != 0:
                                changes_set.add( change ) 
        changes = list( changes_set )
        changes.sort()

        # Setup dimensions and data
        nchanges = len( changes )
        jump_array_data = []
        lmax = self.global_internal_lmax()
        tmp_index.set_angmom([0,0,0])
        tmp_direction = dsl_direction('x')
        ncc = 0
        # Build array of increments/decrements with dimensions 
        # jump( 3*len(changes), ncc(lmax) )
        while tmp_index.angmom()[0] != lmax + 1:
            jump_cc = []
            for change in changes:
                xyz_list = []
                for d in [ 'x', 'y', 'z' ] :
                    tmp_direction.set(d)
                    if tmp_index.angmom()[tmp_direction.index()] + change >= 0:
                        xyz_list.append( cartesian_index_change(\
                                tmp_index.angmom()[:], change, tmp_direction ) )
                    else: 
                        # Set to zero if change not allowed (e.g. if subtracting angmom
                        # would result in negative angular momentum)
                        xyz_list.append( 0 )
                jump_cc.append( xyz_list )
            jump_array_data.append( jump_cc )
            tmp_index.increment()
            ncc += 1
            if ncc % jumps_per_line == 0:
                # Insert special linebreak character recognized by data_array.definition_out()
                # but invalid in C code (will be replaced by a line break on source code).
                jump_array_data.append( '#linebreak#' ) 

        jump_array_dimensions = [ ncc, nchanges, 3 ]
        jump_array_support_data = { 'min change' : min(changes), 'max change' : max(changes) }
        self.global_jump_array().set_dimensions(jump_array_dimensions)
        self.global_jump_array().set_support_data(jump_array_support_data)
        self.global_jump_array().set_data(jump_array_data)

        ### hrr_layer array ###
        la_values_per_line = 1 # determines when special linebreak string appended
        def hrr_maximum_layer_size(la,lb):
            """
            Determine the maximum size of a block of integrals used in the
            Head-Gordon-Pople-type HRR, transferring angular momentum from
            centre b to centre a.
            
                ( a | b ) = ( a - 1i | b + 1i ) - ABi ( a - 1i | b )
        
            The block size is the number of array elements needed to store
            integrals for a particular angular momentum value of a (lq), and 
            for angular momentum lb to la + lb - lq for b.
        
            The size of this block can be used to determine the maximum size
            of the work array needed for the HRR.
            """
            def ncc(l):
                return (l+1)*(l+2)/2 
            def ncc0(l):
                return (l+1)*(l+2)*(l+3)/6
            max_layer_size = 0
            for lq in range(0,la+1):
                layer_size = (ncc0(la+lb-lq)-ncc0(lb-1))*ncc(lq)
                if layer_size > max_layer_size:
                    max_layer_size = layer_size
            return max_layer_size
        max_layer_size_list = []
        # The highest combination of angular momentum for a pair of indexes, after a HRR
        # operation transferring angular momentum from b to a has been completed
        lmaxa = self.global_internal_lmax() - self.generator_options().global_user_lmax()
        lmaxb = self.generator_options().global_user_lmax()
        hrr_layer_array_data = []
        hrr_layer_array_dimensions = [ lmaxa + 1, lmaxb + 1 ]
        for la in range(lmaxa+1):
            la_list = []
            for lb in range(lmaxb+1):
                la_list.append( int( hrr_maximum_layer_size(la,lb) ) )
            hrr_layer_array_data.append(la_list)
            if la % la_values_per_line == 0:
                # Insert special linebreak character recognized by data_array.definition_out()
                # but invalid in C code (will be replaced by a line break on source code).
                hrr_layer_array_data.append( '#linebreak#' ) 

        self.global_hrr_layer_array().set_dimensions( hrr_layer_array_dimensions )
        self.global_hrr_layer_array().set_data( hrr_layer_array_data )

    def detect_special_functions(self):
        """
        Check the base ( zero angular momentum, dsl_integral.base() )expression associated 
        with each dsl_integral in self.dintegral_list() and if a special function, e.g. 
        the Boys function, is required then set the appropriate attributes of generator to 
        ensure that this function is setup correctly.

        A "special function" is defined as a function which requires an additional function
        call prior to the assignment of the base integral expression. Typically this function
        will generate an array of values which will be used to evaluate the base expression for
        a set of auxiliary index values.
        
        This must be called before self.setup_integral_wrappers(), since the integral_wrappers
        may require the special functions to already have been setup.

        Implementation notes:
        * Currently, only the built-in implementation of the Boys function is supported, which
          is detected by checking for the presence of objects of the type 
          self.generator_options().default_boys_f()
        """
        for dintegral in self.dintegral_list():
            # Check for built-in special functions
            # Boys function
            base_expr = dintegral.base()
            default_boys_f = self.generator_options().default_boys_f()
            default_boys_f_list = expr_list_objects( base_expr, default_boys_f )
            if len( default_boys_f_list ) > 0:
                # boys_f object present in base expression, will need to incorporate a call 
                # to the Boys function
                self.set_default_boys_present(True)

            # Check for external special functions 
            pass # not implemented

    def setup_integral_wrappers(self):
        """
        For each dsl_integral in self.dintegral_list(), create an integral_wrapper object, by
        calling self.create_integral_wrapper(). 
        
        This also populates self._wi_list and self._wi_dict (wi is short for "wrapped integral" 
        in this context).
        """
        for dintegral in self.dintegral_list():
            self.create_integral_wrapper(dintegral)

    def create_integral_wrapper(self,integral):
        """
        Create an integral_wrapper object encapsulating a dsl_integral, with associated
        attributes for source code generation and adds this to self._wi_list and self._wi_dict.

        By default the integral_wrapper.add_header_dependency method is called to add 
        the global data header file.

        Detailed setup of each integral_wrapper object is performed by the 
        integral_wrapper.__init__ method.

        integral:   dsl_integral object to be encapsulated in an integral_wrapper object
        """
        assert isinstance(integral,dsl_integral), \
                "integral must be a dsl_integral object"
        # Generate full_name for integral, including suffixes
        if len( self.full_suffix() ) > 0:
            full_name = integral.name() + self.full_suffix()
        else:
            full_name = integral.name()
        source_filename = full_name+".c"
        header_filename = full_name+".h"
        wi = integral_wrapper( integral, full_name, source_filename, header_filename, gen = self)
        wi.add_header_dependency( self.global_data_name()+".h" )
        # Check whether a global_functions dependency is necessary
        global_functions_reqd = False
        for wrr in wi.wrr_list():
            if isinstance(wrr,hrr_wrapper):
                # HRR present, global_functions required
                global_functions_reqd = True
                # Update generator.hrr_present attribute to reflect need for HRR function code
                if self.hrr_present() != True: self.set_hrr_present(True)
                break
            ### Add further checks here for additional global functions ###
            if self.generator_options().contracted() == True:
                # contractions required, global_functions required
                global_functions_reqd = True
                # generator.contract_present should already be True
                assert self.contract_present() == True, \
                        "generator_options.contracted() and generator.contract_present() "+\
                        "should be equal."
            if self.generator_options().spherical_transformed() == True:
                # spherical transformation required, global_functions required
                global_functions_reqd = True
                # generator.sph_present should already be True
                assert self.sph_trans_present() == True, \
                        "generator_options.spherical_transformed() and generator.sph_trans_present() "+\
                        "should be equal."
        if global_functions_reqd == True:
            wi.add_header_dependency( self.global_functions_name()+".h" )
            # Update generator.global_functions_present attribute to reflect need for 
            # global_functions header/source files
            if self.global_functions_present() != True: self.set_global_functions_present(True)
        self._wi_list.append( wi )
        self._wi_dict[wi.integral().name()] = wi
       
    def setup_boys_wrapper(self):
        """
        If self.default_boys_present() == True, create a boys_wrapper object corresponding
        to the default (built-in) Boys function.

        This creates data_array and function_wrapper objects corresponding to the data and 
        code necessary to evaluate the Boys function in integral classes. These are encapsulated 
        in a boys_wrapper object.

        Implementation notes:
          * Only the default boys_f and boys_function classes in boys.py are supported at present.
          * This method relies on the behaviour and interfaces of these classes.
        """

        # Short name for callable_* objects
        boys_function_call       = generator_boys_function.boys_function_call
        boys_function_prototype  = generator_boys_function.boys_function_prototype
        boys_function_source     = generator_boys_function.boys_function_source

        if self.default_boys_present() == True:
            # Boys function code required by at least one integral class
            # NB. The following code relies on the interface and implementation to boys_f and 
            # boys_function classes from boys.py
            
            # [ Current implementation restriction: ]
            # * Full support for external boys function code is not implemented.
            # Setup data_array objects
            boys_function_obj = self.generator_options().default_boys_f().boys_function_obj
            oof_data_array = data_array( boys_function_obj.oof_pointer(),\
                                                [ boys_function_obj.oof_length() ],\
                                                boys_function_obj.oof_data() )
            # Modify smallx_data_array and mediumx_data_array to have special 
            # '#linebreak#' strings after every set of nvals
            smallx_data  = []
            mediumx_data = []
            for raw_data, modified_data in [ \
                    ( boys_function_obj.smallx_boys_data(), smallx_data ),\
                    ( boys_function_obj.mediumx_boys_data(), mediumx_data ) ]:
                for nval_list in raw_data:
                    modified_data.append( nval_list )
                    modified_data.append( '#linebreak#' )
                del modified_data[-1]
            smallx_data_array = data_array( boys_function_obj.smallx_pointer(),\
             [ boys_function_obj.smallx_xpoints() , boys_function_obj.smallx_nvals() ],\
                                                smallx_data )
            mediumx_data_array = data_array( boys_function_obj.mediumx_pointer(),\
             [ boys_function_obj.mediumx_xpoints() , boys_function_obj.mediumx_nvals() ],\
                                                mediumx_data )
            # DEBUG
            #smallx_data_array.declaration_out( printer(), extern = False )
            #smallx_data_array.definition_out( printer() )
            #mediumx_data_array.declaration_out( printer(), extern = False )
            #mediumx_data_array.definition_out( printer() )

            # Setup a const_dict object, to assert that all variables passed in should be
            # const, except work_array.pointer() -- this will be expressed in the
            # function prototype and definition
            const_dict = {}
            for arg in boys_function_obj.boys_dsl_function().args():
                if not arg is boys_function_obj.boys_output_pointer():
                    const_dict[arg] = True
                else:
                    const_dict[arg] = False

            # Setup function_wrapper object (incl. callable_* objects)
            f_dsl       = boys_function_obj.boys_dsl_function()
            f_call      = boys_function_call( f_dsl, gen = self )
            f_prototype = boys_function_prototype( f_dsl, gen = self, const_dict = const_dict )

           
            local_variable_dict = {}
            for k, v in boys_function_obj.boys_local_variable_dict().items():
                local_variable_dict[ k ] = v
    
            ordered_local_variable_list = determine_assignment_order(\
                                  list( local_variable_dict.values() ) )
            local_variable_list = ordered_local_variable_list
            f_source    = boys_function_source( f_dsl, boys_function_obj,\
                                                local_variable_list, local_variable_dict,\
                                                oof_data_array, smallx_data_array, mediumx_data_array,\
                                                gen = self, const_dict = const_dict )
    
            f_wrapper   = function_wrapper( f_dsl, f_call, f_prototype, f_source, gen = self )
            # DEBUG
            #f_call(printer())
            #f_prototype(printer() )
            #f_source(printer() )
    
            self._boys_wrapper_obj = boys_wrapper( f_wrapper, oof_data_array,\
                                                smallx_data_array, mediumx_data_array, gen = self )


    def out(self,use_stdout=False):
        """
        Outputs integral evaluation source code to the appropriate files, which can
        be compiled to create an intception library. 

        This should be called in an input script for intception after the dsl_integral objects
        have been setup, and passed to a generator object.

        The source code is generated based on the attributes of the generator instance
        (e.g. the integrals which were passed to via generator.__init__, and the
        attributes of generator.generator_options() ).

        A printer object is created for each output file (source and header) and the appropriate
        source code output methods are called with these printer objects. These are
        individual methods of generator which call methods of the appropriate objects (e.g.
        integral_wrapper) to output the source code.

        If use_stdout is True, then the entire source code is output to stdout (terminal)
        rather than files (useful for debugging). Extra lines are added in this case to
        mark the beginning and end of individual files.
        """
        ### Start a timer and output status to stdout ###
        stdout_width = self.generator_options().stdout_width()
        if use_stdout == True:
            destination_str = "stdout"
        else:
            destination_str = self.generator_options().output_directory()
        init_str = "Outputting code to "+destination_str+"..." 
        with context_timer(init_str, self.generator_options().stdout_width(),\
                                     self.generator_options().silent()) as t:
            ### Open all files for output ###
            directory_prefix = self.generator_options().output_directory()+'/'
            open_file_list = []
            # Create file objects for data array output
            # Create file objects for Boys function output (if necessary)
            if use_stdout == False:
                global_data_source_file = open(directory_prefix+\
                        self.global_data_name()+'.c','w')
                open_file_list.append( global_data_source_file )
                global_data_header_file = open(directory_prefix+\
                        self.global_data_name()+'.h','w')
                open_file_list.append( global_data_header_file )
                if self.default_boys_present() == True:
                    boys_function_source_file = open(directory_prefix+\
                            self.boys_function_name()+'.c','w')
                    open_file_list.append( boys_function_source_file )
                    boys_function_header_file = open(directory_prefix+\
                            self.boys_function_name()+'.h','w')
                    open_file_list.append( boys_function_header_file )
                if self.global_functions_present() == True:
                    global_functions_source_file = open(directory_prefix+\
                            self.global_functions_name()+'.c','w')
                    open_file_list.append( global_functions_source_file )
                    global_functions_header_file = open(directory_prefix+\
                            self.global_functions_name()+'.h','w')
                    open_file_list.append( global_functions_header_file )
            elif use_stdout == True:
                global_data_source_file = sys.stdout
                global_data_header_file = sys.stdout
                if self.default_boys_present() == True:
                    boys_function_source_file = sys.stdout
                    boys_function_header_file = sys.stdout
                if self.global_functions_present() == True:
                    global_functions_source_file = sys.stdout
                    global_functions_header_file = sys.stdout
            global_data_header_printer = cprinter( converter( None, dsl_direction('x'), None, self ),\
                                                    output_file = global_data_header_file, \
                                                    preamble_comment_str = self.preamble_comment_str() )
            global_data_source_printer = cprinter( converter( None, dsl_direction('x'), None, self ),\
                                                    output_file = global_data_source_file, \
                                                    preamble_comment_str = self.preamble_comment_str() )
            if self.default_boys_present() == True: 
                boys_function_header_printer = cprinter(\
                        converter = converter(None, dsl_direction('x'), None, self ),\
                        output_file = boys_function_header_file, \
                        preamble_comment_str = self.preamble_comment_str() )
                boys_function_source_printer = cprinter(\
                        converter = converter(None, dsl_direction('x'), None, self ),\
                        output_file = boys_function_source_file, \
                        preamble_comment_str = self.preamble_comment_str() )
            if self.global_functions_present() == True: 
                global_functions_header_printer = cprinter(\
                        converter = converter(None, dsl_direction('x'), None, self ),\
                        output_file = global_functions_header_file, \
                        preamble_comment_str = self.preamble_comment_str() )
                global_functions_source_printer = cprinter(\
                        converter = converter(None, dsl_direction('x'), None, self ),\
                        output_file = global_functions_source_file, \
                        preamble_comment_str = self.preamble_comment_str() )

            # Create file objects for integral class output
            for wi in self._wi_list:
                if use_stdout == False:
                    source_file = open(directory_prefix+wi.source_filename(),'w')
                    open_file_list.append( source_file )
                    header_file = open(directory_prefix+wi.header_filename(),'w')
                    open_file_list.append( header_file )
                elif use_stdout == True:
                    source_file = sys.stdout
                    header_file = sys.stdout
                wi.set_source_printer( source_file, \
                        preamble_comment_str = self.preamble_comment_str() )
                wi.set_header_printer( header_file, \
                        preamble_comment_str = self.preamble_comment_str() )

            ### Output data arrays ###
            self.global_data_out( global_data_header_printer, global_data_source_printer, use_stdout )

            ### Output Boys function (if necessary) ###
            if self.default_boys_present() == True:
                self.boys_function_out( boys_function_header_printer, \
                        boys_function_source_printer, use_stdout )

            ### Output global functions (if necessary) ###
            if self.global_functions_present() == True:
                self.global_functions_out( global_functions_header_printer,\
                        global_functions_source_printer, use_stdout )
            
            ### Output source code for main parts of library ###
            pass


            ### Output source code for integral classes ###
            for wi in self._wi_list:
                self.integral_class_out( wi, use_stdout )

            ### Close all open files ###
            for f in open_file_list:
                f.close()
        
        ### Output status and timing information to stdout ###
        if not self.generator_options().silent() :
            stdout_width = self.generator_options().stdout_width()
            print( str( ''.join( [ '-' for i in range(stdout_width) ] ) ) )
            if len( self.message_list() ) > 0:
                # Print messages to user
                print( "Messages:" )
                i = 0
                for msg_tup in self.message_list():
                    # Print over at least 2 lines
                    message = str(i) + ' : ' + msg_tup[0]
                    message_line_list = [ message[j:j+stdout_width ] for j in range(0, len(message), stdout_width ) ] 
                    info    = "[ "+msg_tup[1]+" ]"
                    info_line_list    = [ info[j:j+stdout_width ] for j in range(0, len(info), stdout_width ) ]
                    for line in message_line_list:
                        print( line )
                    for line in info_line_list:
                        print( line )
                    i += 1
                # Output line of hyphens to close program
                print( str( ''.join( [ '-' for i in range(stdout_width) ] ) ) )

    def global_data_out(self,ph,pc,use_stdout=False):
        """
        Outputs header and source files containing global data objects.

        The declaration_out and definition_out methods of the global data data_array
        objects are called to output source code via the printer objects ph and pc.
        
        ph:         printer object that will output to header file.
        pc:         printer object that will output to source file.
        use_stdout: if True, output extra lines to denote start and end of files.
        """
        global_data_filename = self.global_data_name()
        if use_stdout == True: 
            ph.out('### START HEADER FILE ('+global_data_filename+') ###',endl='')
        ### Header file ###
        # cc array
        self.global_cc_array().declaration_out(ph)
        # jump array
        self.global_jump_array().declaration_out(ph)
        if self.hrr_present() == True:
            # hrr_layer array
            self.global_hrr_layer_array().declaration_out(ph)
        if use_stdout == True: 
            ph.out('### END HEADER FILE ('+global_data_filename+') ###',endl='')
            ph.blankline()

        if use_stdout == True: 
            pc.out('### START SOURCE FILE ('+global_data_filename+') ###',endl='')
        ### Source file ###
        # cc array
        self.global_cc_array().definition_out(pc)
        pc.blankline()
        # jump array
        self.global_jump_array().definition_out(pc)
        pc.blankline()
        if self.hrr_present() == True:
            # hrr_layer array
            self.global_hrr_layer_array().definition_out(pc)
        if use_stdout == True: 
            pc.out('### END SOURCE FILE ('+global_data_filename+') ###',endl='')
            pc.blankline()

    def global_functions_out(self,ph,pc,use_stdout=False):
        """
        Outputs header and source files containing global function objects.

        The prototype_out and source_out methods of the global function function_wrapper
        objects are called to output source code via the printer objects ph and pc.
        
        ph:         printer object that will output to header file.
        pc:         printer object that will output to source file.
        use_stdout: if True, output extra lines to denote start and end of files.
        """
        global_functions_filename = self.global_functions_name()
        if use_stdout == True: 
            ph.out('### START HEADER FILE ('+global_functions_filename+') ###',endl='')
        ### Header file ###
        if self.hrr_present() == True:
            self.hrr_function_wrapper_obj().prototype_out(ph)
        if self.contract_present() == True:
            self.contract_function_wrapper_obj().prototype_out(ph)
        if self.sph_trans_present() == True:
            self.sph_trans_function_wrapper_obj().prototype_out(ph)
        ph.blankline()
        if use_stdout == True: 
            ph.out('### END HEADER FILE ('+global_functions_filename+') ###',endl='')
            ph.blankline()

        if use_stdout == True: 
            pc.out('### START SOURCE FILE ('+global_functions_filename+') ###',endl='')
        ### Source file ###
        for hname in self.global_functions_header_dependency_list(): # dependencies between intception files
            pc.include( hname, quotes = True )
        if self.hrr_present() == True:
            self.hrr_function_wrapper_obj().source_out(pc)
        if self.contract_present() == True:
            self.contract_function_wrapper_obj().source_out(pc)
        if self.sph_trans_present() == True:
            self.sph_trans_function_wrapper_obj().source_out(pc)
        if use_stdout == True: 
            pc.out('### END SOURCE FILE ('+global_functions_filename+') ###',endl='')
            pc.blankline()

    def boys_function_out(self,ph,pc,use_stdout=False):
        """
        Outputs header and source files containing Boys function and associated data.
        
        For tabulated Boys function data, the declaration_out and definition_out 
        methods of data_array objects are called to output source code via printer objects
        ph and pc. 

        For the Boys function itself, prototype_out and source_out methods of the corresponding
        function_wrapper object are called to output source code via printer objects ph and pc.
        
        ph:         printer object that will output to header file.
        pc:         printer object that will output to source file.
        use_stdout: if True, output extra lines to denote start and end of files.
        """
        boys_function_filename = self.boys_function_name()
        boys_wrapper_obj = self.boys_wrapper_obj()
        if use_stdout == True: 
            ph.out('### START HEADER FILE ('+boys_function_filename+') ###',endl='')
        ### Header file ###
        # oof data array
        boys_wrapper_obj.oof_data_array().declaration_out(ph, extern = False)
        # smallx data array
        boys_wrapper_obj.smallx_data_array().declaration_out(ph, extern = False)
        # mediumx data array
        boys_wrapper_obj.mediumx_data_array().declaration_out(ph, extern = False)
        # boys function
        boys_wrapper_obj.function_wrapper_obj().prototype_out(ph)
        if use_stdout == True: 
            ph.out('### END HEADER FILE ('+boys_function_filename+') ###',endl='')
            ph.blankline()

        if use_stdout == True: 
            pc.out('### START SOURCE FILE ('+boys_function_filename+') ###',endl='')
        ### Source file ###
        pc.include('math.h')
        pc.blankline()
        # oof data array
        boys_wrapper_obj.oof_data_array().definition_out(pc)
        pc.blankline()
        # smallx data array
        boys_wrapper_obj.smallx_data_array().definition_out(pc)
        pc.blankline()
        # mediumx data array
        boys_wrapper_obj.mediumx_data_array().definition_out(pc)
        pc.blankline()
        # boys function
        boys_wrapper_obj.function_wrapper_obj().source_out(pc)

        if use_stdout == True: 
            pc.out('### END SOURCE FILE ('+boys_function_filename+') ###',endl='')
            pc.blankline()

    def integral_class_out(self,wi,use_stdout=False):
        """
        Outputs all the source code necessary to evaluate integrals of the type
        described in the dsl_integral object encapsulated in the intergal_wrapper object
        wi.

        The "main function" which would be called to evaluate the integral class is output via 
        integral_wrapper.main_function().prototype_out, and 
        integral_wrapper.main_function().source_out
        using the printer objects integral_wrapper.header_printer() and 
        integral_wrapper.source_printer().

        Supporting functions are output using the same printer objects, by calling the 
        prototype_out and source_out methods of function_wrapper objects in 
        integral_wrapper.support_function_list(). These are typically functions for
        sizing the work_array, evaluating RRs and the base expression (zero angular momentum case).

        wi:         integral_wrapper object with a dsl_integral object accessible through
                    wi.integral() which defines the integral class.
        use_stdout: this does not change the printer objects used for output (these should
                    already be set for wi), but does add some extra lines to the output
                    to make it more human-readable (useful for debugging).
        """
        assert isinstance(wi,integral_wrapper),\
                "wi should be of class self.wrapped_integral"
        header_filename = wi.header_filename()
        source_filename = wi.source_filename()
        integral_name = wi.full_name()
        ph = wi.header_printer()
        pc = wi.source_printer() 

        ### HEADER FILE ###
        if use_stdout == True: 
            ph.out('### START HEADER FILE ('+integral_name+') ###',endl='')
        ph.include('stdio.h')
        ph.blankline()
        # Output main integral function prototype
        wi.main_function().prototype_out(ph)
        # Output supporting function prototypes
        for f in wi.support_function_list():
            f.prototype_out(ph)
        if use_stdout == True: 
            ph.out('### END HEADER FILE ('+integral_name+') ###',endl='')
            ph.blankline()

        ### SOURCE FILE ###
        if use_stdout == True: 
            pc.out('### START SOURCE FILE ('+integral_name+') ###',endl='')
        pc.include('stdio.h')
        pc.include('math.h')
        pc.include( header_filename, quotes = True)
        for hname in wi.header_dependency_list(): # dependencies between intception files
            pc.include( hname, quotes = True )
        pc.blankline()
        # Output main integral function
        wi.main_function().source_out(pc)
        pc.blankline()
        # Output supporting functions
        for f in wi.support_function_list():
            f.source_out(pc)
        if use_stdout == True: 
            pc.out('### END SOURCE FILE ('+integral_name+') ###',endl='')
            pc.blankline()

