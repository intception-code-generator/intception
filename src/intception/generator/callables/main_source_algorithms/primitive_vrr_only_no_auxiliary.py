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
Provides the algo_main_source derived class from integral_wrapper_main_function.main_source
which contains the code necessary to generate an integral evaluation code
for general integral classes with the following properties:
    * primitive integrals
    * VRR-only algorithm
    * no auxiliary indexes
"""
from intception.dsl import *
from intception.generator.callables import integral_wrapper_main_function
from intception.generator.callables import integral_wrapper_work_array_size_function
from intception.generator.algo_descriptor import algo_descriptor
from intception.generator.algo_variable_setup import algo_variable_setup
from intception.generator.algo_indexing_setup import algo_indexing_setup
from . import primitive_general

### Create a specific algo_descriptor object for this algorithm ###
algo_desc = algo_descriptor( is_contracted = False, hrr_present = False, aux_index_present = False )

### Create a specific algo_indexing_setup derived class for this algorithm ###
class integral_array_indexing_setup(algo_indexing_setup):
    def __call__(self):
        """
        Setup indexing attributes for integral arrays attached to integral_wrapper
        based on specific algorithm details (sets up integral_wrapper.work_array).
        """
        wintegral = self.wintegral

        # Go through RR list to determine order of indexes in work array 
        work_array_windex_index_list = []
        # Separate vrr_wrapper and hrr_wrapper objects
        for wrr in wintegral.wrr_list():
            assert wrr.rr().rrtype() == 'vrr', \
                "Only VRR RRs allowed in this algorithm (contracted, VRR-only, with auxiliary)"
            if wrr.rr().rrtype() == 'vrr':
                windex = wintegral.windex_index_dict()[ wrr.changing_index().name() ]
                work_array_windex_index_list.append( windex )
        
        # Any additional integral indexes not captured from wrr_list
        for windex in wintegral.windex_index_list():
            assert windex in work_array_windex_index_list,\
                    "Currently, all integral indexes must be incremented in an RR"
        # Add array_index for placement of base function output in work_array
        work_array_index_base = \
                        dsl_scalar( name = 'iwrkb', vartype = 'int' )
        wintegral.work_array().add_special_array_index('base',work_array_index_base )
        
        # Set work_array windex_list to integral indexes (ordered) +
        # auxiliary indexes (unordered)
        wintegral.work_array().set_windex_list( work_array_windex_index_list )

### Create a specific algo_variable_setup derived class for this algorithm ###
# Integral argument list (use general class for primitive integrals)
integral_argument_setup = primitive_general.integral_argument_setup
        
# Automatically assigned local_variables
automatic_local_variable_setup = primitive_general.automatic_local_variable_setup

# Manually assigned local variables
class manual_local_variable_setup(algo_variable_setup):
    def __call__(self,manual_list):
        """
        Setup manual_list of integral local variables (this is called in 
        integral_wrapper.setup_local_variable_list)
        
        manual_list:    a Python list, to be populated with dsl_* objects
                        that will be manually assigned in the main function
                        for an integral class.
        """
        wintegral     = self.wintegral
        algorithm_obj = self.wintegral.algorithm_obj()
        manual_set = set()
        # No HRRs, last wrr_group is a wvrr_group object
        final_wrr_group = algorithm_obj.wvrr_group_list[-1]
        for index in final_wrr_group.end_index_list:
            windex = wintegral.windex_index_dict()[ index.name() ]
            for var in [ windex.work_array_skip(),\
                        windex.work_array_offset(),\
                        windex.array_index() ]:
                manual_set.add( var )

        ### VRR sequence variables ###
        for wvrr_group in algorithm_obj.wvrr_group_list:
            for wrr in wvrr_group.wrr_list:
                # changing/unrolled index
                windex = wintegral.windex_index_dict()[ wrr.changing_index().name() ] 
                manual_set.add( windex.index_loop_max() )
                # loop_index_list indexes
                for index in wrr.loop_index_list():
                    windex = wintegral.windex_index_dict()[ index.name() ]
                    manual_set.add( windex.loop_length() )
        # work_array_offset set at the end of VRR sequence, depends on work_array_length
        for index in final_wrr_group.end_index_list:
            windex = wintegral.windex_index_dict()[ index.name() ]
            manual_set.add( windex.work_array_length() )

        ### Further all algorithm variables ###
        manual_set.add( wintegral.work_array().array_index() )
        # special_array_index_dict contains array_index for base function output, required for all
        # algorithms
        manual_set.add( wintegral.work_array().special_array_index_dict()['base'] )

        # Append unique set elements to list for output
        for var in manual_set:
            manual_list.append(var)

### Create a specific work_array_size_source class for this algorithm ###
class algo_work_array_size_source(primitive_general.primitive_work_array_size_source):
    def serial_out(self,p):
        """
        For integral class using this algorithm, output a simple
        function that returns an integer number for the minimum 
        required size of the work_array.
        
        Primitive integrals:
        The largest the work_array will ever need to be is 
            2 * ncc0(la) * ncc0(lb) * ... ncc0(lX) + m_max + 1
        """

        f = self._f_dsl
        local_variable_list = self._local_variable_list
        # Output code common to all work array size algorithms
        p.funcdef( f, self.const_dict() )

        # Local variable declarations
        for dvar in local_variable_list:
            p.declaration( dvar )

        # Variable assignment lines
        for dvar in local_variable_list:
            if dvar.expr() != None:
                if dvar.is_cartesian() == False:
                    p.out( dvar.assign( dvar.expr() ) )
                elif dvar.is_cartesian() == True:
                    for d in [ 0, 1, 2 ]:
                        p.converter.direction().set(d)
                        p.out( dvar.assign( dvar.expr() ) )

        #### Determine maximum memory requirement for VRR phase ###
        wvrr_group_list = self._wintegral.algorithm_obj().wvrr_group_list
        memreq_expr = self.vrr_memreq_expr(wvrr_group_list[-1] )
        # This is always the maximum memory requirement, so simply return
        # value of expression -- no comparisons required.
        ### Output return statement ###
        p.returnstmt( memreq_expr )

        p.endfuncdef()

class algo_main_source(primitive_general.primitive_algo_main_source):
    def serial_out(self,p):
        """
        Output main integral evaluation code for the case of 
        - primitive
        - VRR-only
        - no-auxiliary index
        integrals requested.
        """
        f = self._f_dsl
        auto_local_variable_list = self._wintegral.auto_local_variable_list()
        manual_local_variable_list = self._wintegral.manual_local_variable_list()

        algorithm = self._wintegral.algorithm_obj()
        final_wrr_group = None # set to the wrr_group for the set of RR operations executed last
        # Output function definition top part
        p.funcdef( f, self.const_dict() )

        ### Local variable declarations ###
        var_to_declare_list = auto_local_variable_list + \
                                manual_local_variable_list
        for dvar in var_to_declare_list:
            p.declaration( dvar )

        ### Automatic variable assignment ###
        # (for variables not depending on primitive exponents)
        for dvar in auto_local_variable_list:
            if hasattr(dvar,'expr') and dvar.expr() != None:
                if dvar.is_cartesian() == False:
                    p.out( dvar.assign( dvar.expr() ) )
                elif dvar.is_cartesian() == True:
                    for d in [ 0, 1, 2 ]:
                        p.converter.direction().set(d)
                        p.out( dvar.assign( dvar.expr() ) )

        wvrr_group_list = algorithm.wvrr_group_list

        # Set final_wrr_group (contains information on the shape of the array of intermediates
        # immediately before these are copied to output array)
        final_wrr_group = wvrr_group_list[-1]

        ### Setup work_array indexing for base class evaluation and VRR sequence ###
        # Set work_array_length and work_array_skip values for sequence of wvrr_groups
        # These will define the shape and indexing of the end result of the sequence of VRRs,
        # which is completely characterised by the last wvrr_group.
        final_wvrr_group = wvrr_group_list[-1]
        # Assign work_array_length expressions for VRR operations
        for index in final_wvrr_group.end_index_list: # incremented indexes
            windex   = self._wintegral.windex_index_dict()[ index.name() ]
            ncc_type = final_wvrr_group.end_ncc_type_dict[ index.name() ]
            l_expr   = final_wvrr_group.end_max_l_expr_dict[ index.name() ]
            assert ncc_type == 'ncc0', \
                    "VRR operations should always result in a ncc_type of 'ncc0'"
            work_array_length_var  = windex.work_array_length()
            work_array_length_expr = ( l_expr + 1 )*( l_expr + 2 )*( l_expr + 3 )/6
            p.out( work_array_length_var.assign( work_array_length_expr ) )

        # Assign work_array_skip expressions for per-Cartesian-component and
        # per-primitive component skips for VRR sequence.
        # Order of indexes does not change, nor does indexing of array during
        # VRR sequence, so these can be set from information in the final_wvrr_group
        skip_index_list = list( reversed( final_wvrr_group.end_order_list[:] ) )
        incremented_index_list = list( reversed( final_wvrr_group.end_index_list[:] ))
        index  = skip_index_list[0]
        windex = self._wintegral.windex_index_dict()[ index.name() ]
        work_array_skip_var  = windex.work_array_skip() # skip for each Cartesian component
        work_array_skip_expr = 1
        if windex.index() in incremented_index_list:
            p.out( work_array_skip_var.assign( work_array_skip_expr ) )
        for index in skip_index_list[1:]:
            if windex.index() in incremented_index_list:
                # Set per-Cartesian component skip
                # ..multiply by work_array_length of previous index (last skip output)
                if work_array_skip_expr == 1:
                    work_array_skip_expr = windex.work_array_length() 
                else:
                    work_array_skip_expr = work_array_skip_expr * windex.work_array_length()
            # Set windex for current index in loop
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            # Output per-Cartesian component skip (if current index is incremented)
            if windex.index() in incremented_index_list:
                work_array_skip_var  = windex.work_array_skip()
                p.out( work_array_skip_var.assign( work_array_skip_expr ) )
                
        # Special base function work_array array_index dsl_scalar object
        work_array_index_base = \
                self._wintegral.work_array().special_array_index_dict()['base']
        # For primitive integrals, the iwrk0 array position is always zero initially
        work_array_index0_initial_expr = 0
        # No HRR or auxiliary index => no need for work_array_index1,
        # work_array_index0 and work_array_index_base are set to zero
        work_array_index0 = self._wintegral.work_array().array_index()
        # Output initial work_array array_index assignments
        p.out( work_array_index0.assign(work_array_index0_initial_expr ) )
        # Output assignment for special base function work_array array_index
        p.out( work_array_index_base.assign( work_array_index0 ) )

        ### VRR sequence (and base evaluation) ###
        self.vrr_sequence_no_auxiliary_out( p, wvrr_group_list,\
            work_array_index0, work_array_index0_initial_expr )

        # Set state of work_array dimensions (offsets) after VRR sequence
        # for integral indexes that have been incremented.
        # The work_array_length variables do not need to be set, since they
        # do not change during the VRR sequence and are set prior to this.
        # wvrr_groups currently always end with 'ncc0' type indexes, and this is 
        # assumed here.
        for index in final_wvrr_group.end_index_list :
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            ncc_type = final_wvrr_group.end_ncc_type_dict[ index.name() ]
            work_array_length_var   = windex.work_array_length()
            work_array_offset_var   = windex.work_array_offset()
            output_array_length_var = windex.output_array_length() 
            # Determine work_array_length_expr and work_array_offset_expr 
            # from ncc_type
            assert ncc_type == 'ncc0',\
                    "VRRs must end with ncc_type == 'ncc0' for all indexes"
            work_array_offset_expr = work_array_length_var-output_array_length_var
            # Output assignments for offset only
            p.out( work_array_offset_var.assign( work_array_offset_expr ) )

        ### Copy correct integral intermediates from intermediate array ###
        final_operation = final_wvrr_group
        self.copy_to_output_array_out(p,final_operation,work_array_index0)
        
        # End of function definition
        p.endfuncdef()

