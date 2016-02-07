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
    * contracted integrals
    * VRR-only algorithm
    * Auxiliary index (2-layer algorithm)
"""
from intception.dsl import *
from intception.dsl_extensions import *
from intception.generator.callables import integral_wrapper_main_function
from intception.generator.callables import integral_wrapper_work_array_size_function
from intception.generator.algo_descriptor import algo_descriptor
from intception.generator.algo_variable_setup import algo_variable_setup
from intception.generator.algo_indexing_setup import algo_indexing_setup
from . import contracted_general

### Create a specific algo_descriptor object for this algorithm ###
algo_desc = algo_descriptor( is_contracted = True, hrr_present = False, aux_index_present = True )

### Create a specific algo_indexing_setup derived class for this algorithm ###
class integral_array_indexing_setup(algo_indexing_setup):
    def __call__(self):
        """
        Setup indexing attributes for integral arrays attached to integral_wrapper
        based on specific algorithm details (sets up integral_wrapper.work_array).
        """
        wintegral = self.wintegral

        # Go through RR list to determine order of indexes in work array (for each auxiliary index
        # combination)
        work_array_windex_index_list = []
        # (see definition in vrr_wrapper class)
        work_array_windex_aux_list = []
        wvrr_list = []
        # Separate vrr_wrapper and hrr_wrapper objects
        hrr_found = False
        for wrr in wintegral.wrr_list():
            assert wrr.rr().rrtype() == 'vrr', \
                "Only VRR RRs allowed in this algorithm (contracted, VRR-only, with auxiliary)"
            if wrr.rr().rrtype() == 'vrr':
                wvrr_list.append(wrr)
                windex = wintegral.windex_index_dict()[ wrr.changing_index().name() ]
                work_array_windex_index_list.append( windex )
        
        # Any additional integral indexes not captured from wrr_list
        for windex in wintegral.windex_index_list():
            assert windex in work_array_windex_index_list,\
                    "Currently, all integral indexes must be incremented in an RR"
        # * Get (unordered) list of auxiliary indexes and add to list
        # * Create new array_index_dict entries for each change value of 
        #   auxiliary index
        all_wvrr_aux_change_dict = {}
        for windex in wintegral.windex_aux_list():
            work_array_windex_aux_list.append( windex )
            all_wvrr_aux_change_dict[ windex.index().name() ] = set()
        for wvrr in wvrr_list:
            for windex in work_array_windex_aux_list:
                for change in wvrr.aux_change_dict()[ windex.index().name() ]:
                    all_wvrr_aux_change_dict[ windex.index().name() ].add( change )
        # Add array_index for placement of base function output in work_array
        # The following code assumes only a single auxiliary index
        work_array_index_base = \
                        dsl_scalar( name = 'iwrkb', vartype = 'int' )
        wintegral.work_array().add_special_array_index('base',work_array_index_base )
        assert len( work_array_windex_aux_list ) <= 1,\
                "Only 0 or 1 auxiliary indexes per integral class supported"
        aux_windex = work_array_windex_aux_list[0]
        aux_change_list = list( all_wvrr_aux_change_dict[ aux_windex.index().name() ] )
        # change = 0 automatically set in array_index_dict for integral_array objects
        if 0 in aux_change_list: aux_change_list.remove(0)
        for change in aux_change_list:
            work_array_index = \
                    dsl_scalar( name = 'iwrk'+str(change), vartype = 'int' )
            work_array_index_start = \
                    dsl_scalar( name = 'iwrk'+str(change)+'_start', vartype = 'int' )
            wintegral.work_array().add_array_index( change, work_array_index,\
                    work_array_index_start )
        # Add special_array_index for contraction "layer":
        # if a iwrk1 is already present due to auxiliary index algorithm, use that,
        # else create a new iwrk1 array_index tuple
        if not 1 in list( wintegral.work_array().array_index_dict().keys() ):
            # if iwrk1 not present, then create iwrk1
            work_array_index = \
                        dsl_scalar( name = 'iwrk1', vartype = 'int' )
            work_array_index_start = \
                    dsl_scalar( name = 'iwrk1_start', vartype = 'int' )
            wintegral.work_array().add_array_index( 1, work_array_index,\
                    work_array_index_start )
        # Add additional work_array array_indexes for accessing layers inside the
        # loops over primitives (add additional entry to end of tuple)
        for k, v in wintegral.work_array().array_index_dict().items() :
            work_array_index_prim = \
                dsl_scalar( name = 'iwrk'+str(k)+'_prim', vartype = 'int' )
            wintegral.work_array().array_index_dict()[ k ] = v + ( work_array_index_prim , )

        # Set work_array windex_list to integral indexes (ordered) +
        # auxiliary indexes (unordered)
        wintegral.work_array().set_windex_list( work_array_windex_index_list +\
                                            work_array_windex_aux_list )
### Create a specific algo_variable_setup derived class for this algorithm ###
# Integral argument list (use general class for contracted integrals)
integral_argument_setup = contracted_general.integral_argument_setup
        
# Automatically assigned local_variables
automatic_local_variable_setup = contracted_general.automatic_local_variable_setup

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
        auxiliary_index = algorithm_obj.wvrr_group_list[0].auxiliary_index
        # Add auxiliary index local variables:
        # Only a single auxiliary index is allowed per integral class, so set these
        # variables per algorithm object
        assert algorithm_obj.m_layers == 2,\
                "only 2-layer auxiliary index algorithm supported"
        windex = wintegral.windex_aux_dict()[ auxiliary_index.name() ]
        manual_set.add( windex.length_dict()['loop min'] ) # nloopm_min
        manual_set.add( windex.loop_length() )     # nloopm
        manual_set.add( windex.array_index() )     # im
        manual_set.add( windex.work_array_skip() ) # iwskipm

        for wvrr_group in algorithm_obj.wvrr_group_list:
            for wrr in wvrr_group.wrr_list:
                # changing/unrolled index
                windex = wintegral.windex_index_dict()[ wrr.changing_index().name() ] 
                # Not required for 1-layer auxiliary index algorithm
                if wvrr_group.m_layers != 1:
                    manual_set.add( windex.index_loop_max() )
        # work_array_offset set at the end of VRR sequence, depends on work_array_length
        for index in final_wrr_group.end_index_list:
            windex = wintegral.windex_index_dict()[ index.name() ]
            manual_set.add( windex.work_array_length() )

        ### Early transform sequence variables ###
        for index in final_wrr_group.end_index_list:
            windex = wintegral.windex_index_dict()[ index.name() ]
            for var in [ \
                        windex.skip_dict()[ 'prim' ],\
                        windex.array_index_dict()['prim'],\
                        windex.array_index_dict()['cont'],\
                       ] :
                manual_set.add( var )
        # Only need iskipX_tmp for indexes which are not the fastest changing index
        for index in final_wrr_group.end_index_list[:-1]:
            windex = wintegral.windex_index_dict()[ index.name() ]
            manual_set.add( windex.skip_dict()['tmp'] )
        # If spherical transformation is done, require nX_tmp and iskipX_tmp for all indexes
        # ... also iwskipX_cont and iskipX_cont_tmp for all indexes
        if self.gen.generator_options().spherical_transformed() == True:
            for index in final_wrr_group.end_order_list:
                windex = wintegral.windex_index_dict()[ index.name() ]
                for var in [\
                            windex.skip_dict()['tmp'], \
                            windex.length_dict()['tmp'], \
                            windex.skip_dict()['contwork'], \
                            windex.skip_dict()['conttmp'] \
                           ]:
                    manual_set.add( var )


        ### Further all algorithm variables ###
        # Auxiliary index present---require multiple array_index and array_index_start variables. 
        assert algorithm_obj.m_layers == 2, "Only 2-layer auxiliary index algorithm supported"
        for array_index_tuple in [ wintegral.work_array().array_index_dict()[0],\
                                   wintegral.work_array().array_index_dict()[1] \
                                 ]:
            for var in array_index_tuple:
                manual_set.add( var )
        # special_array_index_dict contains array_index for base function output, required for all
        # algorithms
        manual_set.add( wintegral.work_array().special_array_index_dict()['base'] )

        # Contractions are required: require multiple array_index and array_index_start variables, 
        # otherwise, only one array_index variable is required.
        for array_index_tuple in [ wintegral.work_array().array_index_dict()[0],\
                                   wintegral.work_array().array_index_dict()[1] \
                                 ]:
            for var in array_index_tuple:
                manual_set.add( var )

        # Append unique set elements to list for output
        for var in manual_set:
            manual_list.append(var)

### Create a specific work_array_size_source class for this algorithm ###
class algo_work_array_size_source(contracted_general.contracted_work_array_size_source):
    def serial_out(self,p):
        """
        For integral class evaluated using only VRRs, simply output the 
        use the work_array length_expr. 

        Contracted integrals:
        We require the largest of
            2 * (ncc0(la) * ncc0(lb) * ... ncc0(lX) * n_contraction + m_max + 1)
        where n_contraction is the product of all nX_prim values, and
            2 * (ncc0(la) * ncc0(lb) * ... ncc0(lX) * n_contraction )
        where n_contraction is the largest product of primitive and
        contraction lengths (see max_n_contraction_out).
        """
        
        f = self._f_dsl
        local_variable_list = self._local_variable_list
        
        # All memreq and n_contraction local variables required
        memreq1 = dsl_scalar(name = 'memreq1', vartype = 'int' )
        memreq2 = dsl_scalar(name = 'memreq2', vartype = 'int' )
        local_variable_list.append( memreq1 )
        local_variable_list.append( memreq2 )
        n_contraction1 = dsl_scalar( name = 'n_contraction1', vartype = 'int' )
        n_contraction2 = dsl_scalar( name = 'n_contraction2', vartype = 'int' )
        local_variable_list.append( n_contraction1 )
        local_variable_list.append( n_contraction2 )

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
        last_wvrr_group = self._wintegral.algorithm_obj().wvrr_group_list[-1]
        # nprima * nprimb * ... * nprimX
        self.n_contraction_out(p,n_contraction1,True) 
        memreq_expr = self.vrr_memreq_expr(last_wvrr_group,n_contraction1)
        # Output work array size required for VRR phase to memreq1
        p.out( memreq1.assign( memreq_expr ) )

        ### Determine maximum memory requirement for early transform phase ###
        self.max_n_contraction_out(p,n_contraction1,n_contraction2)
        memreq_expr = self.early_transform_memreq_expr(last_wvrr_group,n_contraction1)
        # Output work array size required for early transform phase to memreq2
        p.out( memreq2.assign( memreq_expr ) )

        ### Output if block to determine largest memreq value ###
        p.ifblock( dsl_binop( op_gt, memreq2, memreq1 ) )
        p.out( memreq1.assign( memreq2 ) )
        p.endifblock()
        
        ### Output return statement ###
        p.returnstmt( memreq1 )

        p.endfuncdef()


class algo_main_source(contracted_general.contracted_algo_main_source):
    def serial_out(self,p):
        """
        Output main integral evaluation code for the case of 
        - contracted
        - VRR only
        - 1 auxiliary index
        integrals requested.
        """
        f = self._f_dsl
        auto_local_variable_list = self._wintegral.auto_local_variable_list()
        auto_prim_loop_variable_list = self._wintegral.auto_prim_loop_variable_list()
        auto_prim_loop_depends_dict  = self._wintegral.auto_prim_loop_depends_dict()
        manual_local_variable_list = self._wintegral.manual_local_variable_list()

        algorithm = self._wintegral.algorithm_obj()

        # Output function definition top part
        p.funcdef( f, self.const_dict() )

        ### Local variable declarations ###
        var_to_declare_list = auto_local_variable_list +\
                                auto_prim_loop_variable_list +\
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
        early_transform_list = algorithm.early_transform_list

        ### Setup work_array indexing for base class evaluation and VRR sequence ###
        # Set work_array_length and work_array_skip values for sequence of wvrr_groups
        # These will define the shape and indexing of the end result of the sequence of VRRs,
        # which is completely characterised by the last wvrr_group.
        # These values are the same for all groups of primitive exponents for contratced integrals
        # -- they do not need to be evaluated inside the loops over primitives.
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
        # Assume a single auxiliary index, 2-layer algorithm
        assert algorithm.m_layers == 2,\
                "Only the 2-layer auxiliary index algorithm is currently supported."
        assert len( self._wintegral.windex_aux_list() ) == 1,\
                "Maximum of 1 auxiliary index currently supported."    
        # In this algorithm, auxiliary index is the last (slowest changing) index
        skip_index_list.append( final_wvrr_group.auxiliary_index )
        index  = skip_index_list[0]
        # First index cannot be auxiliary, so no need to check windex_aux_dict
        windex = self._wintegral.windex_index_dict()[ index.name() ]
        work_array_skip_var  = windex.work_array_skip() # skip for each Cartesian component
        work_array_skip_expr = 1
        p.out( work_array_skip_var.assign( work_array_skip_expr ) )
        work_array_prim_skip_var  = windex.skip_dict()['prim']
        work_array_prim_skip_expr = windex.work_array_length()
        p.out( work_array_prim_skip_var.assign( work_array_prim_skip_expr ) )
        is_auxiliary = False # first index is always non-auxiliary
        i = 1
        for index in skip_index_list[1:]:
            # Set per-Cartesian component skip
            # ..multiply by work_array_length of previous index (last skip output)
            if work_array_skip_expr == 1:
                work_array_skip_expr = windex.work_array_length() 
            else:
                work_array_skip_expr = work_array_skip_expr * windex.work_array_length()
            # If contracted integrals are required, multiply the work_array_skip_expr by
            # the number or primitives for that index (if previous index is not auxiliary)
            # ..multiply by nX_prim of previous index (last skip output)
            if work_array_skip_expr == 1:
                work_array_skip_expr = windex.length_dict()['prim']
            else:
                work_array_skip_expr = work_array_skip_expr * windex.length_dict()['prim']
            # Determine if current index is auxiliary
            if index.name() in list( self._wintegral.windex_aux_dict().keys() ):
                is_auxiliary = True
            else:
                is_auxiliary = False
                assert index.name() in list( self._wintegral.windex_index_dict().keys() ),\
                        "index must be in windex_index_dict or windex_aux_dict"
            # Set windex for current index in loop
            if is_auxiliary == False:
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                # Output per-Cartesian component skip
                work_array_skip_var  = windex.work_array_skip()
                p.out( work_array_skip_var.assign( work_array_skip_expr ) )
                # Output per-primitive skip
                work_array_prim_skip_var = windex.skip_dict()['prim']
                work_array_prim_skip_expr = work_array_skip_expr * windex.work_array_length()
                p.out( work_array_prim_skip_var.assign( work_array_prim_skip_expr ) )
            elif is_auxiliary == True:
                assert i+1 == len( skip_index_list ),\
                        "auxiliary_index should always be last in skip_index_list"
                windex = self._wintegral.windex_aux_dict()[ index.name() ]
                # Output per-Cartesian component skip
                work_array_skip_var  = windex.work_array_skip()
                p.out( work_array_skip_var.assign( work_array_skip_expr ) )
            i += 1

        ### Setup work_array indexing ###
        # Set iwrk_start values to start of each m-layer/HRR layer
        work_array_index_dict = self._wintegral.work_array().array_index_dict()
        special_work_array_index_dict = self._wintegral.work_array().special_array_index_dict()
        work_array_index0, work_array_index_start0, work_array_index_prim0 = \
            work_array_index_dict[0]
        work_array_index1, work_array_index_start1, work_array_index_prim1 = \
            work_array_index_dict[1]
        work_array_index_base = \
                special_work_array_index_dict['base']
        # Half work array size is where iwrk_start1 should be set
        half_work_array_size_function = self._wintegral._work_array_size_function.f_dsl() / 2
        # For contracted integrals, the iwrk array positions depend on the current primitive 
        # exponent index values
        windex = self._wintegral.windex_index_dict()[ final_wvrr_group.end_order_list[0].name() ]
        work_array_index0_initial_expr = work_array_index_prim0 +\
                windex.array_index_dict()['prim'] * windex.skip_dict()['prim']
        work_array_index1_initial_expr = work_array_index_prim1 +\
                windex.array_index_dict()['prim'] * windex.skip_dict()['prim']
        for index in final_wvrr_group.end_order_list[1:]:
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            work_array_index0_initial_expr = work_array_index0_initial_expr + \
                                        windex.array_index_dict()['prim'] * windex.skip_dict()['prim']
            work_array_index1_initial_expr = work_array_index1_initial_expr + \
                                        windex.array_index_dict()['prim'] * windex.skip_dict()['prim']
        # Simply add work_array_index_start1 to arrive at work_array_index1 start position
        # ...use a 2-layer algorithm, where the work_array is divided into two equal halves
        # If auxiliary index or HRR present, require work_array_index1 to be set
        # Output work_array_index_start values
        p.out( work_array_index_start0.assign( 0 ) )
        p.out( work_array_index_start1.assign( half_work_array_size_function ) )
        # Output assignment for special base function work_array array_index
        # In the 2-layer reduced memory algorithm for auxiliary indexes
        # base function output is stored immediately after the second layer ends
        # and accessed as needed inside loops over auxiliary index values
        assert max( list( work_array_index_dict.keys() ) ) <= 1,\
                "Only increments of 0 and 1 supported in auxiliary index"
        windex = self._wintegral.windex_aux_list()[0]
        p.out( work_array_index_base.assign( \
                work_array_index_start1 + windex.work_array_skip() ) )
        loop_length_min_var = windex.length_dict()['loop min']
        loop_length_var  = windex.loop_length()
        loop_length_min_expr = algorithm.m_min
        loop_length_expr     = algorithm.m_max
        p.out( loop_length_min_var.assign( loop_length_min_expr ) )
        p.out( loop_length_var.assign( loop_length_expr ) )

        assert len( algorithm.early_transform_list ) > 0,\
                "for contracted integrals, the algorithm object must contain contraction "+\
                "objects."

        ### VRR sequence ###
        self.vrr_sequence_with_auxiliary_out( p, wvrr_group_list,\
            auto_prim_loop_variable_list,auto_prim_loop_depends_dict,\
            work_array_index0, work_array_index1, work_array_index_base, \
            work_array_index_start0, work_array_index_start1,\
            work_array_index_prim0, work_array_index_prim1,\
            work_array_index0_initial_expr, work_array_index1_initial_expr )

        # Set state of work_array dimensions (offsets) after VRR sequence
        # for integral indexes that have been incremented.
        # wvrr_groups currently always end with 'ncc0' type indexes, and this is 
        # assumed here.
        # We only need to set work_array_offset values for the indexes in loop_index_list
        # of the first contract operation, since after this, there are no offsets (all
        # Cartesian component ranges are simply 'ncc' type).
        for index in early_transform_list[0].start_index_list :
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            ncc_type = early_transform_list[0].start_ncc_type_dict[ index.name() ]
            work_array_length_var   = windex.work_array_length()
            work_array_offset_var   = windex.work_array_offset()
            output_array_length_var = windex.output_array_length() 
            # Determine work_array_length_expr and work_array_offset_expr 
            # from ncc_type
            assert ncc_type == 'ncc0',\
                    "VRRs must end with ncc_type == 'ncc0' for all indexes"
            l_min_expr = early_transform_list[0].end_min_l_expr_dict[ index.name() ]
            work_array_offset_expr = l_min_expr * (l_min_expr + 1) *(l_min_expr + 2)/6
            # Output assignments for offset only
            p.out( work_array_offset_var.assign( work_array_offset_expr ) )

        ### Early transform sequence ###
        # work_array_index0 is where the final outcome of the VRR sequence is, and must be 
        # swapped to the other work_array_index_start{0,1} "layer" so the result of 
        # contractions can be placed there.
        # ... we need to determine which work_array_index_start{0,1} "layer" corresponds to
        # these array_index values at the end of the VRR sequence, and this can only be known
        # at runtime.
        # Use an if statement -- at the end of the VRR sequence the largest of work_array_index0
        # and work_array_index1 corresponds to work_array_index_start1, and the smallest to
        # work_array_index_start0, so at the beginning of the contraction sequence, these are
        # swapped.
        p.ifblock( dsl_binop( op_gt, work_array_index0, work_array_index1 ) ) # then
        # Set work_array_index for this contraction
        p.out( work_array_index0.assign( work_array_index_start0 ) )
        p.out( work_array_index1.assign( work_array_index_start1 ) )
        # work_array_index start already set for this contraction
        p.elseblock()
        # Set work_array_index for this contraction
        p.out( work_array_index0.assign( work_array_index_start1 ) )
        p.out( work_array_index1.assign( work_array_index_start0 ) )
        # Set work_array_index_start for this contraction
        p.out( work_array_index_start0.assign( work_array_index0 ) )
        p.out( work_array_index_start1.assign( work_array_index1 ) )
        p.endifblock()
        # Use contracted_general routine
        self.early_transform_sequence_out( p, early_transform_list, \
                work_array_index0, work_array_index_start0, \
                work_array_index1, work_array_index_start1, \
                last_transform_swap = False )

        ### Copy correct integral intermediates from intermediate array ###
        final_operation = early_transform_list[-1]
        self.copy_to_output_array_out(p,final_operation,work_array_index_start0)
        
        # End of function definition
        p.endfuncdef()
