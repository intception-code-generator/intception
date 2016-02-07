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
    * VRR and HRR algorithm (2-layers)
    * no auxiliary indexes
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
algo_desc = algo_descriptor( is_contracted = True, hrr_present = True, aux_index_present = False )

### Create a specific algo_indexing_setup derived class for this algorithm ###
class integral_array_indexing_setup(algo_indexing_setup):
    def __call__(self):
        """
        Setup indexing attributes for integral arrays attached to integral_wrapper
        based on specific algorithm details (sets up integral_wrapper.work_array).
        """
        wintegral = self.wintegral

        work_array_windex_index_list = []
        # Go through RR list to determine order of indexes in work array 
        work_array_windex_index_list = []
        hrr_found = False
        for wrr in wintegral.wrr_list():
            assert wrr.rr().rrtype() in [ 'vrr', 'hrr' ], \
                "Only VRR and HRR-type RRs supported currently"
            if wrr.rr().rrtype() == 'vrr':
                assert hrr_found == False, "VRRs may not follow HRRs in algorithm"
                windex = wintegral.windex_index_dict()[ wrr.changing_index().name() ]
                work_array_windex_index_list.append( windex )
            elif wrr.rr().rrtype() == 'hrr':
                windex = wintegral.windex_index_dict()[ wrr.move_to_index().name() ]
                work_array_windex_index_list.append( windex )

        # Any additional integral indexes not captured from wrr_list
        for windex in wintegral.windex_index_list():
            assert windex in work_array_windex_index_list,\
                    "Currently, all integral indexes must be incremented in an RR"
        assert len( wintegral.windex_aux_list() ) == 0, \
                "No auxiliary indexes in this algorithm (contracted, VRR-only, no auxiliary)"
        # Add array_index for placement of base function output in work_array
        # The following code assumes only a single auxiliary index
        work_array_index_base = \
                        dsl_scalar( name = 'iwrkb', vartype = 'int' )
        wintegral.work_array().add_special_array_index('base',work_array_index_base )
        # Add special_array_index for contraction/HRR "layer":
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

        # Set work_array windex_list to integral indexes
        wintegral.work_array().set_windex_list( work_array_windex_index_list )

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
        # HRR present, last wrr_group is a whrr_group object
        # Only a single HRR operation supported currently
        assert len( algorithm_obj.whrr_group_list ) == 1,\
                "Only a single HRR operation per integral class supported."
        final_wrr_group = algorithm_obj.whrr_group_list[-1]
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
                #X manual_set.add( windex.index_loop_max() )
                # loop_index_list indexes
                for index in wrr.loop_index_list():
                    windex = wintegral.windex_index_dict()[ index.name() ]
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


        ### HRR sequence variables ###
        for whrr_group in algorithm_obj.whrr_group_list:
            # move_to_index
            windex = wintegral.windex_index_dict()[ whrr_group.wrr.move_to_index().name() ]
            manual_set.add( windex.index_loop_max() )
            manual_set.add( windex.skip_dict()['conttmp'] )
            manual_set.add( windex.skip_dict()['contwork'] )
            # move_from_index
            windex = wintegral.windex_index_dict()[ whrr_group.wrr.move_from_index().name() ]
            manual_set.add( windex.skip_dict()['conttmp'] )
            manual_set.add( windex.skip_dict()['contwork'] )
            # all indexes looped over
            for index in whrr_group.end_index_list:
                windex = wintegral.windex_index_dict()[ index.name() ]
                skip_hrr_var = windex.skip_dict()['tmp']
                manual_set.add( skip_hrr_var )
            # indexes looped over that are not move_{from,to}_index
            loop_index_list = whrr_group.end_index_list[:]
            loop_index_list.remove( whrr_group.wrr.move_from_index() )
            loop_index_list.remove( whrr_group.wrr.move_to_index() )
            for index in loop_index_list:
                windex = wintegral.windex_index_dict()[ index.name() ]
            # work_array_offset set at end of HRR operation, may depend on work_array_length
            #for index in whrr_group.end_index_list:
            #    if whrr_group.end_ncc_type_dict[ index.name() ] == 'ncc0':
            #        # ncc0 type index -- requires work_array_length to be set
            #        windex = wintegral.windex_index_dict()[ index.name() ]
            #        manual_set.add( windex.work_array_length() )

        ### Further all algorithm variables ###
        # HRR and contractions are required: require multiple array_index and array_index_start variables.
        for array_index, array_index_start in [ wintegral.work_array().array_index_dict()[0],\
                                                wintegral.work_array().array_index_dict()[1] ]:
            manual_set.add( array_index )
            manual_set.add( array_index_start )

        # Append unique set elements to list for output
        for var in manual_set:
            manual_list.append(var)

### Create a specific work_array_size_source class for this algorithm ###
class algo_work_array_size_source(contracted_general.contracted_work_array_size_source):
    def serial_out(self,p):
        """
        For integral class using this algorithm, output a simple
        function that returns an integer number for the minimum 
        required size of the work_array.
        
        Contracted integrals:
        The largest the work_array will ever need to be is the largest of,
            2 * ncc0(lb) * ... ncc0(lX) * n_contraction
        where b, ... X are incremented at the end of the VRR sequence, and
        n_contraction is the largest product of n_prim and n_cont lengths
        that occur in the contraction sequence, and
            2 * hrr_layer[ la ][ lb ] * ... ncc0(lX) * n_contraction
        where angular momentum is moved from index b to index a, and b, ... X
        were previously incremented in the VRR sequence, and where 
        n_contraction is the product of all nX_cont values.
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

        wvrr_group_list = self._wintegral.algorithm_obj().wvrr_group_list
        whrr_group_list = self._wintegral.algorithm_obj().whrr_group_list
        ### Determine maximum memory requirement for early transform phase ###
        # No auxiliary index present
        # VRR to contraction transition will always require more memory than VRR
        # phase (it is 2x the size of VRR memory requirement if no
        # auxiliary index is required).
        self.max_n_contraction_out(p,n_contraction1,n_contraction2)
        memreq_expr = self.early_transform_memreq_expr(wvrr_group_list[-1],n_contraction1)
        # Output work array size required for early transform phase to memreq2
        p.out( memreq1.assign( memreq_expr ) )
        
        ### Determine maximum memory requirement for HRR phase ###
        assert len( whrr_group_list ) == 1,\
            "Only a single HRR per integral class currently supported."
        whrr_group = whrr_group_list[0]
        index_ncc_dict = whrr_group.end_ncc_type_dict
        self.n_contraction_out(p,n_contraction1,False) 
        memreq_expr = self.hrr_memreq_expr( whrr_group, index_ncc_dict, n_contraction1 )
        p.out( memreq2.assign( memreq_expr ) )

        ### Set memreq1 to largest value of memreq1,2 ###
        p.ifblock( dsl_binop( op_gt, memreq2, memreq1 ) )
        p.out( memreq1.assign( memreq2 ) )
        p.endifblock()

        ### Determine memory requirement for late transform phase ###
        assert len( whrr_group_list ) == 1,\
            "Only a single HRR per integral class currently supported."
        whrr_group = whrr_group_list[0]
        memreq_expr = self.late_transform_memreq_expr( whrr_group, n_contraction1 )
        
        p.out( memreq2.assign( memreq_expr ) )

        ### Set memreq1 to largest value of memreq1,2 ###
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
        - VRR + HRR 
        - no auxiliary index
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
        whrr_group_list = algorithm.whrr_group_list
        early_transform_list = algorithm.early_transform_list
        late_transform_list  = algorithm.late_transform_list

        ### Setup work_array indexing for base class evaluation and VRR sequence ###
        # Set work_array_length and work_array_skip values for sequence of wvrr_groups
        # These will define the shape and indexing of the end result of the sequence of VRRs,
        # which is completely characterised by the last wvrr_group.
        # These values are the same for all groups of primitive exponents for contratced integrals
        # -- they do not need to be evaluated inside the loops over primitives.
        final_wvrr_group = wvrr_group_list[-1]
        # Assign work_array_length expressions for VRR operations
        for index in final_wvrr_group.end_order_list: # all integral indexes
            windex   = self._wintegral.windex_index_dict()[ index.name() ]
            try:
                ncc_type = final_wvrr_group.end_ncc_type_dict[ index.name() ]
                assert ncc_type == 'ncc0', \
                        "VRR operations should always result in a ncc_type of 'ncc0'"
                l_expr   = final_wvrr_group.end_max_l_expr_dict[ index.name() ]
                work_array_length_expr = ( l_expr + 1 )*( l_expr + 2 )*( l_expr + 3 )/6
            except KeyError:
                # We expect KeyErrors, where an index is not incremented in the VRR 
                # sequence. In this case, set work_array_length_var to 1
                ncc_type = None
                work_array_length_expr = 1
            work_array_length_var  = windex.work_array_length()
            p.out( work_array_length_var.assign( work_array_length_expr ) )

        # Assign work_array_skip expressions for per-Cartesian-component and
        # per-primitive component skips for VRR sequence.
        # Order of indexes does not change, nor does indexing of array during
        # VRR sequence, so these can be set from information in the final_wvrr_group
        skip_index_list = list( reversed( final_wvrr_group.end_order_list[:] ) )
        index  = skip_index_list[0]
        # First index cannot be auxiliary, so no need to check windex_aux_dict
        windex = self._wintegral.windex_index_dict()[ index.name() ]
        work_array_skip_var  = windex.work_array_skip() # skip for each Cartesian component
        work_array_skip_expr = 1
        p.out( work_array_skip_var.assign( work_array_skip_expr ) )
        work_array_prim_skip_var  = windex.skip_dict()['prim']
        work_array_prim_skip_expr = windex.work_array_length()
        p.out( work_array_prim_skip_var.assign( work_array_prim_skip_expr ) )
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
            # Set windex for current index in loop
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            # Output per-Cartesian component skip (if current index is incremented)
            work_array_skip_var  = windex.work_array_skip()
            p.out( work_array_skip_var.assign( work_array_skip_expr ) )
            # Output per-primitive skip (if contracted integrals are required)
            work_array_prim_skip_var = windex.skip_dict()['prim']
            work_array_prim_skip_expr = work_array_skip_expr * windex.work_array_length()
            p.out( work_array_prim_skip_var.assign( work_array_prim_skip_expr ) )

        ### Setup work_array indexing ###
        # Set iwrk_start values to start of each m-layer/HRR layer
        work_array_index_dict = \
            self._wintegral.work_array().array_index_dict()
        work_array_index0, work_array_index_start0 = \
            self._wintegral.work_array().array_index_dict()[0]
        work_array_index1, work_array_index_start1 = \
            self._wintegral.work_array().array_index_dict()[1]
        # Half work array size is where iwrk_start1 should be set
        half_work_array_size_function = self._wintegral._work_array_size_function.f_dsl() / 2
        # For contracted integrals, the iwrk array positions depend on the current primitive 
        # exponent index values
        windex = self._wintegral.windex_index_dict()[ final_wvrr_group.end_order_list[0].name() ]
        work_array_index0_initial_expr = work_array_index_start1 +\
                windex.array_index_dict()['prim'] * windex.skip_dict()['prim']
        for index in final_wvrr_group.end_order_list[1:]:
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            work_array_index0_initial_expr = work_array_index0_initial_expr + \
                                        windex.array_index_dict()['prim'] * windex.skip_dict()['prim']
        # ...use a 2-layer algorithm, where the work_array is divided into two equal halves
        # If auxiliary index or HRR present, require work_array_index1 to be set
        # Output work_array_index_start values
        p.out( work_array_index_start0.assign( 0 ) )
        p.out( work_array_index_start1.assign( half_work_array_size_function ) )

        ### VRR sequence ###
        self.vrr_sequence_no_auxiliary_out(p,wvrr_group_list,\
            auto_prim_loop_variable_list,auto_prim_loop_depends_dict,\
            work_array_index0, work_array_index1,\
            work_array_index0_initial_expr)
        
        # Set state of work_array dimensions (offsets) after VRR sequence
        # wvrr_groups currently always end with 'ncc0' type indexes, and this is 
        # assumed here.
        # We only need to set work_array_offset values for the indexes in loop_index_list
        # of the first contract operation, since after this, there are no offsets (all
        # Cartesian component ranges are simply 'ncc' type).
        for index in early_transform_list[0].start_order_list :
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            work_array_length_var   = windex.work_array_length()
            work_array_offset_var   = windex.work_array_offset()
            output_array_length_var = windex.output_array_length() 
            if index in early_transform_list[0].start_index_list:
                # Incremented by VRR sequence
                # Determine work_array_length_expr and work_array_offset_expr 
                # from ncc_type
                ncc_type = early_transform_list[0].start_ncc_type_dict[ index.name() ]
                assert ncc_type == 'ncc0',\
                        "VRRs must end with ncc_type == 'ncc0' for all indexes"
                l_min_expr = early_transform_list[0].end_min_l_expr_dict[ index.name() ]
                work_array_offset_expr = l_min_expr * (l_min_expr + 1) *(l_min_expr + 2)/6
            else:
                # Not yet incremented
                work_array_offset_expr = 0
            # Output assignments for offset only
            p.out( work_array_offset_var.assign( work_array_offset_expr ) )

        ### Early transform sequence ###
        # Set work_array_index for transformation sequence 
        # (layers pre-set before VRR)
        p.out( work_array_index0.assign( work_array_index_start0 ) )
        p.out( work_array_index1.assign( work_array_index_start1 ) )
        # Use contracted_general routine
        self.early_transform_sequence_out( p, early_transform_list, \
                work_array_index0, work_array_index_start0, \
                work_array_index1, work_array_index_start1, \
                last_transform_swap = False )

        ### HRRs ###
        if len( late_transform_list ) > 0:
            no_hrr_swap = True
        else:
            no_hrr_swap = False
        self.hrr_sequence_out( p, whrr_group_list, \
                work_array_index0, work_array_index_start0, \
                work_array_index1, work_array_index_start1, \
                no_hrr_swap = no_hrr_swap )

        ### Late transformation sequence ###
        # Use contracted_general routine
        if len( late_transform_list ) > 0:
            self.late_transform_sequence_out( p, late_transform_list, \
                    work_array_index0, work_array_index_start0, \
                    work_array_index1, work_array_index_start1, \
                    last_transform_swap = False )

        ### Copy correct integral intermediates from intermediate array ###
        if len( late_transform_list ) > 0:
            final_operation = late_transform_list[-1]
        else:
            final_operation = whrr_group_list[-1]
        self.copy_to_output_array_out(p,final_operation,work_array_index_start0)

        # End of function definition
        p.endfuncdef()
