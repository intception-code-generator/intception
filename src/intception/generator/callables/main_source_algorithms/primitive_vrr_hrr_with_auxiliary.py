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
    * VRR and HRR algorithm (2-layer algorithm)
    * Auxiliary index (2-layer algorithm)
"""
from intception.dsl import *
from intception.generator.callables import integral_wrapper_main_function
from intception.generator.callables import integral_wrapper_work_array_size_function
from intception.generator.algo_descriptor import algo_descriptor
from intception.generator.algo_variable_setup import algo_variable_setup
from intception.generator.algo_indexing_setup import algo_indexing_setup
from . import primitive_general

### Create a specific algo_descriptor object for this algorithm ###
algo_desc = algo_descriptor( is_contracted = False, hrr_present = True, aux_index_present = True )

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
        whrr_list = []
        # Separate vrr_wrapper and hrr_wrapper objects
        hrr_found = False
        for wrr in wintegral.wrr_list():
            assert wrr.rr().rrtype() in [ 'vrr', 'hrr' ], \
                "Only VRR and HRR-type RRs supported currently"
            if wrr.rr().rrtype() == 'vrr':
                assert hrr_found == False, "VRRs may not follow HRRs in algorithm"
                wvrr_list.append(wrr)
                windex = wintegral.windex_index_dict()[ wrr.changing_index().name() ]
                work_array_windex_index_list.append( windex )
            elif wrr.rr().rrtype() == 'hrr':
                hrr_found = True
                whrr_list.append(wrr)
                windex = wintegral.windex_index_dict()[ wrr.move_to_index().name() ]
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
        # Add special_array_index for HRR "layer":
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

        # Set work_array windex_list to integral indexes (ordered) +
        # auxiliary indexes (unordered)
        wintegral.work_array().set_windex_list( work_array_windex_index_list +\
                                            work_array_windex_aux_list )

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
        # At least 1 HRR, last wrr_group is a whrr_group object
        final_wrr_group = algorithm_obj.whrr_group_list[-1]
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
                manual_set.add( windex.index_loop_max() )
                # loop_index_list indexes
                for index in wrr.loop_index_list():
                    windex = wintegral.windex_index_dict()[ index.name() ]
                    manual_set.add( windex.loop_length() )
        final_wvrr_group = algorithm_obj.wvrr_group_list[-1]
        # work_array_offset set at the end of VRR sequence, depends on work_array_length
        for index in final_wvrr_group.end_index_list:
            windex = wintegral.windex_index_dict()[ index.name() ]
            manual_set.add( windex.work_array_length() )

        ### HRR sequence variables ###
        for whrr_group in algorithm_obj.whrr_group_list:
            # move_to_index
            windex = wintegral.windex_index_dict()[ whrr_group.wrr.move_to_index().name() ]
            manual_set.add( windex.index_loop_max() )
            if len( whrr_group.end_index_list ) > 2:
                # move_to_index loop_length only required if more than 2 indexes
                # as it is used in setting work_array_skip for other indexes
                manual_set.add( windex.loop_length() )
            # move_from_index
            windex = wintegral.windex_index_dict()[ whrr_group.wrr.move_from_index().name() ]
            manual_set.add( windex.loop_length() )
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
                manual_set.add( windex.loop_length() )
            # work_array_offset set at end of HRR operation, may depend on work_array_length
            for index in whrr_group.end_index_list:
                if whrr_group.end_ncc_type_dict[ index.name() ] == 'ncc0':
                    # ncc0 type index -- requires work_array_length to be set
                    windex = wintegral.windex_index_dict()[ index.name() ]
                    manual_set.add( windex.work_array_length() )

        # If spherical transformation is done, require nX_tmp and iskipX_tmp for all indexes
        # ... also iwskipX_cont and iskipX_cont_tmp for all indexes
        # ... also work_array_length for all indexes 
        # [ Not currently supported ]
        if self.gen.generator_options().spherical_transformed() == True:
            assert self.gen.generator_options().spherical_transformed() == False,\
                    "Spherical transformation for primitive integrals not supported."
        #    for index in final_wrr_group.end_order_list:
        #        windex = wintegral.windex_index_dict()[ index.name() ]
        #        for var in [\
        #                    windex.skip_dict()['tmp'], \
        #                    windex.length_dict()['tmp'], \
        #                    windex.work_array_length()
        #                   ]:
        #            manual_set.add( var )
        #    # If the first spherical transformation precedes the HRR and does not
        #    # occur on the fastest changing index, need iX_tmp for that index
        #    # being transformed, also.
        #    if len( algorithm_obj.early_transform_list ) > 0:
        #        t = algorithm_obj.early_transform_list[0]
        #        assert t.optype == 'spherical transform', \
        #                "only spherical_transform operations allowed in "+\
        #                "early_transform_list for primitive integral types."
        #        if not t.index_to_transform is t.end_order_list[-1]:
        #            # Not fastest-changing index -- require special first
        #            # transformation to truncate angular momentum ranges
        #            windex = wintegral.windex_index_dict()[ t.index_to_transform.name() ]
        #            manual_set.add( windex.array_index_dict()['tmp'] )




        ### Further all algorithm variables ###
        # If auxiliary index present, or a HRR is used, or contractions are required
        # require multiple array_index and 
        # array_index_start variables, otherwise, only one array_index variable is required.
        if final_wvrr_group.auxiliary_index != None:
            assert algorithm_obj.m_layers == 2, "Only 2-layer auxiliary index algorithm supported"
        for array_index, array_index_start in [ wintegral.work_array().array_index_dict()[0],\
                                                wintegral.work_array().array_index_dict()[1] ]:
            manual_set.add( array_index )
            manual_set.add( array_index_start )
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
        The largest the work_array will ever need to be is the largest of, and
            2 * ncc0(lb) * ... ncc0(lX) + m_max + 1
        where b, ... X are incremented at the end of the VRR sequence, and
            2 * hrr_layer[ la ][ lb ] * ... ncc0(lX)
        where angular momentum is moved from index b to index a, and b, ... X
        were previously incremented in the VRR sequence.
        """
        f = self._f_dsl
        local_variable_list = self._local_variable_list
        
        # Add local variables to store and compare integers
        memreq1 = dsl_scalar(name = 'memreq1', vartype = 'int' )
        memreq2 = dsl_scalar(name = 'memreq2', vartype = 'int' )
        local_variable_list.append( memreq1 )
        local_variable_list.append( memreq2 )

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
        m_layers       = self._wintegral.algorithm_obj().m_layers
        assert m_layers == 2, "only 2-layer auxiliary index algorithm suppored"
        assert len( whrr_group_list ) == 1,\
            "Only a single HRR per integral class currently supported."
        # Determine maximum memory requirement for VRR phase and VRR to HRR 
        # transition
        # Compare work_array requirement for VRR-only and VRR-to-HRR
        # steps---auxiliary index is required for primitive integrals
        #### Determine maximum memory requirement for VRR phase ###
        memreq_expr = self.vrr_memreq_expr(wvrr_group_list[-1] )
        p.out( memreq1.assign( memreq_expr ) )
        
        ### Determine maximum memory requirement for HRR phase ###
        whrr_group = whrr_group_list[0]
        index_ncc_dict = whrr_group.end_ncc_type_dict
        memreq_expr = self.hrr_memreq_expr( whrr_group, index_ncc_dict )
        p.out( memreq2.assign( memreq_expr ) )

        ### Set memreq1 to largest value of memreq1,2 ###
        p.ifblock( dsl_binop( op_gt, memreq2, memreq1 ) )
        p.out( memreq1.assign( memreq2 ) )
        p.endifblock()

        ### Output return statement ###
        p.returnstmt( memreq1 )

        p.endfuncdef()

class algo_main_source(primitive_general.primitive_algo_main_source):
    def serial_out(self,p):
        """
        Output main integral evaluation code for the case of 
        - primitive
        - VRR and HRR (2-layer algorithm)
        - with auxiliary index (2-layer algorithm)
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
        var_to_declare_list = auto_local_variable_list +\
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
        # [ Spherical transformation not supported currently ]
        assert len( early_transform_list ) == 0,\
                "Spherical transformation for primitive integrals not supported."
        assert len( late_transform_list ) == 0,\
                "Spherical transformation for primitive integrals not supported."

        # Set final_wrr_group (contains information on the shape of the array of intermediates
        # immediately before these are copied to output array)
        if len( whrr_group_list ) > 0:
            final_wrr_group = whrr_group_list[-1]
        else:
            final_wrr_group = wvrr_group_list[-1]

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
        incremented_index_list = list( reversed( final_wvrr_group.end_index_list[:] ))
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
        if windex.index() in incremented_index_list:
            p.out( work_array_skip_var.assign( work_array_skip_expr ) )
        is_auxiliary = False # first index is always non-auxiliary
        for index in skip_index_list[1:]:
            if windex.index() in incremented_index_list:
                # Set per-Cartesian component skip
                # ..multiply by work_array_length of previous index (last skip output)
                if work_array_skip_expr == 1:
                    work_array_skip_expr = windex.work_array_length() 
                else:
                    work_array_skip_expr = work_array_skip_expr * windex.work_array_length()
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
            elif is_auxiliary == True:
                windex = self._wintegral.windex_aux_dict()[ index.name() ]
            # Output per-Cartesian component skip (if current index is incremented)
            if windex.index() in incremented_index_list or is_auxiliary == True:
                work_array_skip_var  = windex.work_array_skip()
                p.out( work_array_skip_var.assign( work_array_skip_expr ) )
                
        # Special base function work_array array_index dsl_scalar object
        work_array_index_base = \
                self._wintegral.work_array().special_array_index_dict()['base']
        # The initial expression for work_array_index0 is the same, regardless of whether an auxiliary index
        # is involved
        # For primitive integrals, the iwrk0 array position is always zero initially
        work_array_index0_initial_expr = 0
        # If an auxiliary index is involved in the VRR sequence, or at least 1 HRR is present,
        # or contractions are required...
        # ...use a 2-layer algorithm, where the work_array is divided into two equal halves
        # If auxiliary index or HRR present, require work_array_index1 to be set
        # Set iwrk_start values to start of each m-layer/HRR layer
        work_array_index_dict = \
            self._wintegral.work_array().array_index_dict()
        work_array_index0, work_array_index_start0 = \
            self._wintegral.work_array().array_index_dict()[0]
        work_array_index1, work_array_index_start1 = \
            self._wintegral.work_array().array_index_dict()[1]
        p.out( work_array_index0.assign( work_array_index0_initial_expr ) )

        # Determine half the size of the work_array (mid point between two "layers")
        half_work_array_size_function = self._wintegral._work_array_size_function.f_dsl() / 2
        p.out( work_array_index1.assign( work_array_index0 + half_work_array_size_function ) )
        # Output assignment for special base function work_array array_index
        # In the 2-layer reduced memory algorithm for auxiliary indexes
        # base function output is stored immediately after the second layer ends
        # and accessed as needed inside loops over auxiliary index values
        assert max( list( work_array_index_dict.keys() ) ) <= 1,\
                "Only increments of 0 and 1 supported in auxiliary index"
        windex = self._wintegral.windex_aux_list()[0]
        p.out( work_array_index_base.assign( \
                work_array_index1 + windex.work_array_skip() ) )
        loop_length_min_var = windex.length_dict()['loop min']
        loop_length_var  = windex.loop_length()
        # Note: The implied auxiliary index (m_layers == 1) algorithm for VRR evaluation
        # has not been implemented for primitive integral types, so loop_length_min_expr == 0 always.
        #loop_length_min_expr = algorithm.m_min
        loop_length_min_expr = 0
        loop_length_expr     = algorithm.m_max
        p.out( loop_length_min_var.assign( loop_length_min_expr ) )
        p.out( loop_length_var.assign( loop_length_expr ) )

        ### VRR sequence (and base evaluation) ###
        self.vrr_sequence_with_auxiliary_out( p, wvrr_group_list,\
            work_array_index0, work_array_index1, work_array_index_base, \
            work_array_index_start0, work_array_index_start1, \
            work_array_index0_initial_expr )

        # Set state of work_array dimensions (offsets) after VRR sequence
        # for integral indexes that have been incremented.
        # The work_array_length variables do not need to be set, since they
        # do not change during the VRR sequence and are set prior to this.
        # wvrr_groups currently always end with 'ncc0' type indexes, and this is 
        # assumed here.
        index_offset_list = final_wvrr_group.end_index_list
        for index in index_offset_list :
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            work_array_length_var   = windex.work_array_length()
            work_array_offset_var   = windex.work_array_offset()
            output_array_length_var = windex.output_array_length() 
            # Incremented by VRR sequence
            # Determine work_array_length_expr and work_array_offset_expr 
            # from ncc_type
            ncc_type = final_wvrr_group.end_ncc_type_dict[ index.name() ]
            assert ncc_type == 'ncc0',\
                    "VRRs must end with ncc_type == 'ncc0' for all indexes"
            l_min_expr = windex.index_value()
            work_array_offset_expr = l_min_expr * (l_min_expr + 1) *(l_min_expr + 2)/6
            # Output assignments for offset only
            p.out( work_array_offset_var.assign( work_array_offset_expr ) )

        ### HRRs ###
        self.hrr_sequence_out(p,whrr_group_list,\
                    work_array_index0, work_array_index_start0, \
                    work_array_index1, work_array_index_start1, \
                    initial_hrr_swap = False )

        # Unswap layer final layer if auxiliary index algorithm of HRR sequence used
        assert algorithm.m_layers == 2,\
                   "Only 2-layer VRR + auxiliary index currently supported"
        p.out( work_array_index0.assign( work_array_index_start0 ) )

        ### Copy correct integral intermediates from intermediate array ###
        final_operation = whrr_group_list[-1]
        self.copy_to_output_array_out(p,final_operation,work_array_index0)
        
        # End of function definition
        p.endfuncdef()

