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
Provides general classes and methods for use in specific primitive algorithm
modules.
"""

from intception.dsl import *
from intception.dsl_extensions import *
from intception.generator.supporting_functions import expr_recursive_search, determine_assignment_order
from intception.generator.algo_variable_setup import algo_variable_setup
from intception.generator.callables import integral_wrapper_main_function
from intception.generator.callables import integral_wrapper_work_array_size_function

class integral_argument_setup(algo_variable_setup):
    def __call__(self,arg_list):
        """
        Setup arg_list of integral function arguments (called in 
        integral_wrapper.setup_integral_argument_list) for primitive
        integrals.
        
        arg_list:       a Python list, to be populated with dsl_* objects
                        that will be the arguments for the main function
                        for an integral class.
        """
        wintegral = self.wintegral
        # For each integral index, add appropriate arguments (these are attachments)
        for windex in wintegral.windex_index_list():
            # Some indexes will have additional attachments (e.g. centre, exponent).
            for attachment in windex.index().attachment_list():
                arg_list.append( attachment )
            # Every index will have a corresponding value (e.g. l for Cartesian Gaussians).
            arg_list.append( windex.index_value() )
            # Every index will have a corresponding length (e.g. ncc(l) for Cartesian Gaussians).
            arg_list.append( windex.output_array_length() )
            # Every index will have a corresponding skip in output array.
            arg_list.append( windex.output_array_skip() )
            # Add pointer to start of spherical transformation coefficients in
            # sph_trans_array (if spherically transformed integrals are requested)
            # [ Not supported currently ]
            assert self.gen.generator_options().spherical_transformed() == False,\
                    "Spherical transformation for primitive integrals not supported."
            #if self.gen.generator_options().spherical_transformed() == True:
            #    arg_list.append( windex.array_index_dict()['sph_trans'] )
       
        # Additional arguments
        for arg in wintegral.integral().arg_list():
            arg_list.append( arg )

        # Add array of spherical transform coefficients (if necessary)
        # [ Not supported currently ]
        assert self.gen.generator_options().spherical_transformed() == False,\
                "Spherical transformation for primitive integrals not supported."
        #if self.gen.generator_options().spherical_transformed() == True:
        #    arg_list.append( wintegral.sph_trans_array().pointer() )

        # Add arguments for work and output array
        arg_list.append( wintegral.work_array().pointer() ) 
        arg_list.append( wintegral.output_array().pointer() )

        # Special cases: 
        # 1 index integral, with translational invariance.
        if len( wintegral.windex_index_list() ) == 1 and wintegral.integral().translational_invariance() == True:
            # If windex has a centre attribute, this can be removed from argument list, since the integral
            # cannot depend on the single centre.
            try:
                windex = wintegral.windex_index_list()[0]
                if windex.index().cen() in arg_list: arg_list.remove( windex.index().cen() )
            except AttributeError:
                # We expect AttributeError where windex.index() does not hav a cen method, and
                # simply skip this step in that case.
                pass
            # Inform user of this using a message.
            message_str = "Integral "+wintegral.integral().name()+" is a 1-index integral with "+\
                          "translational invariance. "+\
                          "The centre of the index has been removed from the integral function "+\
                          "argument list, since this should not affect the value of the integral."
            info_str = str(__name__)+", integral_argument_setup.__call__"
            self.wintegral.generator().add_message( message_str, info_str )


class automatic_local_variable_setup(algo_variable_setup):
    def __call__(self,ordered_list_list,ordered_dict_list):
        """
        Setup ordered lists of integral local variables (this is called in 
        integral_wrapper.setup_local_variable_list) to be automatically assigned.
        
        ordered_list_list:  a Python list of lists, each ordered_list_list[i] is 
                            to be populated with dsl_* objects
                            that will be automatically assigned in the main function
                            for an integral class, ordered such that later elements
                            depend on earlier elements.

        ordered_dict_list:  a Python list of dictionaries, each ordered_dict_list[i]
                            is a dictionary that relates items (keys) in ordered_list_list[i]
                            to some set of values

        Each ordered_list_list[i] has a different use. In this case:
          * ordered_list_list[0] is the list of variables that can be assigned
            before loops over primitive indexes (independent of primitive exponent).
          * ordered_list_list[1] is unused (empty), but kept so a consistent interface
            is presented across primitive and contracted cases.

        Each ordered_dict_list[i] relates to ordered_list_list[i]. In this case:
          * ordered_dict_list[0] is unused (empty dictionary)
          * ordered_dict_list[1] is unused (empty dictionary)
        ... both are kept to present a consistent interface accross primitive and
        contracted cases.
        """
        wintegral = self.wintegral
        assert len(ordered_list_list) == 0, "ordered_list_list must be empty when passed as argument"
        assert len(ordered_dict_list) == 0, "ordered_dict_list must be empty when passed as argument"

        # Append empty lists to ordered_list_list
        ordered_list0   = []
        ordered_list_list.append( ordered_list0 )
        ordered_list_list.append( [] ) # empty, unused
        # Append empty dictionary to ordered_dict_list (for ordered_list_list[0])
        ordered_dict_list.append( {} ) # empty, unused
        ordered_dict_list.append( {} ) # empty, unused

        # Get unordered list of arguments from base expression
        base_expr = wintegral.integral().base()
        unordered_list = expr_list_objects( base_expr, dsl_variable )
        automatic_set  = set( unordered_list )
        # Scan arguments for special functions in base_expr
        # [ Current implementation restriction ]
        # * Only the default Boys function is currently scanned for
        # ... Default Boys function
        default_boys_f = wintegral.generator().generator_options().default_boys_f()
        default_boys_f_list = expr_list_objects( base_expr, default_boys_f )
        assert len( default_boys_f_list ) <= 1, \
                "Only one default_boys_f object per base function expression currently supported"
        for bf in default_boys_f_list:
            bf_dvars = expr_list_objects( bf.x(), dsl_variable )
            for dvar in bf_dvars:
                automatic_set.add( dvar )

        # Scan {h,v}rr_wrapper objects for local variables and add to a set
        for wrr in wintegral.wrr_list():
            assert wrr.rr().rrtype() in ['vrr','hrr'], "Only VRR and HRR type RRs supported currently"
            if wrr.rr().rrtype() == 'vrr':
                rr_expr = wrr.simplified_expr()
            elif wrr.rr().rrtype() == 'hrr':
                rr_expr = wrr.rr().expr()
            rr_unordered_list = expr_list_objects( rr_expr, dsl_variable )
            for v in rr_unordered_list:
                automatic_set.add( v )


        # Recursively search dsl_variable object expr() values for additional
        # dsl_variable objects and add to automatic_set
        tmp_automatic_frozenset = frozenset( automatic_set )
        for var in tmp_automatic_frozenset: 
            l = expr_recursive_search( var ) 
            for v in l:
                automatic_set.add(v)
        
        # Remove any dsl_variables which are main function arguments
        for darg in wintegral.integral_argument_list():
            automatic_set.discard( darg )

        i = 0
        unordered_list = list( automatic_set )
        # Create a list of local variables for the main function for this integral
        # in assignment order, such that all dependencies are satisfied and any
        # dsl_variable object which are components of other objects are replaced
        # by the objects they are a component of.
        tmp_ordered_list = determine_assignment_order( unordered_list )
        for var in tmp_ordered_list:
            # Append to function argument rather than use return, so interface 
            # is consistent with other algo_variable_setup derived classes
            ordered_list_list[0].append( var )

class primitive_work_array_size_source(integral_wrapper_work_array_size_function.work_array_size_source):
    """
    Derived class providing some primitive-integral-specific methods.
    Should be inherited by algorithm-specific work_array_size_source classes in each
    algo_module (e.g. primitive_vrr_only_no_auxiliary).
    """
    def vrr_memreq_expr(self,wvrr_group):
        """
        Returns a dsl_* expression containing the memory required for evaluation
        of a series of VRRs (incrementing the indexes in index_list).
        """
        end_work_array_index_list = wvrr_group.end_index_list
        end_work_array_l_expr_dict = wvrr_group.end_max_l_expr_dict
        for index in end_work_array_index_list:
            # [ Current implementation restriction ]
            assert wvrr_group.end_ncc_type_dict[ index.name() ] == 'ncc0',\
            "Only VRRs which produce ncc0 type ranges of Cartesian components "+\
            "are currently supported."
        index = end_work_array_index_list[0]
        l     = end_work_array_l_expr_dict[ index.name() ]
        ncc0_expr = ( l + 1 ) * ( l + 2 ) * ( l + 3 ) / 6
        memreq_expr = ncc0_expr
        for index in end_work_array_index_list[1:]:
            l     = end_work_array_l_expr_dict[ index.name() ]
            ncc0_expr = ( l + 1 ) * ( l + 2 ) * ( l + 3 ) / 6
            memreq_expr = dsl_binop( op_mul, memreq_expr, ncc0_expr )
        # If an auxiliary index is present, use 2-layer algorithm
        if wvrr_group.auxiliary_index != None:
            assert isinstance(wvrr_group.auxiliary_index,dsl_integer_index),\
                    "Only integer auxiliary indexes supported"
            assert self._wintegral.algorithm_obj().m_layers == 2,\
                    "Only 2-layer auxiliary index algorithm supported"
            m_max = self._wintegral.algorithm_obj().m_max
            # 2-layer auxiliary index algorithm
            memreq_expr = 2 * ( memreq_expr + m_max + 1 )
        return memreq_expr 

    def hrr_memreq_expr(self,whrr_group,index_ncc_dict):
        """
        Returns a dsl_* expression containing the memory required for evaluation
        of a HRR operation.
        """
        assert len(whrr_group.wrr_list) == 1, "Only one HRR allowed per wrr_group"
        whrr = whrr_group.wrr_list[0]
        end_work_array_index_list = whrr_group.end_index_list
        end_work_array_max_l_expr_dict = whrr_group.end_max_l_expr_dict
        end_work_array_min_l_expr_dict = whrr_group.end_min_l_expr_dict
        move_from_index = whrr.move_from_index()
        move_from_end_l = end_work_array_max_l_expr_dict[ move_from_index.name() ]
        move_to_index = whrr.move_to_index()
        move_to_end_l = end_work_array_max_l_expr_dict[ move_to_index.name() ]
        memreq_expr = self.generator().global_hrr_layer_array()[ move_to_end_l ][ move_from_end_l ] 
        other_index_list = end_work_array_index_list[:]
        other_index_list.remove( move_from_index )
        other_index_list.remove( move_to_index )
        for index in other_index_list:
            l     = end_work_array_max_l_expr_dict[ index.name() ]
            l_min = end_work_array_min_l_expr_dict[ index.name() ]
            try:
                assert index_ncc_dict[ index.name() ] in \
                        [ 'ncc', 'ncc0', 'nsph', 'ncc0-ncc0'], \
                        "ncc_type must be one of 'ncc', 'ncc0' and 'ncc0-ncc0'"
                if index_ncc_dict[ index.name() ] == 'ncc0':
                    # Incremented by VRR, use ncc0 (VRRs keep integrals with 
                    # 0 .. l units of angular momentum)
                    ncc_expr = ( l + 1 ) * ( l + 2 ) * ( l + 3 ) / 6
                elif index_ncc_dict[ index.name() ] == 'ncc0-ncc0':
                    # Truncated during previous operation 
                    ncc_expr = ( l + 1 ) * ( l + 2 ) * ( l + 3 ) / 6 - \
                               l_min * (l_min + 1) * (l_min + 2) / 6
                elif index_ncc_dict[ index.name() ] == 'ncc':
                    # Incremented by HRR, or truncated during previous operation
                    ncc_expr = ( l + 1 ) * ( l + 2 ) / 2
                elif index_ncc_dict[ index.name() ] == 'nsph':
                    # [ Not currently supported ]
                    # Spherical transformed and truncated during previous operation
                    assert index_ncc_dict[ index.name() ] != 'nsph',\
                            "Spherical transformation for primitive integrals not supported."
                    ncc_expr = 2 * l + 1
            except KeyError:
                print("An unrecognized index was requested from index_ncc_dict.")
                raise 
            memreq_expr = dsl_binop( op_mul, memreq_expr, ncc_expr )
        memreq_expr = dsl_binop( op_mul, 2, memreq_expr )
        return memreq_expr

    def vrr_to_hrr_memreq_expr(self,whrr_group):
        """
        Returns a dsl_* expression containing the 2x the memory required to store
        the end result of a sequence of VRR operations.

        Since the work_array size function uses 2x the maximum HRR layer size
        as the maximum required memory for a sequence of HRR operations, this
        dsl_* expression ensures that, if the first HRR layer (in which all
        previously incremented indexes go from 0..l and l_move_to = 0) is the
        largest layer, then the maximum work_array size is correct.

        HRRs always follow contraction, so we can assume that all primitive
        lengths have been converted to contracted lengths (i.e. nprim --> ncont)
        """
        assert len(whrr_group.wrr_list) == 1, "Only one HRR allowed per wrr_group"
        whrr = whrr_group.wrr_list[0]
        start_work_array_index_list = whrr_group.start_index_list
        start_work_array_max_l_expr_dict = whrr_group.start_max_l_expr_dict
        index = start_work_array_index_list[0]
        l     = start_work_array_max_l_expr_dict[ index.name() ]
        ncc0_expr = ( l + 1 ) * ( l + 2 ) * ( l + 3 ) / 6
        memreq_expr = ncc0_expr
        for index in start_work_array_index_list[1:]:
            l     = start_work_array_max_l_expr_dict[ index.name() ]
            ncc0_expr = ( l + 1 ) * ( l + 2 ) * ( l + 3 ) / 6
            memreq_expr = dsl_binop( op_mul, memreq_expr, ncc0_expr )
        memreq_expr = 2 * memreq_expr
        return memreq_expr

class primitive_algo_main_source(integral_wrapper_main_function.main_source):
    """
    Derived class providing some primitive-integral-specific methods.
    Should be inherited by algorithm-specific algo_main_source classes in each
    algo_module (e.g. primitive_vrr_only_no_auxiliary).
    """
    def vrr_sequence_no_auxiliary_out(self,p,wvrr_group_list,\
            work_array_index0, work_array_index0_initial_expr ):
        """
        Output code for a sequence of VRRs with no auxiliary indexes (all
        indexes involved are Cartesian Gaussians).

        p:                      cprinter object for code output
        wvrr_group_list:        list of whrr_group objects
                                (see integral_wrapper.algorithm class).
        work_array_index0
        work_array_index_start0
                                dsl_scalar variables which are integer pointers
                                into the work_array (i.e. iwrk0 and
                                iwrk0_start).
        work_array_index0_initial_expr
                                dsl_binop/unop expression which contains the 
                                expression (integer) which represents the start
                                point in the work array for the block of
                                integrals corresponding to work_array_index0
        """

        ### Base class evaluation ###
        # Output call to integral class base() function
        wfunction = self._wintegral.base_function()
        wfunction.call_out(p)

        ### VRR sequence ###
        for wrr_group in wvrr_group_list:
            p.out( work_array_index0.assign( work_array_index0_initial_expr ) )

            # VRR evaluation
            iw = 0 # keep track of number of wrrs that have been output
            for wrr in wrr_group.wrr_list:
                iw += 1
                assert wrr.unrolled_index() is wrr.changing_index(),\
                        "RR must be vertical recurrence relation"
                wfunction = wrr.rr_function()
                changing_windex = self._wintegral.windex_index_dict()[ \
                                            wrr.changing_index().name() ]

                # Set loop length variables for outer loops
                for index in reversed( wrr.loop_index_list() ):
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    loop_length_var  = windex.loop_length()
                    loop_length_expr = wrr_group.loop_length_expr_dict[ windex.index().name() ]
                    p.out( loop_length_var.assign( loop_length_expr ) )

                # Set index_max() value for unrolled loop (inside function)
                lmax_loop_var = changing_windex.index_loop_max()
                lmax_loop_expr = wrr_group.index_loop_max_expr_dict[ changing_windex.index().name() ]
                p.out( lmax_loop_var.assign( lmax_loop_expr ) ) 

                nloops = 0
                # Output loop structure top part
                # ... and set work_array_index() value in last loop before unrolled loop
                # if necessary.
                for index in reversed( wrr.loop_index_list() ):
                    nloops += 1
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    counter     = windex.array_index()
                    loop_start  = 0 
                    loop_length = windex.loop_length()
                    p.forloop( counter.assign(loop_start), \
                            dsl_binop(op_lt, counter , loop_length ), \
                            counter.assign( dsl_binop( op_add, counter, 1 ) ) )

                # Call to rr_function
                # For primitive integrals, use default argument list for VRR function
                wfunction.call_out(p)

                # Update work_array_index()
                work_array_index_dict = \
                    self._wintegral.work_array().array_index_dict()
                if len( wrr.loop_index_list() ) > 0:
                    last_loop_index = list( reversed( wrr.loop_index_list() ) )[-1]
                    windex = self._wintegral.windex_index_dict()[ last_loop_index.name() ]
                    last_loop_skip  = windex.work_array_skip()
                    for change, index_tuple in work_array_index_dict.items():
                        work_array_index = index_tuple[0]
                        p.out( work_array_index.assign( work_array_index + last_loop_skip ) )
                # Output loop structure bottom part
                for loop in range( nloops ) :
                    p.endforloop()

                # Reset work_array array_index values, if they have been modified
                # inside loop structure
                if nloops > 0:
                    p.out( work_array_index0.assign( work_array_index0_initial_expr ) )

    def vrr_sequence_with_auxiliary_out(self, p, wvrr_group_list,\
            work_array_index0, work_array_index1, work_array_index_base, \
            work_array_index_start0, work_array_index_start1, \
            work_array_index0_initial_expr ):
        """
        Output code for a sequence of VRRs with a single auxiliary index, using
        a 2-layer algorithm. This is effective where the integral being assigned
        depends on integrals with the auxiliary index unchanged and incremented.
        For VRRs where the no integral dependencies have an unchanged auxiliary 
        index, a 1-layer algorithm may be more appropriate.

        p:                      cprinter object for code output
        wvrr_group_list:        list of whrr_group objects
                                (see integral_wrapper.algorithm class).
        work_array_index{0,1,base}
        work_array_index_start{0,1}
                                dsl_scalar variables which are integer pointers
                                into the work_array (i.e. iwrk{0,1,b} and
                                iwrk{0,1}_start).
        work_array_index0_initial_expr
                                dsl_binop/unop expression which contains the 
                                expression (integer) which represents the start
                                point in the work array for the block of
                                integrals corresponding to work_array_index0
        """

        ### Base class evaluation ###
        # Output call to integral class base() function
        wfunction = self._wintegral.base_function()
        wfunction.call_out(p)

        ### VRR sequence ###
        ngroups = 0
        for wrr_group in wvrr_group_list:
            # Outside loops over auxiliary indexes and horizontal recurrence relations 
            # [ Current implementation restriction ]
            # This relies on the assumptions:
            # * Only one auxiliary index supported per integral class.
            # * Auxiliary index may only have increments of 0 or 1.
            # * Only one unique VRR supported per integral class.
            # ... i.e. if the VRR expression depends on an aux index, then all
            #          VRRs for integral indexes depend on this.
            nloops_aux = 0
            assert len( self._wintegral.windex_aux_list() ) == 1,\
                    "Only a single auxiliary index is currently supported per integral."
            # Since we are only supporting one auxiliary index per integral class
            # currently, we are not concerned with order. In the future 
            # windex_aux_list() may need to be reordered, or an alternative list
            # used to represent the auxiliary index_wrapper objects.
            waux = self._wintegral.windex_aux_dict()[ wrr_group.auxiliary_index.name() ]
            nloops_aux += 1
            counter     = waux.array_index()
            loop_start  = wrr_group.start_m
            loop_end = wrr_group.end_m
            p.forloop( counter.assign(loop_start), \
                    dsl_binop(op_ge, counter , loop_end ), \
                    counter.assign( dsl_binop( op_sub, counter, 1 ) ) )
            # Set layer start points for next loop over auxiliary index
            p.out( work_array_index_start0.assign( work_array_index0 ) )
            p.out( work_array_index_start1.assign( work_array_index1 ) )
            # Copy base function work_array element to correct position
            waux_array_index = counter
            p.out( self._wintegral.work_array()[ work_array_index0 ].assign(\
                 self._wintegral.work_array()[ work_array_index_base + \
                 waux_array_index ] ) )

            # VRR evaluation
            iw = 0 # keep track of number of wrrs that have been output
            for wrr in wrr_group.wrr_list:
                iw += 1
                assert wrr.unrolled_index() is wrr.changing_index(),\
                        "RR must be vertical recurrence relation"
                wfunction = wrr.rr_function()
                changing_windex = self._wintegral.windex_index_dict()[ \
                                            wrr.changing_index().name() ]

                # Set loop length variables for outer loops
                for index in reversed( wrr.loop_index_list() ):
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    loop_length_var  = windex.loop_length()
                    loop_length_expr = wrr_group.loop_length_expr_dict[ windex.index().name() ]
                    p.out( loop_length_var.assign( loop_length_expr ) )

                # Set index_max() value for unrolled loop (inside function)
                lmax_loop_var = changing_windex.index_loop_max()
                lmax_loop_expr = wrr_group.index_loop_max_expr_dict[ changing_windex.index().name() ]
                p.out( lmax_loop_var.assign( lmax_loop_expr ) ) 

                nloops = 0
                # Output loop structure top part
                # ... and set work_array_index() value in last loop before unrolled loop
                # if necessary.
                for index in reversed( wrr.loop_index_list() ):
                    nloops += 1
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    counter     = windex.array_index()
                    loop_start  = 0 
                    loop_length = windex.loop_length()
                    p.forloop( counter.assign(loop_start), \
                            dsl_binop(op_lt, counter , loop_length ), \
                            counter.assign( dsl_binop( op_add, counter, 1 ) ) )

                # Call to rr_function
                wfunction.call_out(p)

                # Update work_array_index()
                work_array_index_dict = \
                    self._wintegral.work_array().array_index_dict()
                if len( wrr.loop_index_list() ) > 0:
                    last_loop_index = list( reversed( wrr.loop_index_list() ) )[-1]
                    windex = self._wintegral.windex_index_dict()[ last_loop_index.name() ]
                    last_loop_skip  = windex.work_array_skip()
                    for change, index_tuple in work_array_index_dict.items():
                        work_array_index = index_tuple[0]
                        p.out( work_array_index.assign( work_array_index + last_loop_skip ) )
                # Output loop structure bottom part
                for loop in range( nloops ) :
                    p.endforloop()

                # Reset work_array array_index values, if they have been modified
                # inside loop structure
                if nloops > 0:
                    if iw != len( wrr_group.wrr_list ):
                        p.out( work_array_index0.assign( work_array_index_start0 ) )
                        p.out( work_array_index1.assign( work_array_index_start1 ) )
                    else:
                        p.out( work_array_index0.assign( work_array_index0_initial_expr ) )

            # The following code assumes that only increments of 0 and 1 are
            # present and that only a single auxiliary index is present for
            # each integral class
            # Loops over auxiliary index
            assert len( self._wintegral.windex_aux_list() ) == 1,\
                    "Only a single auxiliary index is currently supported per integral."
            assert max( list( work_array_index_dict.keys() ) ) <= 1,\
                    "Only increments of 0 and 1 supported in auxiliary index"
            # src_work_array_index{0,1} and src_work_array_index_start{0,1} are
            # initialized above (in another if wrr_group.auxiliary_index != None
            # block).
            # Swap layers for next loop over auxiliary index
            p.out( work_array_index0.assign( work_array_index_start1 ) )
            p.out( work_array_index1.assign( work_array_index_start0 ) )

            # Output loops structure bottom part (auxiliary indexes)
            for loop in range( nloops_aux ):
                p.endforloop()

            # The following code assumes that only increments of 0 and 1 are
            # present and that only a single auxiliary index is present for
            # each integral class
            # Unswap layers, unless this is the final wvrr_group -- the layers are
            # unswapped at the beginning of the output array copy, so this final
            # unswap is unnecessary
            ngroups += 1
            if ngroups != len( wvrr_group_list ):
                # Unswap layers for next loop structure over auxiliary index
                p.out( work_array_index0.assign( work_array_index_start0 ) )
                p.out( work_array_index1.assign( work_array_index_start1 ) )

    def hrr_sequence_out(self,p,whrr_group_list,\
                    work_array_index0, work_array_index_start0, \
                    work_array_index1, work_array_index_start1,\
                    initial_hrr_swap ):
        """
        Output code for a sequence of HRRs.
        [ Currently limited to 1 per integral class ]

        p:                      cprinter object for code output
        whrr_group_list:        list of whrr_group objects
                                (see integral_wrapper.algorithm class).
        work_array_index{0,1}
        work_array_index_start{0,1}
                                dsl_scalar variables which are integer pointers
                                into the work_array (i.e. iwrk{0,1} and
                                iwrk{0,1}_start).
        initial_hrr_swap:       bool, determines whether code which swaps 
                                work_array_index{0,1} values between 
                                work_array_index_start{0,1} values when
                                first HRR step is executed (i.e. l_move_to -> 1 )
                                is output.
        """
        ### HRRs ###
        for wrr_group in whrr_group_list:
            assert len( wrr_group.wrr_list ) <= 1, "Only one HRR per wrr_group"
            wrr = wrr_group.wrr_list[0]
            wfunction = wrr.rr_function()
            # work_array array_index and array_index_start variables should have been
            # set already, i.e.
            move_to_windex = self._wintegral.windex_index_dict()[ \
                                        wrr.move_to_index().name() ]
            move_from_windex = self._wintegral.windex_index_dict()[ \
                                        wrr.move_from_index().name() ]
            assert move_to_windex.index().is_cartesian() == True,\
                    "move_to_index in HRR must be Cartesian"
            assert move_from_windex.index().is_cartesian() == True,\
                    "move_from_index in HRR must be Cartesian"
            move_to_loop_lmax_var = move_to_windex.index_loop_max()
            move_to_start_min_l = None # move_to_index is not in start_l_expr_dict
            move_to_start_max_l = None # move_to_index is not in start_l_expr_dict
            move_to_end_min_l   = wrr_group.end_min_l_expr_dict[ \
                                            move_to_windex.index().name() ]
            move_to_end_max_l   = wrr_group.end_max_l_expr_dict[ \
                                            move_to_windex.index().name() ]
            move_from_loop_lmax_var = move_from_windex.index_loop_max()
            move_from_start_min_l = wrr_group.start_min_l_expr_dict[ \
                                        move_from_windex.index().name() ]
            move_from_start_max_l = wrr_group.start_max_l_expr_dict[ \
                                        move_from_windex.index().name() ]
            move_from_end_min_l   = wrr_group.end_min_l_expr_dict[ \
                                        move_from_windex.index().name() ]
            move_from_end_max_l   = wrr_group.end_max_l_expr_dict[ \
                                        move_from_windex.index().name() ]

            # Only execute HRR if l_move_to > 0
            # Top part of if statement
            p.ifblock( dsl_binop( op_gt, move_to_end_max_l, 0 ) )

            # Set loop length variables for outer loops
            for index in wrr.loop_index_list():
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                #ncc_type = final_wvrr_group.end_ncc_type_dict[ index.name() ]
                #l_expr   = final_wvrr_group.end_max_l_expr_dict[ index.name() ]
                # Can use end_{max_l_expr,ncc_type}_dict values from HRR group
                # Since these outer loop indexes are truncated during the first
                # step of the HRR
                # process (ncc0 --> ncc) and their max_l_expr is unchanged.
                ncc_type = wrr_group.end_ncc_type_dict[ index.name() ]
                l_max_expr   = wrr_group.end_max_l_expr_dict[ index.name() ]
                l_min_expr   = wrr_group.end_min_l_expr_dict[ index.name() ]
                loop_length_var  =  windex.loop_length() 
                # Determine loop_length_expr based on ncc_type
                if ncc_type == 'ncc':
                    loop_length_expr = \
                            (l_max_expr + 1)*(l_max_expr + 2)/2 
                elif ncc_type == 'ncc0':
                    loop_length_expr = \
                            (l_max_expr + 1)*(l_max_expr + 2)*(l_max_expr + 3)/6
                elif ncc_type == 'ncc0-ncc0':
                    loop_length_expr = \
                            (l_max_expr + 1)*(l_max_expr + 2)*(l_max_expr + 3)/6 -\
                            (l_min_expr)*(l_min_expr + 1)*(l_min_expr + 2)/6
                elif ncc_type == 'nsph':
                    # [ Not currently supported ]
                    # Incremented by HRR, or truncated during previous operation
                    assert index_ncc_dict[ index.name() ] != 'nsph',\
                            "Spherical transformation for primitive integrals not supported."
                    loop_length_expr = 2 * l_max_expr + 1
                else:
                    raise Exception('ncc_type unrecognized')
                p.out( loop_length_var.assign( loop_length_expr ) )

            # Output loop_length and work_array skip variables
            # Variables with lower angmom in move_to_index take values from previous iteration 
            # higher angmom layer, or from last operation in VRR sequence
            # Output dummy work_array_skip variable for move_to_index (this is not 
            # needed in the first loop over move_to_index, as l_move_to = 0, but
            # needs to be set to avoid using unpredictable behaviour from uninitialized
            # value)
            work_array_skip_var = move_to_windex.work_array_skip() 
            p.out( work_array_skip_var.assign(0) )
            for index in wrr_group.end_index_list:
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                # move_to_index layer
                work_array_skip_var   = windex.work_array_skip()
                # move_to_index - 1 layer
                hrr_skip_var   = windex.skip_dict()['tmp']
                # hrr_length_var never needs to be set
                p.out( hrr_skip_var.assign( work_array_skip_var ) )
            # Output work_array_skip variable for move_from_index:
            assert wrr_group.end_index_list[-1] is move_from_windex.index(),\
                    "last index in end_index_list mist be move_from_index"
            work_array_skip_var = move_from_windex.work_array_skip() 
            p.out( work_array_skip_var.assign(1) )

            if initial_hrr_swap == True:
                # Initial work_array array_index assignments: 
                # If no auxiliary index in VRR sequence, iwrk0 and iwrk1 need swapping, 
                # since the result of the VRR sequence. If the VRR sequence involves
                # an auxiliary index and a 2-layer algorithm has been used,
                # then this swap will already have been done.
                p.out( work_array_index0.assign( work_array_index_start1 ) )
                p.out( work_array_index1.assign( work_array_index_start0 ) )

            ## Special first angmom transfer ##
            # l_move_to = 1
            # First angular momentum transfer (from l_move_to = 0 to l_move_to = 1)
            # is special, since the indexing of the endpoint of the VRR sequence 
            # will not necessarily have move_to_index and move_from_index as the
            # two fastest changing indexes.
            # Additionally, the loops size of other integral indexes (previously
            # incremented) may be truncated (ncc0 length to ncc)

            # Set l_move_from
            counter    = move_to_loop_lmax_var
            p.out( counter.assign( 1 ) )

            # Set layer start points for this iteration over move_to_index
            p.out( work_array_index_start0.assign( work_array_index0 ) )
            p.out( work_array_index_start1.assign( work_array_index1 ) )

            # Output loop_length
            # move_to_index (ncc_type = 'ncc' always)
            # this is always the second fastest changing index in the HRR sequence
            # -- do not need to set work_array_length is only two indexes
            if len( wrr_group.end_index_list ) > 2:
                # l_move_to = 1, ncc(l_move_to) = 3
                loop_length_var  = move_to_windex.loop_length() 
                p.out( loop_length_var.assign( 3 ) )
            # move_from_index (special ncc expression, ncc0(l_mt+l_mf-l_mf_loop) - ncc0(l_mf) )
            l_expr1 = move_from_start_max_l
            l_expr2 = move_from_end_max_l
            loop_length_var  = move_from_windex.loop_length()
            loop_length_expr = (l_expr1)*(l_expr1+1)*(l_expr1+2)/6 -\
                                l_expr2*(l_expr2+1)*(l_expr2+2)/6
            p.out( loop_length_var.assign( loop_length_expr ) )
            previous_windex = wrr_group.end_index_list[-1]
            skip_expr = move_from_windex.loop_length()
            for index in reversed( wrr_group.end_index_list[:-1] ):
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                work_array_skip_var = windex.work_array_skip()
                # Output iwskip value
                p.out( work_array_skip_var.assign( skip_expr ) )
                # Update skip_expr
                skip_expr = skip_expr * windex.loop_length()

            # Output loops over indexes in wrr.loop_index_list (not move_{to,from}_index)
            nloops = 0
            for index in wrr.loop_index_list():
                nloops += 1
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                counter     = windex.array_index()
                loop_start  = 0
                loop_length = windex.loop_length()
                # contains information about work_array size
                p.forloop( counter.assign(loop_start),\
                        dsl_binop(op_lt, counter , loop_length ),\
                        counter.assign( dsl_binop( op_add, counter, 1 ) ) )
            # For the l_move_to = 0 "layer", the indexing is different to 
            # subsequent layers: move_from_index and move_to_index may not be the
            # two fastest changing indexes. Therefore, set iwrk1 
            # by multiplying skip values, rather than incrementing each inner loop
            if len( wrr.loop_index_list() ) > 0:
                index = wrr.loop_index_list()[0]
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                index_expr0 = windex.array_index() * windex.work_array_skip()
                index_expr1 = ( windex.array_index() + windex.work_array_offset() )*\
                        windex.skip_dict()['tmp']
                for index in wrr.loop_index_list()[1:]:
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    index_expr0 = index_expr0 + windex.array_index() * windex.work_array_skip()
                    index_expr1 = index_expr1 + \
                            ( windex.array_index() + windex.work_array_offset() ) *\
                            windex.skip_dict()['tmp']
                p.out( work_array_index0.assign(\
                        work_array_index_start0 + index_expr0 ) )
                p.out( work_array_index1.assign(\
                        work_array_index_start1 + index_expr1 ) )

            # Move work_array_index1 to start of l_move_from integrals 
            # The integrals from 0..l_move_from can be ignored in this VRR sequence
            # to HRR sequence transition
            l_expr      = move_from_end_max_l
            offset_expr = move_from_windex.skip_dict()['tmp'] *\
                    l_expr * (l_expr + 1) * (l_expr + 2) / 6 
            p.out( work_array_index1.assign( work_array_index1 + offset_expr ) )

            # Call to rr_function
            wfunction.call_out(p)

            # Output loop structure bottom part for loop_index_list indexes
            for loop in range( nloops ) :
                p.endforloop()

            # Swap iskipX_hrr values
            for index in wrr_group.end_index_list:
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                # move_to_index layer
                work_array_skip_var   = windex.work_array_skip()
                # move_to_index - 1 layer
                hrr_skip_var   = windex.skip_dict()['tmp'] 
                # hrr_length_var never needs to be set
                p.out( hrr_skip_var.assign( work_array_skip_var ) )

            # Swap layers for next loop over move_to_index
            p.out( work_array_index0.assign( work_array_index_start1 ) )
            p.out( work_array_index1.assign( work_array_index_start0 ) )
            # For subsequent angmom transfers, start with lmax_loop = 2 for 
            # move_from_index, since the first angmom transfer was just done
            loop_start = 2

            ## Standard angmom transfers ##
            # angular momentum ranges of previously incremented indexes and
            # move_from index are truncated, either by previous angmom transfer
            # for this wrr_group, or a previous HRR or contraction operation.
            nloops = 0
            # Output loop over l for move_to_index
            nloops += 1
            counter    = move_to_loop_lmax_var 
            loop_length = move_to_end_max_l 
            p.forloop( counter.assign(loop_start),\
                    dsl_binop(op_le, counter , loop_length ),\
                    counter.assign( dsl_binop( op_add, counter, 1 ) ) )

            # Set layer start points for this iteration over move_to_index
            p.out( work_array_index_start0.assign( work_array_index0 ) )
            p.out( work_array_index_start1.assign( work_array_index1 ) )

            # Output loop_length
            # move_to_index (ncc_type = 'ncc' always)
            # this is always the second fastest changing index in the HRR sequence
            # -- do not need to set work_array_length is only two indexes
            if len( wrr_group.end_index_list ) > 2:
                l_expr = move_to_loop_lmax_var 
                loop_length_var  = move_to_windex.loop_length()
                loop_length_expr = (l_expr+1)*(l_expr+2)/2
                p.out( loop_length_var.assign( loop_length_expr ) )
            # move_from_index (special ncc expression, ncc0(l_mt+l_mf-l_mf_loop) - ncc0(l_mf) )
            l_expr1 = move_from_start_max_l - move_to_loop_lmax_var
            l_expr2 = move_from_end_max_l
            loop_length_var  = move_from_windex.loop_length()
            loop_length_expr = (l_expr1+1)*(l_expr1+2)*(l_expr1+3)/6 -\
                    l_expr2*(l_expr2+1)*(l_expr2+2)/6 
            p.out( loop_length_var.assign( loop_length_expr ) )
            previous_windex = wrr_group.end_index_list[-1]
            skip_expr = move_from_windex.loop_length()
            for index in reversed( wrr_group.end_index_list[:-1] ):
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                work_array_skip_var = windex.work_array_skip()
                # Output iwskip value
                p.out( work_array_skip_var.assign( skip_expr ) )
                # Update skip_expr
                skip_expr = skip_expr * windex.loop_length()

            # Output loops over indexes in wrr.loop_index_list (not move_{to,from}_index)
            for index in wrr.loop_index_list():
                nloops += 1
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                counter     = windex.array_index()
                loop_start  = 0
                loop_length = windex.loop_length()
                # contains information about work_array size
                p.forloop( counter.assign(loop_start),\
                        dsl_binop(op_lt, counter , loop_length ),\
                        counter.assign( dsl_binop( op_add, counter, 1 ) ) )

            # Set work_array_index pointers
            if len( wrr.loop_index_list() ) > 0:
                index = wrr.loop_index_list()[0]
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                index_expr0 = windex.array_index() * windex.work_array_skip()
                index_expr1 = windex.array_index() * windex.skip_dict()['tmp']
                for index in wrr.loop_index_list()[1:]:
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    index_expr0 = index_expr0 + windex.array_index() * windex.work_array_skip()
                    index_expr1 = index_expr1 + windex.array_index() *\
                            windex.skip_dict()['tmp']
                p.out( work_array_index0.assign(\
                        work_array_index_start0 + index_expr0 ) )
                p.out( work_array_index1.assign(\
                        work_array_index_start1 + index_expr1 ) )
           
            # Call to rr_function
            wfunction.call_out(p)

            # Output loop structure bottom part for loop_index_list indexes
            for loop in range( nloops-1 ) :
                p.endforloop()

            # Swap iskipX_hrr values
            for index in wrr_group.end_index_list:
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                # move_to_index layer
                work_array_skip_var   = windex.work_array_skip()
                # move_to_index - 1 layer
                hrr_skip_var   = windex.skip_dict()['tmp']
                # hrr_length_var never needs to be set
                p.out( hrr_skip_var.assign( work_array_skip_var ) )

            # Swap layers for next loop over move_to_index
            p.out( work_array_index0.assign( work_array_index_start1 ) )
            p.out( work_array_index1.assign( work_array_index_start0 ) )

            # Close final loop over move_to_index l values
            p.endforloop()

            # Set state of work_array dimensions (lengths, offsets) after HRR 
            for index in wrr_group.end_index_list :
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                ncc_type = wrr_group.end_ncc_type_dict[ windex.index().name() ]
                l_expr   = wrr_group.end_max_l_expr_dict[ windex.index().name() ]
                work_array_length_var   = windex.work_array_length()
                work_array_offset_var   = windex.work_array_offset()
                output_array_length_var = windex.output_array_length() 
                # Determine work_array_length_expr and work_array_offset_expr 
                # from ncc_type
                assert ncc_type in [ 'ncc', 'nsph' ],\
                        "ncc_type must be [ 'ncc', 'nsph' ],  for output phase"
                if ncc_type == 'ncc':
                    #work_array_length_expr = (l_expr + 1)*(l_expr + 2)/2
                    work_array_offset_expr = 0
                    work_array_offset_expr = work_array_offset_expr
                    # Offset is zero, no need to set length
                    p.out( work_array_offset_var.assign( work_array_offset_expr ) )
                elif ncc_type == 'nsph':
                    # [ Not currently supported ]
                    assert index_ncc_dict[ index.name() ] != 'nsph',\
                            "Spherical transformation for primitive integrals not supported."
                    work_array_length_expr = 2 * l_expr + 1
                    work_array_offset_expr = 0
                    # Offset is zero, no need to set length
                    p.out( work_array_offset_var.assign( work_array_offset_expr ) )


            # Only execute HRR if l_move_to > 0
            # If l_move_to == 0, update move_from_index length and offset
            p.elseblock()
            work_array_offset_var = move_to_windex.work_array_offset() 
            work_array_skip_var   = move_to_windex.work_array_skip()
            # Output assignment for work_array_offset, where l_move_to == 0
            # (no need to set work_array_length, as nothing depends on this)
            p.out( work_array_offset_var.assign(0) )
            # Set work_array_skip for move_to_index = 0, since this needs 
            # initializing
            p.out( work_array_skip_var.assign(0) )

            # End of if..else..block
            p.endifblock()

    def copy_to_output_array_out(self,p,final_operation,work_array_index0):
        """
        Output code for copying integrals from the work_array to the output_array
        following completion of all operations.

        p:                      cprinter object for code output
        final_operation:        object of wvrr_group, whrr_group, contraction,
                                spherical_transform etc (see 
                                integral_wrapper.algorithm class) -- last operation
                                output in integral main function.
        work_array_index_start0:
                                dsl_scalar variable which is an integer pointers
                                into the work_array at the start of the integral
                                block for output
        """
        nloops = 0
        output_windex_loop_list  = []
        for index in final_operation.end_index_list:
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            output_windex_loop_list.append( windex )

        # Output loops structure top part
        for windex in output_windex_loop_list:
            nloops += 1
            counter     = windex.array_index() 
            loop_start  = 0
            loop_length = windex.output_array_length()
            p.forloop( counter.assign(loop_start),\
                    dsl_binop(op_lt, counter , loop_length ),\
                    counter.assign( dsl_binop( op_add, counter, 1 ) ) )

        # Generate correct DSL expressions for output and work array indexing
        #work_array_index_expr = src_work_array_index0
        work_array_index_expr = work_array_index0
        output_array_index_expr = None           
        for i in range(0, len(output_windex_loop_list) ):
            windex = output_windex_loop_list[i]
            work_expr = dsl_binop( op_mul, windex.array_index() + \
                            windex.work_array_offset(), windex.work_array_skip() ) 
            output_expr = dsl_binop( op_mul, windex.array_index(), \
                            windex.output_array_skip() )
            work_array_index_expr = dsl_binop( op_add, work_array_index_expr, \
                            work_expr )
            if output_array_index_expr == None:
                output_array_index_expr = output_expr
            else:
                output_array_index_expr = dsl_binop( op_add, \
                        output_array_index_expr, output_expr )
        output_array_ref = self._wintegral.output_array()[ output_array_index_expr ]
        work_array_ref = self._wintegral.work_array()[ work_array_index_expr ]
        # Assignment of work_array elements to output array
        p.out( output_array_ref.assign( work_array_ref ) )

        # Output loop structure bottom part
        for loop in range( nloops ) :
            p.endforloop()
