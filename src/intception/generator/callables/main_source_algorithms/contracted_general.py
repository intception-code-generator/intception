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
Provides general classes and methods for use in specific contracted algorithm
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
        integral_wrapper.setup_integral_argument_list) for contracted
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
            # Contracted integrals are required---add nprim and ncont lengths
            assert isinstance( windex.index(), dsl_cartesian_gaussian ),\
                    "only dsl_cartesian_gaussian type indexes are supported for contracted integrals."
            arg_list.append( windex.length_dict()['prim'] )
            arg_list.append( windex.length_dict()['cont'] )
            arg_list.append( windex.skip_dict()['cont'] )
            # Add pointer to start of contraction coefficients in contract_array
            arg_list.append( windex.contract_array().pointer() )
            # Add pointer to start of spherical transformation coefficients in
            # sph_trans_array (if spherically transformed integrals are requested)
            if self.gen.generator_options().spherical_transformed() == True:
                arg_list.append( windex.sph_trans_array().pointer() )
        
        # Additional arguments
        for arg in wintegral.integral().arg_list():
            arg_list.append( arg )

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
          * ordered_list_list[1] is the list of variables that depend on 
            primitive exponents, and must be assigned inside loops over primitive
            exponents.

        Each ordered_dict_list[i] relates to ordered_list_list[i]. In this case:
          * ordered_dict_list[0] is unused (empty dictionary)
          * ordered_dict_list[1] has as keys items from ordered_dict_list[1] and as
            values sets of per-primitive variables that the key depends on, used for
            determining where in a set of loops over primitives the key dsl_* object
            should be assigned in output code.
        """
        wintegral = self.wintegral
        assert len(ordered_list_list) == 0, "ordered_list_list must be empty when passed as argument"
        assert len(ordered_dict_list) == 0, "ordered_dict_list must be empty when passed as argument"

        # Scan base expr for local variables and add to an unordered list
        base_expr = wintegral.integral().base()
        unordered_list = expr_list_objects( base_expr, dsl_variable )
        automatic_set0  = set( unordered_list )

        # Append empty lists to ordered_list_list
        ordered_list0   = []
        ordered_list1  = []
        ordered_list_list.append( ordered_list0 )
        ordered_list_list.append( ordered_list1 )
        # Append empty dictionary to ordered_dict_list (for ordered_list_list[0])
        ordered_dict0 = {}
        ordered_dict_list.append( ordered_dict0 )
        
        # Create list of automatically assigned sets of variables
        # Variables in later lists may depend on variables in previous list, but
        # not vice versa, e.g. automatic_set_list[1] may depend on automatic_set_list[0]
        # This allows separate blocks of interdependent variables to be assigned
        automatic_set_list = [ automatic_set0 ] 
        # Also create a list of lists of ordered variables, corresponding to the 
        # automatic_set_list, but where variables have been ordered according to
        # dependencies
        # ( This is ordered_list0 and ordered_list1, created above )

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
                automatic_set0.add( dvar )

        # Scan {h,v}rr_wrapper objects for local variables and add to a set
        for wrr in wintegral.wrr_list():
            assert wrr.rr().rrtype() in ['vrr','hrr'], "Only VRR and HRR type RRs supported currently"
            if wrr.rr().rrtype() == 'vrr':
                rr_expr = wrr.simplified_expr()
            elif wrr.rr().rrtype() == 'hrr':
                rr_expr = wrr.rr().expr()
            rr_unordered_list = expr_list_objects( rr_expr, dsl_variable )
            for v in rr_unordered_list:
                automatic_set0.add( v )


        # Recursively search dsl_variable object expr() values for additional
        # dsl_variable objects and add to automatic_set0
        tmp_automatic_frozenset = frozenset( automatic_set0 )
        for var in tmp_automatic_frozenset: 
            l = expr_recursive_search( var ) 
            for v in l:
                automatic_set0.add(v)

        # Split automatically assigned variables into two sets, one containing variables
        # which depend (directly or indirectly) on primitive exponents, and one containing
        # variables which do not depend on exponents. 
        # This should not cause any issues of depencies not being fulfilled provided that
        # the variables depending on the primitive exponents are assigned after those 
        # that do not depend on the exponents.
        # (The "tree" derived variables depending on exponents can only depend on the 
        # variables outside of that "tree" -- the reverse is not true.)
        tmp_automatic_frozenset = frozenset( automatic_set0 )
        automatic_set1_depends_dict = {} # relates members of automatic_set1 (keys) to
                                         # a the variables (exponents) they depend on
                                         # (recursively)
        automatic_set1 = set() # variables depending on exponents
        automatic_set_list.append( automatic_set1 )
        ordered_list_list.append( ordered_list1 )
        for var in tmp_automatic_frozenset:
            l = expr_recursive_search( var )
            for v in wintegral.contraction_info_obj().contraction_convert_dict.keys():
                if v in l:
                    automatic_set1.add( var )
                    if var in automatic_set1_depends_dict.keys():
                        automatic_set1_depends_dict[ var ].add( v )
                    else:
                        automatic_set1_depends_dict[ var ] = set([ v ])

                    automatic_set0.discard( var )

        # Append automatic_set1_depends_dict to ordered_dict_list (as ordered_dict_list[1])
        # for output
        ordered_dict_list.append( automatic_set1_depends_dict ) 
        
        # Automatic set list must have length 2, since we only deal with two sets of
        # interdependent automatically-assigned variables
        assert len(automatic_set_list ) == 2, "len(automatic_set_list) must be 2"

        # Remove any dsl_variables which are main function arguments
        for aset in automatic_set_list:
            for darg in wintegral.integral_argument_list():
                aset.discard( darg )

        i = 0
        for aset in automatic_set_list:
            unordered_list = list( aset )
            # Create a list of local variables for the main function for this integral
            # in assignment order, such that all dependencies are satisfied and any
            # dsl_variable object which are components of other objects are replaced
            # by the objects they are a component of.
            tmp_list = determine_assignment_order( unordered_list )
            # Add ordered list entries to corresponding list for output
            for var in tmp_list:
                ordered_list_list[i].append( var )
            i += 1


class contracted_work_array_size_source(integral_wrapper_work_array_size_function.work_array_size_source):
    """
    Derived class providing some contracted-integral-specific methods.
    Should be inherited by algorithm-specific work_array_size_source classes in each
    algo_module (e.g. contracted_vrr_only_no_auxiliary).
    """
    def vrr_memreq_expr(self,wvrr_group,n_contraction):
        """
        Returns a dsl_* expression containing the memory required for evaluation
        of a series of VRRs (incrementing the indexes in index_list).

        If an auxiliary index is present, this asks for more memory than strictly
        necessary, but the indexing of this memory is simpler than the minimum memory 
        requirement:
            ncc0(a) * ncc0(b) * ... ncc0(X) * ( 1 + n_contraction) + m_max + 1
        since the memory asked for is easily divisible by two:
            2 * ( ncc0(a) * ncc0(b) * ... ncc0(X) * n_contraction + m_max + 1 )
        and is (at most) slightly larger than the memory required for the 
        contractions that follow:
            2 * ( ncc0(a) * ncc0(b) * ... ncc0(X) * n_contraction )

        If the final wvrr_group has m_layers == 1, then the memory requirement for
        the VRR sequence could be less than the value calculated by the expression
        returned by this function, since the implied auxiliary index algorithm
        requires
            ( ncc0(a) * ncc0(b) * ... ncc0(X) * n_contraction ).
        However, since the subsequent early_transform step always requests
            2 * ( ncc0(a) * ncc0(b) * ... ncc0(X) * n_contraction1
        where n_contraction1 is the largest product of nprim and ncont values involved
        in the early_transform sequence, the amount of additional memory requested 
        by the VRR sequence would be at most 2 * ( m_max + 1 ). 
        Since requesting a (possibly) slightly larger block of memory simplifies the
        indexing of integrals, we do not change the behaviour of the work_array_size
        function where the final wvrr_group has m_layers == 1.
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
        # Contracted--- multiply by maximum product of nprim and ncont
        # values
        memreq_expr = dsl_binop( op_mul, memreq_expr, n_contraction )
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

    def hrr_memreq_expr(self,whrr_group,index_ncc_dict,n_contraction):
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
                        [ 'ncc', 'ncc0', 'ncc0-ncc0', 'nsph'], \
                        "ncc_type must be one of 'ncc', 'ncc0' and 'ncc0-ncc0'"
                if index_ncc_dict[ index.name() ] == 'ncc0':
                    # Incremented by VRR, use ncc0 (VRRs keep integrals with 
                    # 0 .. l units of angular momentum)
                    ncomp_expr = ( l + 1 ) * ( l + 2 ) * ( l + 3 ) / 6
                elif index_ncc_dict[ index.name() ] == 'ncc0-ncc0':
                    # Truncated during previous contraction/HRR operation 
                    ncomp_expr = ( l + 1 ) * ( l + 2 ) * ( l + 3 ) / 6 - \
                               l_min * (l_min + 1) * (l_min + 2) / 6
                elif index_ncc_dict[ index.name() ] == 'ncc':
                    # Incremented by HRR, or truncated during previous operation
                    ncomp_expr = ( l + 1 ) * ( l + 2 ) / 2
                elif index_ncc_dict[ index.name() ] == 'nsph':
                    # Incremented by HRR, or truncated during previous operation
                    ncomp_expr = 2 * l + 1
            except KeyError:
                print("An unrecognized index was requested from index_ncc_dict.")
                raise 
            memreq_expr = dsl_binop( op_mul, memreq_expr, ncomp_expr )
        memreq_expr = dsl_binop( op_mul, 2, memreq_expr )
        memreq_expr = dsl_binop( op_mul, memreq_expr, n_contraction )
        return memreq_expr


    def early_transform_memreq_expr(self,wvrr_group,n_contraction):
        """
        Returns a dsl_* expression containing the memory required for evaluation
        for the initial transition from VRR phase to early transform (pre-HRR) phase.

        This requests double the size of the product of all Cartesian component
        indexes at the end of the sequence of VRR operations multiplied by 
        the largest combination of primitive and contraction length indexes 
        (n_contraction), i.e.
            memreq = 2 * ncc0(a) * ncc0(b) * .. ncc0(X) * n_contraction

        This will always be greater than the memory requirement of any 
        single contraction operation, since n_contraction is the largest
        combination of nprim and ncont values, and the number of Cartesian
        components does not change during contraction.

        We do not need to consider the size of any subsequent spherical
        transformation steps, since these will always require less memory than the
        initial contraction.

        Implementation note:
        This will in general request slightly more memory than is really 
        necessary, since we truncate the number of angular momentum components 
        in some indexes in the first contraction step.

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
        memreq_expr = 2 * memreq_expr * n_contraction
        return memreq_expr

    def late_transform_memreq_expr(self,whrr_group,n_contraction):
        """
        Returns a dsl_* expression containing the memory required for evaluation
        for the initial transition from VRR phase to late transform (post-HRR) phase.

        The n_contraction object should have previously been assigned (in code
        output) to the product of contraction length indexes.

        This requests double the size of the product of all Cartesian /spherical 
        component indexes at the end of the sequence of HRR operations multiplied by 
        the contraction length indexes (n_contraction) -- all contractions should
        be completed by this point, i.e.
            memreq = 2 * ncomp(a) * ncomp(b) * .. ncomp(X) * n_contraction
        where ncomp is ncc, or nsph, depending on which indexes have been
        spherically transformed prior to the HRR sequence.

        This will always be greater than the memory requirement of any 
        subsequent spherical transformation, since nsph <= ncc.
        """
        end_work_array_index_list = whrr_group.end_index_list
        end_work_array_l_expr_dict = whrr_group.end_max_l_expr_dict
        end_work_array_ncc_type_dict = whrr_group.end_ncc_type_dict
        for index in end_work_array_index_list:
            # [ Current implementation restriction ]
            assert end_work_array_ncc_type_dict[ index.name() ] in ['ncc','nsph'],\
            "After HRR, angular momentum ranges should be truncated such that "+\
            "ncc_type == 'ncc' or 'nsph'."
        index = end_work_array_index_list[0]
        l     = end_work_array_l_expr_dict[ index.name() ]
        if end_work_array_ncc_type_dict[ index.name() ] == 'ncc':
            ncc_expr = ( l + 1 ) * ( l + 2 ) / 2
        elif end_work_array_ncc_type_dict[ index.name() ] == 'nsph':
            ncc_expr = 2 * l + 1
        memreq_expr = ncc_expr
        for index in end_work_array_index_list[1:]:
            l     = end_work_array_l_expr_dict[ index.name() ]
            if end_work_array_ncc_type_dict[ index.name() ] == 'ncc':
                ncc_expr = ( l + 1 ) * ( l + 2 ) / 2
            elif end_work_array_ncc_type_dict[ index.name() ] == 'nsph':
                ncc_expr = 2 * l + 1
            memreq_expr = dsl_binop( op_mul, memreq_expr, ncc_expr )
        memreq_expr = 2 * memreq_expr * n_contraction
        return memreq_expr

    def max_n_contraction_out(self,p,n_contraction1,n_contraction2):
        """
        Outputs code for determining the maximum required product of 
        primitive and contraction lengths, based on order of 
        contraction of indexes.

        The result is that n_contraction1 is set to the largest out of 
            na_prim * nb_prim * ... * nX_prim,
            na_cont * nb_prim * ... * nX_prim,
            na_cont * nb_cont * ... * nX_prim,
            ...
            na_cont * nb_cont * ... * nX_cont,
        where the order of contraction is a,b,c,...X
        """
        windex_contract_order_list = []
        for t in self._wintegral.algorithm_obj().early_transform_list:
            if t.optype == 'contraction':
                windex = self._wintegral.windex_index_dict()[\
                         t.index_to_contract.name() ]
                windex_contract_order_list.append( windex )
        expr_to_compare_list = []
        windex = windex_contract_order_list[0]
        # all nprim dimensions
        expr = windex.length_dict()['prim']
        for windex in windex_contract_order_list[1:]:
            expr = expr * windex.length_dict()['prim']
        expr_to_compare_list.append( expr )
        # contractions
        for windex in windex_contract_order_list:
            n_prim = windex.length_dict()['prim']
            n_cont = windex.length_dict()['cont']
            expr = expr_find_and_replace( expr, n_prim, n_cont )
            expr_to_compare_list.append( expr )
        # Loop over comparisons to determine largest expr size
        p.out( n_contraction1.assign( expr_to_compare_list[0] ) )
        for expr in expr_to_compare_list[1:]:
            p.out( n_contraction2.assign( expr ) )
            p.ifblock( dsl_binop( op_lt, n_contraction1, n_contraction2 ) )
            p.out( n_contraction1.assign( n_contraction2 ) )
            p.endifblock()

    def n_contraction_out(self,p,n_contraction,all_prim):
        """
        Outputs code assigning n_contraction to the product of 
        primitive/contraction lengths.

        all_prim:   if True, output product of primitive lengths
                    if False, output product of contraction lengths
        """
        if all_prim == True:    keystr = 'prim'
        elif all_prim == False: keystr = 'cont'
        windex_contract_order_list = []
        for t in self._wintegral.algorithm_obj().early_transform_list:
            if t.optype == 'contraction':
                windex = self._wintegral.windex_index_dict()[\
                         t.index_to_contract.name() ]
                windex_contract_order_list.append( windex )
        windex = windex_contract_order_list[0]
        expr = windex.length_dict()[keystr]
        for windex in windex_contract_order_list[1:]:
            expr = expr * windex.length_dict()[keystr]
        p.out( n_contraction.assign( expr ) )

class contracted_algo_main_source(integral_wrapper_main_function.main_source):
    """
    Derived class providing some contracted-integral-specific methods.
    Should be inherited by algorithm-specific algo_main_source classes in each
    algo_module (e.g. contracted_vrr_only_no_auxiliary).
    """
    def vrr_sequence_no_auxiliary_out(self,p,wvrr_group_list,\
            auto_prim_loop_variable_list,auto_prim_loop_depends_dict,\
            work_array_index0, work_array_index1,\
            work_array_index0_initial_expr ):
        """
        Output code for a sequence of VRRs with no auxiliary indexes (all
        indexes involved are Cartesian Gaussians), with outer loops over 
        primitive exponents.

        p:                      cprinter object for code output
        wvrr_group_list:        list of whrr_group objects
                                (see integral_wrapper.algorithm class).
        auto_prim_loop_variable_list
                                list of dsl_variable objects to be assigned
                                inside loops over primitive indexes, ordered
                                such that output in this order will satisy
                                interdependencies in assignments
        auto_prim_loop_depends_dict
                                dictionary, with auto_prim_loop_variable_list
                                items as keys and values are sets of per-primitive
                                dsl_variable objects that the keys depend on,
                                (i.e. exponents of Cartesian Gaussians)
        work_array_index{0,1}
                                dsl_scalar variables which are integer pointers
                                into the work_array (e.g. iwrk{0,1})
        work_array_index0_initial_expr
                                dsl_binop/unop expression which contains the 
                                expression (integer) which represents the start
                                point in the work array for the block of
                                integrals corresponding to work_array_index0
        """
        ### Loops over primitives ###
        # If contracted integrals are required, insert loops over primitive exponents
        nloops_prim = 0
        prim_depend_set = set() # set of exponent dsl_* objects that have been
                                # "set" inside the loop structure
        iauto_prim = 0
        unassigned_list = auto_prim_loop_variable_list[:]
        final_wvrr_group = wvrr_group_list[-1]
        for index in final_wvrr_group.end_order_list:
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            prim_array_index = windex.array_index_dict()[ 'prim' ]
            nprim            = windex.length_dict()[ 'prim' ]
            start     = prim_array_index.assign( 0 )
            condition = dsl_binop( op_lt, prim_array_index, nprim )
            iteration = prim_array_index.assign( prim_array_index + 1 )
            p.forloop( start, condition, iteration )
            prim_depend_set.add( windex.index().exp() )

            ### Automatic variable assignment  ###
            # (for variables depending on primitive exponents )
            # The following loop structure allows auto_prim_loop_variable_list,
            # (which contains a list of automatically assigned variables, in an
            # order where sequential assignment satisfies interdependencies)
            # to be iterated through in such a way that variables are output
            # when their dependency on a loop index is satisifed.
            # This avoids all variables being inefficiently assigned for each
            # cycle of the innermost loop.
            # This is based on the fact that any variable which recursively 
            # depends on a set of per-primitive exponents
            # cannot depend on a variable which depends on a superset of these
            # e.g. a variable depending on { xa, xb }, cannot depend on any
            # variable depending on {xa, xb, xc}, though it can depend on 
            # variables depending on {xa} or {xb}.
            tmp_unassigned_list = unassigned_list[:]
            for dvar in tmp_unassigned_list:
                dvar_depend_set = auto_prim_loop_depends_dict[ dvar ]
                # If loops over primitive exponent indexes have been output, then
                # prim_depend_set contains the exponents. 
                # If dvar depends (recursively, via expr attributes) on
                # primitive exponents, then these feature in dvar_depend_set
                # If prim_depend_set is a superset of dvar_depend_set, then all
                # loops necessary to evaluate dvar have been output, and dvar
                # can be assigned and removed from unassigned_list.
                # If dvar_depend_set is a superset of prim_depend_set, then
                # further loops over primitive exponent indexes must be output.
                if prim_depend_set >= dvar_depend_set:
                    assert hasattr(dvar,'expr') and dvar.expr() != None,\
                    "all autoassign variables must have expr attribute != None."
                    if dvar.is_cartesian() == False:
                        p.out( dvar.assign( dvar.expr() ) )
                    elif dvar.is_cartesian() == True:
                        for d in [ 0, 1, 2 ]:
                            p.converter.direction().set(d)
                            p.out( dvar.assign( dvar.expr() ) )
                    unassigned_list.remove( dvar )
                    iauto_prim += 1
            nloops_prim += 1

        assert iauto_prim == len( auto_prim_loop_variable_list ),\
            "Number of variables output is inconsistent with length of list: "+\
            "all automatically assigned variables must be output."

        p.out( work_array_index0.assign( work_array_index0_initial_expr ) )
        # Use work_array_index1 to store work_array_index0_initial_expr while
        # work_array_index0 is incremented
        p.out( work_array_index1.assign( work_array_index0 ) )

        ### Base class evaluation ###
        # Output call to integral class base() function
        wfunction = self._wintegral.base_function()
        # Contracted integrals are required --- modify the arguments used for
        # calling VRR functions so that correct elements of exponent array are
        # passed
        tmp_arg_list = self._wintegral.contraction_info_obj().contraction_convert_list(\
                        wfunction.f_dsl().args() )
        # Replace work_array_index_base with work_array_index0, since no auxiliary indexes
        # are present
        work_array_index_base = \
                self._wintegral.work_array().special_array_index_dict()['base']
        arg_list = []
        for arg in tmp_arg_list:
            if arg is work_array_index_base:
                arg_list.append( work_array_index0 )
            else:
                arg_list.append( arg )
        # Call to rr_function with new_arg_list
        wfunction.call_out(p,arg_list)

        ### VRR sequence ###
        ngroups = 0
        for wrr_group in wvrr_group_list:
            # VRR evaluation
            iw = 0 # keep track of number of wrrs that have been output
            for wrr in wrr_group.wrr_list:
                iw += 1
                assert wrr.unrolled_index() is wrr.changing_index(),\
                        "RR must be vertical recurrence relation"
                wfunction = wrr.rr_function()
                changing_windex = self._wintegral.windex_index_dict()[ \
                                            wrr.changing_index().name() ]

                nloops = 0
                # Output loop structure top part
                # ... and set work_array_index() value in last loop before unrolled loop.
                # TODO 
                # This indexing could be done more efficiently, using multiple counters
                # and addition, rather than multiplication, but for testing purposes,
                # the simpler multiplication will be used.
                loop_index_list = list( reversed( wrr.loop_index_list() ) )
                work_array_index_expr = work_array_index1
                if len( loop_index_list ) > 1:
                    for index in loop_index_list[:-1]:
                        nloops += 1
                        windex = self._wintegral.windex_index_dict()[ index.name() ]
                        counter     = windex.array_index()
                        loop_start  = 0 
                        loop_length = windex.work_array_length()
                        p.forloop( counter.assign(loop_start), \
                                dsl_binop(op_lt, counter , loop_length ), \
                                counter.assign( dsl_binop( op_add, counter, 1 ) ) )
                        # Accumulate indexing expressions in work_array_index_expr
                        work_array_index_expr = work_array_index_expr + \
                                windex.array_index() * windex.work_array_skip()
                    
                    # ... set work_array_index0 value in last loop before unrolled loop.
                    p.out( work_array_index0.assign( work_array_index_expr ) )
                
                if len( loop_index_list ) > 0:
                    # Output inner loop over calls to vrr function
                    index = loop_index_list[-1]
                    nloops += 1
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    counter     = windex.array_index()
                    loop_start  = 0 
                    loop_length = windex.work_array_length()
                    p.forloop( counter.assign(loop_start), \
                        dsl_binop(op_lt, counter , loop_length ), \
                        counter.assign( dsl_binop( op_add, counter, 1 ) ) )

                # Call to rr_function
                # Contracted integrals are required --- modify the arguments used for
                # calling VRR functions so that correct elements of exponent array are
                # passed
                tmp_arg_list = self._wintegral.contraction_info_obj().contraction_convert_list(\
                            wfunction.f_dsl().args() )
                # Replace lmax_loop_var with l_expr from wrr_group, since this is VRR only with
                # no auxiliary index
                lmax_loop_var  = changing_windex.index_loop_max()
                lmax_loop_expr = wrr_group.end_max_l_expr_dict[ changing_windex.index().name() ]
                arg_list = []
                for arg in tmp_arg_list:
                    if arg is lmax_loop_var:
                        #arg_list.append( changing_windex.index_value() )
                        arg_list.append( lmax_loop_expr )
                    else:
                        arg_list.append( arg )
                # Call to rr_function with new_arg_list
                wfunction.call_out(p,arg_list)

                # Update work_array_index()
                if len( loop_index_list ) > 0:
                    last_loop_index = list( reversed( wrr.loop_index_list() ) )[-1]
                    windex = self._wintegral.windex_index_dict()[ last_loop_index.name() ]
                    last_loop_skip  = windex.work_array_skip()
                    p.out( work_array_index0.assign( work_array_index0 + last_loop_skip ) )
                # Output loop structure bottom part
                for loop in range( nloops ) :
                    p.endforloop()

                # If this is not the last VRR in the wvrr_group_list, then reset
                # work_array_index0 back to initial_expr (using work_array_index1)
                if iw != len( wrr_group.wrr_list ):
                    p.out( work_array_index0.assign( work_array_index1 ) )

        ### Close loops over primitives ###
        # If contracted integrals are required, close loops over primitives
        if self.generator().generator_options().contracted() == True:
            for loop in range( nloops_prim ):
                p.endforloop()

    def vrr_sequence_with_auxiliary_out(self, p, wvrr_group_list,\
            auto_prim_loop_variable_list,auto_prim_loop_depends_dict,\
            work_array_index0, work_array_index1, work_array_index_base, \
            work_array_index_start0, work_array_index_start1, \
            work_array_index_prim0, work_array_index_prim1,\
            work_array_index0_initial_expr, work_array_index1_initial_expr ):
        """
        Output code for a sequence of VRRs with no auxiliary indexes (all
        indexes involved are Cartesian Gaussians), with outer loops over 
        primitive exponents.

        p:                      cprinter object for code output
        wvrr_group_list:        list of whrr_group objects
                                (see integral_wrapper.algorithm class).
        auto_prim_loop_variable_list
                                list of dsl_variable objects to be assigned
                                inside loops over primitive indexes, ordered
                                such that output in this order will satisy
                                interdependencies in assignments
        auto_prim_loop_depends_dict
                                dictionary, with auto_prim_loop_variable_list
                                items as keys and values are sets of per-primitive
                                dsl_variable objects that the keys depend on,
                                (i.e. exponents of Cartesian Gaussians)
        work_array_index{0,1,base}
        work_array_index_start{0,1}
        work_array_index_prim{0,1}
                                dsl_scalar variables which are integer pointers
                                into the work_array (e.g. iwrk{0,1,b})
        work_array_index{0,1}_initial_expr
                                dsl_binop/unop expression which contains the 
                                expression (integer) which represents the start
                                point in the work array for the block of
                                integrals corresponding to work_array_index{0,1}
        """

        ### Loops over primitives ###
        # If contracted integrals are required, insert loops over primitive exponents
        nloops_prim = 0
        prim_depend_set = set() # set of exponent dsl_* objects that have been
                                # "set" inside the loop structure
        iauto_prim = 0
        unassigned_list = auto_prim_loop_variable_list[:]
        final_wvrr_group = wvrr_group_list[-1]
        for index in final_wvrr_group.end_order_list:
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            prim_array_index = windex.array_index_dict()[ 'prim' ]
            nprim            = windex.length_dict()[ 'prim' ]
            start     = prim_array_index.assign( 0 )
            condition = dsl_binop( op_lt, prim_array_index, nprim )
            iteration = prim_array_index.assign( prim_array_index + 1 )
            p.forloop( start, condition, iteration )
            prim_depend_set.add( windex.index().exp() )

            ### Automatic variable assignment  ###
            # (for variables depending on primitive exponents )
            # The following loop structure allows auto_prim_loop_variable_list,
            # (which contains a list of automatically assigned variables, in an
            # order where sequential assignment satisfies interdependencies)
            # to be iterated through in such a way that variables are output
            # when their dependency on a loop index is satisifed.
            # This avoids all variables being inefficiently assigned for each
            # cycle of the innermost loop.
            # This is based on the fact that any variable which recursively 
            # depends on a set of per-primitive exponents
            # cannot depend on a variable which depends on a superset of these
            # e.g. a variable depending on { xa, xb }, cannot depend on any
            # variable depending on {xa, xb, xc}, though it can depend on 
            # variables depending on {xa} or {xb}.
            tmp_unassigned_list = unassigned_list[:]
            for dvar in tmp_unassigned_list:
                dvar_depend_set = auto_prim_loop_depends_dict[ dvar ]
                # If loops over primitive exponent indexes have been output, then
                # prim_depend_set contains the exponents. 
                # If dvar depends (recursively, via expr attributes) on
                # primitive exponents, then these feature in dvar_depend_set
                # If prim_depend_set is a superset of dvar_depend_set, then all
                # loops necessary to evaluate dvar have been output, and dvar
                # can be assigned and removed from unassigned_list.
                # If dvar_depend_set is a superset of prim_depend_set, then
                # further loops over primitive exponent indexes must be output.
                if prim_depend_set >= dvar_depend_set:
                    assert hasattr(dvar,'expr') and dvar.expr() != None,\
                    "all autoassign variables must have expr attribute != None."
                    if dvar.is_cartesian() == False:
                        p.out( dvar.assign( dvar.expr() ) )
                    elif dvar.is_cartesian() == True:
                        for d in [ 0, 1, 2 ]:
                            p.converter.direction().set(d)
                            p.out( dvar.assign( dvar.expr() ) )
                    unassigned_list.remove( dvar )
                    iauto_prim += 1
            nloops_prim += 1

        assert iauto_prim == len( auto_prim_loop_variable_list ),\
            "Number of variables output is inconsistent with length of list: "+\
            "all automatically assigned variables must be output."

        # Set iwrk0_prim and iwrk1_prim to iwrk0_start and iwrk1_start so that
        # the final result of the VRR sequence is always placed in the same "layer"
        p.out( work_array_index_prim0.assign( work_array_index_start0 ) )
        p.out( work_array_index_prim1.assign( work_array_index_start1 ) )
        # Set iwrk0 and iwrk1 to initial values (per-primitive)
        p.out( work_array_index0.assign( work_array_index0_initial_expr ) )
        p.out( work_array_index1.assign( work_array_index1_initial_expr ) )


        ### Base class evaluation ###
        # Output call to integral class base() function
        wfunction = self._wintegral.base_function()
        wfunction.call_out(p)

        ### VRR sequence ###
        ngroups = 0
        for wrr_group in wvrr_group_list:
            # [ Current implementation restriction ]
            # This relies on the assumptions:
            # * Only one auxiliary index supported per integral class.
            # * Auxiliary index may only have increments of 0 or 1.
            # * Only one unique VRR supported per integral class.
            # ... i.e. if the VRR expression depends on an aux index, then all
            #          VRRs for integral indexes depend on this.
            
            assert wrr_group.m_layers in [1,2],\
                    "Only single of double-layer algorithm supported for VRRs with auxiliary indexes."
            if wrr_group.m_layers == 1:
                # Single-layer, with auxiliary
                # Implementation restriction -- must be the last wvrr_group, with a single wrr
                assert len( wrr_group.wrr_list ) == 1,\
                        "Only a single wrr can feature in a m_layer == 1 wvrr_group operation currently."
                assert wrr_group is wvrr_group_list[-1],\
                        "Currently, only a single m_layer == 1 wvrr_group operation is supported, and "+\
                        "this must be the final wvrr_group."

                # VRR evaluation
                iw = 0 # keep track of number of wrrs that have been output
                for wrr in wrr_group.wrr_list:
                    iw += 1
                    assert wrr.unrolled_index() is wrr.changing_index(),\
                            "RR must be vertical recurrence relation"
                    wfunction = wrr.rr_function()
                    changing_windex = self._wintegral.windex_index_dict()[ \
                                                wrr.changing_index().name() ]

                    nloops = 0
                    # Output loop structure top part
                    # ... and set work_array_index() value in last loop before unrolled loop.
                    # TODO 
                    # This indexing could be done more efficiently, using multiple counters
                    # and addition, rather than multiplication, but for testing purposes,
                    # the simpler multiplication will be used.
                    loop_index_list = list( reversed( wrr.loop_index_list() ) )
                    work_array_index_expr = work_array_index1
                    if len( loop_index_list ) > 1:
                        for index in loop_index_list[:-1]:
                            nloops += 1
                            windex = self._wintegral.windex_index_dict()[ index.name() ]
                            counter     = windex.array_index()
                            loop_start  = 0 
                            loop_length = windex.work_array_length()
                            p.forloop( counter.assign(loop_start), \
                                    dsl_binop(op_lt, counter , loop_length ), \
                                    counter.assign( dsl_binop( op_add, counter, 1 ) ) )
                            # Accumulate indexing expressions in work_array_index_expr
                            work_array_index_expr = work_array_index_expr + \
                                    windex.array_index() * windex.work_array_skip()
                        
                        # ... set work_array_index0 value in last loop before unrolled loop.
                        p.out( work_array_index0.assign( work_array_index_expr ) )
                    
                    if len( loop_index_list ) > 0:
                        # Output inner loop over calls to vrr function
                        index = loop_index_list[-1]
                        nloops += 1
                        windex = self._wintegral.windex_index_dict()[ index.name() ]
                        counter     = windex.array_index()
                        loop_start  = 0 
                        loop_length = windex.work_array_length()
                        p.forloop( counter.assign(loop_start), \
                            dsl_binop(op_lt, counter , loop_length ), \
                            counter.assign( dsl_binop( op_add, counter, 1 ) ) )

                    # Call to rr_function
                    # Contracted integrals are required --- modify the arguments used for
                    # calling VRR functions so that correct elements of exponent array are
                    # passed
                    tmp_arg_list = self._wintegral.contraction_info_obj().contraction_convert_list(\
                                wfunction.f_dsl().args() )
                    # Replace lmax_loop_var with l_expr from wrr_group, since this is VRR only with
                    # no auxiliary index
                    lmax_loop_var  = changing_windex.index_loop_max()
                    lmax_loop_expr = wrr_group.end_max_l_expr_dict[ changing_windex.index().name() ]
                    arg_list = []
                    for arg in tmp_arg_list:
                        if arg is lmax_loop_var:
                            #arg_list.append( changing_windex.index_value() )
                            arg_list.append( lmax_loop_expr )
                        else:
                            arg_list.append( arg )
                    # Call to rr_function with new_arg_list
                    wfunction.call_out(p,arg_list)

                    # Update work_array_index()
                    if len( loop_index_list ) > 0:
                        last_loop_index = list( reversed( wrr.loop_index_list() ) )[-1]
                        windex = self._wintegral.windex_index_dict()[ last_loop_index.name() ]
                        last_loop_skip  = windex.work_array_skip()
                        p.out( work_array_index0.assign( work_array_index0 + last_loop_skip ) )
                    # Output loop structure bottom part
                    for loop in range( nloops ) :
                        p.endforloop()

                # Set prim index back to starting point
                p.out( work_array_index0.assign( work_array_index_prim0 ) )

            elif wrr_group.m_layers == 2:
                # Double-layer, with auxiliary
                # Outside loops over auxiliary indexes
                nloops_aux = 0
                if wrr_group.auxiliary_index != None:
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
                    p.out( work_array_index_prim0.assign( work_array_index0 ) )
                    p.out( work_array_index_prim1.assign( work_array_index1 ) )
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

                    ## Set loop length variables for outer loops
                    #for index in reversed( wrr.loop_index_list() ):
                    #    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    #    loop_length_var  = windex.loop_length()
                    #    loop_length_expr = wrr_group.loop_length_expr_dict[ windex.index().name() ]
                    #    p.out( loop_length_var.assign( loop_length_expr ) )

                    # Set index_max() value for unrolled loop (inside function)
                    lmax_loop_var = changing_windex.index_loop_max()
                    lmax_loop_expr = wrr_group.index_loop_max_expr_dict[ changing_windex.index().name() ]
                    p.out( lmax_loop_var.assign( lmax_loop_expr ) ) 

                    nloops = 0
                    # Output loop structure top part
                    # ... and set work_array_index() value in last loop before unrolled loop.
                    # TODO 
                    # This indexing could be done more efficiently, using multiple counters
                    # and addition, rather than multiplication, but for testing purposes,
                    # the simpler multiplication will be used.
                    loop_index_list = list( reversed( wrr.loop_index_list() ) )
                    work_array_index0_expr = work_array_index_prim0
                    work_array_index1_expr = work_array_index_prim1
                    if len( loop_index_list ) > 1:
                        for index in loop_index_list[:-1]:
                            nloops += 1
                            windex = self._wintegral.windex_index_dict()[ index.name() ]
                            counter     = windex.array_index()
                            loop_start  = 0 
                            loop_length = windex.work_array_length()
                            p.forloop( counter.assign(loop_start), \
                                    dsl_binop(op_lt, counter , loop_length ), \
                                    counter.assign( dsl_binop( op_add, counter, 1 ) ) )
                            # Accumulate indexing expressions in work_array_index_expr
                            work_array_index0_expr = work_array_index0_expr + \
                                    windex.array_index() * windex.work_array_skip()
                            work_array_index1_expr = work_array_index1_expr + \
                                    windex.array_index() * windex.work_array_skip()
                        
                        # ... set work_array_index0 value in last loop before unrolled loop.
                        p.out( work_array_index0.assign( work_array_index0_expr) )
                        p.out( work_array_index1.assign( work_array_index1_expr) )
                    
                    if len( loop_index_list ) > 0:
                        # Output inner loop over calls to vrr function
                        index = loop_index_list[-1]
                        nloops += 1
                        windex = self._wintegral.windex_index_dict()[ index.name() ]
                        counter     = windex.array_index()
                        loop_start  = 0 
                        loop_length = windex.work_array_length()
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
                        # Do not reset if this is the last in wrr_list, since they 
                        # will be set during swapping 
                        if iw != len( wrr_group.wrr_list ):
                            p.out( work_array_index0.assign( work_array_index_prim0 ) )
                            p.out( work_array_index1.assign( work_array_index_prim1 ) )

                # The following code assumes that only increments of 0 and 1 are
                # present and that only a single auxiliary index is present for
                # each integral class
                # Loops over auxiliary index
                assert len( self._wintegral.windex_aux_list() ) == 1,\
                        "Only a single auxiliary index is currently supported per integral."
                assert max( list( work_array_index_dict.keys() ) ) <= 1,\
                        "Only increments of 0 and 1 supported in auxiliary index"
                # Swap prim layer pointers for next loop over auxiliary index
                p.out( work_array_index0.assign( work_array_index_prim1 ) )
                p.out( work_array_index1.assign( work_array_index_prim0 ) )
                # Output loops structure bottom part (auxiliary indexes)
                for loop in range( nloops_aux ):
                    p.endforloop()

                # The following code assumes that only increments of 0 and 1 are
                # present and that only a single auxiliary index is present for
                # each integral class
                # Unswap prim layers for next loop structure over auxiliary index
                p.out( work_array_index0.assign( work_array_index_prim0 ) )
                p.out( work_array_index1.assign( work_array_index_prim1 ) )
                ngroups += 1

        ### Close loops over primitives ###
        # If contracted integrals are required, close loops over primitives
        if self.generator().generator_options().contracted() == True:
            for loop in range( nloops_prim ):
                p.endforloop()


    def early_transform_sequence_out(self,p,early_transform_list,\
                    work_array_index0, work_array_index_start0, \
                    work_array_index1, work_array_index_start1, \
                    last_transform_swap):
        """
        Output code for a sequence of contractions and spherical transformations
        which precede a HRR (if present).

        p:                      cprinter object for code output
        early_transform_list:   list of contraction or spherical_transform objects
                                (see integral_wrapper.algorithm class).
        work_array_index{0,1}
        work_array_index_start{0,1}
                                dsl_scalar variables which are integer pointers
                                into the work_array (i.e. iwrk{0,1} and
                                iwrk{0,1}_start).
        last_transform_swap:    bool, determines whether code which swaps 
                                work_array_index{0,1} values between 
                                work_array_index_start{0,1} values is output
                                before exiting the method.

        Implementation restrictions:
        * Contractions should always precede spherical transformations, so that
          all contractions are completed before spherical transformation occurs.
        """
        # Early transformation sequence occurs before a HRR, if present
        # Early transformation sequence contains spherical transformations and 
        # contractions
        ### Early transform sequence ###
        contracted_index_list = []
        sph_transed_index_list = []
        contraction_list = []
        sph_trans_list   = []
        # Separate contractions and spherical transforms and check that contractions
        # always precede spherical transformations
        sph_trans_found = False
        for t in early_transform_list:
            if t.optype == 'contraction':
                assert sph_trans_found == False, \
                        "spherical transformations cannot precede contractions"
                contraction_list.append( t )
            elif t.optype == 'spherical transform':
                sph_trans_list.append( t )
        assert len( contraction_list ) == len( self._wintegral.windex_index_list() ),\
                "Number of contraction operations must be the same as number of integral indexes."
        # First transformation is special --- angular momentum ranges of integrals are truncated
        # This is always a contraction
        t = early_transform_list[0]
        assert t in contraction_list, \
                "The first transformation should always be a contraction."
        contracted_index_list.append(t.index_to_contract)
        for index in t.loop_index_list:
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            iskip     = windex.work_array_skip()
            tmp_iskip = windex.skip_dict()['tmp']
            p.out( tmp_iskip.assign( iskip ) )
        # Reset work_array_length variables (ncc_type is 'ncc', 'ncc0-ncc0') for 
        # all indexes
        for index in t.end_index_list:
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            nw       = windex.work_array_length()
            l        = t.end_max_l_expr_dict[ index.name() ]
            if t in contraction_list:
                ncc_type = t.end_ncc_type_dict[ index.name() ]
                assert ncc_type in [ 'ncc', 'ncc0-ncc0' ],\
                        "After first contraction, all indexes should have "+\
                        "ncc_type == 'ncc' or 'ncc0-ncc0'."
                if ncc_type == 'ncc':
                    l_expr   = (l+1)*(l+2)/2
                elif ncc_type == 'ncc0-ncc0':
                    l_expr   = (l+1)*(l+2)*(l+3)/6 - windex.work_array_offset()
            p.out( nw.assign( l_expr ) )
        # Set work_array_skip variables (work_array_index0 layer) for all indexes
        index = t.end_order_list[-1]
        windex = self._wintegral.windex_index_dict()[ index.name() ]
        work_array_skip_var  = windex.work_array_skip() # skip for each Cartesian component
        work_array_skip_expr = 1
        p.out( work_array_skip_var.assign( work_array_skip_expr ) )
        for index in reversed( t.end_order_list[:-1] ):
            if work_array_skip_expr == 1:
                work_array_skip_expr = windex.work_array_length() 
            else:
                work_array_skip_expr = work_array_skip_expr * windex.work_array_length()
            if windex.index() in contracted_index_list:
                work_array_skip_expr = work_array_skip_expr * windex.length_dict()['cont']
            else:
                work_array_skip_expr = work_array_skip_expr * windex.length_dict()['prim']
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            work_array_skip_var = windex.work_array_skip()
            p.out( work_array_skip_var.assign( work_array_skip_expr ) )

        # Output loops over indexes in t.loop_index_list
        # Since we are truncating the angular momentum range, we must explicitly
        # loop over the Cartesian/spherical component and primitive indexes
        # For first transformation, we can assume that no indexes have yet been
        # contracted.
        for k, v in t.start_cont_dict.items():
            assert v == 'prim', \
                    "at start of first transformation in early transform list, "+\
                    "all indexes should be uncontracted (primitive)."
        cont_str = 'prim'
        nloops = 0
        for index in t.loop_index_list:
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            counter_nprim     = windex.array_index_dict()[cont_str]
            loop_start_nprim  = 0 
            loop_length_nprim = windex.length_dict()[cont_str]
            counter_ncomp       = windex.array_index()
            loop_start_ncomp    = 0
            loop_length_ncomp   = windex.work_array_length() 
            # Loop over primitives
            nloops += 1
            p.forloop( counter_nprim.assign(loop_start_nprim),\
                    dsl_binop(op_lt, counter_nprim , loop_length_nprim ),\
                    counter_nprim.assign( dsl_binop( op_add, counter_nprim, 1 ) ) )
            # Loop over components
            nloops += 1
            p.forloop( counter_ncomp.assign(loop_start_ncomp),\
                    dsl_binop(op_lt, counter_ncomp , loop_length_ncomp ),\
                    counter_ncomp.assign( dsl_binop( op_add, counter_ncomp, 1 ) ) )

        # First transformation -- special, since angular momentum ranges are
        # truncated and indexing must be done by multiplication to "cut"
        # out unecessary blocks of integrals
        if len( t.loop_index_list ) > 0:
            index = t.loop_index_list[0]
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            windex_to_contract = self._wintegral.windex_index_dict()[ t.index_to_contract.name() ]
            index_to_contract_offset_var = windex_to_contract.work_array_offset()
            index_expr0 = work_array_index_start0 + \
                          windex.array_index() * windex.work_array_skip() + \
                          windex.array_index_dict()[cont_str] * \
                          windex.work_array_skip() * windex.work_array_length()
            index_expr1 = work_array_index_start1 + index_to_contract_offset_var  +\
                          ( windex.array_index() + windex.work_array_offset() )*\
                          windex.skip_dict()['tmp'] + \
                          windex.array_index_dict()[cont_str] * windex.skip_dict()[cont_str]
            for index in t.loop_index_list[1:]:
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                index_expr0 = index_expr0 + \
                              windex.array_index() * windex.work_array_skip() + \
                              windex.array_index_dict()[cont_str] * \
                              windex.work_array_skip() * windex.work_array_length()
                index_expr1 = index_expr1 + \
                              ( windex.array_index() + windex.work_array_offset() )*\
                              windex.skip_dict()['tmp'] + \
                              windex.array_index_dict()[cont_str] * windex.skip_dict()[cont_str]
        else:
            # No additional indexes to loop over
            windex_to_contract = self._wintegral.windex_index_dict()[ t.index_to_contract.name() ]
            index_to_contract_offset_var = windex_to_contract.work_array_offset()
            index_expr0 = work_array_index_start0
            index_expr1 = work_array_index_start1 + index_to_contract_offset_var

        # Update work_array_index variables
        p.out( work_array_index0.assign( index_expr0 ) )
        p.out( work_array_index1.assign( index_expr1 ) )

        # Output call to contract function
        t.wfunction.call_out(p)

        # Close loops over primitive and Cartesian component indexes
        for loop in range( nloops ) :
            p.endforloop()

        if len( early_transform_list ) > 0 :
            # Swap layers
            p.out( work_array_index0.assign( work_array_index_start1 ) )
            p.out( work_array_index1.assign( work_array_index_start0 ) )

        # After first transformation, ncc_type should be 'ncc' or 'ncc0-ncc0'
        # for all indexes
        for index in t.end_index_list :
            windex   = self._wintegral.windex_index_dict()[ index.name() ]
            ncc_type = t.end_ncc_type_dict[ windex.index().name() ]
            # Determine work_array_length_expr and work_array_offset_expr 
            # from ncc_type
            assert ncc_type in [ 'ncc', 'ncc0-ncc0' ],\
                    "ncc_type must be 'ncc' or 'ncc0-ncc0' after first transformation"

        # Subsequent transformations (contractions, spherical transformations)
        i = 1
        # Implementation restriction:
        # * Contractions should always precede spherical transformations.
        for t in early_transform_list[1:]:
            if t in sph_trans_list:
                assert len( contracted_index_list ) == len( contraction_list ),\
                        "All indexes should be contracted before spherical "+\
                        "transformation occurs."
                if self.generator().generator_options().spherical_transform_conditionals() == True:
                    # Output if block, so spherical transform not done if l < 2
                    windex = self._wintegral.windex_index_dict()[ t.index_to_transform.name() ]
                    l      = windex.index_value()
                    p.ifblock( dsl_binop( op_ge, l, 2 ) )
            # Set layer start points for this contraction
            p.out( work_array_index_start0.assign( work_array_index0 ) )
            p.out( work_array_index_start1.assign( work_array_index1 ) )
            # Set tmp skip variables (for work_array_index1 layer) for indexes being looped over
            for index in t.loop_index_list:
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                iskip     = windex.work_array_skip()
                tmp_iskip = windex.skip_dict()['tmp']
                p.out( tmp_iskip.assign( iskip ) )
            if t in contraction_list:
                # For contraction, the number of components does not change after the first
                # contraction (which truncates angular momentum ranges.
                contracted_index_list.append(t.index_to_contract)
            elif t in sph_trans_list:
                # For spherical transformation, the number of components changes, so we must set
                # the work_array_length variable
                sph_transed_index_list.append(t.index_to_transform)
                windex    = self._wintegral.windex_index_dict()[ t.index_to_transform.name() ]
                nw        = windex.work_array_length()
                nw_tmp    = windex.length_dict()['tmp']
                iwskip    = windex.work_array_skip()
                tmp_iskip = windex.skip_dict()['tmp']
                l         = t.end_max_l_expr_dict[ t.index_to_transform.name() ]
                assert t.end_ncc_type_dict[ t.index_to_transform.name() ] == 'nsph',\
                        "After sph_trans, transformed index should have ncc_type == 'nsph'."
                l_expr   = 2*l + 1
                # Assign nw_tmp to previous nw value
                p.out( nw_tmp.assign( nw ) )
                # Assign nw to new number of spherical component
                p.out( nw.assign( l_expr ) )
                # Assign iskipX_tmp for index_to_transform
                p.out( tmp_iskip.assign( iwskip ) )
            # Set work_array_skip variables (work_array_index0 layer) for all indexes
            index = t.end_order_list[-1]
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            cont_str = t.end_cont_dict[ windex.index().name() ]
            work_array_skip_var  = windex.work_array_skip() # skip for each Cartesian component
            work_array_skip_expr = 1
            p.out( work_array_skip_var.assign( work_array_skip_expr ) )
            for index in reversed( t.end_order_list[:-1] ):
                if work_array_skip_expr == 1:
                    work_array_skip_expr = windex.work_array_length() 
                else:
                    work_array_skip_expr = work_array_skip_expr * windex.work_array_length()
                work_array_skip_expr = work_array_skip_expr * windex.length_dict()[cont_str]
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                cont_str = t.end_cont_dict[ windex.index().name() ]
                work_array_skip_var = windex.work_array_skip()
                p.out( work_array_skip_var.assign( work_array_skip_expr ) )

            # If spherical transform, set iwskipX_cont and iskipX_cont_tmp for index_to_transform
            if t in sph_trans_list:
                windex = self._wintegral.windex_index_dict()[ t.index_to_transform.name() ]
                p.out( windex.skip_dict()[ 'contwork' ].assign( windex.work_array_skip() * windex.work_array_length() ) )
                p.out( windex.skip_dict()[ 'conttmp' ].assign( windex.skip_dict()['tmp'] * windex.length_dict()['tmp'] ) )

            # Output loops over indexes in t.loop_index_list
            nloops = 0
            if len( t.loop_index_list ) > 0:
                for index in t.loop_index_list:
                    nloops += 1
                    cont_str = t.end_cont_dict[ index.name() ]
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    counter     = windex.array_index()
                    loop_start  = 0
                    loop_length = windex.length_dict()[cont_str] * windex.work_array_length()
                    # contains information about work_array size
                    p.forloop( counter.assign(loop_start),\
                            dsl_binop(op_lt, counter , loop_length ),\
                            counter.assign( dsl_binop( op_add, counter, 1 ) ) )
            if t in sph_trans_list:
                # Add additional loop over contractions of windex_to_transform
                nloops += 1
                windex = self._wintegral.windex_index_dict()[ t.index_to_transform.name() ]
                cont_str = t.end_cont_dict[ index.name() ]
                assert cont_str == 'cont', "cont_str must be 'cont' for index_to_transform"
                counter     = windex.array_index_dict()['cont']
                loop_start  = 0
                loop_length = windex.length_dict()[cont_str]
                p.forloop( counter.assign(loop_start),\
                        dsl_binop(op_lt, counter , loop_length ),\
                        counter.assign( dsl_binop( op_add, counter, 1 ) ) )

            # Output call to contract function
            t.wfunction.call_out(p)

            # Update work_array_index0 and work_array_index1
            if nloops > 0:
                # Subsequent contractions -- can simply iteratively accumulate
                # indexing, since the integrals are continuous (no gaps to 
                # discard).
                if t in contraction_list:
                    assert len(t.loop_index_list) > 0,\
                            "loop_index_list for a contraction must be > 0 for work_array "+\
                            "index incrementing code to be output."
                    index = t.loop_index_list[-1]
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    # Add iwskip, iwskip_tmp
                    index_expr0 = work_array_index0 + windex.work_array_skip()
                    index_expr1 = work_array_index1 + windex.skip_dict()['tmp']
                elif t in sph_trans_list:
                    index = t.index_to_transform
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    # Add iwskip_cont, iwskip_cont_tmp
                    index_expr0 = work_array_index0 + windex.skip_dict()['contwork']
                    index_expr1 = work_array_index1 + windex.skip_dict()['conttmp']
                # Update work_array_index variables
                p.out( work_array_index0.assign( index_expr0 ) )
                p.out( work_array_index1.assign( index_expr1 ) )

            # Output bottom part of loops over indexes in c.loop_index_list
            for loop in range( nloops ) :
                p.endforloop()

            if i+1 != len( early_transform_list ) or last_transform_swap == True:
                # Swap layers
                p.out( work_array_index0.assign( work_array_index_start1 ) )
                p.out( work_array_index1.assign( work_array_index_start0 ) )

            if t in sph_trans_list and self.generator().generator_options().spherical_transform_conditionals() == True:
                # Output end of if block, so spherical transform not done if l < 2
                p.endifblock()
            i += 1

        # After final transformation, ncc_type should be 'ncc',  'nsph' or  'ncc0-ncc0' for all indexes
        for index in t.end_index_list :
            windex   = self._wintegral.windex_index_dict()[ index.name() ]
            ncc_type = t.end_ncc_type_dict[ windex.index().name() ]
            # Determine work_array_length_expr and work_array_offset_expr 
            # from ncc_type
            assert ncc_type in [ 'ncc', 'nsph', 'ncc0-ncc0' ],\
                    "ncc_type must be 'ncc', 'nsph', or 'ncc0-ncc0' for early contraction sequence"

        # At the end of the early transform sequence, all indexes should be
        # contracted
        if len( early_transform_list ) > 0:
            t = early_transform_list[-1]
            for index in t.end_order_list:
                assert t.end_cont_dict[ index.name() ] == 'cont',\
                        "All indexes must be contracted at end of early transform sequence."

    def late_transform_sequence_out(self,p,late_transform_list,\
                    work_array_index0, work_array_index_start0, \
                    work_array_index1, work_array_index_start1, \
                    last_transform_swap):
        """
        Output code for a sequence of spherical transformations, following 
        a HRR sequence.

        p:                      cprinter object for code output
        late_transform_list:    list of contraction or spherical_transform objects
                                (see integral_wrapper.algorithm class).
        work_array_index{0,1}
        work_array_index_start{0,1}
                                dsl_scalar variables which are integer pointers
                                into the work_array (i.e. iwrk{0,1} and
                                iwrk{0,1}_start).

        last_transform_swap:    bool, determines whether code which swaps 
                                work_array_index{0,1} values between 
                                work_array_index_start{0,1} values is output
                                before exiting the method.
        """
        # Late transformation sequence occurs after a HRR, if present
        # Late transformation sequence only contains spherical transformations
        sph_transed_index_list = []
        sph_trans_list   = []
        # Separate contractions and spherical transforms
        for t in late_transform_list:
            assert t.optype == "spherical transform",\
                "only spherical transformations may occur in the late contraction sequence (post-HRR)"
            sph_trans_list.append( t )
        
        # Spherical transformations
        i = 0
        for t in late_transform_list:
            if self.generator().generator_options().spherical_transform_conditionals() == True:
                # Output if block, so spherical transform not done if l < 2
                windex = self._wintegral.windex_index_dict()[ t.index_to_transform.name() ]
                l      = windex.index_value()
                p.ifblock( dsl_binop( op_ge, l, 2 ) )
            # Set layer start points for this contraction
            p.out( work_array_index_start0.assign( work_array_index0 ) )
            p.out( work_array_index_start1.assign( work_array_index1 ) )
            # Set tmp skip variables (for work_array_index1 layer) for indexes being looped over
            for index in t.loop_index_list:
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                iskip     = windex.work_array_skip()
                tmp_iskip = windex.skip_dict()['tmp']
                p.out( tmp_iskip.assign( iskip ) )
            # For spherical transformation, the number of components changes, so we must set
            # the work_array_length variable
            sph_transed_index_list.append(t.index_to_transform)
            windex    = self._wintegral.windex_index_dict()[ t.index_to_transform.name() ]
            nw        = windex.work_array_length()
            nw_tmp    = windex.length_dict()['tmp']
            iwskip    = windex.work_array_skip()
            tmp_iskip = windex.skip_dict()['tmp']
            l         = t.end_max_l_expr_dict[ t.index_to_transform.name() ]
            assert t.end_ncc_type_dict[ t.index_to_transform.name() ] == 'nsph',\
                    "After sph_trans, transformed index should have ncc_type == 'nsph'."
            l_expr   = 2*l + 1
            # Assign nw_tmp to previous nw value
            p.out( nw_tmp.assign( nw ) )
            # Assign nw to new number of spherical component
            p.out( nw.assign( l_expr ) )
            # Assign iskipX_tmp for index_to_transform
            p.out( tmp_iskip.assign( iwskip ) )
            # Set work_array_skip variables (work_array_index0 layer) for all indexes
            index = t.end_order_list[-1]
            windex = self._wintegral.windex_index_dict()[ index.name() ]
            cont_str = t.end_cont_dict[ windex.index().name() ]
            work_array_skip_var  = windex.work_array_skip() # skip for each Cartesian component
            work_array_skip_expr = 1
            p.out( work_array_skip_var.assign( work_array_skip_expr ) )
            for index in reversed( t.end_order_list[:-1] ):
                if work_array_skip_expr == 1:
                    work_array_skip_expr = windex.work_array_length() 
                else:
                    work_array_skip_expr = work_array_skip_expr * windex.work_array_length()
                work_array_skip_expr = work_array_skip_expr * windex.length_dict()[cont_str]
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                cont_str = t.end_cont_dict[ windex.index().name() ]
                work_array_skip_var = windex.work_array_skip()
                p.out( work_array_skip_var.assign( work_array_skip_expr ) )

            # If spherical transform, set iwskipX_cont and iskipX_cont_tmp for index_to_transform
            windex = self._wintegral.windex_index_dict()[ t.index_to_transform.name() ]
            p.out( windex.skip_dict()[ 'contwork' ].assign( windex.work_array_skip() * windex.work_array_length() ) )
            p.out( windex.skip_dict()[ 'conttmp' ].assign( windex.skip_dict()['tmp'] * windex.length_dict()['tmp'] ) )

            # Output loops over indexes in t.loop_index_list
            nloops = 0
            if len( t.loop_index_list ) > 0:
                for index in t.loop_index_list:
                    nloops += 1
                    cont_str = t.end_cont_dict[ index.name() ]
                    windex = self._wintegral.windex_index_dict()[ index.name() ]
                    counter     = windex.array_index()
                    loop_start  = 0
                    loop_length = windex.length_dict()[cont_str] * windex.work_array_length()
                    # contains information about work_array size
                    p.forloop( counter.assign(loop_start),\
                            dsl_binop(op_lt, counter , loop_length ),\
                            counter.assign( dsl_binop( op_add, counter, 1 ) ) )
            # Add additional loop over contractions of windex_to_transform
            nloops += 1
            windex = self._wintegral.windex_index_dict()[ t.index_to_transform.name() ]
            cont_str = t.end_cont_dict[ index.name() ]
            assert cont_str == 'cont', "cont_str must be 'cont' for index_to_transform"
            counter     = windex.array_index_dict()['cont']
            loop_start  = 0
            loop_length = windex.length_dict()[cont_str]
            p.forloop( counter.assign(loop_start),\
                    dsl_binop(op_lt, counter , loop_length ),\
                    counter.assign( dsl_binop( op_add, counter, 1 ) ) )

            # Output call to contract function
            t.wfunction.call_out(p)

            # Update work_array_index0 and work_array_index1
            if nloops > 0:
                # Subsequent contractions -- can simply iteratively accumulate
                # indexing, since the integrals are continuous (no gaps to 
                # discard).
                index = t.index_to_transform
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                # Add iwskip_cont, iwskip_cont_tmp
                index_expr0 = work_array_index0 + windex.skip_dict()['contwork']
                index_expr1 = work_array_index1 + windex.skip_dict()['conttmp']
                # Update work_array_index variables
                p.out( work_array_index0.assign( index_expr0 ) )
                p.out( work_array_index1.assign( index_expr1 ) )

            # Output bottom part of loops over indexes in c.loop_index_list
            for loop in range( nloops ) :
                p.endforloop()

            if i+1 != len( late_transform_list ) or last_transform_swap == True:
                # Swap layers
                p.out( work_array_index0.assign( work_array_index_start1 ) )
                p.out( work_array_index1.assign( work_array_index_start0 ) )

            if self.generator().generator_options().spherical_transform_conditionals() == True:
                # Output end of if block, so spherical transform not done if l < 2
                p.endifblock()
            i += 1

        # After final transformation, ncc_type should be 'nsph'
        for index in late_transform_list[-1].end_index_list :
            windex   = self._wintegral.windex_index_dict()[ index.name() ]
            ncc_type = t.end_ncc_type_dict[ windex.index().name() ]
            # Determine work_array_length_expr and work_array_offset_expr 
            # from ncc_typea
            assert ncc_type in [ 'nsph' ],\
                    "ncc_type must be 'nsph' for late contraction sequence"

    def hrr_sequence_out(self,p,whrr_group_list,\
                    work_array_index0, work_array_index_start0, \
                    work_array_index1, work_array_index_start1, \
                    no_hrr_swap):
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

        no_hrr_swap:            bool, determines whether code which swaps 
                                work_array_index{0,1} values between 
                                work_array_index_start{0,1} values when
                                no HRR is executed (i.e. l_move_to = 0 )
                                is output before exiting the method.
        """

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

            # Output work_array skip variables
            # Variables with lower angmom in move_to_index take values from previous iteration 
            # higher angmom layer, or from last operation in VRR sequence
            # Output dummy work_array_skip variable for move_to_index (this is not 
            # needed in the first loop over move_to_index, as l_move_to = 0, but
            # needs to be set to avoid using unpredictable behaviour from uninitialized
            # value)
            # move_to_index work_array per-contraction skip
            work_array_skip_var = move_to_windex.work_array_skip() 
            work_array_length_var = move_to_windex.work_array_length()
            hrr_cont_skip_var = move_to_windex.skip_dict()['conttmp']
            p.out( hrr_cont_skip_var.assign( work_array_skip_var * work_array_length_var ) )
            # move_from_index work_array per-contraction skip
            work_array_skip_var = move_from_windex.work_array_skip() 
            work_array_length_var = move_from_windex.work_array_length()
            hrr_cont_skip_var = move_from_windex.skip_dict()['conttmp']
            p.out( hrr_cont_skip_var.assign( work_array_skip_var * work_array_length_var ) )
            #p.out( work_array_skip_var.assign(0) )
            for index in wrr_group.end_index_list:
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                # move_to_index layer
                work_array_skip_var   = windex.work_array_skip()
                # move_to_index - 1 layer
                hrr_skip_var      = windex.skip_dict()['tmp']
                # hrr_length_var never needs to be set
                p.out( hrr_skip_var.assign( work_array_skip_var ) )
            # Output work_array_skip variable for move_from_index:
            assert wrr_group.end_index_list[-1] is move_from_windex.index(),\
                    "last index in end_index_list must be move_from_index"
            # If the above assertion is true, then the following assignment is unecessary!
            #work_array_skip_var = move_from_windex.work_array_skip() 
            #p.out( work_array_skip_var.assign(1) )

            # Initial work_array array_index assignments: 
            # If no auxiliary index in VRR sequence, iwrk0 and iwrk1 need swapping, 
            # since the result of the VRR sequence. If the VRR sequence involves
            # an auxiliary index and a 2-layer algorithm has been used,
            # then this swap will already have been done.
            # No auxiliary index present, iwrk0 and iwrk1 need swapping
            p.out( work_array_index0.assign( work_array_index_start1 ) )
            p.out( work_array_index1.assign( work_array_index_start0 ) )

            ## Standard angmom transfers ##
            loop_start = 1
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
            l_expr = move_to_loop_lmax_var 
            loop_length_var  = move_to_windex.work_array_length()
            loop_length_expr = (l_expr+1)*(l_expr+2)/2
            p.out( loop_length_var.assign( loop_length_expr ) )
            # move_from_index (special ncc expression, ncc0(l_mt+l_mf-l_mf_loop) - ncc0(l_mf) )
            l_expr1 = move_from_start_max_l - move_to_loop_lmax_var
            l_expr2 = move_from_end_max_l
            loop_length_var  = move_from_windex.work_array_length()
            offset_var       = move_from_windex.work_array_offset() # set in contraction sequence
            #loop_length_expr = (l_expr1+1)*(l_expr1+2)*(l_expr1+3)/6 -\
            #        l_expr2*(l_expr2+1)*(l_expr2+2)/6 
            loop_length_expr = (l_expr1+1)*(l_expr1+2)*(l_expr1+3)/6 - offset_var
            p.out( loop_length_var.assign( loop_length_expr ) )
            previous_windex = wrr_group.end_index_list[-1]
            skip_expr = loop_length_var * move_from_windex.length_dict()['cont']
            for index in reversed( wrr_group.end_index_list[:-1] ):
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                work_array_skip_var = windex.work_array_skip()
                loop_length_var  = windex.work_array_length()
                # Output iwskip value
                p.out( work_array_skip_var.assign( skip_expr ) )
                # Update skip_expr
                skip_expr = skip_expr * loop_length_var * windex.length_dict()['cont']
            # Only move_from and move_to_index need iwskip_cont to be set
            for windex in [ move_from_windex, move_to_windex ]:
                work_array_skip_var = windex.work_array_skip()
                loop_length_var     = windex.work_array_length()
                work_array_skip_cont_var = windex.skip_dict()['contwork']
                # Output iwskip_cont value
                p.out( work_array_skip_cont_var.assign( \
                        work_array_skip_var * loop_length_var ) )

            # Output loops over indexes in wrr.loop_index_list (not move_{to,from}_index)
            for index in wrr.loop_index_list():
                nloops += 1
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                counter     = windex.array_index()
                loop_start  = 0
                loop_length = windex.work_array_length() * windex.length_dict()['cont']
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

            # Swap iskipX_tmp values
            for index in wrr_group.end_index_list:
                windex = self._wintegral.windex_index_dict()[ index.name() ]
                # move_to_index layer
                work_array_skip_var   = windex.work_array_skip()
                # move_to_index - 1 layer
                hrr_skip_var   = windex.skip_dict()['tmp']
                # hrr_length_var never needs to be set
                p.out( hrr_skip_var.assign( work_array_skip_var ) )
            # Swap iwskipX_cont_tmp values for move_from_ and move_to_index
            for windex in [ move_from_windex, move_to_windex ]:
                work_array_skip_cont_tmp_var = windex.skip_dict()['conttmp']
                work_array_skip_cont_var = windex.skip_dict()['contwork']
                # Output iwskip_cont value
                p.out( work_array_skip_cont_tmp_var.assign( work_array_skip_cont_var ) )


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
                        "ncc_type must be 'ncc' or 'nsph' after HRR"

            # Only execute HRR if l_move_to > 0, otherwise just swap iwrk0, iwrk1
            if no_hrr_swap == True:
                p.elseblock()
                # Swap layers for subsequent operation sequence
                p.out( work_array_index0.assign( work_array_index_start1 ) )
                p.out( work_array_index1.assign( work_array_index_start0 ) )
            p.endifblock()

    def copy_to_output_array_out(self,p,final_operation,work_array_index_start0):
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
        # Since we allow iskip_cont to be freely set, we must loop over
        # contraction and cartesian component indexes separately
        for windex in output_windex_loop_list :
            counter_ncc       = windex.array_index() 
            loop_start_ncc    = 0
            loop_length_ncc   = windex.output_array_length()
            counter_ncont     = windex.array_index_dict()['cont']
            loop_start_ncont  = 0
            loop_length_ncont = windex.length_dict()['cont']
            # Loop over contraction index
            nloops += 1
            p.forloop( counter_ncont.assign(loop_start_ncont),\
                    dsl_binop(op_lt, counter_ncont , loop_length_ncont ),\
                    counter_ncont.assign( dsl_binop( op_add, counter_ncont, 1 ) ) )
            # Loop over Cartesian component index
            nloops += 1
            p.forloop( counter_ncc.assign(loop_start_ncc),\
                    dsl_binop(op_lt, counter_ncc , loop_length_ncc ),\
                    counter_ncc.assign( dsl_binop( op_add, counter_ncc, 1 ) ) )

        # Generate correct DSL expressions for output and work array indexing
        #work_array_index_expr = src_work_array_index0
        work_array_index_expr = work_array_index_start0
        output_array_index_expr = None           
        for i in range(0, len(output_windex_loop_list) ):
            windex = output_windex_loop_list[i]
            work_expr = ( windex.array_index() +\
                        windex.array_index_dict()['cont'] * windex.output_array_length() ) *\
                        windex.work_array_skip()
            output_expr = windex.array_index() * windex.output_array_skip() +\
                          windex.array_index_dict()['cont'] * windex.skip_dict()['cont']
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

