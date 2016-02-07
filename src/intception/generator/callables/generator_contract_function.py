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
Provides the contract_call, contract_prototype and contract_source 
classes which are derived from callable_* classes in intception.generator.wrappers.

The contract_source class contains the code necessary to generate a self-contained
function for performing a single contraction of a primitive Cartesian index.
"""

from intception.dsl import *
from intception.generator.wrappers import function_wrapper

class contract_call(function_wrapper.callable_function_call):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function calls for a contraction.
    """
    pass

class contract_prototype(function_wrapper.callable_prototype):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function prototypes for a contraction.
    """
    pass

class noblas_contract_source(function_wrapper.callable_source):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate executable source code for a contraction.
    """
    def __init__(self,f_dsl, \
                 work_array, contract_array, \
                 n_pre, n_prim, n_cont, iskip_prim, \
                 gen, const_dict = {} ):
        function_wrapper.callable_source.__init__(self, f_dsl, gen, const_dict = const_dict)
        self._work_array     = work_array     # general_array object
        self._contract_array = contract_array # general_array object
        self._dn_pre         = n_pre          # dsl_scalar
        self._dn_prim        = n_prim         # dsl_scalar
        self._dn_cont        = n_cont         # dsl_scalar
        self._diskip_prim    = iskip_prim     # dsl_scalar

    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            self.serial_out(p)

    def serial_out(self,p):
        """
        Output serial variant of contract function using printer object p.
        """
        f = self._f_dsl

        # Function arguments
        work_array     = self._work_array
        diwrk0         = self._work_array.array_index()
        diwrk1         = self._work_array.special_array_index_dict()['contract']
        contract_array = self._contract_array
        dn_pre         = self._dn_pre
        dn_prim        = self._dn_prim
        dn_cont        = self._dn_cont
        diskip_prim    = self._diskip_prim

        ### Naive implementation for debugging ###
        # Local variables
        diprim = dsl_scalar('iprim',vartype='int')
        dicont = dsl_scalar('icont',vartype='int')
        dipre  = dsl_scalar('ipre',vartype='int')
        local_variable_list = [ diprim, dicont, dipre ]

        # Start function definition
        p.funcdef( f, self.const_dict() )

        p.out("// Naive implementation for testing")

        # Local variable declarations
        for dvar in local_variable_list:
            p.declaration( dvar )

        # Loop structure top part
        nloops = 0
        p.forloop( dicont.assign(0), dsl_binop( op_lt, dicont, dn_cont ), dicont.assign(dicont+1) )
        nloops += 1
        p.forloop( dipre.assign(0), dsl_binop( op_lt, dipre, dn_pre ), dipre.assign(dipre+1) )
        nloops += 1
        # Zero the work_array element before acculumating contractions
        p.out( work_array[ diwrk0 + dipre + dicont * dn_pre ].assign(0.0) )
        p.forloop( diprim.assign(0), dsl_binop( op_lt, diprim, dn_prim ), diprim.assign(diprim+1) )
        nloops += 1

        # Accumulate contractions in work_array element
        p.out( work_array[ diwrk0 + dipre + dicont * dn_pre ].assign(\
                dsl_binop( op_add, work_array[ diwrk0 + dipre + dicont * dn_pre ], \
                dsl_binop( op_mul, work_array[ diwrk1 + dipre + diprim * diskip_prim ], \
                contract_array[ diprim + dicont * dn_prim ] ) ) ) )

        # Loop structure bottom part
        for iloop in range(nloops):
            p.endforloop()

        p.returnstmt()
        # Close function definition
        p.endfuncdef()
        ###

    def vector_out(self,p):
        f = self._f_dsl
        p.funcdef( f )
        p.out("NOTHING TO SEE HERE",endl='')
        p.endfuncdef( f )

class cblas_contract_source(function_wrapper.callable_source):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate executable source code for a contraction.
    """
    def __init__(self,f_dsl, \
                 work_array, contract_array, \
                 n_pre, n_prim, n_cont, iskip_prim, \
                 gen, const_dict = {} ):
        function_wrapper.callable_source.__init__(self, f_dsl, gen, const_dict = const_dict)
        self._work_array     = work_array     # general_array object
        self._contract_array = contract_array # general_array object
        self._dn_pre         = n_pre          # dsl_scalar
        self._dn_prim        = n_prim         # dsl_scalar
        self._dn_cont        = n_cont         # dsl_scalar
        self._diskip_prim    = iskip_prim     # dsl_scalar

    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            self.serial_out(p)

    def serial_out(self,p):
        """
        Output serial variant of contract function using printer object p.
        """
       
        f = self._f_dsl

        # Function arguments
        work_array     = self._work_array
        diwrk0         = self._work_array.array_index()
        diwrk1         = self._work_array.special_array_index_dict()['contract']
        contract_array = self._contract_array
        dn_pre         = self._dn_pre
        dn_prim        = self._dn_prim
        dn_cont        = self._dn_cont
        diskip_prim    = self._diskip_prim

        # Local variables
        diprim = dsl_scalar('iprim',vartype='int')
        dicont = dsl_scalar('icont',vartype='int')
        dipre  = dsl_scalar('ipre',vartype='int')
        local_variable_list = [ diprim, dicont, dipre ]

        # Preset parameters for determining when to call cblas (set in generator_options object)
        cblas_dgemm_min_lengths = self.generator().generator_options().cblas_dgemm_min_lengths()
        for i in cblas_dgemm_min_lengths: assert isinstance( i, int ),\
                "elements of cblas_dgemm_min_lengths must be a Python integer"
        cblas_min_n_pre  = cblas_dgemm_min_lengths[0] # m
        cblas_min_n_cont = cblas_dgemm_min_lengths[1] # n
        cblas_min_n_prim = cblas_dgemm_min_lengths[2] # k

        # Function for outputting naive loop implementation of matrix multiply
        def naive_loops_out(n_cont_loop,n_pre_loop,n_prim_loop):
            nloops = 0
            work_array_expr0 = diwrk0
            work_array_expr1 = diwrk1
            contract_array_expr = dsl_zero() 
            # Loop structure top part
            if n_cont_loop == True:
                work_array_expr0 = work_array_expr0 + dicont * dn_pre
                contract_array_expr = dicont * dn_prim
                p.forloop( dicont.assign(0), dsl_binop( op_lt, dicont, dn_cont ), dicont.assign(dicont+1) )
                nloops += 1
            if n_pre_loop == True:
                work_array_expr0 = work_array_expr0 + dipre
                work_array_expr1 = work_array_expr1 + dipre 
                p.forloop( dipre.assign(0), dsl_binop( op_lt, dipre, dn_pre ), dipre.assign(dipre+1) )
                nloops += 1
            if n_prim_loop == True:
                # Zero the work_array element before acculumating contractions
                work_array_expr1 = work_array_expr1 + diprim * diskip_prim
                if isinstance(contract_array_expr,dsl_zero):
                    contract_array_expr = diprim
                else:
                    contract_array_expr = contract_array_expr + diprim
                p.out( work_array[ work_array_expr0 ].assign(0.0) )
                p.forloop( diprim.assign(0), dsl_binop( op_lt, diprim, dn_prim ), diprim.assign(diprim+1) )
                nloops += 1
                # Accumulate in work_array element
                mm_expr = dsl_binop( op_add, work_array[ work_array_expr0 ], \
                            dsl_binop( op_mul, work_array[ work_array_expr1 ], \
                            contract_array[ contract_array_expr ] ) )
            else:
                # Single assignment to work_array element
                mm_expr = dsl_binop( op_mul, work_array[ work_array_expr1 ], \
                            contract_array[ contract_array_expr ] )

            # Output work_array assignment
            p.out( work_array[ work_array_expr0 ].assign( mm_expr ) )

            # Loop structure bottom part
            for iloop in range(nloops):
                p.endforloop()

        ### CBLAS implementation  ###
        # cblas_dgemm dsl_function object
        #void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
        #                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
        #                 const int K, const double alpha, const double  *A,
        #                 const int lda, const double  *B, const int ldb,
        #                 const double beta, double  *C, const int ldc)
        dorder  = dsl_value(name='Order', value='CblasColMajor', vartype='enum')
        dtransA = dsl_value(name='TransA', value='CblasNoTrans', vartype='enum')
        dtransB = dsl_value(name='TransB', value='CblasNoTrans', vartype='enum')
        dm      = dn_pre
        dn      = dn_cont
        dk      = dn_prim
        dalpha  = dsl_value(name='alpha',value = 1.0, vartype = 'double')
        dA      = dsl_unop(op_ref, work_array.pointer()[ diwrk1 ] )
        dlda    = diskip_prim
        dB      = contract_array.pointer()
        dldb    = dn_prim
        dbeta   = dsl_value(name='alpha',value = 0.0, vartype = 'double')
        dC      = dsl_unop(op_ref, work_array.pointer()[ diwrk0 ] )
        dldc    = dn_pre
        arg_list = [ dorder, dtransA, dtransB, dm, dn, dk, dalpha, dA, dlda, dB, dldb, dbeta, dC, dldc ]
        f_dsl = dsl_function(name='cblas_dgemm', args=arg_list, vartype='void')
        # Create callable objects
        cblas_dgemm_call = function_wrapper.callable_function_call(f_dsl,self.generator() )
        cblas_dgemm_prototype = None # unneeded
        cblas_dgemm_source    = None # unneeded
        cblas_dgemm_function_wrapper_obj = function_wrapper( f_dsl, cblas_dgemm_call, \
                                                             cblas_dgemm_prototype, cblas_dgemm_source,\
                                                             self.generator() )

        # Start function definition
        p.funcdef( f, const_dict )

        # Local variable declarations
        for dvar in local_variable_list:
            p.declaration( dvar )

        p.out("// Requires CBLAS to be installed on the system and correctly linked during compilation")

        # Output special case for very small matrices
        p.ifblock( dsl_binop( op_eq, dn_cont, 1 ) ) 
        p.ifblock( dsl_binop( op_eq, dn_pre, 1 ) )
        p.ifblock( dsl_binop( op_eq, dn_prim, 1 ) )
        # n_cont == n_pre == n_prim == 1
        naive_loops_out( False, False, False )
        p.returnstmt()
        p.elifblock( dsl_binop( op_lt, dn_prim, cblas_min_n_prim ) )
        # n_cont == n_pre == 1, n_prim < cblas_min_k
        naive_loops_out( False, False, True )
        p.returnstmt()
        p.endifblock()
        p.elifblock( dsl_binop( op_lt, dn_pre, cblas_min_n_pre ) )
        p.ifblock( dsl_binop( op_eq, dn_prim, 1 ) )
        # n_cont == 1, n_pre < cblas_min_m, n_prim == 1 
        naive_loops_out( False, True, False )
        p.returnstmt()
        p.elifblock( dsl_binop( op_lt, dn_prim, cblas_min_n_prim ) )
        # n_cont == 1, n_pre < cblas_min_m, n_prim < cblas_min_k
        naive_loops_out( False, True, True )
        p.returnstmt()
        p.endifblock()
        p.endifblock()
        p.elifblock( dsl_binop( op_lt, dn_cont, cblas_min_n_cont ) )
        p.ifblock( dsl_binop( op_eq, dn_pre, 1 ) )
        p.ifblock( dsl_binop( op_eq, dn_prim, 1 ) )
        # n_cont < c_blas_min_n, n_pre == n_prim == 1
        naive_loops_out( True, False, False )
        p.returnstmt()
        p.elifblock( dsl_binop( op_lt, dn_prim, cblas_min_n_prim ) )
        # n_cont < c_blas_min_n, n_pre == 1, n_prim < cblas_min_k
        naive_loops_out( True, False, True )
        p.returnstmt()
        p.endifblock()
        p.elifblock( dsl_binop( op_lt, dn_pre, cblas_min_n_pre ) )
        p.ifblock( dsl_binop( op_eq, dn_prim, 1 ) )
        # n_cont < cblas_min_n, n_pre < cblas_min_m, n_prim == 1 
        naive_loops_out( True, True, False )
        p.returnstmt()
        p.elifblock( dsl_binop( op_lt, dn_prim, cblas_min_n_prim ) )
        # n_cont < cblas_min_n, n_pre < cblas_min_m, n_prim < cblas_min_k
        naive_loops_out( True, True, True )
        p.returnstmt()
        p.endifblock()
        p.endifblock()
        p.endifblock()

        # Output call to cblas_dgemm
        cblas_dgemm_function_wrapper_obj.call_out(p)
        
        p.returnstmt()
        # Close function definition
        p.endfuncdef()
        ###

    def vector_out(self,p):
        f = self._f_dsl
        p.funcdef( f )
        p.out("NOTHING TO SEE HERE",endl='')
        p.endfuncdef( f )

