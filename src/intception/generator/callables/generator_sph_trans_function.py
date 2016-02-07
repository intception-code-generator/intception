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
Provides the sph_trans_call, sph_trans_prototype and sph_trans_source 
classes which are derived from callable_* classes in intception.generator.wrappers.

The sph_trans_source class contains the code necessary to generate a self-contained
function for performing a single spherical transformation of a primitive Cartesian index.
"""

from intception.dsl import *
from intception.generator.wrappers import function_wrapper

class sph_trans_call(function_wrapper.callable_function_call):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function calls for a spherical transformation.
    """
    pass

class sph_trans_prototype(function_wrapper.callable_prototype):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function prototypes for a spherical transformation.
    """
    pass

class noblas_sph_trans_source(function_wrapper.callable_source):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate executable source code for a spherical transformation.
    """
    def __init__(self,f_dsl, \
                 work_array, sph_trans_array, \
                 n_pre, n_cart, n_sph, iskip_cart, \
                 gen, const_dict = {} ):
        function_wrapper.callable_source.__init__(self, f_dsl, gen, const_dict = const_dict)
        self._work_array     = work_array     # general_array object
        self._sph_trans_array = sph_trans_array # general_array object
        self._dn_pre         = n_pre          # dsl_scalar
        self._dn_cart        = n_cart         # dsl_scalar
        self._dn_sph         = n_sph          # dsl_scalar
        self._diskip_cart    = iskip_cart     # dsl_scalar

    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            self.serial_out(p)

    def serial_out(self,p):
        """
        Output serial variant of sph_trans function using printer object p.
        """
        f = self._f_dsl

        # Function arguments
        work_array     = self._work_array
        diwrk0         = self._work_array.array_index()
        diwrk1         = self._work_array.special_array_index_dict()['sph_trans']
        sph_trans_array = self._sph_trans_array
        dn_pre         = self._dn_pre
        dn_cart        = self._dn_cart
        dn_sph         = self._dn_sph
        diskip_cart    = self._diskip_cart

        ### Naive implementation for debugging ###
        # Local variables
        dicart = dsl_scalar('icart',vartype='int')
        disph = dsl_scalar('isph',vartype='int')
        dipre  = dsl_scalar('ipre',vartype='int')
        local_variable_list = [ dicart, disph, dipre ]

        # Start function definition
        p.funcdef( f, self.const_dict() )

        p.out("// Naive implementation for testing")

        # Local variable declarations
        for dvar in local_variable_list:
            p.declaration( dvar )

        # Loop structure top part
        nloops = 0
        p.forloop( disph.assign(0), dsl_binop( op_lt, disph, dn_sph ), disph.assign(disph+1) )
        nloops += 1
        p.forloop( dipre.assign(0), dsl_binop( op_lt, dipre, dn_pre ), dipre.assign(dipre+1) )
        nloops += 1
        # Zero the work_array element before acculumating spherical components
        p.out( work_array[ diwrk0 + dipre + disph * dn_pre ].assign(0.0) )
        p.forloop( dicart.assign(0), dsl_binop( op_lt, dicart, dn_cart ), dicart.assign(dicart+1) )
        nloops += 1

        # Accumulate spherical components in work_array element
        p.out( work_array[ diwrk0 + dipre + disph * dn_pre ].assign(\
                dsl_binop( op_add, work_array[ diwrk0 + dipre + disph * dn_pre ], \
                dsl_binop( op_mul, work_array[ diwrk1 + dipre + dicart * diskip_cart ], \
                sph_trans_array[ dicart + disph * dn_cart ] ) ) ) )

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

class cblas_sph_trans_source(function_wrapper.callable_source):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate executable source code for a spherical transformation.
    """
    def __init__(self,f_dsl, \
                 work_array, sph_trans_array, \
                 n_pre, n_cart, n_sph, iskip_cart, \
                 gen, const_dict = {} ):
        function_wrapper.callable_source.__init__(self, f_dsl, gen, const_dict = const_dict)
        self._work_array     = work_array     # general_array object
        self._sph_trans_array = sph_trans_array # general_array object
        self._dn_pre         = n_pre          # dsl_scalar
        self._dn_cart        = n_cart         # dsl_scalar
        self._dn_sph        = n_sph         # dsl_scalar
        self._diskip_cart    = iskip_cart     # dsl_scalar

    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            self.serial_out(p)

    def serial_out(self,p):
        """
        Output serial variant of sphract function using printer object p.
        """
       
        f = self._f_dsl

        # Function arguments
        work_array     = self._work_array
        diwrk0         = self._work_array.array_index()
        diwrk1         = self._work_array.special_array_index_dict()['sph_trans']
        sph_trans_array = self._sph_trans_array
        dn_pre         = self._dn_pre
        dn_cart        = self._dn_cart
        dn_sph         = self._dn_sph
        diskip_cart    = self._diskip_cart

        # Local variables
        dicart = dsl_scalar('icart',vartype='int')
        disph  = dsl_scalar('isph',vartype='int')
        dipre  = dsl_scalar('ipre',vartype='int')
        local_variable_list = [ dicart, disph, dipre ]

        # Preset parameters for determining when to call cblas (set in generator_options object)
        cblas_dgemm_min_lengths = self.generator().generator_options().cblas_dgemm_min_lengths()
        for i in cblas_dgemm_min_lengths: assert isinstance( i, int ),\
                "elements of cblas_dgemm_min_lengths must be a Python integer"
        cblas_min_n_pre  = cblas_dgemm_min_lengths[0] # m
        cblas_min_n_sph  = cblas_dgemm_min_lengths[1] # n
        cblas_min_n_cart = cblas_dgemm_min_lengths[2] # k

        # Function for outputting naive loop implementation of matrix multiply
        def naive_loops_out(n_sph_loop,n_pre_loop,n_cart_loop):
            nloops = 0
            work_array_expr0 = diwrk0
            work_array_expr1 = diwrk1
            sph_trans_array_expr = dsl_zero()
            # Loop structure top part
            if n_sph_loop == True:
                work_array_expr0 = work_array_expr0 + disph * dn_pre
                sph_trans_array_expr = disph * dn_cart
                p.forloop( disph.assign(0), dsl_binop( op_lt, disph, dn_sph ), disph.assign(disph+1) )
                nloops += 1
            if n_pre_loop == True:
                work_array_expr0 = work_array_expr0 + dipre
                work_array_expr1 = work_array_expr1 + dipre 
                p.forloop( dipre.assign(0), dsl_binop( op_lt, dipre, dn_pre ), dipre.assign(dipre+1) )
                nloops += 1
            if n_cart_loop == True:
                # Zero the work_array element before acculumating spherical components
                work_array_expr1 = work_array_expr1 + dicart * diskip_cart
                if isinstance(sph_trans_array_expr,dsl_zero):
                    sph_trans_array_expr = dicart
                else:
                    sph_trans_array_expr = sph_trans_array_expr + dicart
                p.out( work_array[ work_array_expr0 ].assign(0.0) )
                p.forloop( dicart.assign(0), dsl_binop( op_lt, dicart, dn_cart ), dicart.assign(dicart+1) )
                nloops += 1
                # Accumulate in work_array element
                mm_expr = dsl_binop( op_add, work_array[ work_array_expr0 ], \
                            dsl_binop( op_mul, work_array[ work_array_expr1 ], \
                            sph_trans_array[ sph_trans_array_expr ] ) )
            else:
                # Single assignment to work_array element
                mm_expr = dsl_binop( op_mul, work_array[ work_array_expr1 ], \
                            sph_trans_array[ sph_trans_array_expr ] )

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
        dn      = dn_sph
        dk      = dn_cart
        dalpha  = dsl_value(name='alpha',value = 1.0, vartype = 'double')
        dA      = dsl_unop(op_ref, work_array.pointer()[ diwrk1 ] )
        dlda    = diskip_cart
        dB      = sph_trans_array.pointer()
        dldb    = dn_cart
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
        p.funcdef( f, self.const_dict() )

        # Local variable declarations
        for dvar in local_variable_list:
            p.declaration( dvar )

        p.out("// Requires CBLAS to be installed on the system and correctly linked during compilation")

        # Output special case for very small matrices
        p.ifblock( dsl_binop( op_eq, dn_sph, 1 ) ) 
        p.ifblock( dsl_binop( op_eq, dn_pre, 1 ) )
        p.ifblock( dsl_binop( op_eq, dn_cart, 1 ) )
        # n_sph == n_pre == n_cart == 1
        naive_loops_out( False, False, False )
        p.returnstmt()
        p.elifblock( dsl_binop( op_lt, dn_cart, cblas_min_n_cart ) )
        # n_sph == n_pre == 1, n_cart < cblas_min_k
        naive_loops_out( False, False, True )
        p.returnstmt()
        p.endifblock()
        p.elifblock( dsl_binop( op_lt, dn_pre, cblas_min_n_pre ) )
        p.ifblock( dsl_binop( op_eq, dn_cart, 1 ) )
        # n_sph == 1, n_pre < cblas_min_m, n_cart == 1 
        naive_loops_out( False, True, False )
        p.returnstmt()
        p.elifblock( dsl_binop( op_lt, dn_cart, cblas_min_n_cart ) )
        # n_sph == 1, n_pre < cblas_min_m, n_cart < cblas_min_k
        naive_loops_out( False, True, True )
        p.returnstmt()
        p.endifblock()
        p.endifblock()
        p.elifblock( dsl_binop( op_lt, dn_sph, cblas_min_n_sph ) )
        p.ifblock( dsl_binop( op_eq, dn_pre, 1 ) )
        p.ifblock( dsl_binop( op_eq, dn_cart, 1 ) )
        # n_sph < c_blas_min_n, n_pre == n_cart == 1
        naive_loops_out( True, False, False )
        p.returnstmt()
        p.elifblock( dsl_binop( op_lt, dn_cart, cblas_min_n_cart ) )
        # n_sph < c_blas_min_n, n_pre == 1, n_cart < cblas_min_k
        naive_loops_out( True, False, True )
        p.returnstmt()
        p.endifblock()
        p.elifblock( dsl_binop( op_lt, dn_pre, cblas_min_n_pre ) )
        p.ifblock( dsl_binop( op_eq, dn_cart, 1 ) )
        # n_sph < cblas_min_n, n_pre < cblas_min_m, n_cart == 1 
        naive_loops_out( True, True, False )
        p.returnstmt()
        p.elifblock( dsl_binop( op_lt, dn_cart, cblas_min_n_cart ) )
        # n_sph < cblas_min_n, n_pre < cblas_min_m, n_cart < cblas_min_k
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

