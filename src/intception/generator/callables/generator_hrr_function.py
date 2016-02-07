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
Provides the hrr_call, hrr_prototype and hrr_source 
classes which are derived from callable_* classes in intception.generator.wrappers.

The hrr_source class contains the code necessary to generate a self-contained
function for performing a single HRR step, shifting one unit of angular momentum
from one Cartesian index to another.
"""

from intception.dsl import *
from intception.generator.wrappers import function_wrapper

class hrr_call(function_wrapper.callable_function_call):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function calls specific to a HRR.
    """
    pass

class hrr_prototype(function_wrapper.callable_prototype):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate function prototypes specific to a HRR.
    """
    pass

class primitive_hrr_source(function_wrapper.callable_source):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate executable source code specific to a HRR.
    """
    def __init__(self,f_dsl, Rfromto, lmax_move_from, \
                            l_move_from, loop_l_move_to,\
                            iskip_move_from0, iskip_move_to0, \
                            iskip_move_from1, iskip_move_to1, \
                            work_array, gen, const_dict = {} ):
        function_wrapper.callable_source.__init__(self, f_dsl, gen, const_dict = const_dict)
        self._dRfromto = Rfromto
        self._dlmax_move_from = lmax_move_from
        self._dl_move_from = l_move_from
        self._dloop_l_move_to = loop_l_move_to
        self._diskip_move_from0 = iskip_move_from0
        self._diskip_move_to0 = iskip_move_to0
        self._diskip_move_from1 = iskip_move_from1
        self._diskip_move_to1 = iskip_move_to1
        self._work_array = work_array

    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            self.serial_out(p)

    def serial_out(self,p):
        f = self._f_dsl
        direction = p.converter.direction()
        work_array       = self._work_array
        dRfromto         = self._dRfromto
        dlmax_move_from   =self._dlmax_move_from
        dl_move_from     = self._dl_move_from
        dloop_l_move_to  = self._dloop_l_move_to
        diskip_move_from0= self._diskip_move_from0
        diskip_move_to0  = self._diskip_move_to0
        diskip_move_from1= self._diskip_move_from1
        diskip_move_to1  = self._diskip_move_to1
        diwrk0            = self._work_array.array_index()
        diwrk1            = self._work_array.special_array_index_dict()[ 'hrr' ]

        # Create a local variable list of dsl_* objects
        dmx_from = dsl_scalar(name = 'mx_from', vartype = 'int' )
        dmy_from = dsl_scalar(name = 'my_from', vartype = 'int' )
        dmz_from = dsl_scalar(name = 'mz_from', vartype = 'int' )
        dmx_to   = dsl_scalar(name = 'mx_to', vartype = 'int' )
        dmy_to   = dsl_scalar(name = 'my_to', vartype = 'int' )
        dmz_to   = dsl_scalar(name = 'mz_to', vartype = 'int' )
        di_from  = dsl_scalar(name = 'i_from', vartype = 'int' )
        di_to    = dsl_scalar(name = 'i_to', vartype = 'int' )
        dip1_from  = dsl_scalar(name = 'ip1_from', vartype = 'int' )
        dim1_to    = dsl_scalar(name = 'im1_to', vartype = 'int' )
        dl_mf    = dsl_scalar( name = 'l_mf', vartype = 'int' )
        dioff_mf = dsl_scalar( name = 'ioff_mf', vartype = 'int' )
        local_variable_list = [ dmx_from, dmy_from, dmz_from, \
                                dmx_to,   dmy_to,   dmz_to, \
                                di_from,     di_to, \
                                dip1_from,   dim1_to, \
                                dl_mf, dioff_mf ]
        p.funcdef( f, self.const_dict() )
        # Local variable declarations
        for dvar in local_variable_list:
            p.declaration( dvar )

        # Define internal function for loops over move_from_index
        def move_from_loops(direction):
            """
            Outputs loops over move_from_index.

            Nested function has read access to all variables of enclosing function.
            """
            assert isinstance(direction,dsl_direction),\
                    "direction must be dsl_direction"
            # Loop over angular momentum level in move_from_index
            p.forloop( dl_mf.assign(dl_move_from),\
                        dsl_binop( op_le, dl_mf, dlmax_move_from ),\
                        dl_mf.assign( dl_mf + 1 ) ) # l_mf
            # Loop over Cartesian components of move_from_index
            p.forloop( dmx_from.assign( dl_mf ), \
                        dsl_binop(op_ge,dmx_from,0), \
                        dmx_from.assign( dmx_from - 1 ) ) # mx_from
            p.forloop( dmy_from.assign( dl_mf - dmx_from ),\
                        dsl_binop(op_ge,dmy_from,0),\
                        dmy_from.assign( dmy_from - 1 ) ) # my_from
            p.out( dmz_from.assign( dl_mf - dmx_from - dmy_from ) )
            if direction.index() == 0:
                p.out( dip1_from.assign( \
                 dmz_from + ( dmy_from + dmz_from )*( dmy_from + dmz_from + 1)/2 +\
                 ( dmx_from + dmy_from + dmz_from + 1 )*( dmx_from + dmy_from + dmz_from + 2)*\
                 ( dmx_from + dmy_from + dmz_from + 3 )/6 - dioff_mf ) )
            elif direction.index() == 1:
                p.out( dip1_from.assign( \
                 dmz_from + ( dmy_from + dmz_from + 1)*( dmy_from + dmz_from + 2)/2 +\
                 ( dmx_from + dmy_from + dmz_from + 1)*( dmx_from + dmy_from + dmz_from + 2)*\
                 ( dmx_from + dmy_from + dmz_from + 3 )/6 - dioff_mf ) )
            elif direction.index() == 2:
                p.out( dip1_from.assign( \
                 dmz_from + 1 + ( dmy_from + dmz_from + 1 )*( dmy_from + dmz_from + 2)/2 +\
                 ( dmx_from + dmy_from + dmz_from + 1)*( dmx_from + dmy_from + dmz_from + 2)*\
                 ( dmx_from + dmy_from + dmz_from + 3 )/6 - dioff_mf ) )
            # Output HRR expr
            # Preconvert work_array index expressions, since these will not be automatically
            # processed by cprinter (not part of dsl_* expr directly accessible to cprinter)
            #src_i0_i0   = p.converter.expr_full_convert( \
            #        diwrk0 + di_to * diskip_move_to0 + di_from * diskip_move_from0 )
            #src_im1_ip1 = p.converter.expr_full_convert( \
            #        diwrk1 + dim1_to * diskip_move_to1 + dip1_from * diskip_move_from1) 
            #src_im1_i1  = p.converter.expr_full_convert( \
            #        diwrk1 + dim1_to * diskip_move_to1 + di_from * diskip_move_from1 )
            #hrr_expr =  work_array[ src_im1_ip1 ] + dRfromto * work_array[ src_im1_i1 ]
            #p.out( work_array[ src_i0_i0 ].assign( hrr_expr ) )
            # Automatic conversion of index expressions has been implemented, so pre-conversion
            # unnecessary
            i0_i0   = diwrk0 + di_to * diskip_move_to0 + di_from * diskip_move_from0
            im1_ip1 = diwrk1 + dim1_to * diskip_move_to1 + dip1_from * diskip_move_from1
            im1_i1  = diwrk1 + dim1_to * diskip_move_to1 + di_from * diskip_move_from1
            hrr_expr =  work_array[ im1_ip1 ] + dRfromto * work_array[ im1_i1 ]
            p.out( work_array[ i0_i0 ].assign( hrr_expr ) )
            p.out( di_from.assign( di_from + 1 ) )
            # Close loops over move_from_index Cartesian components
            p.endforloop() # my_from 
            p.endforloop() # mx_from
            # Close loop over move_from_index angular momentum level
            p.endforloop() # l_mf

        p.out( dioff_mf.assign( \
                dl_move_from * ( dl_move_from + 1 ) * (dl_move_from + 2) / 6 ) )
        # Loop over Cartesian components of move_to_index where mx != 0 (d = 0 )
        p.out( di_to.assign( 0 ) )
        direction.set('x')
        p.forloop( dmx_to.assign( dloop_l_move_to ),\
                    dsl_binop(op_ge,dmx_to,1),\
                    dmx_to.assign( dmx_to - 1 ) ) # mx_to
        p.forloop( dmy_to.assign( dloop_l_move_to - dmx_to ),\
                    dsl_binop(op_ge,dmy_to,0),\
                    dmy_to.assign( dmy_to - 1 ) ) # my_to
        p.out( dmz_to.assign( dloop_l_move_to - dmx_to - dmy_to ) )
        p.out( dim1_to.assign( dmz_to + ( dmy_to + dmz_to )*(dmy_to + dmz_to + 1)/2 ) )
        p.out( di_from.assign( 0 ) )
        #### move_from_index loops ###
        move_from_loops(direction)
        p.out( di_to.assign( di_to + 1 ) )
        # Close loops over move_to_index Cartesian components
        p.endforloop() # my_to
        p.endforloop() # mx_to
        
        # Loop over Cartesian components of move_to_index where mx == 0, my !=0
        # (d = 1)
        direction.set('y')
        p.out( dmx_to.assign(0) )
        p.forloop( dmy_to.assign( dloop_l_move_to - dmx_to ), \
                    dsl_binop(op_ge,dmy_to,1), \
                    dmy_to.assign( dmy_to - 1 ) ) # my_to
        p.out( dmz_to.assign( dloop_l_move_to - dmy_to ) )
        p.out( dim1_to.assign( dmz_to + ( dmy_to + dmz_to - 1 )*(dmy_to + dmz_to )/2 ) )
        p.out( di_from.assign( 0 ) )
        #### move_from_index loops ###
        move_from_loops(direction)
        p.out( di_to.assign( di_to + 1 ) )
        # Close loops over move_to_index Cartesian components
        p.endforloop() # my_to

        # Loop over Cartesian components of move_to_index where mx == 0, my == 0
        # (d = 2)
        direction.set('z')
        p.out( dmy_to.assign(0) )
        p.out( dmz_to.assign( dloop_l_move_to ) )
        p.out( dim1_to.assign( dmz_to - 1 + ( dmz_to - 1 )*( dmz_to )/2 ) )
        p.out( di_from.assign( 0 ) )
        move_from_loops(direction)
        # Close function definition
        p.endfuncdef()

    def vector_out(self,p):
        f = self._f_dsl
        p.funcdef( f )
        p.out("NOTHING TO SEE HERE",endl='')
        p.endfuncdef( f )

class contracted_hrr_source(function_wrapper.callable_source):
    """
    Instances are callable objects carrying with them all the information
    necessary to generate executable source code specific to a HRR.
    """
    def __init__(self,f_dsl, Rfromto, lmax_move_from, \
                            l_move_from, loop_l_move_to,\
                            iskip_move_from0, iskip_move_to0, \
                            iskip_move_from1, iskip_move_to1, \
                            ncont_move_from,  ncont_move_to, \
                            iskip_cont_move_from0, iskip_cont_move_to0, \
                            iskip_cont_move_from1, iskip_cont_move_to1, \
                            work_array, gen, const_dict ):
        function_wrapper.callable_source.__init__(self, f_dsl, gen, const_dict)
        self._dRfromto = Rfromto
        self._dlmax_move_from = lmax_move_from
        self._dl_move_from = l_move_from
        self._dloop_l_move_to = loop_l_move_to
        self._diskip_move_from0 = iskip_move_from0
        self._diskip_move_to0 = iskip_move_to0
        self._diskip_move_from1 = iskip_move_from1
        self._diskip_move_to1 = iskip_move_to1
        self._dncont_move_from = ncont_move_from
        self._dncont_move_to   = ncont_move_to
        self._diskip_cont_move_from0 = iskip_cont_move_from0
        self._diskip_cont_move_to0 = iskip_cont_move_to0
        self._diskip_cont_move_from1 = iskip_cont_move_from1
        self._diskip_cont_move_to1 = iskip_cont_move_to1
        self._work_array = work_array

    def __call__(self,p):
        #assert isinstance(p,printer), "p must be a printer object"
        if self.generator().generator_options().parallelism().vectorized() == True:
            raise Exception("vector output not implemented yet!")
            self.vector_out(p)
        elif self.generator().generator_options().parallelism().vectorized() == False:
            self.serial_out(p)

    def serial_out(self,p):
        f = self._f_dsl
        direction = p.converter.direction()
        work_array       = self._work_array
        dRfromto         = self._dRfromto
        dlmax_move_from   =self._dlmax_move_from
        dl_move_from     = self._dl_move_from
        dloop_l_move_to  = self._dloop_l_move_to
        diskip_move_from0= self._diskip_move_from0
        diskip_move_to0  = self._diskip_move_to0
        diskip_move_from1= self._diskip_move_from1
        diskip_move_to1  = self._diskip_move_to1
        dncont_move_from = self._dncont_move_from
        dncont_move_to   = self._dncont_move_to
        diskip_cont_move_from0= self._diskip_cont_move_from0
        diskip_cont_move_to0  = self._diskip_cont_move_to0
        diskip_cont_move_from1= self._diskip_cont_move_from1
        diskip_cont_move_to1  = self._diskip_cont_move_to1
        diwrk0            = self._work_array.array_index()
        diwrk1            = self._work_array.special_array_index_dict()[ 'hrr' ]

        # Create a local variable list of dsl_* objects
        dmx_from = dsl_scalar(name = 'mx_from', vartype = 'int' )
        dmy_from = dsl_scalar(name = 'my_from', vartype = 'int' )
        dmz_from = dsl_scalar(name = 'mz_from', vartype = 'int' )
        dmx_to   = dsl_scalar(name = 'mx_to', vartype = 'int' )
        dmy_to   = dsl_scalar(name = 'my_to', vartype = 'int' )
        dmz_to   = dsl_scalar(name = 'mz_to', vartype = 'int' )
        di_from  = dsl_scalar(name = 'i_from', vartype = 'int' )
        di_to    = dsl_scalar(name = 'i_to', vartype = 'int' )
        dicont_from  = dsl_scalar(name = 'icont_from', vartype = 'int' )
        dicont_to    = dsl_scalar(name = 'icont_to', vartype = 'int' )
        dip1_from  = dsl_scalar(name = 'ip1_from', vartype = 'int' )
        dim1_to    = dsl_scalar(name = 'im1_to', vartype = 'int' )
        dl_mf    = dsl_scalar( name = 'l_mf', vartype = 'int' )
        dioff_mf = dsl_scalar( name = 'ioff_mf', vartype = 'int' )
        local_variable_list = [ dmx_from, dmy_from, dmz_from, \
                                dmx_to,   dmy_to,   dmz_to, \
                                di_from,     di_to, \
                                dicont_from,     dicont_to, \
                                dip1_from,   dim1_to, \
                                dl_mf, dioff_mf ]
        p.funcdef( f, self.const_dict() )
        # Local variable declarations
        for dvar in local_variable_list:
            p.declaration( dvar )

        # Define internal function for loops over move_from_index
        def move_from_loops(direction):
            """
            Outputs loops over move_from_index.

            Nested function has read access to all variables of enclosing function.
            """
            assert isinstance(direction,dsl_direction),\
                    "direction must be dsl_direction"
            # Loop over contraction in move_from_index
            p.forloop( dicont_from.assign(0),\
                        dsl_binop( op_lt, dicont_from, dncont_move_from ),\
                        dicont_from.assign( dicont_from + 1 ) )
            # Set counter to zero
            p.out( di_from.assign( 0 ) )
            # Loop over angular momentum level in move_from_index
            p.forloop( dl_mf.assign(dl_move_from),\
                        dsl_binop( op_le, dl_mf, dlmax_move_from ),\
                        dl_mf.assign( dl_mf + 1 ) ) # l_mf
            # Loop over Cartesian components of move_from_index
            p.forloop( dmx_from.assign( dl_mf ), \
                        dsl_binop(op_ge,dmx_from,0), \
                        dmx_from.assign( dmx_from - 1 ) ) # mx_from
            p.forloop( dmy_from.assign( dl_mf - dmx_from ),\
                        dsl_binop(op_ge,dmy_from,0),\
                        dmy_from.assign( dmy_from - 1 ) ) # my_from
            p.out( dmz_from.assign( dl_mf - dmx_from - dmy_from ) )
            if direction.index() == 0:
                p.out( dip1_from.assign( \
                 dmz_from + ( dmy_from + dmz_from )*( dmy_from + dmz_from + 1)/2 +\
                 ( dmx_from + dmy_from + dmz_from + 1 )*( dmx_from + dmy_from + dmz_from + 2)*\
                 ( dmx_from + dmy_from + dmz_from + 3 )/6 - dioff_mf ) )
            elif direction.index() == 1:
                p.out( dip1_from.assign( \
                 dmz_from + ( dmy_from + dmz_from + 1)*( dmy_from + dmz_from + 2)/2 +\
                 ( dmx_from + dmy_from + dmz_from + 1)*( dmx_from + dmy_from + dmz_from + 2)*\
                 ( dmx_from + dmy_from + dmz_from + 3 )/6 - dioff_mf ) )
            elif direction.index() == 2:
                p.out( dip1_from.assign( \
                 dmz_from + 1 + ( dmy_from + dmz_from + 1 )*( dmy_from + dmz_from + 2)/2 +\
                 ( dmx_from + dmy_from + dmz_from + 1)*( dmx_from + dmy_from + dmz_from + 2)*\
                 ( dmx_from + dmy_from + dmz_from + 3 )/6 - dioff_mf ) )
            # Output HRR expr
            # Preconvert work_array index expressions, since these will not be automatically
            # processed by cprinter (not part of dsl_* expr directly accessible to cprinter)
            #src_i0_i0   = p.converter.expr_full_convert( \
            #        diwrk0 + di_to * diskip_move_to0 + di_from * diskip_move_from0 )
            #src_im1_ip1 = p.converter.expr_full_convert( \
            #        diwrk1 + dim1_to * diskip_move_to1 + dip1_from * diskip_move_from1) 
            #src_im1_i1  = p.converter.expr_full_convert( \
            #        diwrk1 + dim1_to * diskip_move_to1 + di_from * diskip_move_from1 )
            #hrr_expr =  work_array[ src_im1_ip1 ] + dRfromto * work_array[ src_im1_i1 ]
            #p.out( work_array[ src_i0_i0 ].assign( hrr_expr ) )
            # Automatic conversion of index expressions has been implemented, so pre-conversion
            # unnecessary
            i0_i0   = diwrk0 + di_to * diskip_move_to0 + di_from * diskip_move_from0 + \
                        dicont_from * diskip_cont_move_from0 + dicont_to * diskip_cont_move_to0
            im1_ip1 = diwrk1 + dim1_to * diskip_move_to1 + dip1_from * diskip_move_from1 +\
                        dicont_from * diskip_cont_move_from1 + dicont_to * diskip_cont_move_to1
            im1_i1  = diwrk1 + dim1_to * diskip_move_to1 + di_from * diskip_move_from1 +\
                        dicont_from * diskip_cont_move_from1 + dicont_to * diskip_cont_move_to1
            hrr_expr =  work_array[ im1_ip1 ] + dRfromto * work_array[ im1_i1 ]
            p.out( work_array[ i0_i0 ].assign( hrr_expr ) )
            p.out( di_from.assign( di_from + 1 ) )
            # Close loops over move_from_index Cartesian components
            p.endforloop() # my_from 
            p.endforloop() # mx_from
            # Close loop over move_from_index angular momentum level
            p.endforloop() # l_mf
            # Close loop over contraction
            p.endforloop() # icont_from

        p.out( dioff_mf.assign( \
                dl_move_from * ( dl_move_from + 1 ) * (dl_move_from + 2) / 6 ) )
        # Loop over contractions of move_to_index 
        p.forloop( dicont_to.assign(0),\
                    dsl_binop( op_lt, dicont_to, dncont_move_to ),\
                    dicont_to.assign( dicont_to + 1 ) )
        # Loop over Cartesian components of move_to_index where mx != 0 (d = 0 )
        p.out( di_to.assign( 0 ) )
        direction.set('x')
        p.forloop( dmx_to.assign( dloop_l_move_to ),\
                    dsl_binop(op_ge,dmx_to,1),\
                    dmx_to.assign( dmx_to - 1 ) ) # mx_to
        p.forloop( dmy_to.assign( dloop_l_move_to - dmx_to ),\
                    dsl_binop(op_ge,dmy_to,0),\
                    dmy_to.assign( dmy_to - 1 ) ) # my_to
        p.out( dmz_to.assign( dloop_l_move_to - dmx_to - dmy_to ) )
        p.out( dim1_to.assign( dmz_to + ( dmy_to + dmz_to )*(dmy_to + dmz_to + 1)/2 ) )
        #### move_from_index loops ###
        move_from_loops(direction)
        p.out( di_to.assign( di_to + 1 ) )
        # Close loops over move_to_index Cartesian components
        p.endforloop() # my_to
        p.endforloop() # mx_to
        
        # Loop over Cartesian components of move_to_index where mx == 0, my !=0
        # (d = 1)
        direction.set('y')
        p.out( dmx_to.assign(0) )
        p.forloop( dmy_to.assign( dloop_l_move_to - dmx_to ), \
                    dsl_binop(op_ge,dmy_to,1), \
                    dmy_to.assign( dmy_to - 1 ) ) # my_to
        p.out( dmz_to.assign( dloop_l_move_to - dmy_to ) )
        p.out( dim1_to.assign( dmz_to + ( dmy_to + dmz_to - 1 )*(dmy_to + dmz_to )/2 ) )
        #### move_from_index loops ###
        move_from_loops(direction)
        p.out( di_to.assign( di_to + 1 ) )
        # Close loops over move_to_index Cartesian components
        p.endforloop() # my_to

        # Loop over Cartesian components of move_to_index where mx == 0, my == 0
        # (d = 2)
        direction.set('z')
        p.out( dmy_to.assign(0) )
        p.out( dmz_to.assign( dloop_l_move_to ) )
        p.out( dim1_to.assign( dmz_to - 1 + ( dmz_to - 1 )*( dmz_to )/2 ) )
        move_from_loops(direction)

        # Close loop over contractions of move_to_index 
        p.endforloop() # icont_to

        # Close function definition
        p.endfuncdef()

    def vector_out(self,p):
        f = self._f_dsl
        p.funcdef( f )
        p.out("NOTHING TO SEE HERE",endl='')
        p.endfuncdef( f )

