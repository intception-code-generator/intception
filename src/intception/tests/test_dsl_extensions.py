#!/usr/bin/env python3
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
import unittest
import math
from intception.dsl import *
from intception.dsl_extensions import *

class TestObjects(unittest.TestCase):
    """Tests for individual objects in the dsl_extensions module.
    Unit tests are only written for non-trivial methods, i.e.
    not methods which just return a member variable."""
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_dsl_rr(self):
        ga = dsl_cartesian_gaussian('a')
        gb = dsl_cartesian_gaussian('b')
        name = 'test_rr'
        rrtype = 'vrr'
        changing_index = ga
        expr = ga + gb
        test_rr = dsl_rr(name,rrtype,changing_index,expr)
        # Test __init__()
        self.assertRaises(AssertionError,dsl_rr,'name','xxx',ga,ga+gb)
        # Test __str__()
        self.assertEqual( src( expr ), str( test_rr ) )
#       print( test_rr.debug_str() )

    def test_dsl_integral(self):
        # Test __init__(), parse_args(), set_allowed_types() together
        # Create necessary DSL objects to pass to dsl_integral instance
        ga = dsl_cartesian_gaussian('a')
        gb = dsl_cartesian_gaussian('b')
        gc = dsl_cartesian_gaussian('c')
        m  = dsl_integer_index('m')
        C  = gc.cen()
        dsc = 'Two-index nuclear attraction integrals over primitive Cartesian Gaussians'
        nm  = 'nuclear_2idx'
        integral1 = dsl_integral(ga,gb,m,C,name=nm,desc=dsc)
#        print( integral1.debug_str() )
        self.assertIn( ga, integral1._index_list )
        self.assertIn( gb, integral1._index_list )
        self.assertIn( m, integral1._aux_list )
        self.assertIn( C, integral1._arg_list )
        self.assertEqual( 'intception_integral_'+nm+'_ab_m_cC', integral1._name )
        self.assertEqual( dsc, integral1._desc )
        # Test __init__(), then set_allowed_types() 
        # (with implementation restrictions)
        dsc = 'Three-index overlap integral over primitive Cartesian Gaussians'
        integral2 = dsl_integral(ga,gb,gc,name='integral2',desc=dsc)
        integral2._indextypes = [ dsl_function ]
        self.assertRaises( AssertionError, integral2.set_allowed_types )
        integral2 = dsl_integral(ga,gb,gc,name='integral2',desc=dsc)
        integral2._auxtypes = [ dsl_value ]
        self.assertRaises( AssertionError, integral2.set_allowed_types )
        integral2 = dsl_integral(ga,gb,gc,name='integral2',desc=dsc)
        integral2._argtypes = [ dsl_integral ]
        self.assertRaises( AssertionError, integral2.set_allowed_types )
        # ... test with copy_from set
        self.assertRaises( AssertionError, dsl_integral, ga, gb, gc, m, copy_from = integral2 )
        self.assertRaises( AssertionError, dsl_integral, dsl_cartesian_gaussian('x'),\
                copy_from = integral2 )
        integral2 = dsl_integral( ga,gb,gc,m,C, name=nm, desc=dsc )
        # Copy and modify attributes from integral2
        integral3 = dsl_integral(ga, gb+1, gc-1, copy_from = integral2 )
        self.assertIs( integral3.index_binop( ga ).op(), op_add )
        self.assertEqual( integral3.index_binop( ga ).right(), 0 )
        self.assertIs( integral3.index_binop( gb ).op(), op_add )
        self.assertEqual( integral3.index_binop( gb ).right(), 1 )
        self.assertIs( integral3.index_binop( gc ).op(), op_sub )
        self.assertEqual( integral3.index_binop( gc ).right(), 1 )
        self.assertIs( integral3.aux_binop( m ).op(), op_add )
        self.assertEqual( integral3.aux_binop( m ).right(), 0 )
        self.assertIn( C, integral3.arg_list() )
        # Implicitly copy non-keyword arguments from integral3
        integral4 = dsl_integral( copy_from = integral3 ) 
        self.assertIs( integral4.index_binop( ga ).op(), op_add )
        self.assertEqual( integral4.index_binop( ga ).right(), 0 )
        self.assertIs( integral4.index_binop( gb ).op(), op_add )
        self.assertEqual( integral4.index_binop( gb ).right(), 1 )
        self.assertIs( integral4.index_binop( gc ).op(), op_sub )
        self.assertEqual( integral4.index_binop( gc ).right(), 1 )
        self.assertIs( integral4.aux_binop( m ).op(), op_add )
        self.assertEqual( integral4.aux_binop( m ).right(), 0 )
        self.assertIn( C, integral4.arg_list() )
        # Modify ga and gb indexes from integral4 while copying gc unchanged
        integral5 = dsl_integral( ga+1, gb-2, m+3, copy_from = integral4 ) 
        self.assertIs( integral5.index_binop( ga ).op(), op_add )
        self.assertEqual( integral5.index_binop( ga ).right(), 1 )
        self.assertIs( integral5.index_binop( gb ).op(), op_sub )
        self.assertEqual( integral5.index_binop( gb ).right(), 2 )
        self.assertIs( integral5.index_binop( gc ).op(), op_sub )
        self.assertEqual( integral5.index_binop( gc ).right(), 1 )
        self.assertIs( integral5.aux_binop( m ).op(), op_add )
        self.assertEqual( integral5.aux_binop( m ).right(), 3 )
        self.assertIn( C, integral5.arg_list() )

        # Test parse_args from class
        arg_tuples_list = [ \
                            ( ga, gb, gc, m ),\
                            ( ga + 1, gb - 1 , gc, m ),\
                            ( ga + 0, gb - 3 , gc + 3, m + 1 ),\
                            ( ga + 9, m + 2 ),\
                            ( gc, m + 2 ) \
                          ]
        index_binop_list = [ \
                            [ dsl_binop(op_add,ga,0), dsl_binop(op_add,gb,0), dsl_binop(op_add,gc,0) ],\
                            [ dsl_binop(op_add,ga,1), dsl_binop(op_sub,gb,1), dsl_binop(op_add,gc,0) ],\
                            [ dsl_binop(op_add,ga,0), dsl_binop(op_sub,gb,3), dsl_binop(op_add,gc,3) ],\
                            [ dsl_binop(op_add,ga,9) ],\
                            [ dsl_binop(op_add,gc,0) ] \
                          ]
        aux_binop_list   = [ \
                           [ dsl_binop(op_add,m,0) ],\
                           [ dsl_binop(op_add,m,0) ],\
                           [ dsl_binop(op_add,m,1) ],\
                           [ dsl_binop(op_add,m,2) ],\
                           [ dsl_binop(op_add,m,2) ],\
                          ]
        integral_list = []                         
        for a in arg_tuples_list:
            integral_list.append( dsl_integral(name='misc') )
            dsl_integral.parse_args(integral_list[-1],a)
        for i in integral_list:
            ii = integral_list.index(i)
            for local_list, obj_list, obj_dict in zip( [ index_binop_list[ii], aux_binop_list[ii] ],\
                                                       [ i._index_list,    i._aux_list ],\
                                                       [ i._index_dict,    i._aux_dict ] ):
                for binop in local_list:
                    self.assertIn( binop.left(), obj_list )
                    obj_binop = obj_dict[ binop.left().name() ]
                    self.assertTrue( isinstance( obj_binop, dsl_binop ) )
                    self.assertIs( binop.op(), obj_binop.op() )
                    self.assertEqual( binop.left().name(), obj_binop.left().name() )
                    self.assertEqual( binop.right(), obj_binop.right() )
        # Check implementation restrictions are applied correctly
        # No decrements in auxiliaries
        self.assertRaises(AssertionError,dsl_integral.parse_args,dsl_integral(name='misc'),( ga, m - 1 ) )
        # Only one auxiliary allowed
        self.assertRaises(AssertionError,dsl_integral.parse_args,dsl_integral(name='misc'),\
                ( dsl_integer_index('x'),dsl_integer_index('y'),dsl_integer_index('z') ) )

        # Test generate_unique_name()
        test_integral_list = [ \
                               dsl_integral(ga,gb,gc,name='test1'),\
                               dsl_integral(gc,m,C,name='test2'),\
                               dsl_integral(ga,m,gb,C,name='test3'),\
                               dsl_integral(C,m,dsl_scalar(name = 'x'),ga,gb,name='test4')
                             ]
        test_integral_names = [ \
                               'intception_integral_test1_abc',\
                               'intception_integral_test2_c_m_cC',\
                               'intception_integral_test3_ab_m_cC',\
                               'intception_integral_test4_ab_m_cCx'
                              ]
        for integral, name in zip( test_integral_list, test_integral_names ):
            self.assertEqual( integral.name(), name )

        # Test set_base() method
        test_integral = dsl_integral(name='test_integal')
        xp   = dsl_scalar('xp')
        RAB2 = dsl_scalar('RAB2')
        pi   = dsl_scalar('pi')
        base_expr = ( xp )**(3.0/2.0) * ( pi * xp )**(3.0/2.0) * exp( -xp*RAB2)
        test_integral.set_base( base_expr )
        self.assertIs( test_integral.base(), base_expr )
        # ... cannot contain dsl_integral() objects
        self.assertRaises(AssertionError,test_integral.set_base, dsl_integral(name='test') * xp )
        # ... base expr can only be an object of a class which inherits from dsl_base 
        self.assertRaises(AssertionError,test_integral.set_base, 'spam' )

        # Test add_*() methods...
        # ... no args for __init__() --> no add_* methods called
        test_integral = dsl_integral(name='test')
        list_of_lists = [ test_integral._index_list, test_integral._aux_list,\
                          test_integral._arg_list ]
        for l in list_of_lists:
            self.assertEqual( len( l ), 0 )

        # ... add_index()
        test_integral = dsl_integral(name='test_integral')
        index_list  = [ ga, gb, gc, dsl_cartesian_gaussian('d') ]
        index_change_list = [ 0, 1, 1, 2 ]
        index_op_list     = [op_add, op_sub, op_add, op_sub ]
        for i, c, o in zip( index_list, index_change_list, index_op_list ):
            if c == 0:
                test_integral.add_index(i)
            elif o is op_add:
                test_integral.add_index(i,c)
            else:
                test_integral.add_index(i,c,o)
            self.assertIn( i, test_integral._index_list )
            self.assertIs( i, test_integral._index_dict[ i.name() ].left() )
            self.assertEqual( c, test_integral._index_dict[ i.name() ].right() )
            self.assertIs( o, test_integral._index_dict[ i.name() ].op() )
        self.assertEqual( test_integral._index_list, index_list )
        # test assertions
        test_integral = dsl_integral(indextypes=[dsl_cartesian_gaussian],name='test_integral')
        test_integral.add_index(ga)
        # only allowed types
        self.assertRaises(AssertionError,test_integral.add_index, dsl_integer_index('x'))
        # no duplicates
        self.assertRaises(AssertionError,test_integral.add_index, ga )
        
        # ... add_auxiliary()
        integral_list = [ dsl_integral(name='test1'), dsl_integral(name='test2'), dsl_integral(name='test3') ]
        aux_list = [ m, dsl_integer_index('n'), dsl_integer_index('o') ]
        aux_change_list = [ 0, 1, 2 ]
        aux_op_list = [ op_add, op_add, op_add ]
        for integral, a, c, o in zip( integral_list, aux_list, aux_change_list, aux_op_list ):
            if c == 0:
                integral.add_auxiliary(a)
            elif o is op_add:
                integral.add_auxiliary(a,c)
            else:
                integral.add_auxiliary(a,c,o)
            self.assertIn( a, integral._aux_list )
            self.assertIs( a, integral._aux_dict[ a.name() ].left() )
            self.assertEqual( c, integral._aux_dict[ a.name() ].right() )
            self.assertIs( o, integral._aux_dict[ a.name() ].op() )
        # test assertions
        test_integral = dsl_integral(auxtypes=[dsl_integer_index],name='test_integral')
        test_integral.add_auxiliary(m)
        # only allowed types
        self.assertRaises(AssertionError,test_integral.add_auxiliary, dsl_cartesian_gaussian('x') )
        # no duplicates
        self.assertRaises(AssertionError,test_integral.add_auxiliary, m )

        # add_argument()
        test_integral = dsl_integral(name='test_integral')
        argument_list  = [ C, dsl_scalar(name = 'x'), dsl_position(name = 'D'), dsl_scalar(name = 'y',vartype='int')]
        for a in argument_list:
            test_integral.add_argument(a)
            self.assertIn( a, test_integral._arg_list )
            self.assertIs( a, test_integral._arg_dict[ a.name() ] )
        self.assertEqual( test_integral._arg_list, argument_list )
        # test assertions
        test_integral = dsl_integral(argtypes=[dsl_scalar,dsl_position],name='test_integral')
        test_integral.add_argument(C)
        # only allowed types
        self.assertRaises(AssertionError,test_integral.add_argument, dsl_integer_index('x'))
        # no duplicates
        self.assertRaises(AssertionError,test_integral.add_index, C )

        # _add_rr(), add_vrr(), add_hrr()
        xa = ga.exp()
        A  = ga.cen()
        xb = gb.exp()
        B  = gb.cen()
        xc = gc.exp()
        C  = gc.cen()
        # Simple artificial example testing add_vrr and permuting of indexes
        test_integral = dsl_integral(ga,gb,name='test_integral')
        vrr_expr = dsl_binop(op_add, ga * test_integral.int(ga-2,gb), xa * test_integral.int(ga-1,gb) )
        # ga * (a-2,b) + xa * (a-1,b)
        # g1    int1     dp1    int2
        test_integral.add_vrr('testvrr',ga, expr = vrr_expr )
        integral_vrr = test_integral.rr('testvrr',ga).expr()
        g1   = integral_vrr.left().left()
        int1 = integral_vrr.left().right()
        dp1  = integral_vrr.right().left()
        int2 = integral_vrr.right().right()
        self.assertIs( integral_vrr, vrr_expr )
        self.assertIs( g1, ga )
        self.assertIs( dp1, xa )
        self.assertTrue( isinstance( int1, dsl_integral ) )
        self.assertEqual( int1.index_binop( ga ).op(), op_sub )
        self.assertEqual( int1.index_binop( ga ).right(), 2 )
        self.assertEqual( int1.index_binop( gb ).op(), op_add )
        self.assertEqual( int1.index_binop( gb ).right(), 0 )
        self.assertEqual( int2.index_binop( ga ).op(), op_sub )
        self.assertEqual( int2.index_binop( ga ).right(), 1 )
        self.assertEqual( int2.index_binop( gb ).op(), op_add )
        self.assertEqual( int2.index_binop( gb ).right(), 0 )
        self.assertIs( test_integral.rr('testvrr',ga).expr().right().left(), xa )
        ### Permutation of indexes in VRRs is disabled ###
        ### Explicit definition of each VRR expr is required ###
        ## gb * (a,b-2) + xb * (a,b-1)
        ## g1    int1     dp1    int2
        #test_integral.add_vrr('testvrr',gb) # vrr_expr should be permuted
        #integral_vrr = test_integral.rr('testvrr',gb).expr()
        #g1   = integral_vrr.left().left()
        #int1 = integral_vrr.left().right()
        #dp1  = integral_vrr.right().left()
        #int2 = integral_vrr.right().right()
        #self.assertIsNot( integral_vrr, vrr_expr )
        #self.assertIs( g1, gb )
        #self.assertIs( dp1, xb )
        #self.assertTrue( isinstance( int1, dsl_integral ) )
        #self.assertEqual( int1.index_binop( ga ).op(), op_add )
        #self.assertEqual( int1.index_binop( ga ).right(), 0 )
        #self.assertEqual( int1.index_binop( gb ).op(), op_sub )
        #self.assertEqual( int1.index_binop( gb ).right(), 2 )
        #self.assertEqual( int2.index_binop( ga ).op(), op_add )
        #self.assertEqual( int2.index_binop( ga ).right(), 0 )
        #self.assertEqual( int2.index_binop( gb ).op(), op_sub )
        #self.assertEqual( int2.index_binop( gb ).right(), 1 )
        ## ... adding a vrr with index that is not initially defined with integral class
        #self.assertRaises(AssertionError, test_integral.add_vrr, 'testvrr', gc )
        ###
        # ... adding a vrr with index that is not initially defined with integral class
        self.assertRaises(AssertionError, test_integral.add_vrr, 'testvrr', gc, ( gc - 1) * int1 )

        # More typical example testing add_vrr()
        default_namer.reset_counters()
        abc = dsl_integral(ga,gb,gc,name='test')
        o_o_xabc  = dsl_scalar( expr = 1.0 / (xa+xb+xc) )
        o_o_2xabc = dsl_scalar( expr = 0.5 * o_o_xabc )
        cAcB = dsl_position( expr = A - B )
        cAcC = dsl_position( expr = A - C )
        cBcC = dsl_position( expr = B - C )
        GA = dsl_position( expr = xb * o_o_xabc * -cAcB +\
                            xc * o_o_xabc * -cAcC )
        GB = dsl_position( expr = xa * o_o_xabc * cAcB +\
                            xc * o_o_xabc * -cBcC )
        GC = dsl_position( expr = xb * o_o_xabc * cBcC +\
                            xa * o_o_xabc * cAcC )
        # Truncated 3-index overlap intergal vrr
        # GA * (a-1,b,c) + o_o_2xabc * ( ga-1 ) * (a-2,b,c)
        # dpP    int1        dp          g1         int2
        vrra_expr = dsl_binop( op_add, GA * abc.int(ga-1,gb,gc), \
                       dsl_binop( op_mul, o_o_2xabc, \
                           dsl_binop( op_mul, (ga - 1), abc.int(ga-2,gb,gc) ) \
                           ) \
                       )
        vrrb_expr = dsl_binop( op_add, GB * abc.int(ga,gb-1,gc), \
                       dsl_binop( op_mul, o_o_2xabc, \
                           dsl_binop( op_mul, (gb - 1), abc.int(gb,gb-2,gc) ) \
                           ) \
                       )
        vrrc_expr = dsl_binop( op_add, GC * abc.int(ga,gb,gc-1), \
                       dsl_binop( op_mul, o_o_2xabc, \
                           dsl_binop( op_mul, (gc - 1), abc.int(ga,gb,gc-2) ) \
                           ) \
                       )
        abc.add_vrr('vrr', ga, expr = vrra_expr)
        abc.add_vrr('vrr', gb, expr = vrrb_expr)
        abc.add_vrr('vrr', gc, expr = vrrc_expr)
        vrra = abc.rr('vrr',ga)
        vrrb = abc.rr('vrr',gb)
        vrrc = abc.rr('vrr',gc)
         
        ### Permutation of indexes in VRRs is disabled ###
        ### Explicit definition of each VRR expr is required ###
        #abc.add_vrr('vrr', gb ) # should give permuted RR expression
        #vrrb = abc.rr('vrr',gb)
        #abc.add_vrr('vrr', gc ) # should give permuted RR expression
        #vrrc = abc.rr('vrr',gc)
        ###
        vrr_list = [ vrra, vrrb, vrrc ]
        vrr_data = [ 
                    { # vrra
                       'int1' : { 'a' : [ op_sub, 1 ], 'b' : [op_add,0], 'c' : [op_add,0] },\
                       'int2' : { 'a' : [ op_sub, 2 ], 'b' : [op_add,0], 'c' : [op_add,0] },\
                       'g1'   : ga 
                    },\
                    { # vrrb
                       'int1' : { 'a' : [ op_add, 0 ], 'b' : [op_sub,1], 'c' : [op_add,0] },\
                       'int2' : { 'a' : [ op_add, 0 ], 'b' : [op_sub,2], 'c' : [op_add,0] },\
                       'g1'   : gb 
                    },\
                    { # vrrc
                       'int1' : { 'a' : [ op_add, 0 ], 'b' : [op_add,0], 'c' : [op_sub,1] },\
                       'int2' : { 'a' : [ op_add, 0 ], 'b' : [op_add,0], 'c' : [op_sub,2] },\
                       'g1'   : gc 
                    }
                ]
        for i in range( len(vrr_list) ):
            v = vrr_list[i] 
            d = vrr_data[i]
            int1 = v.expr().left().right() 
            int2 = v.expr().right().right().right()
            g1   = v.expr().right().right().left().left()
            self.assertIs( d['g1'], g1 )
            for g in [ ga, gb, gc ]:
                self.assertIs( int1.index_binop( g ).op(), d['int1'][ g.name() ][0] )
                self.assertEqual( int1.index_binop( g ).right(), d['int1'][ g.name() ][1] )
                self.assertIs( int2.index_binop( g ).op(), d['int2'][ g.name() ][0] )
                self.assertEqual( int2.index_binop( g ).right(), d['int2'][ g.name() ][1] )

        # Simple example testing add_hrr()
        test_integral = dsl_integral( ga, gb, gc, name='test_integral' )
        test_integral.add_vrr('vrr',ga, ga * test_integral.int( ga-1, gb, gc ) )
        test_integral.add_hrr('hrr',ga, gb )
        # HRR expression automatically set to
        # ( a+1, b-1 ) + AB * ( a, b-1 )
        #     int1       dpP    int2
        hrr = test_integral.rr('hrr',gb)
        hrr_data = { # hrr
                       'int1' : { 'a' : [ op_add, 1 ], 'b' : [op_sub, 1 ], 'c' : [ op_add, 0 ] },\
                       'int2' : { 'a' : [ op_add, 0 ], 'b' : [op_sub, 1 ], 'c' : [ op_add, 0 ] },\
                       'dpP'  : { 'left' : A, 'right' : B } \
                    }
        int1 = hrr.expr().left()
        int2 = hrr.expr().right().right()
        dpP  = hrr.expr().right().left()
        self.assertTrue( isinstance( hrr.expr(), dsl_binop ) )
        self.assertTrue( isinstance( int1, dsl_integral ) )
        self.assertTrue( isinstance( int2, dsl_integral ) )
        self.assertTrue( isinstance( dpP, dsl_position ) )
        for g in [ ga, gb, gc ]:
                self.assertIs( int1.index_binop( g ).op(), hrr_data['int1'][ g.name() ][0] )
                self.assertEqual( int1.index_binop( g ).right(), hrr_data['int1'][ g.name() ][1] )
                self.assertIs( int2.index_binop( g ).op(), hrr_data['int2'][ g.name() ][0] )
                self.assertEqual( int2.index_binop( g ).right(), hrr_data['int2'][ g.name() ][1] )
        self.assertIs( dpP.expr().left(), hrr_data['dpP']['left'] )
        self.assertIs( dpP.expr().right(), hrr_data['dpP']['right'] )
        # ... adding a HRR without a previously added VRR incrementing move_from
        test_integral = dsl_integral(ga,gb,name='test_integral')
        self.assertRaises(AssertionError,test_integral.add_hrr,'hrr',ga,gb)

        
class TestFunctions(unittest.TestCase):
    """Tests for functions not attached to a particular object.
    These will generally require instances of DSL objects and
    objects from dsl_extensions and are thus integration tests."""
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_expr_permute(self):
        # Set some basic DSL index and variable objects
        ga = dsl_cartesian_gaussian('a')
        xa = ga.exp()
        A  = ga.cen()
        gb = dsl_cartesian_gaussian('b')
        xb = gb.exp()
        B  = gb.cen()
        gc = dsl_cartesian_gaussian('c')
        xc = gc.exp()
        C  = gc.cen()
        m  = dsl_integer_index('m')

        # Test some basic permutations
        xp = dsl_scalar( expr = xa + xb )
        self.assertTrue( isinstance( xp.expr(), dsl_binop ) )
        self.assertIs( xp.expr().op(), op_add )
        self.assertIs( xp.expr().left(), xa )
        self.assertIs( xp.expr().right(), xb )
        xp_permuted = expr_permute( xp, ga, gb )
        #self.assertTrue( isinstance( xp.expr(), dsl_binop ) )
        #self.assertIs( xp_permuted.expr().op(), op_add )
        #self.assertIs( xp_permuted.expr().left(), xb )
        #self.assertIs( xp_permuted.expr().right(), xa )
        # New basic_symmetry_check will prevent permutation of xa and xb in xp.expr()
        self.assertTrue( isinstance( xp.expr(), dsl_binop ) )
        self.assertIs( xp_permuted, xp )
        self.assertIs( xp_permuted.expr().op(), op_add )
        self.assertIs( xp_permuted.expr().left(), xa )
        self.assertIs( xp_permuted.expr().right(), xb )

        o_o_xabc = dsl_scalar( expr = xp + xc )
        self.assertTrue( isinstance( o_o_xabc.expr(), dsl_binop ) )
        self.assertIs( o_o_xabc.expr().left(), xp )
        self.assertIs( o_o_xabc.expr().right(), xc )
        self.assertIs( o_o_xabc.expr().left().expr().left(), xa )
        self.assertIs( o_o_xabc.expr().left().expr().right(), xb )
        o_o_xabc_permuted = expr_permute( o_o_xabc, gb, ga )
        self.assertTrue( isinstance( o_o_xabc_permuted.expr(), dsl_binop ) )
        #self.assertIsNot( o_o_xabc_permuted.expr().left(), xp ) # new dsl_scalar created
        #self.assertIs( o_o_xabc_permuted.expr().right(), xc )
        #self.assertIs( o_o_xabc_permuted.expr().left().expr().left(), xb )
        #self.assertIs( o_o_xabc_permuted.expr().left().expr().right(), xa )
        # New basic_symmetry check will prevent permutation of xa and xb in xp.expr()
        self.assertIs( o_o_xabc_permuted.expr().left(), xp )
        self.assertIs( o_o_xabc_permuted.expr().right(), xc )
        self.assertIs( o_o_xabc_permuted.expr().left().expr().left(), xa )
        self.assertIs( o_o_xabc_permuted.expr().left().expr().right(), xb )
        o_o_xabc_permuted = expr_permute( o_o_xabc, gc, ga )
        self.assertTrue( isinstance( o_o_xabc_permuted.expr(), dsl_binop ) )
        self.assertIsNot( o_o_xabc_permuted.expr().left(), xp ) # new dsl_scalar created
        self.assertIs( o_o_xabc_permuted.expr().right(), xa )
        self.assertIs( o_o_xabc_permuted.expr().left().expr().left(), xc )
        self.assertIs( o_o_xabc_permuted.expr().left().expr().right(), xb )
        self.assertTrue( isinstance( o_o_xabc.expr(), dsl_binop ) ) # original unchanged
        self.assertIs( o_o_xabc.expr().left(), xp )
        self.assertIs( o_o_xabc.expr().right(), xc )
        self.assertIs( o_o_xabc.expr().left().expr().left(), xa )
        self.assertIs( o_o_xabc.expr().left().expr().right(), xb )

        AB = dsl_position( expr = A - B )
        self.assertTrue( isinstance( AB.expr(), dsl_binop ) )
        self.assertIs( AB.expr().left(), A )
        self.assertIs( AB.expr().right(), B )
        AB_permuted = expr_permute( AB, ga, gb )
        self.assertTrue( isinstance( AB_permuted.expr(), dsl_binop ) )
        self.assertIs( AB_permuted.expr().left(), B )
        self.assertIs( AB_permuted.expr().right(), A )
        AB_permuted = expr_permute( AB, ga, gc )
        self.assertTrue( isinstance( AB_permuted.expr(), dsl_binop ) )
        self.assertIs( AB_permuted.expr().left(), C )
        self.assertIs( AB_permuted.expr().right(), B )

        test_expr = exp( xc * AB )
        self.assertTrue( isinstance( test_expr, dsl_unop ) )
        self.assertTrue( isinstance( test_expr.arg(), dsl_binop ) )
        self.assertIs( test_expr.arg().left(), xc )
        self.assertIs( test_expr.arg().right(), AB )
        self.assertIs( test_expr.arg().right().expr().left(), A )
        self.assertIs( test_expr.arg().right().expr().right(), B )
        test_expr_permuted = expr_permute( test_expr, ga, gb )
        self.assertTrue( isinstance( test_expr_permuted, dsl_unop ) )
        self.assertTrue( isinstance( test_expr_permuted.arg(), dsl_binop ) )
        self.assertIs( test_expr_permuted.arg().left(), xc )
        self.assertIsNot( test_expr_permuted.arg().right(), AB ) # new dsl_position created
        self.assertIs( test_expr_permuted.arg().right().expr().left(), B )
        self.assertIs( test_expr_permuted.arg().right().expr().right(), A )
        test_expr_permuted = expr_permute( test_expr, ga, gc )
        self.assertTrue( isinstance( test_expr_permuted, dsl_unop ) )
        self.assertTrue( isinstance( test_expr_permuted.arg(), dsl_binop ) )
        self.assertIs( test_expr_permuted.arg().left(), xa )
        self.assertIsNot( test_expr_permuted.arg().right(), AB ) # new dsl_position created
        self.assertIs( test_expr_permuted.arg().right().expr().left(), C )
        self.assertIs( test_expr_permuted.arg().right().expr().right(), B )

        # Test exceptions / assertions
        AB = dsl_position( expr = A - B )
        gd = dsl_cartesian_gaussian('d') 
        gd._attachment_dict = {} # ... should result in KeyError on permutation
        self.assertRaises(KeyError, expr_permute, AB, ga, gd )

        test_expr = exp( xc * AB )
        self.assertRaises(AssertionError, expr_permute, test_expr, A, B ) # only permute dsl_index

        test_expr = xa * dsl_integral( ga, gb, m, name = 'test' )
        self.assertRaises(AssertionError, expr_permute, test_expr, ga, gc ) # gc not in dsl_integral
        self.assertRaises(AssertionError, expr_permute, test_expr, ga, m ) # ga, m not same type

        # Test permutation of expressions containing dsl_integral objects
        # (Requires expr_list_objects from dsl_extensions.py)
        default_namer.reset_counters()
        abc  = dsl_integral( ga, gb, gc, m, name = 'abc' )
        xp   = dsl_scalar( expr = xa + xb )
        AB   = dsl_position( expr = A - B )

        # AB * (a-1,b,c) + xp * (a-1,b,c-1)
        # dpP    int1      dp       int2
        expr_0 = dsl_binop( op_add, dsl_binop( op_mul, AB, abc.int(ga-1,gb,gc) ), \
                             dsl_binop( op_mul, xp, abc.int(ga-2,gb,gc-1,m+1) ) )
        # ... permuted expressions
        expr_1 = expr_permute( expr_0, ga, gb )
        expr_2 = expr_permute( expr_0, gb, ga )
        expr_3 = expr_permute( expr_0, gc, ga )
        expr_4 = expr_permute( expr_0, gb, gc )
        expr_list = [ expr_0, expr_1, expr_2, expr_3, expr_4 ]
        # ... check outcome of permutations against data in dictionary
        expr_dict = { 
                    '0' : { 'dpP' : { 'left' : A, 'right' : B },\
                       'int1': { 'a' : [ op_sub, 1 ], \
                                 'b' : [ op_add, 0 ], \
                                 'c' : [ op_add, 0 ], \
                                 'm' : [ op_add, 0 ] \
                               },\
                       'dp'  : { 'left' : xa, 'right' : xb },\
                        'int2': { 'a' : [ op_sub, 2 ], \
                                  'b' : [ op_add, 0 ], \
                                  'c' : [ op_sub, 1 ], \
                                  'm' : [ op_add, 1 ] \
                               }\
                          },\
                    '1' : { 'dpP' : { 'left' : B, 'right' : A },\
                       'int1': { 'a' : [ op_add, 0 ], \
                                 'b' : [ op_sub, 1 ], \
                                 'c' : [ op_add, 0 ], \
                                 'm' : [ op_add, 0 ] \
                               },
#                       'dp'  : { 'left' : xb, 'right' : xa },\
                       'dp'  : { 'left' : xa, 'right' : xb },\
                        'int2': { 'a' : [ op_add, 0 ], \
                                  'b' : [ op_sub, 2 ], \
                                  'c' : [ op_sub, 1 ], \
                                  'm' : [ op_add, 1 ] \
                               }\
                          },\
                    '2' : { 'dpP' : { 'left' : B, 'right' : A },\
                       'int1': { 'a' : [ op_add, 0 ], \
                                 'b' : [ op_sub, 1 ], \
                                 'c' : [ op_add, 0 ], \
                                 'm' : [ op_add, 0 ] \
                               },
#                       'dp'  : { 'left' : xb, 'right' : xa },\
                       'dp'  : { 'left' : xa, 'right' : xb },\
                        'int2': { 'a' : [ op_add, 0 ], \
                                  'b' : [ op_sub, 2 ], \
                                  'c' : [ op_sub, 1 ], \
                                  'm' : [ op_add, 1 ] \
                               }\
                          },\
                    '3' : { 'dpP' : { 'left' : C, 'right' : B },\
                       'int1': { 'a' : [ op_add, 0 ], \
                                 'b' : [ op_add, 0 ], \
                                 'c' : [ op_sub, 1 ], \
                                 'm' : [ op_add, 0 ] \
                               },
                       'dp'  : { 'left' : xc, 'right' : xb },\
                        'int2': { 'a' : [ op_sub, 1 ], \
                                  'b' : [ op_add, 0 ], \
                                  'c' : [ op_sub, 2 ], \
                                  'm' : [ op_add, 1 ] \
                               }\
                          },\
                    '4' : { 'dpP' : { 'left' : A, 'right' : C },\
                       'int1': { 'a' : [ op_sub, 1 ], \
                                 'b' : [ op_add, 0 ], \
                                 'c' : [ op_add, 0 ], \
                                 'm' : [ op_add, 0 ] \
                               },
                       'dp'  : { 'left' : xa, 'right' : xc },\
                        'int2': { 'a' : [ op_sub, 2 ], \
                                  'b' : [ op_sub, 1 ], \
                                  'c' : [ op_add, 0 ], \
                                  'm' : [ op_add, 1 ] \
                               }\
                          }\
                    }
        for i in range( len( expr_list ) ):
            expr = expr_list[i]
            expr_data = expr_dict[ str(i) ]
            dpP  = expr.left().left()
            int1 = expr.left().right()
            dp   = expr.right().left()
            int2 = expr.right().right()
            for v, vname in [ (dpP, 'dpP'), (dp, 'dp' ) ]:
                self.assertIs( v.expr().left(), expr_data[vname]['left'] )
                self.assertIs( v.expr().right(), expr_data[vname]['right'] )
            for integral, intname in [ (int1, 'int1'), (int2, 'int2') ]:
                for g in [ ga, gb, gc ]:
                    self.assertIs( integral.index_binop( g ).op(), expr_data[intname][ g.name() ][0] )
                    self.assertIs( integral.index_binop( g ).right(), expr_data[intname][ g.name() ][1] )

        # Check for swapping auxiliary indexes
        x = dsl_scalar( name = 'x' )
        m = dsl_integer_index(name = 'm')
        n = dsl_integer_index(name = 'n')
#        amn = dsl_integral( ga, m, n )
#        expr_0 = exp( x * amn( a-1, m+1, n ) )
#        expr_1 = expr_permute( expr_0, m, n )
#        expr_list = [ expr_0, expr_1 ]
#        expr_data = [ { 'm' : [ op_add, 1 ], 'n' : [ op_add, 0 ] },\
#                      { 'm' : [ op_add, 0 ], 'n' : [ op_add, 1 ] } ]
#        for i in len( expr_list ):
#            x = expr_list[i]
#            x_data = expr_data[i]
#            self.assertTrue( isinstance( x, dsl_unop ) )
#            self.assertTrue( isinstance( x.arg(), dsl_binop ) )
#            self.assertIs( x.arg().left(), x )
#            self.assertTrue( isinstance( x.arg().right(), dsl_integral ) )
#            for a in [ m, n ]:
#                self.assertIs( x.arg().right().aux_binop( a ).op(), x_data[ a.name() ][0] )
#                self.assertIs( x.arg().right().aux_binop( a ).right(), x_data[ a.name() ][1] )
        # 02/09/2014: Only one auxiliary currently allowed - check for AssertionError
        self.assertRaises( AssertionError, dsl_integral, ga, m, n, x )

    def test_expr_basic_symmetry_check(self):
        # Testing for the situation where expr_basic_symmetry_check is to be applied --
        # for the checking of symmetry in expr set for dsl_variable objects with respect
        # to permutation of integral indexes.

        # Set some basic DSL index and variable objects
        ga = dsl_cartesian_gaussian('a')
        xa = ga.exp()
        A  = ga.cen()
        gb = dsl_cartesian_gaussian('b')
        xb = gb.exp()
        B  = gb.cen()
        gc = dsl_cartesian_gaussian('c')
        xc = gc.exp()
        C  = gc.cen()
        m  = dsl_integer_index('m')

        # Some simple dsl_variable objects with expr set
        xp = dsl_scalar( expr = xa + xb )
        xq = dsl_scalar( expr = xb + xc )
        o_o_xabc1 = dsl_scalar( expr = 1.0 / ( xp + xc ) )
        o_o_xabc2 = dsl_scalar( expr = 1.0 / ( xa + xq ) )
        o_o_xp  = dsl_scalar( expr = 1.0 / xp )
        RAB  = dsl_position( expr = A - B )
        ABC = dsl_position( expr = A * B * C )
        AB  = dsl_position( expr = A * B )
        # Permutation symmetric
        self.assertTrue( expr_basic_symmetry_check( xp, ga, gb ) )
        self.assertTrue( expr_basic_symmetry_check( xp, gb, ga ) )
        self.assertTrue( expr_basic_symmetry_check( o_o_xabc1, ga, gb ) )
        self.assertTrue( expr_basic_symmetry_check( o_o_xabc1, gb, ga ) )
        self.assertTrue( expr_basic_symmetry_check( o_o_xabc2, gb, gc ) )
        self.assertTrue( expr_basic_symmetry_check( o_o_xabc2, gc, gb ) )
        self.assertTrue( expr_basic_symmetry_check( o_o_xp, ga, gb ) )
        self.assertTrue( expr_basic_symmetry_check( ABC, ga, gb ) )
        self.assertTrue( expr_basic_symmetry_check( ABC, gb, ga ) )
        # Not permutation symmetric
        self.assertFalse( expr_basic_symmetry_check( xp, ga, gc ) )
        self.assertFalse( expr_basic_symmetry_check( xp, gc, ga ) )
        # A - B is not fully associative
        self.assertFalse( expr_basic_symmetry_check( RAB, ga, gb ) )
        self.assertFalse( expr_basic_symmetry_check( RAB, gb, ga ) )
        self.assertFalse( expr_basic_symmetry_check( RAB, ga, gc ) )
        self.assertFalse( expr_basic_symmetry_check( RAB, gc, ga ) )
        self.assertFalse( expr_basic_symmetry_check( o_o_xp, gc, gb ) )
        self.assertFalse( expr_basic_symmetry_check( o_o_xp, gb, gc ) )
        # Since xp.expr(), as part of o_o_xabc1 is permuted, this is False
        self.assertFalse( expr_basic_symmetry_check( o_o_xabc1, gb, gc ) )
        self.assertFalse( expr_basic_symmetry_check( o_o_xabc1, gc, gb ) )
        # Since xq.expr(), as part of o_o_xabc2 is permuted, this is False
        self.assertFalse( expr_basic_symmetry_check( o_o_xabc2, ga, gb ) )
        self.assertFalse( expr_basic_symmetry_check( o_o_xabc2, gb, ga ) )

        # More complicated expressions
        s1 = dsl_scalar()
        p1 = dsl_position()
        v1 = dsl_value(name = 'testvalue')
        f1 = dsl_function(name = 'testfunction')
        z1 = dsl_zero()
        g1 = dsl_cartesian_gaussian( name = 'g1' )
        # additional dsl_variable and dsl_value objects involved
        self.assertTrue( expr_basic_symmetry_check( xa + xb + s1, ga, gb ) )
        self.assertTrue( expr_basic_symmetry_check( A + B + s1, ga, gb ) )
        self.assertTrue( expr_basic_symmetry_check( A + B - ( s1 + z1 ), ga, gb ) )
        # 
        self.assertFalse( expr_basic_symmetry_check( A + p1 - ( xa - xb ), ga, gb ) )
        self.assertFalse( expr_basic_symmetry_check( A - B  - ( f1 + s1 ), ga, gb ) )
        self.assertFalse( expr_basic_symmetry_check( v1* A  -  p1 * B , ga, gb ) )
        # dsl_unop objects involved
        self.assertTrue( expr_basic_symmetry_check( g1 + s1 - exp( z1 ), ga, gb ) ) 
        self.assertTrue( expr_basic_symmetry_check( g1 + s1 - exp( xb + xa ), ga, gb ) )
        #
        self.assertFalse( expr_basic_symmetry_check( g1 + s1 - exp( xb - xa ), ga, gb ) )
        self.assertFalse( expr_basic_symmetry_check( g1 + s1 - exp( xb + xa ), gc, gb ) )
        self.assertFalse( expr_basic_symmetry_check( -B + A , ga, gb ) )
        # Some combinations of dsl_variable objects with expr set
        self.assertTrue( expr_basic_symmetry_check( xp + xc, ga, gb ) )
        self.assertTrue( expr_basic_symmetry_check( xp + AB, ga, gb ) )
        self.assertTrue( expr_basic_symmetry_check( AB + (A*B)*C, ga, gb ) )
        #
        self.assertFalse( expr_basic_symmetry_check( xp + xq, gb, gc ) )
        self.assertFalse( expr_basic_symmetry_check( xp + RAB, ga, gb ) )
        self.assertFalse( expr_basic_symmetry_check( AB + A*(B*C), ga, gb ) )

        # Test AssertionError from presence of dsl_integral object
        int1 = dsl_integral( ga, gb, gc, name='test' )
        bad_scalar = dsl_scalar( expr = xa + xb + xc + int1 )
        self.assertRaises(AssertionError,expr_basic_symmetry_check,int1, ga, gb)
        self.assertRaises(AssertionError,expr_basic_symmetry_check,xa+xb+int1, ga, gb)
        self.assertRaises(AssertionError,expr_basic_symmetry_check,bad_scalar, ga, gb)

    def test_expr_deepish_copy(self):
        # Create some DSL objects for use in expressions
        v1 = dsl_value(name = 'v1', value = 42 )
        s1 = dsl_scalar(name = 's1' )
        p1 = dsl_position( name = 'p1' )
        z1 = dsl_zero()
        f1 = dsl_function('f1')
        in1= dsl_integer_index('in1')
        g1 = dsl_cartesian_gaussian(name = 'g1')
        s2 = dsl_scalar(name = 's2' )
        p2 = dsl_position( name = 'p2' )
        g2 = dsl_cartesian_gaussian(name = 'g2')
        int1 = dsl_integral( g1, g2, name = 'test' )
        int2 = dsl_integral( g1, in1, name = 'test' )
        # Create and test some DSL expressions
        expr =  s1 + s2
        copied_expr = expr_deepish_copy( expr )
        self.assertIsNot( copied_expr, expr )
        self.assertIs( copied_expr.left(), s1 )
        self.assertIs( copied_expr.right(), s2 )
        #
        expr = dsl_binop( op_add, p1 * int1, p2 * int2 )
        copied_expr = expr_deepish_copy( expr )
        self.assertIsNot( copied_expr, expr )
        self.assertIsNot( copied_expr.left(), expr.left() )
        self.assertIsNot( copied_expr.right(), expr.right() )
        self.assertIs( copied_expr.left().left(), expr.left().left() )
        self.assertIs( copied_expr.left().right(), expr.left().right() )
        self.assertIs( copied_expr.right().left(), expr.right().left() )
        self.assertIs( copied_expr.right().right(), expr.right().right() )
        #
        expr = dsl_binop( op_mul, z1, s1 * v1 )
        copied_expr = expr_deepish_copy( expr )
        self.assertIsNot( copied_expr, expr )
        self.assertIs( copied_expr.left(), expr.left() )
        self.assertIsNot( copied_expr.right(), expr.right() )
        self.assertIs( copied_expr.right().left(), expr.right().left() )
        self.assertIs( copied_expr.right().right(), expr.right().right() )
        #
        expr = dsl_unop( op_exp, s1 * f1 )
        copied_expr = expr_deepish_copy( expr )
        self.assertIsNot( copied_expr, expr )
        self.assertIsNot( copied_expr.arg(), expr.arg() )
        self.assertIs( copied_expr.arg().left(), expr.arg().left() )
        self.assertIs( copied_expr.arg().right(), expr.arg().right() )
        #
        expr = dsl_binop( op_div, dsl_unop( op_int, s2 ), dsl_unop( op_exp, int2 ) )
        copied_expr = expr_deepish_copy( expr )
        self.assertIsNot( copied_expr, expr )
        self.assertIsNot( copied_expr.left(), expr.left() )
        self.assertIsNot( copied_expr.right(), expr.right() )
        self.assertIs( copied_expr.left().arg(), expr.left().arg() )
        self.assertIs( copied_expr.right().arg(), expr.right().arg() )
                        
    def test_expr_find_and_replace(self):
        # Create some DSL objects for use in expressions
        v1 = dsl_value(name = 'v1', value = 42 )
        s1 = dsl_scalar(name = 's1' )
        p1 = dsl_position( name = 'p1' )
        z1 = dsl_zero()
        f1 = dsl_function('f1')
        in1= dsl_integer_index('in1')
        g1 = dsl_cartesian_gaussian(name = 'g1')
        s2 = dsl_scalar(name = 's2' )
        p2 = dsl_position( name = 'p2' )
        g2 = dsl_cartesian_gaussian(name = 'g2')
        int1 = dsl_integral( g1, g2, name = 'int1' )
        int2 = dsl_integral( g1, in1, name = 'int2' )
        # Create and test some DSL expressions
        expr =  s1 + s2
        modified_expr = expr_find_and_replace( expr, s1, dsl_zero() )
        self.assertTrue( isinstance( modified_expr, dsl_binop ) )
        self.assertIs( modified_expr.op(), op_add )
        self.assertTrue( isinstance( modified_expr.left(), dsl_zero ) )
        self.assertIs( modified_expr.right(), s2 )
        #
        expr = dsl_binop( op_add, p1 * int1, p2 * int2 )
        modified_expr = expr_find_and_replace( expr, int1, v1 )
        self.assertTrue( isinstance( modified_expr, dsl_binop ) )
        self.assertIs( modified_expr.op(), op_add )
        self.assertTrue( isinstance( modified_expr.left(), dsl_binop ) )
        self.assertTrue( isinstance( modified_expr.right(), dsl_binop ) )
        self.assertIs( modified_expr.left().left(), p1 )
        self.assertIs( modified_expr.left().right(), v1 )
        self.assertIs( modified_expr.right().left(), p2 )
        self.assertIs( modified_expr.right().right(), int2 )
        #
        expr = dsl_binop( op_mul, z1, s1 * v1 )
        modified_expr = expr_find_and_replace( expr, z1, in1 )
        self.assertIs( modified_expr.op(), op_mul )
        self.assertTrue( isinstance( modified_expr, dsl_binop ) )
        self.assertTrue( isinstance( modified_expr.right(), dsl_binop ) )
        self.assertIs( modified_expr.left(), in1 )
        self.assertIs( modified_expr.right().op(), op_mul )
        self.assertIs( modified_expr.right().left(), s1 )
        self.assertIs( modified_expr.right().right(), v1 )
        #
        expr = dsl_unop( op_exp, s1 * f1 )
        modified_expr = expr_find_and_replace( expr, s1, f1 )
        self.assertIs( modified_expr.op(), op_exp )
        self.assertTrue( isinstance(modified_expr, dsl_unop ) )
        self.assertTrue( isinstance(modified_expr.arg(), dsl_binop ) )
        self.assertIs( modified_expr.arg().left(), f1 )
        self.assertIs( modified_expr.arg().right(), f1 )
        #
        expr = dsl_binop( op_div, dsl_unop( op_int, s2 ), dsl_unop( op_exp, int2 ) )
        modified_expr = expr_find_and_replace( expr, s2, z1 )
        self.assertIs( modified_expr.op(), op_div )
        self.assertTrue( isinstance( modified_expr, dsl_binop ) )
        self.assertIs( modified_expr.left().op(), op_int )
        self.assertTrue( isinstance( modified_expr.left(), dsl_unop ) )
        self.assertIs( modified_expr.left().arg(), z1 )
        self.assertIs( modified_expr.right().op(), op_exp )
        self.assertTrue( isinstance( modified_expr.right(), dsl_unop ) )
        self.assertIs( modified_expr.right().arg(), int2 )

    def test_expr_simplify(self):
        # Create some DSL objects for use in expressions
        v1 = dsl_value(name = 'v1', value = 42 )
        s1 = dsl_scalar(name = 's1' )
        p1 = dsl_position( name = 'p1' )
        z1 = dsl_zero()
        f1 = dsl_function('f1')
        in1= dsl_integer_index('in1')
        g1 = dsl_cartesian_gaussian(name = 'g1')
        s2 = dsl_scalar(name = 's2' )
        p2 = dsl_position( name = 'p2' )
        g2 = dsl_cartesian_gaussian(name = 'g2')
        int1 = dsl_integral( g1, g2, name = 'int1' )
        int2 = dsl_integral( g1, in1, name = 'int2' )
        # Create and test some DSL expressions
        expr =  s1 + s2 * z1
        simplified_expr = expr_simplify( expr )
        self.assertIs( simplified_expr, s1 )
        # check that original expr is unchanged
        self.assertIs( expr.op(), op_add )
        self.assertTrue( isinstance( expr, dsl_binop ) )
        self.assertIs( expr.left(), s1 )
        self.assertTrue( isinstance( expr.right(), dsl_binop ) )
        self.assertIs( expr.right().op(), op_mul )
        self.assertIs( expr.right().left(), s2 )
        self.assertIs( expr.right().right(), z1 ) 
        #
        expr =  s1 * z1 - s2
        simplified_expr = expr_simplify( expr )
        self.assertTrue( isinstance( simplified_expr, dsl_unop ) )
        self.assertIs( simplified_expr.op(), op_neg )
        self.assertIs( simplified_expr.arg(), s2 )
        # Check that original expr is unchanged
        self.assertIs( expr.op(), op_sub )
        self.assertTrue( isinstance( expr, dsl_binop ) )
        self.assertTrue( isinstance( expr.left(), dsl_binop ) )
        self.assertIs( expr.left().op(), op_mul ) 
        self.assertIs( expr.left().left(), s1 )
        self.assertIs( expr.left().right(), z1 )
        self.assertIs( expr.right(), s2 )
        #
        expr = dsl_binop( op_add, p1 * int1, dsl_zero() * int2 )
        simplified_expr = expr_simplify( expr )
        self.assertTrue( isinstance( simplified_expr, dsl_binop ) )
        self.assertIs( simplified_expr.op(), op_mul )
        self.assertIs( simplified_expr.left(), p1 )
        self.assertIs( simplified_expr.right(), int1 )
        #
        expr = dsl_binop( op_mul, z1, s1 * v1 )
        simplified_expr = expr_simplify( expr )
        self.assertTrue( isinstance( simplified_expr, dsl_zero ) )
        #
        expr = dsl_unop( op_exp, z1 * f1 ) 
        simplified_expr = expr_simplify( expr )
        self.assertIs( simplified_expr.op(), op_exp )
        self.assertTrue( isinstance(simplified_expr, dsl_unop ) )
        self.assertTrue( isinstance(simplified_expr.arg(), dsl_zero ) )
        #
        expr = dsl_binop( op_div, dsl_unop( op_int, s2 * dsl_zero() - 1 ), dsl_unop( op_exp, int2 * dsl_zero() ) )
        simplified_expr = expr_simplify( expr )
        self.assertIs( simplified_expr.op(), op_div )
        self.assertTrue( isinstance( simplified_expr, dsl_binop ) )
        self.assertIs( simplified_expr.left().op(), op_int )
        self.assertEqual( simplified_expr.left().arg(), -1 )
        self.assertTrue( isinstance( simplified_expr.left(), dsl_unop ) )
        self.assertIs( simplified_expr.right().op(), op_exp )
        self.assertTrue( isinstance( simplified_expr.right(), dsl_unop ) )
        self.assertTrue( isinstance( simplified_expr.right().arg(), dsl_zero ) )
        # Check that dsl_binop consisting of number.Real objects and dsl_integer_index objects
        # are correctly simplified
        in1.set_value(3)
        expr = dsl_binop( op_sub, 0 * in1, dsl_zero() + int1 )
        simplified_expr = expr_simplify( expr )
        self.assertTrue( isinstance( simplified_expr, dsl_unop ) )
        self.assertIs( simplified_expr.op(), op_neg )
        self.assertIs( simplified_expr.arg(), int1 )
        #
        in1.set_value(2)
        expr = dsl_binop( op_sub, in1 * 3.3, dsl_zero() * int2 - int1 )
        simplified_expr = expr_simplify( expr )
        self.assertTrue( isinstance( simplified_expr, dsl_binop ) )
        self.assertIs( simplified_expr.op(), op_sub )
        self.assertEqual( simplified_expr.left(), 6.6 )
        self.assertTrue( isinstance( simplified_expr.right(), dsl_unop ) )
        self.assertIs( simplified_expr.right().arg(), int1 )

    def test_expr_list_objects(self):
        # Create some DSL objects for use in expressions
        v1 = dsl_value(name = 'v1', value = 42 )
        s1 = dsl_scalar(name = 's1' )
        p1 = dsl_position( name = 'p1' )
        z1 = dsl_zero()
        f1 = dsl_function('f1')
        in1= dsl_integer_index('in1')
        g1 = dsl_cartesian_gaussian(name = 'g1')
        s2 = dsl_scalar(name = 's2' )
        p2 = dsl_position( name = 'p2' )
        g2 = dsl_cartesian_gaussian(name = 'g2')
        int1 = dsl_integral( g1, g2, name = 'int1' )
        int2 = dsl_integral( g1, in1, name = 'int2' )
        # Create and test some DSL expressions
        expr =  s1 + s2 + s1
        scalar_list = expr_list_objects( expr, dsl_scalar )
        position_list = expr_list_objects( expr, dsl_position )
        integral_list = expr_list_objects( expr, dsl_integral )
        self.assertEqual( len(scalar_list), 2 )
        self.assertIn( s1, scalar_list )
        self.assertIn( s2, scalar_list )
        self.assertEqual( len(position_list), 0 )
        self.assertEqual( len(integral_list), 0 )
        #
        expr = dsl_binop( op_add, p1 * int1, p2 * p2 * int2 )
        scalar_list   = expr_list_objects( expr, dsl_scalar )
        position_list = expr_list_objects( expr, dsl_position )
        integral_list = expr_list_objects( expr, dsl_integral )
        self.assertEqual( len(scalar_list), 0 )
        self.assertEqual( len(position_list), 2 )
        self.assertIn( p1, position_list )
        self.assertIn( p2, position_list )
        self.assertEqual( len(integral_list), 2 )
        self.assertIn( int1, integral_list )
        self.assertIn( int2, integral_list )
        #
        expr = dsl_binop( op_mul, z1, s1 * v1 )
        scalar_list = expr_list_objects( expr, dsl_scalar )
        position_list = expr_list_objects( expr, dsl_position )
        integral_list = expr_list_objects( expr, dsl_integral )
        zero_list = expr_list_objects( expr, dsl_zero )
        value_list = expr_list_objects( expr, dsl_value )
        self.assertEqual( len(scalar_list), 1 )
        self.assertIn( s1, scalar_list )
        self.assertEqual( len(position_list), 0 )
        self.assertEqual( len(integral_list), 0 )
        self.assertEqual( len(zero_list), 1 )
        self.assertIn( z1, zero_list )
        self.assertEqual( len(value_list), 1 )
        self.assertIn( v1, value_list )
        #
        expr = dsl_unop( op_exp, s1 * f1 )
        scalar_list = expr_list_objects( expr, dsl_scalar )
        function_list  = expr_list_objects( expr, dsl_function )
        self.assertEqual( len(scalar_list), 1 )
        self.assertIn( s1, scalar_list )
        self.assertEqual( len(function_list), 1 )
        self.assertIn( f1, function_list )
        # Test with parent classes (should list all derived classes)
        expr = dsl_binop( op_div, dsl_unop( op_int, in1 * s2 + int2 ), dsl_unop( op_exp, int1 * int2 * g2 ) )
        variable_list = expr_list_objects( expr, dsl_variable ) # dsl_scalar, dsl_position
        index_list = expr_list_objects( expr, dsl_index ) # dsl_integer_index, dsl_cartesian_gaussian
        integral_list = expr_list_objects( expr, dsl_integral )
        self.assertEqual( len( variable_list ), 1 )
        self.assertIn( s2, variable_list )
        self.assertEqual( len( index_list ), 2 )
        self.assertIn( in1, index_list )
        self.assertIn( g2, index_list )
        self.assertEqual( len( integral_list ), 2 )
        self.assertIn( int1, integral_list )
        self.assertIn( int2, integral_list )
