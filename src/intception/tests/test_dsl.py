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
# Note that the src() function from printer.py is imported via intception.dsl.
# This is important, as when it is passed an object, it first tries to call a
# src method, and if this is not available, a __str__ method is called.

class TestDSLObjects(unittest.TestCase):
    """Tests for individual objects which define the DSL.
    Unit tests are only written for non-trivial methods, i.e.
    not methods which just return a member variable."""
    def setUp(self):
        # Reset counters for default autonamer(s)
        default_namer.reset_counters()
    
    def tearDown(self):
        pass

    def test_dsl_op(self):
        # Test __eq__ vs is
        op1 = dsl_op('+','+','infix',prec=6,assoc='full')
        op2 = dsl_op('+','+','infix',prec=6,assoc='full')
        op3 = dsl_op('abs','abs','function')
        self.assertEqual(op1,op2)
        self.assertEqual(op1,op1)
        self.assertEqual(op2,op2)
        self.assertNotEqual(op1,op3)
        self.assertNotEqual(op2,op3)
        self.assertIsNot(op1,op2)
        del op1, op2, op3

    def test_dsl_base(self):
        # Test that overloaded operators function actually return
        # expected dsl_binop and dsl_unop objects
        b1 = dsl_base()
        b2 = dsl_base()
        i1 = 1
        dp1 = 1.0
        # __add__(), __radd__()
        op = op_add
        tests  = [ b1 + b2, b1 + i1, i1 + b2, b1 + dp1, dp1 + b2 ]
        lefts  = [ b1, b1, i1, b1, dp1 ]
        rights = [ b2, i1, b2, dp1,  b2 ]
        for ii in range( len(tests) ):
            self.assertTrue(isinstance(tests[ii],dsl_binop))
            self.assertIs(tests[ii].op(), op )
            self.assertIs(tests[ii].left(), lefts[ii] )
            self.assertIs(tests[ii].right(), rights[ii] )
        # __sub__(), __rsub__()
        op = op_sub
        tests  = [ b1 - b2, b1 - i1, i1 - b2, b1 - dp1, dp1 - b2 ]
        lefts  = [ b1, b1, i1, b1, dp1 ]
        rights = [ b2, i1, b2, dp1,  b2 ]
        for ii in range( len(tests) ):
            self.assertTrue(isinstance(tests[ii],dsl_binop))
            self.assertIs(tests[ii].op(), op )
            self.assertIs(tests[ii].left(), lefts[ii] )
            self.assertIs(tests[ii].right(), rights[ii] )
        # __mul__(), __rmul__()
        op = op_mul
        tests  = [ b1 * b2, b1 * i1, i1 * b2, b1 * dp1, dp1 * b2 ]
        lefts  = [ b1, b1, i1, b1, dp1 ]
        rights = [ b2, i1, b2, dp1,  b2 ]
        for ii in range( len(tests) ):
            self.assertTrue(isinstance(tests[ii],dsl_binop))
            self.assertIs(tests[ii].op(), op )
            self.assertIs(tests[ii].left(), lefts[ii] )
            self.assertIs(tests[ii].right(), rights[ii] )
        # __truediv__(), __rtruediv__()
        op = op_div
        tests  = [ b1 / b2, b1 / i1, i1 / b2, b1 / dp1, dp1 / b2 ]
        lefts  = [ b1, b1, i1, b1, dp1 ]
        rights = [ b2, i1, b2, dp1,  b2 ]
        for ii in range( len(tests) ):
            self.assertTrue(isinstance(tests[ii],dsl_binop))
            self.assertIs(tests[ii].op(), op )
            self.assertIs(tests[ii].left(), lefts[ii] )
            self.assertIs(tests[ii].right(), rights[ii] )
        # __floordiv__(), __rfloordiv__()
        op = op_floor
        tests  = [ b1 // b2, b1 // i1, i1 // b2, b1 // dp1, dp1 // b2 ]
        lefts  = [ b1, b1, i1, b1, dp1 ]
        rights = [ b2, i1, b2, dp1,  b2 ]
        for ii in range( len(tests) ):
            self.assertTrue(isinstance(tests[ii],dsl_binop))
            self.assertIs(tests[ii].op(), op )
            self.assertIs(tests[ii].left(), lefts[ii] )
            self.assertIs(tests[ii].right(), rights[ii] )
        # __pow__(), __rpow__()
        op = op_pow
        tests  = [ b1 ** b2, b1 ** i1, i1 ** b2, b1 ** dp1, dp1 ** b2 ]
        lefts  = [ b1, b1, i1, b1, dp1 ]
        rights = [ b2, i1, b2, dp1,  b2 ]
        for ii in range( len(tests) ):
            self.assertTrue(isinstance(tests[ii],dsl_binop))
            self.assertIs(tests[ii].op(), op )
            self.assertIs(tests[ii].left(), lefts[ii] )
            self.assertIs(tests[ii].right(), rights[ii] )
        # assign()
        op = op_assign
        tests  = [ b1.assign(b2), b1.assign(i1), b1.assign(dp1) ]
        lefts  = [ b1, b1, b1 ]
        rights = [ b2 ,i1, dp1 ]
        for ii in range( len(tests) ):
            self.assertTrue(isinstance(tests[ii],dsl_binop))
            self.assertIs(tests[ii].op(), op )
            self.assertIs(tests[ii].left(), lefts[ii] )
            self.assertIs(tests[ii].right(), rights[ii] )
        # __neg__()
        test = -b1
        self.assertTrue(isinstance(test,dsl_unop))
        self.assertIs(test.op(),op_neg)
        self.assertIs(test.arg(),b1)
        # __pos__()
        test = +b1
        self.assertTrue(isinstance(test,dsl_unop))
        self.assertIs(test.op(),op_pos)
        self.assertIs(test.arg(),b1)
        # __abs__()
        test = abs(b1)
        self.assertTrue(isinstance(test,dsl_unop))
        self.assertIs(test.op(),op_abs)
        self.assertIs(test.arg(),b1)

    def test_dsl_unop(self): # inherits from dsl_base
        op1 = dsl_op('u+','+','prefix')
        op2 = dsl_op('u-','-','prefix')
        arg1 = dsl_scalar(name = 'arg1')
        arg2 = dsl_scalar(name = 'arg2')
        # Test src() 
        unop1 = dsl_unop( op1, arg1 )
        unop2 = dsl_unop( op2, unop1 )
        unop3 = dsl_unop( op2, arg2 )
        self.assertEqual( src( unop1 ), '+arg1' )
        self.assertEqual( src( unop2 ), '-( +arg1 )' )
        self.assertEqual( src( unop3 ), '-arg2' )

    def test_dsl_binop(self): # inherits from dsl_base
        op1 = dsl_op('*','*','infix',prec=5,assoc='full')
        op2 = dsl_op('+','+','infix',prec=6,assoc='full')
        op3 = dsl_op('-','-','infix',prec=6,assoc='left')
        op4 = dsl_op('**','pow','function')
        arg1 = dsl_scalar(name = 'arg1')
        arg2 = dsl_scalar(name = 'arg2')
        arg3 = dsl_scalar(name = 'arg3')
        arg4 = dsl_scalar(name = 'arg4')
        # Test src()
        binop1 = dsl_binop( op2, arg1, arg2 )
        self.assertEqual( src( binop1 ), 'arg1 + arg2' )
        binop2 = dsl_binop( op1, arg3, binop1 )
        self.assertEqual( src( binop2 ), 'arg3 * ( arg1 + arg2 )' )
        binop3 = dsl_binop( op4, binop2, arg3 )
        self.assertEqual( src( binop3 ), 'pow( arg3 * ( arg1 + arg2 ), arg3 )' )
        binop4 = dsl_binop( op3, arg1, arg2 )
        self.assertEqual( src( binop4 ), 'arg1 - arg2' )
        binop5 = dsl_binop( op2, arg3, binop4 )
        self.assertEqual( src( binop5 ), 'arg3 + arg1 - arg2' )
        binop6 = dsl_binop( op2, binop4, arg3 )
        self.assertEqual( src( binop6 ), 'arg1 - arg2 + arg3' )
        binop7 = dsl_binop( op3, arg1, binop1 )
        self.assertEqual( src( binop7 ), 'arg1 - ( arg1 + arg2 )' )
        binop8 = dsl_binop( op3, binop1, arg1 )
        self.assertEqual( src( binop8 ), 'arg1 + arg2 - arg1' )

    def test_dsl_index(self):
        i1 = dsl_index( constant = False, cartesian = False )
        i2 = dsl_index( constant = True, cartesian = False )
        # Test add_attachment(), attachment()
        var1 = dsl_scalar(name = 'var1')
        var2 = dsl_scalar(name = 'var2', constant = True)
        pos1 = dsl_scalar(name = 'pos1')
        pos2 = dsl_scalar(name = 'pos2', constant = True)
        i1.add_attachment('testvar1',var1)
        i2.add_attachment('testvar2',var2)
        i1.add_attachment('testpos1',pos1)
        i2.add_attachment('testpos2',pos2)
        self.assertIs( i1.attachment('testvar1'), var1 )
        self.assertIs( i2.attachment('testvar2'), var2 )
        self.assertIs( i1.attachment('testpos1'), pos1 )
        self.assertIs( i2.attachment('testpos2'), pos2 )
        self.assertRaises( AssertionError, i1.add_attachment, 'dummy1', var2 )
        self.assertRaises( AssertionError, i1.add_attachment, 'dummy2', pos2 )
        self.assertRaises( AssertionError, i2.add_attachment, 'dummy1', var1 )
        self.assertRaises( AssertionError, i2.add_attachment, 'dummy2', pos1 )

        # Test jump()
        d = dsl_direction(0)
        for ii in range(-100,100):
            d.set( ii % 3 )
            # Default method is jump = change
            self.assertEqual(i1.jump(change = ii),ii)
            # Direction is irrelevant
            self.assertEqual(i1.jump(change = ii), i1.jump(change = ii,direction = d) )

    def test_dsl_value(self): # inherits from dsl_base
        # Test src() (calls __str__() )
        value1 = dsl_value(name = 'value', value = 42 )
        self.assertEqual( src(value1), "42" )
        # Test assign() (NotImplemented)
        self.assertEqual(value1.assign(dsl_scalar(name = 'x')),NotImplemented)
        self.assertEqual(value1.assign(),NotImplemented)

    def test_dsl_variable(self): # inherits from dsl_base
        # Test __init__()
        v1 = dsl_variable(name = 'v1', length = 1, expr = None, vartype = 'int', prefix = 'i', cartesian = False, constant = False, component = None )
        self.assertEqual( v1.name(), 'v1' )
        self.assertEqual( v1.autonamed(), False )
        self.assertEqual( v1.type(), 'int' )
        self.assertEqual( v1.prefix(), 'i' )
        test_namer1 = namer()
        test_namer2 = namer()
        v2 = dsl_variable( name = None, length = 1, expr = v1**2.1, vartype = 'double', prefix = 'dp', cartesian = False, constant = False, component = None, namerobj = test_namer1 )
        self.assertEqual( v2.name(), test_namer2(v2) ) 
        self.assertEqual( v2.autonamed(), True )
        self.assertEqual( v2.type(), 'double')
        self.assertEqual( v2.prefix(), 'dp')
        # Test src() (calls __str__() )
        self.assertEqual( v1.name(), src( v1 ) )
        self.assertEqual( v2.name(), src( v2 ) )
    
    def test_dsl_zero(self):
        # Test src() (calls __str__() )
        self.assertEqual(src(dsl_zero()),"0")
        # Test assign() (NotImplemented)
        self.assertEqual(dsl_zero().assign(dsl_scalar(name = 'x')),NotImplemented)
        self.assertEqual(dsl_zero().assign(),NotImplemented)


    def test_dsl_scalar(self): # inherits from dsl_variable
        # Test autonaming
        for counter in range(1,10):
            i = dsl_scalar(vartype='int')
            self.assertEqual(i.name(),prefix_select(i.type())+'_'+src(counter))
        for counter in range(1,10):
            dp = dsl_scalar(vartype='double')
            self.assertEqual(dp.name(),prefix_select(dp.type())+'_'+src(counter))
        for counter in range(10,20):
            i = dsl_scalar(vartype='int')
            self.assertEqual(i.name(),prefix_select(i.type())+'_'+src(counter))

        # Test assign()
        var1 = dsl_scalar(name = 'var1')
        var2 = dsl_scalar(expr = var1**2 )
        test1 = var1.assign(expr = var2+3)
        test2 = var2.assign()
        self.assertTrue(isinstance(test1,dsl_binop))
        self.assertIs(test1.op(),op_assign)
        self.assertIs(test1.left(),var1)
        self.assertTrue(isinstance(test1.right(),dsl_binop))
        self.assertTrue(isinstance(test2,dsl_binop))
        self.assertIs(test2.op(),op_assign)
        self.assertIs(test2.left(),var2)
        self.assertTrue(isinstance(test2.right(),dsl_binop))

        
    def test_dsl_position(self): # inherits from dsl_variable
        pos1 = dsl_position(name = 'pos1')
        pos2 = dsl_position(expr = pos1+42 )
        # Test assign()
        test1 = pos1.assign(expr = pos2*21)
        test2 = pos2.assign()
        self.assertTrue(isinstance(test1,dsl_binop))
        self.assertIs(test1.op(),op_assign)
        self.assertIs(test1.left(),pos1)
        self.assertTrue(isinstance(test1.right(),dsl_binop))
        self.assertTrue(isinstance(test2,dsl_binop))
        self.assertIs(test2.op(),op_assign)
        self.assertIs(test2.left(),pos2)
        self.assertTrue(isinstance(test2.right(),dsl_binop))

        # Test __getitem__()
        self.assertTrue(isinstance(pos1[0],dsl_scalar))
        self.assertEqual(pos1[0].name(),pos1.name()+'[0]')
        self.assertTrue(isinstance(pos2[2],dsl_scalar))
        self.assertEqual(pos2[2].name(),pos2.name()+'[2]')

    def test_dsl_function(self): # inherits from dsl_base
        args = [ dsl_scalar(name = 'xa'), dsl_position(name = 'A'), dsl_scalar(name = 'la',vartype='int'),
                 dsl_scalar(name = 'na',vartype='int'),\
                 dsl_scalar(name = 'xb'), dsl_position(name = 'B'), dsl_scalar(name = 'lb',vartype='int'),
                 dsl_scalar(name = 'nb',vartype='int'),\
                 dsl_scalar(name = 'Z',vartype='double *') ]
        local_vars = [ dsl_scalar(name = 'i',vartype='int'),dsl_scalar(name = 'xi') ]
        function1 = dsl_function('function', args, local_vars, vartype = 'void' )
        # Test src() (calls __str__() )
        self.assertEqual( src(function1), 'function( xa, A, la, na, xb, B, lb, nb, Z )' )

        # NB call() and prototype() are deprecated
        # Test call()
        #self.assertEqual( function1.call(), 'function( xa, A, la, na, xb, B, lb, nb, Z )' ) 
        #self.assertEqual( function1.call(args = [ dsl_scalar(name = 'test1'), dsl_scalar(name = 'test2') ] ),\
        #                                     'function( test1, test2 )' )
        # Test prototype()
        #self.assertEqual( function1.prototype(), \
        # 'void function( double xa, double A[3], int la, int na, double xb, double B[3], int lb, int nb, double * Z )' )

        # Test assign() (NotImplemented)
        self.assertEqual(function1.assign(dsl_scalar(name = 'x')),NotImplemented)
        self.assertEqual(function1.assign(),NotImplemented)


    def test_dsl_integer_index(self):
        index1 = dsl_integer_index('index1', value=42)
        index2 = dsl_integer_index('index2')
        var1 = dsl_scalar(name = 'var1')
        var2 = dsl_scalar(name = 'var2')
        # Test __init__()
        self.assertRaises( AssertionError,dsl_integer_index, 'index', 1.0 )
        self.assertRaises( AssertionError,dsl_integer_index, 'index', dsl_scalar(name = 'x') )
        self.assertEqual(index1.value(), 42 )
        self.assertEqual(index2.value(), 0 )

        # Test add_var(), var(), var_key(), var_list()
        index1.add_var(var1,'testvar1')
        index1.add_var(var2,'testvar2')
        self.assertIs( index1.var('testvar1'), var1 )
        self.assertIs( index1.var('testvar2'), var2 )
        self.assertEqual( index1.var_key('var1'), 'testvar1' )
        self.assertEqual( index1.var_key('var2'), 'testvar2' )
        self.assertIn( var1, index1.var_list() )
        self.assertIn( var2, index1.var_list() )

        # Test increment(), decrement()
        index1.set_value(+42)
        index2.set_value(-42)
        start1 = 42
        start2 =-42
        i1 = start1
        i2 = start2
        imax = 1000
        for i in range(imax):
            self.assertEqual( i1,index1.value() )
            self.assertEqual( i2,index2.value() )
            index1.increment()
            index2.increment()
            i1 +=1
            i2 +=1
        for i in range(imax):
            self.assertEqual( i1,index1.value() )
            self.assertEqual( i2,index2.value() )
            index1.decrement()
            index2.decrement()
            i1 -=1
            i2 -=1
        self.assertEqual( start1, index1.value() )
        self.assertEqual( start2, index2.value() )

        # Test assign() (NotImplemented)
        self.assertEqual(index1.assign(dsl_scalar(name = 'x')),NotImplemented)
        self.assertEqual(index1.assign(),NotImplemented)

    def test_dsl_cartesian_gaussian(self): # inherits from dsl_base
        g1 = dsl_cartesian_gaussian('g1')
        angmom2 = [ 2, 3, 5 ]
        centre2 = dsl_position(name = 'C2')
        exponent2 = dsl_scalar(name = 'x2')
        g2 = dsl_cartesian_gaussian('g2', angmom2, centre2, exponent2 )
        # Test __init__()
        self.assertEqual( g1.angmom(), [ 0, 0, 0 ] )
        self.assertEqual( g2.angmom(), [ 2, 3, 5 ] )
        self.assertEqual( g1.cen().name(), 'cG1' )
        self.assertIs( g1.attachment('centre'), g1.cen() )
        self.assertEqual( g2.cen().name(), 'C2' )
        self.assertIs( g2.attachment('centre'), centre2 )
        self.assertEqual( g1.exp().name(), 'xg1' )
        self.assertIs( g1.attachment('exponent'), g1.exp()  )
        self.assertEqual( g2.exp().name(), 'x2' )
        self.assertIs( g2.attachment('exponent'), exponent2 )
#        self.assertEqual( g1.var('l').name(), 'lg1' )
#        self.assertEqual( g1.var('i').name(), 'ig1' )
#        self.assertEqual( g1.var('n').name(), 'ng1' )
#        self.assertEqual( g1.var('nt').name(), 'ntg1' )
#        self.assertEqual( g1.var('iskip').name(), 'iskipg1' )
#        self.assertEqual( g1.var('iskipt').name(), 'iskiptg1' )
#        self.assertEqual( g1.var('ioff').name(), 'ioffg1' )

        # Test __getitem__()
        self.assertTrue( isinstance( g1[1], dsl_integer_index ) )
        self.assertEqual( src( g1[0] ), '0' )
        self.assertEqual( src( g1[1] ), '0' )
        self.assertEqual( src( g1[2] ), '0' )
        self.assertEqual( src( g2[0] ), '2' )
        self.assertEqual( src( g2[1] ), '3' )
        self.assertEqual( src( g2[2] ), '5' )
        
        # Test add_var(), var(), var_key(), var_list()
#        var1 = dsl_scalar(name = 'var1')
#        var2 = dsl_scalar(name = 'var2')
#        g1.add_var(var1,'testvar1')
#        g1.add_var(var2,'testvar2')
#        self.assertIs( g1.var('testvar1'), var1 )
#        self.assertIs( g1.var('testvar2'), var2 )
#        self.assertEqual( g1.var_key('var1'), 'testvar1' )
#        self.assertEqual( g1.var_key('var2'), 'testvar2' )
#        self.assertIn( var1, g1.var_list() )
#        self.assertIn( var2, g1.var_list() )

        # Test add_position(), position(), position_key(), position_list()
#        pos1 = dsl_scalar(name = 'pos1')
#        pos2 = dsl_scalar(name = 'pos2')
#        g1.add_position(pos1,'testpos1')
#        g1.add_position(pos2,'testpos2')
#        self.assertIs( g1.position('testpos1'), pos1 )
#        self.assertIs( g1.position('testpos2'), pos2 )
#        self.assertEqual( g1.position_key('pos1'), 'testpos1' )
#        self.assertEqual( g1.position_key('pos2'), 'testpos2' )
#        self.assertIn( pos1, g1.position_list() )
#        self.assertIn( pos2, g1.position_list() )


        # Test increment(), decrement(), index()
        g1.set_angmom([0,0,0])
        lmax = 10
        angmom_list = []
        for l in range(0,lmax+1):
            for mx in range(l, -1, -1):
                for my in range(l-mx, -1, -1):
                    mz = l - mx - my
                    index = figurate(l,3) + figurate(my+mz,2) + figurate(mz,1) + 1
                    self.assertEqual( g1.angmom(), [mx,my,mz] )
                    self.assertEqual( g1.index(), index )
                    angmom_list.append( g1.angmom()[:] )
                    g1.increment()
            
        for angmom in reversed( angmom_list ):
            g1.decrement()
            self.assertEqual( g1.angmom(), angmom )
            self.assertEqual( g1.index(), index )
            index -= 1
        self.assertEqual( g1.angmom(), [0,0,0] )
        self.assertEqual( g1.index(), 1 ) 
                     
        # Test assign() (NotImplemented)
        self.assertEqual(g1.assign(dsl_scalar(name = 'x')),NotImplemented)
        self.assertEqual(g1.assign(),NotImplemented)

        # Test jump()
        # Overloaded from dsl_index for Cartesian indexes
        # Calls cartesian_index_change using angmom from self
        g1.set_angmom([0,0,0])
        min_change = -10
        max_change = 10
        lmax       = 6 
        d = dsl_direction('x')
        for component in ['x','y','z']:
            d.set(component)
            while g1.angmom()[2] != lmax:
                a1 = g1.angmom()
                l1 = sum( a1 )
                if a1[ d.index() ] + min_change < 0:
                    min_change0 = -a1[ d.index() ]
                else:
                    min_change0 = min_change
                for change in range( min_change0, max_change):
                    index1 = figurate(l1,3) + figurate(a1[1]+a1[2],2) + figurate(a1[2],1) + 1
                    a2 = a1[:]
                    a2[ d.index() ] += change
                    l2 = sum( a2 )
                    index2 = figurate(l2,3) + figurate(a2[1]+a2[2],2) + figurate(a2[2],1) + 1
                    self.assertEqual( index2 - index1, g1.jump(change,d) )
                g1.increment()

    def test_dsl_direction(self): # inherits from dsl_base
        # Test __init__(), set()
        d0 = dsl_direction(0)
        d1 = dsl_direction(1)
        d2 = dsl_direction(2)
        dx = dsl_direction('x')
        dy = dsl_direction('y')
        dz = dsl_direction('z')
        self.assertEqual( d0.label(), 'x' )
        self.assertEqual( d0.index(), 0 )
        self.assertEqual( d1.label(), 'y' )
        self.assertEqual( d1.index(), 1 )
        self.assertEqual( d2.label(), 'z' )
        self.assertEqual( d2.index(), 2 )
        self.assertEqual( dx.label(), 'x' )
        self.assertEqual( dx.index(), 0 )
        self.assertEqual( dy.label(), 'y' )
        self.assertEqual( dy.index(), 1 )
        self.assertEqual( dz.label(), 'z' )
        self.assertEqual( dz.index(), 2 )
        d = dsl_direction(0)
        self.assertEqual( d.index(), 0 )
        self.assertEqual( d.label(), 'x' )
        d.set(1) 
        self.assertEqual( d.index(), 1 )
        self.assertEqual( d.label(), 'y' )
        d.set(2)
        self.assertEqual( d.index(), 2 )
        self.assertEqual( d.label(), 'z' )
        d.set('x') 
        self.assertEqual( d.index(), 0 )
        self.assertEqual( d.label(), 'x' )
        d.set('y') 
        self.assertEqual( d.index(), 1 )
        self.assertEqual( d.label(), 'y' )
        d.set('z') 
        self.assertEqual( d.index(), 2 )
        self.assertEqual( d.label(), 'z' )
        var = dsl_scalar(name = 'var')
        d.set( var )
        self.assertIs( d.index() , var )
        self.assertIs( d.label() , var )
        # Check that correct exceptions are raised for nonsensical arguments
        self.assertRaises( AssertionError, d.set, 100 )
        self.assertRaises( AssertionError, d.set, 'r' )
        self.assertRaises( Exception, d.set, dsl_base() )
        


    def test_dsl_pointer(self): # inherits from dsl_base
        a1 = dsl_pointer('a',vartype='double')
        b1 = dsl_pointer('b',vartype='int')
        var = dsl_scalar(name = 'var',vartype='int')
        val = dsl_value('val', 101, vartype='int' )
        # Test vartype() and pointertype()
        self.assertEqual( a1.vartype(), 'double' )
        self.assertEqual( b1.vartype(), 'int' )
        #self.assertEqual( a1.pointertype(), 'double *' )
        #self.assertEqual( b1.pointertype(), 'int *' )

        # Test __getitem__()
        self.assertTrue( isinstance( a1[3], dsl_pointer ) )
        self.assertEqual( src( a1[3] ), 'a[ 3 ]' )
        self.assertTrue( isinstance( a1[var][2], dsl_pointer ) )
        self.assertEqual( src( a1[var][val] ), 'a[ var ][ 101 ]' )
        self.assertTrue( isinstance( b1[a1], dsl_pointer ) )
        self.assertEqual( src( b1[a1] ), 'b[ a ]' )

        # Test assign()
        test1 = a1.assign(var*val)
        test2 = a1[10].assign(val*2)
        self.assertTrue(isinstance(test1,dsl_binop))
        self.assertIs(test1.op(),op_assign)
        self.assertIs(test1.left(),a1)
        self.assertTrue(isinstance(test1.right(),dsl_binop))
        self.assertTrue( src(test1), 'a = var * 101' )
        self.assertTrue(isinstance(test2,dsl_binop))
        self.assertIs(test2.op(),op_assign)
        self.assertTrue(isinstance(test2.left(),dsl_pointer) )
        self.assertTrue(isinstance(test2.right(),dsl_binop))
        self.assertTrue( src(test2), 'a[ 10 ] = 101 * 2' )

class TestDSLExtraFunctions(unittest.TestCase):
    """Tests additional functions which support the DSL objects."""
    def binomial_factorial(self,n,k):
        """Calculate binomial coefficient using math.factorial, for testing against binomial
        coefficients generated by other means."""
        if n >= k: 
            return int( round( math.factorial(n) / ( math.factorial(k) * math.factorial(n - k) ) ) )
        else:
            return 0

    def setUp(self):
        pass

    def text_exp(self):
        """Test that exp returns a dsl_unop which is converted to the correct C source
        code for exp evaluation."""
        arg = dsl_scalar(name = 'arg')
        self.assertEqual( src( exp(arg) ),\
                src( dsl_unop( op_exp, arg ) ) )  

    def test_nint(self):
        """Test that nint returns a dsl_unop which is converted to the correct C source
        code for int( round( ) ) evaluation."""
        arg = dsl_scalar(name = 'arg')
        self.assertEqual( src( nint(arg) ),\
                src( dsl_unop( op_int, dsl_unop(op_round, arg ) ) ) )

    def test_norm(self):
        """Test that norm behaves as expected, returning dsl_binops that correctly
        evaluate to the Euclidean norm, |A| = ( A[0]**2 + A[1]**2 + ... + A[n]**2 )**0.5 in
        C source code form."""
        pos = dsl_position(name = 'pos')
        self.assertEqual( src( norm(pos) ) ,\
                src( eval( '( pos[0]**2.0 + pos[1]**2.0 + pos[2]**2.0 )**0.5' ) ) )

    def test_binomial(self):
        """Test that binomial outputs correct values for, n = 0..nmax, k = 0, n
        by comparing against evaluation of math.factorial."""
        nmax = 20
        for n in range(0,nmax):
            for k in range(0,n):
                self.assertEqual( binomial(n,k), self.binomial_factorial(n,k) )

    def test_figurate(self):
        """Test that figurate outputs correct values for, n = 0..nmax, d = 1..dmax
        by comparing against evaluation using explicit formulae for linear, triangular 
        and tetrahedral numbers, and also evaluation using math.factorial."""
        nmax = 20
        dmax = 20
        for n in range(0,nmax):
            # Linear numbers
            self.assertEqual( n,figurate(n,1) )
            # Triangular numbers
            self.assertEqual( int( round( n*(n+1)/2 ) ), figurate(n,2) )
            # Tetrahedral numbers
            self.assertEqual( int( round( n*(n+1)*(n+2)/6 ) ), figurate(n,3) )
            # General
            for d in range(1,dmax):
                self.assertEqual( self.binomial_factorial(n+d-1,d), figurate(n,d) )

    def test_cartesian_index_change(self):
        """Test values output by cartesian_index_change against an alternative set 
        of values calculated using explicit formula (Eq. 4.3, Andy May's thesis, p.90),
        as opposed to the stepping +/- 1 used in the dsl module."""
        min_change = -10
        max_change = 10
        lmax       = 6 
        d = dsl_direction('x')
        for l1 in range(0, lmax+1):
           for component in ['x','y','z']:
               d.set(component)
               for mx in range(l1, -1, -1):
                   for my in range(l1-mx, -1, -1):
                       mz = l1 - mx - my
                       a1 = [ mx, my, mz ]
                       if a1[ d.index() ] + min_change < 0:
                           min_change0 = -a1[ d.index() ]
                       else:
                           min_change0 = min_change
                       for change in range( min_change0, max_change):
                           index1 = figurate(l1,3) + figurate(a1[1]+a1[2],2) + figurate(a1[2],1) + 1
                           a2 = [ mx, my, mz ]
                           a2[ d.index() ] += change
                           l2 = sum( a2 )
                           index2 = figurate(l2,3) + figurate(a2[1]+a2[2],2) + figurate(a2[2],1) + 1
                           self.assertEqual( index2 - index1, cartesian_index_change(a1,change,d) )

    def tearDown(self):
        pass
    
class TestAllDSL(unittest.TestCase):
    """Tests which check that the DSL behaves as expected, i.e. 
    expressions are correctly parsed and represented by binop / unop objects.
    These are in fact integration tests for the entire DSL, rather than unit tests."""
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_simple_expressions(self):
        # Simple expressions with dsl_scalar, dsl_value, dsl_position
        var_a = dsl_scalar(name = 'var_a')
        var_b = dsl_scalar(name = 'var_b')
        var_c = dsl_scalar(name = 'var_c')
        var_d = dsl_scalar(name = 'var_d')
        pos_a = dsl_position(name = 'pos_a')
        pos_b = dsl_position(name = 'pos_b')
        pos_c = dsl_position(name = 'pos_c')
        pos_d = dsl_position(name = 'pos_d')
        val_a = dsl_value(name = 'val_a',value = 1,vartype='int')
        val_b = dsl_value(name = 'val_b',value = 2,vartype='int')
        val_c = dsl_value(name = 'val_c',value = 3.0,vartype='double')
        val_d = dsl_value(name = 'val_d',value = 4.0,vartype='double')
        # ...operator precedence and associativity testing for binary operations gives
        #    expected expressions as strings with expected bracket placement
        expr_list = [ \
                  ( var_a - pos_a[0] ) - pos_b[0],\
                  var_a - ( pos_a[0] - pos_b[0] ),\
                  ( var_a - pos_a[0] ) - pos_b[0],\
                  ( var_a * var_b ) + pos_a[0],\
                  var_a * ( var_b + pos_a[0] ),\
                  var_a * ( val_a * pos_a[0] ),\
                  ( var_a / val_b ) * pos_a[1],\
                  var_a / (var_c * pos_a[1]),\
                  var_a**var_b**var_c/val_d \
                    ]
        str_list = [ \
                  'var_a - pos_a[0] - pos_b[0]',\
                  'var_a - ( pos_a[0] - pos_b[0] )',\
                  'var_a - pos_a[0] - pos_b[0]',\
                  'var_a * var_b + pos_a[0]',\
                  'var_a * ( var_b + pos_a[0] )',\
                  'var_a * 1 * pos_a[0]',\
                  'var_a / 2 * pos_a[1]',\
                  'var_a / ( var_c * pos_a[1] )',\
                  'pow( var_a, pow( var_b, var_c ) ) / 4.0'\
                  ]
        for i in range( len( expr_list ) ):
            self.assertEqual( src( expr_list[i] ) , str_list[i] ) 
        # ...combinations of dsl_binop and dsl_unop objects
        expr_list = [ \
                    exp( var_a ** val_d ) / ( var_a + var_b * exp( -val_b * pos_b[0]**val_c ) ),\
                    -(val_a + val_b) + -exp( pos_a[0]**2.0 + pos_a[1]**2.0 + pos_a[2]**2.0 ) / abs( var_a ),\
                    val_a // val_b * abs( pos_d[2] * ( pos_a[0] + pos_b[0] + pos_c[0] ) ) \
                    ]
        str_list = [ \
                    'exp( pow( var_a, 4.0 ) ) / ( var_a + var_b * exp( ( -2 ) * pow( pos_b[0], 3.0 ) ) )',\
                    '( -( 1 + 2 ) ) + ( -exp( pow( pos_a[0], 2.0 ) + pow( pos_a[1], 2.0 ) + pow( pos_a[2], 2.0 ) ) ) / abs( var_a )',\
                    'floor( 1, 2 ) * abs( pos_d[2] * ( pos_a[0] + pos_b[0] + pos_c[0] ) )' \
                    ]
        for i in range( len( expr_list ) ):
            self.assertEqual( src( expr_list[i] ) , str_list[i] ) 
        # ...check that explicit specification of expressions using dsl_binop, dsl_unop objects produces
        #    expected outcome
        explicit_list = [\
                dsl_binop( op_add, var_a, dsl_unop( op_exp, dsl_binop( op_pow, var_d, val_d ) ) ),\
                dsl_binop( op_sub, var_a, dsl_binop( op_add, var_b, val_b ) ),\
                dsl_binop( op_add, dsl_binop( op_sub, var_a, var_b), val_b ),\
                dsl_binop( op_div, var_a, dsl_binop( op_div, var_b, var_c ) ),\
                dsl_binop( op_div, dsl_binop( op_div, var_a, var_b ), var_c ),\
                dsl_binop( op_mul, pos_a[0], dsl_binop( op_sub, pos_a[1], pos_a[2] ) ),\
                dsl_binop( op_sub, dsl_binop( op_mul, pos_a[0], pos_a[1]), pos_a[2] ),\
                dsl_binop( op_add, var_a, dsl_unop( op_neg, dsl_binop( op_mul, var_a, var_b ) ) ),\
                dsl_unop( op_neg, dsl_binop( op_pow, var_c, var_d ) ),\
                ]
        overload_list = [
                var_a + exp( var_d**val_d ),\
                var_a - ( var_b + val_b ),\
                var_a - var_b + val_b,\
                var_a / (var_b / var_c ),\
                ( var_a / var_b ) / var_c,\
                pos_a[0] * ( pos_a[1] - pos_a[2] ),\
                pos_a[0] * pos_a[1] - pos_a[2],\
                var_a + -( var_a * var_b ),\
                -var_c**var_d\
                ]
        str_list = [
                'var_a + exp( pow( var_d, 4.0 ) )',\
                'var_a - ( var_b + 2 )',\
                'var_a - var_b + 2',\
                'var_a / ( var_b / var_c )',\
                'var_a / var_b / var_c',\
                'pos_a[0] * ( pos_a[1] - pos_a[2] )',\
                'pos_a[0] * pos_a[1] - pos_a[2]', \
                'var_a + ( -( var_a * var_b ) )',\
                '-pow( var_c, var_d )'\
                ]
        for i in range( len( explicit_list ) ):
            self.assertEqual( src( explicit_list[i] ) , src( overload_list[i] ) ) 
            self.assertEqual( src( explicit_list[i] ) , str_list[i] ) 
        
    def test_complicated_expressions(self):
        # Expressions including dsl_function, dsl_integer_index, dsl_cartesian_gaussian, dsl_array
        # Define some dsl_scalar, dsl_position, and dsl_value objects
        var_a = dsl_scalar(name = 'var_a',vartype='int')
        var_b = dsl_scalar(name = 'var_b',vartype='double')
        pos_a = dsl_position(name = 'pos_a')
        pos_b = dsl_position(name = 'pos_b')
        val_a = dsl_value(name = 'val_a',value = 1,vartype='int')
        val_b = dsl_value(name = 'val_b',value = 2.0,vartype='double')
        # Define some dsl_function,
        function_1 = dsl_function('function_1', args = [ var_a, pos_a ],\
                                  local_vars = [ pos_b, var_b ],\
                                  vartype = 'void')
        function_2 = dsl_function('function_2',vartype='double')
        # ...dsl_integer_index,
        index_1 = dsl_integer_index('index_1',value=100)
        index_2 = dsl_integer_index('index_2')
        # ...dsl_cartesian_gaussian,
        g1 = dsl_cartesian_gaussian('g1',angmom = [ 2, 0, 0 ], centre = pos_a, exponent = var_b )
        g2 = dsl_cartesian_gaussian('g2')
        # ... and dsl_pointer objects.
        pointer_1 = dsl_pointer('pointer_1', vartype = 'int' )
        pointer_2 = dsl_pointer('pointer_2')
        # Define some dsl_scalars and dsl_positions based on others
        x = dsl_scalar(expr=var_a**2.0,vartype='double')
        y = dsl_scalar(expr=var_b + val_a, vartype='int')
        X = dsl_position(expr=pos_a + pos_b)
        Y = dsl_position(expr=pos_a * pos_b)
        # Now make some expressions:
        expr_list = [
                var_b * function_2,\
                pos_a * exp( -var_b**2.0 + var_b/10.0 ) * function_1,\
                pos_a - pos_b * var_b + val_b * pointer_1[ var_a ][ 0 ][ 0 ],\
                pointer_2[ var_a ] ** val_b + pointer_2[ var_a + 1 ] ** val_b,\
                exp( X[0]**2.0 + X[1]**2.0 + X[2]**2.0 ) * function_1,\
                2.0 * dsl_unop( op_dbl, -var_a * var_b ) * ( X - Y ),\
                X * val_b - Y * var_b + x - dsl_unop( op_dbl, y ),\
                exp( norm( Y )**2.0 )*pointer_1[ 0 ][ 0 ][ 0 ] - exp( norm( X )**2.0 )*pointer_1[0][0][0],\
                nint( var_a / var_b ) + nint( x / 1.874 ),\
                pos_a[0]*pos_b[0]*pointer_2[ var_a + val_a ]*function_1\
                ]
        str_list = [
                'var_b * function_2()',\
                'pos_a * exp( ( -pow( var_b, 2.0 ) ) + var_b / 10.0 ) * function_1( var_a, pos_a )',\
                'pos_a - pos_b * var_b + 2.0 * pointer_1[ var_a ][ 0 ][ 0 ]',\
                'pow( pointer_2[ var_a ], 2.0 ) + pow( pointer_2[ var_a + 1 ], 2.0 )',\
                'exp( pow( '+X.name()+'[0], 2.0 ) + pow( '+X.name()+'[1], 2.0 ) + pow( '+X.name()+'[2], 2.0 ) ) * function_1( var_a, pos_a )',\
                '2.0 * ( (double) ( ( -var_a ) * var_b ) ) * ( '+X.name()+' - '+Y.name()+' )',\
                ''+X.name()+' * 2.0 - '+Y.name()+' * var_b + '+x.name()+' - ( (double) '+y.name()+' )',\
                'exp( pow( pow( pow( '+Y.name()+'[0], 2.0 ) + pow( '+Y.name()+'[1], 2.0 ) + pow( '+Y.name()+'[2], 2.0 ), 0.5 ), 2.0 ) ) * pointer_1[ 0 ][ 0 ][ 0 ] - exp( pow( pow( pow( '+X.name()+'[0], 2.0 ) + pow( '+X.name()+'[1], 2.0 ) + pow( '+X.name()+'[2], 2.0 ), 0.5 ), 2.0 ) ) * pointer_1[ 0 ][ 0 ][ 0 ]',\
                 '( (int) round( var_a / var_b ) ) + ( (int) round( '+x.name()+' / 1.874 ) )',\
                 'pos_a[0] * pos_b[0] * pointer_2[ var_a + 1 ] * function_1( var_a, pos_a )'\
                ]
        for i in range( len( expr_list ) ):
            self.assertEqual( src( expr_list[i] ), str_list[i] )

    def test_typical_expressions(self):
        # Now test with some expressions using the kind of DSL objects that would
        # be used in an intception input file
        ga = dsl_cartesian_gaussian('a')
        A  = ga.cen() # name = 'cA'
        xa = ga.exp()
        gb = dsl_cartesian_gaussian('b')
        B  = gb.cen( )# name = 'cB'
        xb = gb.exp()
        gc = dsl_cartesian_gaussian('c')
        C  = gc.cen() # name = 'cC'
        xc = gc.exp()

        # (0|0), zero angular momentum 2-index overlap integral
        pi = dsl_scalar( name = 'pi', expr = math.pi )
        xp = dsl_scalar( expr = xa + xb )
        o_o_xp = dsl_scalar( expr = 1.0/ xp )
        xaxb_o_xp = dsl_scalar( expr = xa * xb * o_o_xp )
        RAB2 = dsl_scalar( expr = (A[0]-B[0])**2.0+(A[1]-B[1])**2.0+(A[2]-B[2])**2.0)
        expr = ( pi / xp )**(3/2) * exp( -xaxb_o_xp*RAB2)
        expr_str =\
              'pow( pi / '+xp.name()+', 1.5 ) * exp( ( -'+xaxb_o_xp.name()+' ) * '+RAB2.name()+' )'
        self.assertEqual( src( expr ), expr_str )
        
        # |P-C|^2, used in various expressions for zero angular momentum integrals
        RPC2 = dsl_scalar( expr = ( ( (xa * (A[0]-C[0]) + xb * ( B[0]-C[0] ))* o_o_xp )**2.0 +\
                            ( (xa * (A[1]-C[1]) + xb * ( B[1]-C[1] ))* o_o_xp )**2.0 +\
                            ( (xa * (A[2]-C[2]) + xb * ( B[2]-C[2] ))* o_o_xp )**2.0 ) )
        expr = RPC2.expr() 
        expr_str = 'pow( ( '+xa.name()+' * ( cA[0] - cC[0] ) + '+xb.name()+' * ( cB[0] - cC[0] ) ) * '+o_o_xp.name()+', 2.0 ) + '\
                   'pow( ( '+xa.name()+' * ( cA[1] - cC[1] ) + '+xb.name()+' * ( cB[1] - cC[1] ) ) * '+o_o_xp.name()+', 2.0 ) + '\
                   'pow( ( '+xa.name()+' * ( cA[2] - cC[2] ) + '+xb.name()+' * ( cB[2] - cC[2] ) ) * '+o_o_xp.name()+', 2.0 )' 
        self.assertEqual( src( expr ), expr_str )

        # (0|0|0), zero-angular momentum 3-index overlap integral
        xabc = dsl_scalar( expr = xp + xc )
        o_o_xabc = dsl_scalar( expr = 1.0 / xabc )
        
        expr = ( xp * o_o_xabc )**(3.0/2.0) *\
               ( pi / xp )**(3.0/2.0) * exp( -xaxb_o_xp*RAB2) *\
               exp( - xp * xc * o_o_xabc *RPC2 ) 
        expr_str = \
              'pow( '+xp.name()+' * '+o_o_xabc.name()+', 1.5 ) * pow( pi / '+xp.name()+', 1.5 )'\
              ' * exp( ( -'+xaxb_o_xp.name()+' ) * '+RAB2.name()+' ) * exp( ( -'+xp.name()+' ) * xc * '+o_o_xabc.name()+' * '+RPC2.name()+' )'
        self.assertEqual( src( expr ), expr_str )

        # (0|1/rC|0), zero angular momentum 2-index nuclear attraction integral
        # ... with mock class representing boys_f (see boys.py for real class)
        class boys_f(dsl_base):
            def __init__(self,m,x,varname='boys_f'):
                self._m = m # should be dsl_integer_index
                self._x = x # should be dsl_scalar, or dsl_binop
                self._varname = varname

            def __str__(self):
                out = []
                out.append(self._varname)
                out.append('[ ')
                out.append( src(self._m) )
                out.append(' ]')
                return ''.join(out)
        m = dsl_integer_index(name = 'm',value=34)
        x = dsl_scalar(name = 'x',vartype='double')
        expr = 2.0 * pi / xp * exp( - xaxb_o_xp * RAB2 ) * boys_f( m, x )
        expr_str = '2.0 * pi / '+xp.name()+' * exp( ( -'+xaxb_o_xp.name()+' ) * '+RAB2.name()+' ) * boys_f[ 34 ]'
        self.assertEqual( src( expr ), expr_str )

