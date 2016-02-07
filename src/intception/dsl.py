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
import sys
import copy
import numbers
from intception import classtools
from intception.namer import namer
from intception.printer import src

class dsl_op:
    """Object representing operators on DSL objects."""
    def __init__(self,pyname,outname,syntax,prec = 0,assoc = None,optype = None):
        """Initialize dsl_op object.
        - pyname is the string used to represent the operator in Python
        - outname is the string which the operator is converted to in output
        - syntax is a string defining the syntax of the operator, e.g.
          'infix', 'prefix', 'function' etc.
        - prec is a number defining operator precedence -- lower numbers have
          higher precedence when the expression is parsed
        - assoc defines the associativity of the operator, e.g.
          'left', 'right' or 'full'
        - optype is not implemented.
        """
        self._pyname  = pyname
        self._outname = outname
        self._syntax  = syntax
        # prec only need be set if syntax != 'function' since precedence is
        # always unambiguous in function notation.
        # assoc only need to be set if syntax == 'infix'
        # since function syntax is always unambiguous regarding operation
        # and unary operations do not have this concept, and we only allow
        # postfix and prefix syntax to be used for unary operations
        self._prec    = prec # operator precedence level, 0 if function
        self._assoc   = assoc # set 'left' if left associative, 'right' if
                              # right associative, anything else (or leave as
                              # None) for fully associative
        self._optype  = optype # Not implemented yet, but type-checking for
                               # a DSL operator may be needed in future

    def __str__(self):
        if self._syntax == 'infix' :
            return ' '+self._outname+' '
        else :
            return self._outname
    
    def __eq__(self,other):
        """Equality test checks that the namespaces of two instances are the same"""
        if isinstance(other,dsl_op):
            return self.__dict__ == other.__dict__
        else:
            return NotImplemented
            
    def syntax(self):
        return self._syntax

    def prec(self):
        return self._prec

    def assoc(self):
        return self._assoc

    def optype(self):
        return self._optype

    def pyname(self):
        return self._pyname

#op_dict = {
#            # Arithmetic operators
#            '+' : [ '+', 'infix' ],            # addition
#            '-' : [ '-', 'infix' ],            # subtraction 
#            '*' : [ '*', 'infix' ],            # multiplication
#            '/' : [ '/', 'infix' ],            # division
#            '**': [ 'pow', 'function' ],       # exponentiation
#            '%' : [ 'mod', 'function' ],       # modulo
#            '//' : [ 'floor', 'function' ],    # floor division
#            # Boolean operators
#            'not'  : [ '!', 'prefix' ],        # not
#            'and'  : [ '&&', 'infix' ],        # and 
#            'or'   : [ '||', 'infix' ],        # or
#            '=='   : [ '==', 'infix' ],        # equal to
#            '!='   : [ '!=', 'infix' ],        # not equal to
#            '<'    : [ '<', 'infix' ],         # less than
#            '<='   : [ '<=', 'infix' ],        # less than or equal to
#            '>'    : [ '>', 'infix' ],         # greater than
#            '>='   : [ '>=', 'infix' ]         # greater than or equal to
#    }

# C operator precedence:
# http://en.cppreference.com/w/cpp/language/operator_precedence
# 3     type casting
# 5     * / % 
# 6     + -

op_int = dsl_op('int','(int) ','prefix',prec=3) # NB. (int) in C truncates the double, not rounds
op_dbl = dsl_op('double','(double) ','prefix',prec=3)

op_mul = dsl_op('*','*','infix',prec=5,assoc='full')
op_div = dsl_op('/','/','infix',prec=5,assoc='left')
op_mod = dsl_op('%','%','infix',prec=5,assoc='left')
op_add = dsl_op('+','+','infix',prec=6,assoc='full')
op_sub = dsl_op('-','-','infix',prec=6,assoc='left')
op_floor = dsl_op('//','floor','function')
op_pow = dsl_op('**','pow','function')
op_neg = dsl_op('u-','-','prefix')
op_pos = dsl_op('u+','+','prefix')
op_abs = dsl_op('abs','abs','function')
op_exp = dsl_op('exp','exp','function')
op_max = dsl_op('max','max','function')
op_min = dsl_op('min','min','function')
op_round = dsl_op('round','round','function')
op_sqrt = dsl_op('sqrt','sqrt','function')

# Not overloaded in dsl_base
op_not = dsl_op('not ','!','prefix')
op_and = dsl_op(' and ','&&','infix',prec=13,assoc='left')
op_or  = dsl_op(' or ','||','infix',prec=14,assoc='left')

op_eq = dsl_op('==','==','infix',prec=9,assoc='left')
op_ne = dsl_op('!=','!=','infix',prec=9,assoc='left')
op_lt = dsl_op('<','<','infix',prec=8,assoc='left')
op_le = dsl_op('<=','<=','infix',prec=8,assoc='left')
op_gt = dsl_op('>','>','infix',prec=8,assoc='left')
op_ge = dsl_op('>=','>=','infix',prec=8,assoc='left')

# Special operators for (de)referencing pointers in C
op_deref = dsl_op(None,'*','prefix') # dereference a pointer
op_ref   = dsl_op(None,'&','prefix') # obtain address of (reference) value

# Assignment operator is special. It cannot be overloaded in Python. 
# It should only ever be used through an assign( ) method of a class
# and will evaluate directly to a string output, rather than a binop.
# No precedence or associativity need be specified because of this.
op_assign = dsl_op('=','=','infix',prec=10)

op_dict = {
            ### BINARY OPERATIONS ###
            # Assignment
            '=' : op_assign,   # assignment (special, see above)
            # Arithmetic operators
            '+' : op_add,      # addition
            '-' : op_sub,      # subtraction 
            '*' : op_mul,      # multiplication
            '/' : op_div,      # division
            '**': op_pow,      # exponentiation
            '%' : op_mod,      # modulo
            '//' : op_floor,   # floor division (cmath)
            # Boolean operators
            'not'  : [ '!', 'prefix' ],        # not
            'and'  : [ '&&', 'infix' ],        # and 
            'or'   : [ '||', 'infix' ],        # or
            '=='   : [ '==', 'infix' ],        # equal to
            '!='   : [ '!=', 'infix' ],        # not equal to
            '<'    : [ '<', 'infix' ],         # less than
            '<='   : [ '<=', 'infix' ],        # less than or equal to
            '>'    : [ '>', 'infix' ],         # greater than
            '>='   : [ '>=', 'infix' ],        # greater than or equal to
            ### UNARY OPERATIONS ###
            'u-' : op_neg,     # unary negation
            'u+' : op_pos,     # unary positive
            'abs': op_abs,     # absolute value (cmath)
            'exp': op_exp,     # exponential
            'max': op_max,     # maximum
            'min': op_min,     # minimum
            'sqrt': op_sqrt    # take a square root
    }

# Prefix dictionary (for variable naming), related C types to variable short 
# name prefixes

prefix_dict = { 
        'int'       : 'i',
        'int *'     : 'iP',
        'double'    : 'dp',
        'double *'  : 'dpP',
        'void'      : 'v',
        'void *'    : 'vP'
}


def prefix_select(vartype):
    assert isinstance(vartype,str), 'vartype must be a string'
    try:
        return prefix_dict[vartype]
    except KeyError:
        raise Exception('vartype '+str(vartype)+' not recognized')

# Instance of namer that will be used in naming variables
default_namer = namer()

# Functions and operators that are not default Python class members and
# cannot be overloaded in the normal way
def exp(arg):
    """Returns a dsl_unop object which corresponds to C source code:
        exp(arg)
    """
    return dsl_unop(op_exp,arg)

def sqrt(arg):
    """Returns a dsl_unop object which corresponds to C source code:
        sqrt(arg)
    """
    return dsl_unop(op_sqrt,arg)

def norm(vector):
    """Euclidian norm of an array. 
    Returns a dsl_binop object that corrsponds to C source code:
        pow( pow( vector[0], 2.0 ) + pow(vector[1], 2.0) + pow(vector[2], 2.0), 0.5 )
    """
    if not isinstance(vector,dsl_position):
        raise Exception('Only dsl_position objects may be arguments')
    squares = None
    for i in range(len(vector)):
        if squares == None:
            squares = dsl_binop(op_pow,vector[i],2.0)
        else:
            squares = dsl_binop(op_add,squares,dsl_binop(op_pow,vector[i],2.0) )
    return dsl_binop( op_pow, squares, 0.5 )

def nint(arg): 
    """Equivalent to Fortra nint(real arg), i.e. rounds double to nearest 
    integer. Rounds away from zero for halfway cases. 
    Returns a dsl_unop object which corresponds to C source code:
        (int) round(arg)
    """
    return dsl_unop(op_int, dsl_unop(op_round, arg) )

# Useful mathematical functions for indexing of Cartesian Gaussians
def binomial(n,k):
    """Returns an integer binomial coefficient ( n, k ), using an iterative algorithm."""
    if n < k:
        return 0
    elif n == k:
        return 1
    elif k == 0:
        return 1
    else:
        # Fast iterative algorithm for determining binomial coefficients
        # See: http://blog.plover.com/math/choose.html
        r = 1
        for d in range(1, k+1):
            r = r * n
            n = n - 1
            r = int( round( r / d ) )
        return r

def figurate(n,d):
    """Return the integer figurate number f^n_d (see Andy May's thesis, p.89),
    using a fast iterative algorithm for determining binomial coefficients."""
    return binomial(n+d-1,d)

class NegativeAngmomError(Exception):
     """Custom exception class for cases where the change in a Cartesian
     index would result in a negative (nonsense) negative angular momentum,
     and thus negative array reference in an array of Cartesian components"""
     def __init__(self,message,value):
         self._message = message
         self._value = value

     def message(self):
         return self._message 

     def value(self):
         return self._value

def cartesian_index_change(start_angmom,change,direction):
    """Return the index change for a 3D Cartesian angular momentum vector in an
    array where the ordering is based on the figurate numbers (see Andy May's 
    thesis, p.89-90 for details of indexing). The change in angular momentum
    is defined by the non-zero integer 'change' and dsl_direction instance
    'direction'."""

    def index_step(angmom,change,direction_index):
        l = sum( angmom )
        if change > 0:
            inc = figurate(l+1,2)
            if direction_index == 1: # y
                inc += angmom[1] + angmom[2] + 1
            elif direction_index == 2: # z
                inc += angmom[1] + angmom[2] + 2
            return inc
        elif change < 0:
            dec = -figurate(l,2)
            if direction_index == 1: # y
                dec -= angmom[1] + angmom[2]
            elif direction_index == 2: # z
                dec -= angmom[1] + angmom[2] + 1
            return dec
        else:
            raise Exception('change should be a non-zero integer')
    if change > 0:
        angmom_step = 1
    elif change < 0:
        angmom_step = -1
    if start_angmom[direction.index()] + change < 0:
        raise NegativeAngmomError('start_angmom[direction] + change < 0',\
                start_angmom[direction.index()] + change)
        #pass
    index_change = 0
    angmom = start_angmom[:]
    for i in range( abs(change) ):
        index_change = index_change +\
                index_step(angmom,change,direction.index() )
        angmom[direction.index()] = angmom[direction.index()] + angmom_step
    return index_change

    def index_cartesian(angmom):
        """Provides the integer index of an angular momentum vector in a list
        of angular momentum vectors ordered in the standard way"""
        l = sum(angmom)
        j = angmom[1]
        k = angmom[2]
        return int( (l+2)*(l+1)*l/6 + (j+k+1)*(j+k)/2 + k + 1 )

# Core DSL objects
class dsl_base(classtools.classtools):
    """Base object defining operator overloads as dsl_binop or dsl_unop objects.
    This is intended to be a parent class, and instances of it should not be created.

    Inherits from classtools.classtools so that useful generic methods can be used
    for all dsl_* objects which inherit from dsl_base."""
    # Binary operations built into Python
    def __add__(self,other):
        return dsl_binop(op_dict['+'], self, other ) 
    def __radd__(self,other):
        return dsl_binop(op_dict['+'], other, self ) 
    def __sub__(self,other):
        return dsl_binop(op_dict['-'], self, other ) 
    def __rsub__(self,other):
        return dsl_binop(op_dict['-'], other, self ) 
    def __mul__(self,other):
        return dsl_binop(op_dict['*'], self, other ) 
    def __rmul__(self,other):
        return dsl_binop(op_dict['*'], other, self ) 
    def __truediv__(self,other):
        return dsl_binop(op_dict['/'], self, other ) 
    def __rtruediv__(self,other):
        return dsl_binop(op_dict['/'], other, self ) 
    def __floordiv__(self,other):
        return dsl_binop(op_dict['//'], self, other ) 
    def __rfloordiv__(self,other):
        return dsl_binop(op_dict['//'], other, self ) 
    def __pow__(self,other):
        return dsl_binop(op_dict['**'], self, other ) 
    def __rpow__(self,other):
        return dsl_binop(op_dict['**'], other, self ) 
    # Non-overloaded binary operations
    def assign(self,expr):
        return dsl_binop(op_assign,self,expr)
    # Unary operations built into Python
    def __neg__(self):
        return dsl_unop(op_dict['u-'], self )
    def __pos__(self):
        return dsl_unop(op_dict['u+'], self )
    def __abs__(self):
        return dsl_unop(op_dict['abs'], self )

class dsl_unop(dsl_base):
    """Class representing a unary operation in a DSL expression.
    The unary operation has an associated operator, op, which should a dsl_op object.
    There is one argument, which should be an object of one of the dsl_* classes which
    inherit from dsl_base.
    Examples of unary operations include negation, abs(), exp().
    dsl_unop objects inherit the overloaded operators from dsl_base.
    The overloaded operators in dsl_base which correspond to a unary operation return
    a dsl_unop object. For example, when an object that inherits from dsl_base is
    negated, (e.g. -dsl_scalar('x')) the return value for that expression is
        dsl_unop(op_neg,dsl_scalar('x'))
    """
    def __init__(self, op, arg):
        self._op  = op
        self._arg =  arg
    def __str__(self):
        return self.src()
    def src(self):
        """
        Returns a string representing the operation for source code output.
        The format of the string is dependent on the syntax of the dsl_op.
        
        Note:
        dsl_binop and dsl_unop are the only dsl_* classes with src() methods.
        These are necessary because generator classes may manipulate trees of
        binary and unary operations, whilst leaving the dsl_binop and dsl_unop
        objects intact. In this case, when src() is called on the expression, 
        the dsl_binop and dsl_unop objects must call src() on their arguments.
        """
        brackets = False
        out      = []
        if isinstance(self.arg(),dsl_binop) :
            if self.arg().op().syntax() == 'infix' :
                brackets = True
        elif isinstance(self.arg(),dsl_unop) :
            if self.arg().op().syntax() == 'prefix' or\
               self.arg().op().syntax() == 'suffix' :
                brackets = True
        if self.op().syntax() == 'prefix' :
            out.append( src( self.op() ) )
            if brackets : out.append( '( ' )
            out.append( src( self.arg() ) )
            if brackets : out.append( ' )' )
        elif self.op().syntax() == 'postfix' :
            if brackets : out.append( '( ' )
            out.append( src( self.arg() ) )
            if brackets : out.append( ' )' )
            out.append( src( self.op() ) )
        elif self.op().syntax() == 'function' :
            out.append( src( self.op() ) )
            out.append( '( ' )
            out.append( src( self.arg() ) )
            out.append( ' )' )
        else:
            raise Exception('Operator syntax '+self.op().syntax()
                    +' not understood for unary operator')
        return ''.join(out)
    def op(self):
        return self._op
    def arg(self):
        return self._arg

class dsl_binop(dsl_base):
    """Class representing a binary operation in a DSL expression.
    The binary operation has an associated operator, op, which should a dsl_op object.
    There are two arguments, which should be objects of one of the dsl_* classes which
    inherit from dsl_base.
    Examples of unary operations include addition, subtraction, multiplication, exponentiation.
    dsl_binop objects inherit the overloaded operators from dsl_base.
    The overloaded operators in dsl_base which correspond to a binary operation return
    a dsl_binop object. For example, when an objects that inherit from dsl_base are
    added, (e.g. dsl_scalar('x') + dsl_scalar('y') ) the return value for that expression is
        dsl_binop( op_add, dsl_scalar('x'), dsl_scalar('y') )
    """
    def __init__(self, op, left, right):
        self._op    = op
        self._left  =  left
        self._right = right
    def __str__(self):
        return self.src()
    def src(self):
        """Returns a string representing the operation for source code output.
        The format of the string is dependent on the syntax of the dsl_op,
        the precedence and associativity of the op associated with this dsl_binop, 
        as well as the precedence and associativity of dsl_op objects associated with
        the left or right arguments (if these are dsl_binop objects). Whether brackets
        are added is determined by operator precedence and associativity as defined
        for the dsl_op objects.

        Note:
        dsl_binop and dsl_unop are the only dsl_* classes with src() methods.
        These are necessary because generator classes may manipulate trees of
        binary and unary operations, whilst leaving the dsl_binop and dsl_unop
        objects intact. In this case, when src() is called on the expression, 
        the dsl_binop and dsl_unop objects must call src() on their arguments.
        """

        left_brackets  = False
        right_brackets = False
        out = []
        # Only include brackets where absolutely necessary:
        if isinstance(self.left(),dsl_binop) and\
                self.left().op().syntax() == 'infix':
           if self.left().op().prec() > self.op().prec() :
               left_brackets = True
           if self.left().op().prec() == self.op().prec() and\
                   self.op().assoc() == 'right' :
               left_brackets = True
        if isinstance(self.right(),dsl_binop) and\
                self._right.op().syntax() == 'infix':
           if self._right.op().prec() > self.op().prec() :
                right_brackets = True
           if self._right.op().prec() == self.op().prec() and\
                   self.op().assoc() == 'left' :
                right_brackets = True
        if isinstance(self.left(),dsl_unop) and\
                ( self._left.op().syntax() == 'prefix' or\
                  self._left.op().syntax() == 'postfix' ) :
                    left_brackets = True
        if isinstance(self.right(),dsl_unop) and\
                ( self._right.op().syntax() == 'prefix' or\
                  self._right.op().syntax() == 'postfix' ) :
                    right_brackets = True
        if isinstance(self.left(),dsl_unop) and\
                ( self._left.op().syntax() == 'prefix' or\
                  self._left.op().syntax() == 'postfix' ) :
                    left_brackets = True
        if isinstance(self.right(),dsl_unop) and\
                ( self._right.op().syntax() == 'prefix' or\
                  self._right.op().syntax() == 'postfix' ) :
                    right_brackets = True
        if isinstance(self.left(),numbers.Real):
            if self.left() < 0:
                left_brackets = True # unary negative part of number
        if isinstance(self.right(),numbers.Real):
            if self.right() < 0:
                right_brackets = True # unary negative part of number
        # Include brackets always
#        if isinstance(self._left,dsl_binop) : left_brackets = True
#        if isinstance(self._right,dsl_binop): right_brackets = True
        if self.op().syntax() == 'infix' :
            if left_brackets : out.append('( ')
            out.append( src( self.left() ) )
            if left_brackets : out.append(' )')
            out.append( src( self.op() ) )
            if right_brackets : out.append('( ')
            out.append( src( self.right() ) )
            if right_brackets : out.append(' )')
        elif self.op().syntax() == 'function' :
            out.append( src( self.op() ) )
            out.append( '( ' )
            out.append( src( self.left() ) )
            out.append( ', ' )
            out.append( src( self.right() ) ) 
            out.append( ' )' )
        else :
            raise Exception('Operator syntax '+self.op().syntax()
                    +' not understood for binary operator')
        return ''.join(out)
    def op(self):
        return self._op
    def left(self):
        return self._left
    def right(self):
        return self._right

class dsl_index(dsl_base):
    """
    Class that represents objects which can behave as integral indexes.
    Any object inheriting from dsl_index should call
        dsl_index.__init__(self) 
    in it's own __init__() method to ensure correct initialization.
    """
    def __init__(self, cartesian, constant):

        self._attachment_dict = {}
        self._attachment_list = []
        self._is_cartesian = cartesian # a Cartesian index will have a src() method that
                                       # depends on Cartesian direction
        # The boolean argument constant is used to determine, whether the attachments will be created 
        # as constants.Note that within Intception, "constant" means that the variable should not 
        # change within the
        self._constant_attachments = constant

    def add_attachment(self,key,value):
        """
        Attaches a DSL object to the dsl_cartesian_gaussian instance. 
        key   : a string that generically refers to the object (value), which will be common
                among other dsl_cartesian_gaussian objects with a corresponding object
                attached.
        value : the attached object, which is typically related to the index.
        By genercially referring to an attached object, these can be permuted along with
        indexes in DSL expressions, e.g. when swapping dsl_cartesian_gaussians ga and gb,
        and occurences of xa, xb can be swapped, and any derived variables containing these
        can be replaced with variables where these indices have been swapped.

        Note that add_attachment checks whether an attachment has a is_constant attribute, and
        if so, checks that the value of this agrees with self.constant_attachments().
        """
        assert isinstance( key, str ), 'key must be a string'
        try:
            assert value.is_constant() == self.constant_attachments(), \
                    "The is_constant attribute of value should agree with the "+\
                    "constant_attachments attribue of the dsl_index object."
        # Expect attribute errors for vaues without an is_constant attribute
        except AttributeError: pass
        self._attachment_dict[ key ] = value
        self._attachment_list.append( value )

    def attachment(self,key):
        """Returns an object that has been attached."""
        return self._attachment_dict[key]

    def attachment_dict(self):
        """Returns the dictionary of attachments."""
        return self._attachment_dict
        
    def attachment_list(self):
        """Returns an unordered list of attachments."""
        return self._attachment_list

    def jump(self,change,direction=None):
        """
        Calculates the change in array index (in a 1D array) for a given 
        increment/decrement (change) and Cartesian direction (if applicable) 
        of increment/decrement.

        For an integer index, direction is irrelevant (this is the default method
        of the parent class dsl_index).
        """
        return change

    def is_cartesian(self):
        return self._is_cartesian

    def constant_attachments(self):
        return self._constant_attachments

class dsl_value(dsl_base):
    """Class representing literal values in output source code.
    dsl_value objects inherit the overloaded operators from dsl_base and are
    when converted to a string (during source code output) are represented
    by their value, rather than name."""
    def __init__(self, name, value = None, vartype = 'double'):
        # Behaves like a dsl_scalar, but src() returns the value of the 
        # variable. Use this for literal constants that you want to manipulate
        # in expressions like other DSL objects.
        # Recommend storing self._value as a string if it is a double precision
        # value to avoid any floating point error.
        self._name   = name
        self._value  = value
        self._vartype= vartype
    def __str__(self):
        """
        Returns a string showing only the value of the object.
        """
        return str(self._value)
    def set_value(self,expr):
        self._value = value
    def value(self):
        return self._value
    def type(self):
        return self._vartype
    def name(self):
        return self._name
    def assign(self,expr=None):
        """Override method in parent class."""
        return NotImplemented 

class dsl_pointer(dsl_base):
    """
    Class representing a pointer in the output source code.
    The objects carry a name and type, but no other information about the data
    pointed to. This allows flexibility when dealing with variable length and
    dynamically allocated arrays.
    """
    def __init__(self, name, vartype = 'int', constant = False):
        """name is the string representing the name of the object in source code.
        vartype is the type of variable in output source code (default int).
        constant is a boolean which indicates whether the object should be considered
        a constant (note that at present, this means that the data pointed to by a 
        pointer is constant, not the pointer itself)."""
        self._name   = name
        self._vartype= vartype
        # note that 
        self._is_constant = constant
    def __str__(self):
        return self._name
    def vartype(self):
        return self._vartype
    def type(self):
        # type() defaults to vartype
        return self.vartype()
    def name(self):
        return self._name
    def __len__(self):
        # dsl_pointer objects are always of length 1
        return 1
    def __getitem__(self,key):
        return dsl_pointer( src( self.name() )+'[ '+src(key)+' ]' ) 
    def assign(self,expr = None):
        """Returns a dsl_binop corresponding to an assignment operation of the 
        object to some valid DSL expression."""
        assert expr != None, 'a valid expression must be provided for assigment'
        return dsl_base.assign(self,expr)
    def is_constant(self):
        return self._is_constant

class dsl_variable(dsl_base):
    """
    Class that represents objects that act as variables in source code.

    Any object inheriting from dsl_index should call
        dsl_variable.__init__(self ...) 
    in it's own __init__() method to ensure correct initialization.

    If name == None, then an instance of the namer class should be called
    to automatically name the object. The default namer is default_namer
    which is instantiated in this file.
    """
    def __init__(self, name, length, expr, vartype, prefix, cartesian, component, constant, namerobj = default_namer):
        self._length = length
        self._vartype= vartype
        self._prefix = prefix
        self._expr = expr
        if name == None:
            self._name = namerobj(self)
            self._autonamed = True
        else:
            self._name = name
            self._autonamed = False
        self._is_cartesian = cartesian
        # for automatic dependency detection, vector components need to refer back to the
        # dsl_position which they are a component of
        self._component_of = component
        # To allow for constant variable types, constant should be set to True or False
        # In the context of Intception, a constant is a variable which should not be modified
        # inside of the Intception generated code -- it may be modified or set outside of this.
        self._is_constant     = constant
    def __str__(self):
        return self._name
    def __len__(self):
        return self._length
    def set_expr(self,expr):
        self._expr = expr
    def expr(self):
        return self._expr
    def type(self):
        return self._vartype
    def name(self):
        return self._name
    def prefix(self):
        return self._prefix
    def autonamed(self):
        return self._autonamed
    def is_cartesian(self):
        return self._is_cartesian
    def component_of(self):
        return self._component_of
    def is_constant(self):
        return self._is_constant

class dsl_zero(dsl_base):
    """Class representing the number zero. Instances of this class are useful in DSL expressions
    for indicating the presence of zero, which can be treated separately to othe objects and used
    to simplify expressions, e.g. we know that 
        dsl_binop(op_mul, dsl_zero(), dsl_scalar() ) 
    can be removed from an expression, since it represents zero multiplied by a variable."""
    def __init__(self):
        pass

    def __str__(self):
        return str(0)

    def assign(self,expr=None):
        """Override method in parent class."""
        return NotImplemented 

class dsl_function(dsl_base):
    """Class representing function (callable units) in output source code.
    These carry with them a list of arguments and local variables, but no
    information about the instructions contained in the callable unit.
    A string containing a call to the function (based on the items in self._args)
    is generated when the src() method is called.
    """
    def __init__(self, name, args = [], local_vars = [],vartype = 'double' ):
        self._name = name
        self._args = args # list of arguments
        self._local_vars = local_vars # list of local variables
        self._vartype  = vartype
    def __str__(self):
        out = []
        out.append( self._name )
        args = self._args
        if len(args) > 0:
            out.append( '( ' )
            for arg in args:
                out.append( src( arg ) )
                out.append( ', ' )
            del out[-1]
            out.append( ' )' )
        else:
            out.append( '()' )
        return ''.join( out )
    def type(self):
        return self._vartype
    def name(self):
        return self._name
    def args(self):
        return self._args
    def set_args(self,args):
        self._args = args
    def local_vars(self):
        return self._local_vars
    def set_local_vars(self,local_vars):
        self._local_vars = local_vars
    def assign(self,expr=None):
        """Override method in parent class."""
        return NotImplemented 

class dsl_scalar(dsl_variable):
    """Class representing scalar variables in output source code. For example, a dsl_scalar object might
    represent a single double precision number, or single integer. 
    dsl_scalar objects inherit the overloaded operators from dsl_base."""
    def __init__(self, name = None, expr = None, vartype = 'double', component_of = None, constant = False):
        prefix = prefix_select(vartype)
        dsl_variable.__init__(self,name,1,expr,vartype,prefix,cartesian = False,component = component_of, constant = constant)
        self._autoassign = True # when outputting source code, variable
        # will be assigned to self._expr immediately after local variable
        # declarations. Set to False if variable depends on other assignments
        # or will be assigned later in the code.
    def autoassign(self):
        return self._autoassign
    def component_of(self):
        return self._component_of
    def assign(self,expr = None):
        """Returns a dsl_binop corresponding to an assignment operation of the 
        object to some valid DSL expression.
        If no expr is given, then an expr set previously is used."""
        if expr == None:
            assert self._expr != None, 'no expr given and none previously set.'
            expr = self._expr
        return dsl_base.assign(self,expr)
    def set_autoassign(self,autoassign):
        self._autoassign = autoassign # True or False

class dsl_position(dsl_variable):
    """Class representing double precision 3-element arrays in output source code.
    These represent 3D Cartesian vector positions.
    dsl_position objects inherit the overloaded operators from dsl_base."""
    def __init__(self, name = None, expr = None, vartype='double', constant = False):
        prefix = prefix_select(vartype)+'3'
        dsl_variable.__init__(self,name,3,expr,vartype,prefix,cartesian = True,component = None, constant = constant )
        self._autoassign = True # when outputting source code, variable
        # will be assigned to self._expr immediately after local variable
        # declarations. Set to False if variable depends on other assignments
        # or will be assigned later in the code.
    def __str__(self):
        return self._name
    def set_expr(self,expr):
        self._expr = expr
    def expr(self):
        return self._expr
    def type(self):
        return self._vartype
    def name(self):
        return self._name
    def __getitem__(self,key):
        return dsl_scalar( src( self.name() )+'['+src(key)+']', component_of = self) 
    def __len__(self):
        return self._length
    def assign(self, expr = None):
        """Returns a dsl_binop corresponding to an assignment operation of the 
        object to some valid DSL expression.
        If no expr is given, then an expr set previously is used."""
        if expr == None:
            assert self._expr != None, 'no expr given and none previously set.'
            expr = self._expr
        return dsl_base.assign(self,expr)
    def set_autoassign(self,autoassign):
        self._autoassign = autoassign # True or False
    def autoassign(self):
        return self._autoassign

# Derived DSL objects
class dsl_integer_index(dsl_index):
    """Class representing an integer index in output source code. For example, a dsl_integer_index object might
    represent an auxiliary index which can be incremented or decremented in units of 1.
    src() returns the value of the index, rather than the name (as in dsl_scalar).
    dsl_integer_index objects inherit the overloaded operators from dsl_base."""
    def __init__(self,name,value = None, constant = False):
        """name is the string name associated with the index.
        value is the initial integer value associated with the index (if none is set, this is 0).
        A dictionary of variables associated with the dsl_integer_index object is also initialized.
        This is useful for attaching variables relevant to the index, for example
        a dsl_scalar called 'i' + self._name is attached as a counter in an array related to this
        variable and a dsl_scalar called 'iskipt' + self._name is attached, and refers to the 
        skip in array elements for each increment in the dsl_integer_index's value.
        The boolean argument constant is used to determine, whether the attachments should be 
        set to constants (i.e. is_constant == True) or not.
        Note that within Intception, "constant" means that the variable should not change within the
        context of the generated code -- it may be modified or set externally."""
        dsl_index.__init__(self, cartesian = False, constant = constant )
        if value != None:
            assert isinstance(value,int), 'value must be an integer value'
        self._name   = name
        self._vars   = {}
        self._var_keys = {}
        if value == None:
            self._value = 0
        else :
            self._value = value # some integer value
    def __str__(self):
        """Return a string containing the value of the index."""
        return str( self._value )
    def add_var(self,var,key):
        """Add a variable to the dictionary of variables. Also add the key to dictionary
        of variable keys. Referring to a variable by a key allows variables to be generically
        named for multiple dsl_integer_index objects and accessed as such, e.g.
        a counter dsl_scalar('i'+self._name,vartype='int') can be stored with the key 'i', then
        for multiple dsl_integer_index objects, a counter of this type can be created and referred to
        generically by 'i', rather than 'ia', 'ib', etc."""
        self._vars[key] = var
        self._var_keys[var.name()] = key
    def var(self,key):
        """Returns the variable associated with key."""
        return self._vars[key]
    def var_key(self,name):
        """Returns the key associated with a variable name."""
        return self._var_keys[ name ] 
    def var_list(self):
        """Returns an (unordered) list of the variables stored in the self._vars dictionary."""
        return list( self._vars.values() )
    def name(self):
        return self._name
    def value(self):
        return self._value
    def index(self):
        return self._value
    def set_value(self,value):
        assert isinstance(value,int), 'value must be an integer value'
        self._value = value
    def increment(self):
        self._value += 1
    def decrement(self):
        self._value -= 1
    def assign(self,expr=None):
        """Override method in parent class."""
        return NotImplemented 


class dsl_cartesian_gaussian(dsl_index):
    """Class representing an Cartesian Gaussian index in output source code. 
    For example, a dsl_cartesian_gaussian object might
    represent a primitive Gaussian index 
    dsl_cartesian_gaussian objects have an associated angular momentum [ mx, my, mz ] which 
    can be incremented or decremented in units of angular momentum ([1,0,0], [0,1,0], [0,0,1])
    and can be indexed in a vector of Cartesian components.
    src() returns the name of the index, as in dsl_scalar.
    dsl_cartesian_gaussian objects inherit the overloaded operators from dsl_base."""
    def __init__(self,name, angmom = None, centre = None, exponent = None, constant = False):
        """name is the string name associated with the Cartesian Gaussian.
        angmom is a 3-element list of integers specifying the angular momentum in the x, y, z
        directions, e.g. [ 1, 0, 0 ].
        Dictionaries of variables and positions associated with the dsl_cartesian_gaussian object 
        are also initialized.
        These are particularly useful when Cartesian indices are permuted in expressions,
        since any positions or variables that are associated with the indices can be identified
        from these dictionaries and permuted accordingly.
        Some positions and variables are inserted into these dictionaries by default, including
        the centre, exponent, and some dsl_scalar objects useful in the generated code, such as
        associated counters, lengths, array skips and total angular momentum.
        The boolean argument constant is used to determine, whether the default variable attachments 
        (exponent, centre) will be created as constants. This does not affect the ability of the
        dsl_cartesian_gaussian object itself to have different angular momentum values associated 
        with it.
        Note that within Intception, "constant" means that the variable should not change within the
        context of the generated code -- it may be modified or set externally."""
        # Initialize variables from parent class
        dsl_index.__init__(self, cartesian = True, constant = constant) 
        self._is_cartesian = True
        self._name   = name
        if angmom == None:
            self._angmom = [ 0,0,0 ]
        else :
            self._angmom = angmom # some 3-index list
        if centre == None :
            centre = dsl_position('c'+self._name.upper(), constant = constant )
        elif not isinstance(centre,dsl_position):
            raise Exception('centre must be a dsl_position object')
        self._centre = centre
        if exponent == None :
            exponent = dsl_scalar('x'+self._name, constant = self.constant_attachments() )
        elif not isinstance(exponent,dsl_scalar):
            raise Exception('exponent must be a dsl_scalar object')
        self._exponent = exponent
        self.add_attachment( 'exponent', self._exponent )
        self.add_attachment( 'centre', self._centre )
    def __str__(self):
        return self._name
    def __getitem__(self,key):
        return dsl_integer_index( src( self.name() ), self._angmom[key] ) 
    def cen(self):
        return self._centre
    def exp(self):
        return self._exponent
    def name(self):
        return self._name
    def angmom(self):
        return self._angmom
    def set_angmom(self,angmom):
        assert sum( angmom ) >= 0, 'total angular momentum cannot be < 0'
        self._angmom = angmom
    def increment(self):
        """Increments to the next level of angular momentum, using only the 
        current level of angular momentum as input"""
        l = sum(self._angmom)
        assert l >= 0, 'total angular momentum cannot be < 0'
        if self._angmom[2] == l:
            l = l + 1
            self._angmom = [ l , 0 , 0 ]
        elif self._angmom[1] == 0:
            self._angmom[0] -= 1
            self._angmom[1]  = l - self._angmom[0]
        else : 
            self._angmom[1] -= 1
        self._angmom[2] = l - self._angmom[1] - self._angmom[0]
    def decrement(self):
        """Increments to the previous level of angular momentum, using only the 
        current level of angular momentum as input"""
        l = sum(self._angmom)
        assert l > 0, 'total angular momentum cannot be < 0'
        if self._angmom[0] == l:
            l = l - 1
            self._angmom = [ 0 , 0 , l ]
        elif self._angmom[1] == l-self._angmom[0]:
            self._angmom[0] += 1
            self._angmom[1]  = 0
        else : 
            self._angmom[1] += 1
        self._angmom[2] = l - self._angmom[1] - self._angmom[0]

    def jump(self,change,direction):
        """
        Calculates the change in array index (in a 1D array) for a given 
        increment/decrement (change) and Cartesian direction (if applicable) 
        of increment/decrement.

        For a Cartesian index, direction is required.
        """
        assert isinstance( direction, dsl_direction ),\
                'direction must be dsl_direction'
        return cartesian_index_change(self._angmom,change,direction)

    def index(self,change=None,direction=None):
        """Provides the integer index of an angular momentum vector in a list
        of angular momentum vectors ordered in the standard way"""
        l = sum(self._angmom)
        j = self._angmom[1]
        k = self._angmom[2]
        if change == None and direction == None:
            return int( round( (l+2)*(l+1)*l/6 + (j+k+1)*(j+k)/2 + k + 1 ) )
        else:
            assert isinstance(change,int), 'change must be an integer'
            assert isinstance(direction,dsl_direction),\
                    'direction must be an instance of dsl_direction'
            return self.index() + self.jump( change, direction )
    def assign(self,expr=None):
        """Override method in parent class."""
        return NotImplemented 
    def constant_attachments(self):
        return self._constant_attachments



# Additional DSL objects
class dsl_direction:
    """Class that represents a Cartesian direction (x, y, z)."""
    def __init__(self,d):
        self.set(d)

    def index(self):
        return self._index

    def label(self):
        return self._label

    def set(self,d):
        """Sets the direction of the object (x, y, z).
        d can be a number (0, 1, 2) or a letter, 'x', 'y', 'z'.
        It can also be represented by a dsl_scalar object (so that
        the direction can be set at run time rather than code
        generation time)."""
        dlabel_to_index = { 'x' : 0, 'y' : 1, 'z' : 2 }
        dindex_to_label = { 0 : 'x', 1 : 'y', 2 : 'z' }
        if isinstance(d,str):
            assert d in [ 'x', 'y', 'z' ], 'd must label a Cartesian direction (x, y, z)'
            self._index = dlabel_to_index[d]
            self._label = d
        elif isinstance(d,int):
            assert d in [ 0, 1, 2 ], 'd must label a Cartesian direction (0, 1, 2)'
            self._index = d
            self._label = dindex_to_label[d]
        elif isinstance(d,dsl_scalar):
            # Direction is set during run time rather than at code-generation time.
            self._index = d
            self._label = d
        else:
            raise Exception("direction must be str 'x','y','z', int 0,1,2 or a dsl_scalar object")



