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
import types # needed to check for builtin_function_or_method 
from intception import classtools
from intception.dsl import *
from intception.supporting_functions import isinstance_list

# Functions that directly operate on DSL expressions and DSL objects
def expr_permute(expr,a,b) : 
    """
    Permute indices a and b in dsl_integral, dsl_cartesian_gaussian and
    dsl_variable objects in an RR expression. The returned expression is
    a copy of expr, but with the dsl_index objects a and b swapped.
    Additionally, any dsl_variable (or derived class) objects which are
    attached to the index objects (i.e. returned by dsl_index.attachment_list)
    are also swapped in the expression. 
    dsl_variable objects for which expr() returns a DSL expression containing 
    other dsl_variable objects will be replaced with new dsl_variable objects
    with the expr() permuted in the same way. 
    This is done recursively, so a dsl_variable.expr() is treated in the same
    way as the DSL expression containing the dsl_variable.
    If expr contains dsl_integral objects, then these must have both a and b
    associated with them. 
    For dsl_integral objects, the increments/decrements in a and b are swapped, so
    if a = ga, b = gb and the expression contains
        dsl_integral( ga+1, gb-1 )  
    This dsl_integral object is replaced by a new dsl_integral
        dsl_integral( ga-1, gb+1 )
    """
    ### Begin expr_permute code ###
    assert type(a) == type(b), 'type(a) != type (b)'
    assert isinstance(a,dsl_index), 'a is not a dsl_index object'
    assert isinstance(b,dsl_index), 'b is not a dsl_index object'
    if expr is a:   
        expr = b
    elif expr is b: 
        expr = a
    elif isinstance(expr,dsl_binop):
        expr = dsl_binop( expr.op(), expr_permute( expr.left(),a,b),\
                              expr_permute( expr.right(),a,b ) )
#        print('out:', type(expr.left() ), type(expr.right() ) )
    elif isinstance(expr,dsl_unop):
        expr = dsl_unop( expr.op(), expr_permute( expr.arg(),a,b) )
    elif isinstance(expr,dsl_integral):
        new_int = expr.int()
        if isinstance_list( a,expr.indextypes() ) and\
           isinstance_list( b,expr.indextypes() ):
               if a in expr.index_list() and b in expr.index_list():
                    a_binop = expr.index_binop( a )
                    b_binop = expr.index_binop( b )
                    new_int.add_index(a,change=b_binop.right(),op=b_binop.op(),update=True)
                    new_int.add_index(b,change=a_binop.right(),op=a_binop.op(),update=True)
               if a in expr.index_list(): 
                   assert b in expr.index_list(), 'dsl_integral instance only has '+\
                        a.name()+'index. Cannot permute.'
               elif b in expr.index_list(): 
                   assert a in expr.index_list(), 'dsl_integral instance only has '+\
                        b.name()+'index. Cannot permute.'
        elif isinstance_list( a,expr.auxtypes() ) and\
             isinstance_list( b,expr.auxtypes() ):
               if a in expr.aux_list() and b in expr.aux_list():
                    a_binop = expr.aux_binop( a )
                    b_binop = expr.aux_binop( b )
                    new_int.add_auxiliary(a,change=b_binop.right(),op=b_binop.op(),\
                            update=True)
                    new_int.add_auxiliary(b,change=a_binop.right(),op=a_binop.op(),\
                            update=True)
        expr = new_int
    elif isinstance(expr,dsl_variable):
        # If dsl_variable objects are directly attached to dsl_index objects a,b
        # then simply swap the objects in the DSL expression.
        if expr in a.attachment_list() :
            for k, v in a.attachment_dict().items():
                if expr is v: 
                    key = k 
                    break
            expr = b.attachment( key )
        elif expr in b.attachment_list() :
            for k, v in b.attachment_dict().items():
                if expr is v: 
                    key = k 
                    break
            expr = a.attachment( key )
        # If dsl_variable objects are not directly attached, then do a recursive
        # search of the variable dependencies to see if they carry objects attached
        # to dsl_index object a,b. If they do, create a new dsl_variable object
        # of the appropriate derived type and permute the expression of that.
        elif expr.expr() != None:
            # Permutation symmetry check
            if expr_basic_symmetry_check( expr.expr(), a, b ) == False:
                new_variable = expr.__class__(expr = expr_permute( expr.expr(), a, b ), vartype=expr.type() )
                expr = new_variable
    return expr

def expr_basic_symmetry_check(expr, a, b):
    """
    Basic check for symmetry in permuting dsl_index objects a and b in an expression 
    of dsl_* objects.

    Returns True if the only permutations that occur involve swapping arguments of a 
    dsl_binop object with a fully associative operator. 

    This is used to prevent new variables being created during permuting of dsl_index
    objects (using expr_permute), where the new variable would clearly be the same 
    as the variable it is derived from. 
    This will not pick up cases of symmetry where permutations are not simple swaps of
    arguments of a dsl_binop, so is not a full check for permutational symmetry.

    If a dsl_integral object is in expr, then an AssertionError will occur as the
    symmetry of dsl_integral objects with respect to index permutation is not considered.

    Implementation notes:
    * Ideally this would be replaced with a subroutine that is able to completely determine
      whether a permuted expression is mathematically equivalent to the non-permuted 
      expression. Then, all redundant creation of new dsl_variable objects could be avoided.
    * This basic check will miss various symmetries, but will at least reduce the number of
      redundant variables created.
    """
    symmetry = True
    def check(expr, a, b):
        nonlocal symmetry
        assert not isinstance(expr,dsl_integral),\
                "expr_basic_symmetry_check cannot process dsl_integral objects"
        if isinstance( expr, dsl_binop ):
            l = expr.left()
            r = expr.right()
            if expr.op().assoc() == 'full':
                # Operator not fully associative
                if (l is a and r is b) or (r is a and l is b):
                    # fully associative and permutation of index
                    return
                elif isinstance(l,dsl_variable) and isinstance(r,dsl_variable):
                    l_key = None
                    r_key = None
                    if l in a.attachment_list() and r in b.attachment_list():
                        for k, v in a.attachment_dict().items():
                            if l is v: 
                                l_key = k 
                                break
                        for k, v in b.attachment_dict().items():
                            if r is v: 
                                r_key = k 
                                break
                    elif l in b.attachment_list() and r in a.attachment_list():
                        for k, v in b.attachment_dict().items():
                            if l is v: 
                                l_key = k 
                                break
                        for k, v in a.attachment_dict().items():
                            if r is v: 
                                r_key = k 
                                break
                    if l_key != None and r_key != None and l_key == r_key:
                        # fully associative and permutation of attached variables
                        return
            check(l, a, b )
            check(r, a, b )
        elif isinstance( expr, dsl_unop ):
            check( expr.arg(), a, b )
            return
        elif (expr is a) or (expr is b):
            # Assume not a symmetric permutation, since can only detect swaps of
            # dsl_binop arguments
            symmetry = False
            return
        elif isinstance( expr, dsl_variable ):
            if ( expr in a.attachment_list() ) or \
                    ( expr in b.attachment_list() ):
                # Assume not a symmetric permutation, since can only detect swaps of
                # dsl_binop arguments
                symmetry = False
                return
            elif expr.expr() != None:
                check( expr.expr(), a, b )
                return
    check( expr, a, b )
    return symmetry


def expr_deepish_copy(expr):
    """Creates a copy of an expression tree of dsl_binop objects, but 
    does not copy the tips of the tree so that other dsl objects are the
    same as those in the original expression."""
    if isinstance(expr,dsl_binop):
        expr = dsl_binop( expr.op(), expr_deepish_copy( expr.left() ),\
                              expr_deepish_copy( expr.right() ) )
    elif isinstance(expr,dsl_unop):
        expr = dsl_unop( expr.op(), expr_deepish_copy( expr.arg() ) )
    return expr

def expr_find_and_replace(expr,pattern,replacement):
    """Descend and expression tree of dsl_binop objects, and replace
    occurences of a specific dsl object with something else."""
    if isinstance(expr,dsl_binop):
        expr = dsl_binop( expr.op(), 
               expr_find_and_replace( expr.left(), pattern,replacement ),\
               expr_find_and_replace( expr.right(), pattern, replacement ) )
    elif isinstance(expr,dsl_unop):
        expr = dsl_unop( expr.op(), expr_find_and_replace( expr.arg(), pattern, replacement ) )
    elif expr is pattern:
        expr = replacement
    return expr

def expr_simplify(expr):
    """Simplify an expression made of DSL objects. Remove terms that are
    clearly zero. Returns a simplifies expression consisting of DSL objects."""
    def rough_parse(expr):
        if isinstance(expr,dsl_binop):
            expr = dsl_binop( expr.op(), rough_parse( expr.left() ),\
                                  rough_parse( expr.right() ) )
            if isinstance(expr.left(),dsl_zero):
                if expr.op() == op_mul:
                    expr = dsl_zero()
                elif expr.op() == op_div:
                    expr = dsl_zero()
            elif isinstance(expr.right(),dsl_zero):
                if expr.op() == op_mul:
                    expr = dsl_zero()
                elif expr.op() == op_div:
                    raise Exception('Error, division by zero detected')
#            elif ( isinstance(expr.left(),dsl_integer_index) or \
#                   isinstance(expr.left(),numbers.Real) ) and \
#                 ( isinstance(expr.right(),dsl_integer_index) or \
#                   isinstance(expr.right(),numbers.Real) ):
            else:
                try:
                    eval_expr = eval( src( expr.left() ) + expr.op().pyname() + \
                                       src( expr.right() ) ) 
                    if eval_expr == 0:
                        expr = dsl_zero()
                    elif isinstance( eval_expr, numbers.Real ):
                        expr = eval_expr
                except (NegativeAngmomError,NameError,SyntaxError):
                    # Expected where expr.left() or expr.right() do not evaluate to
                    # numbers.Real
                    pass
                except TypeError:
                    print( "A TypeError was encountered. This may be due to a collision with a "+\
                           "Python built-in method or function." )
                    for v in [ expr.left(), expr.right() ]:
                        if isinstance( eval( src( v ) ), types.BuiltinFunctionType ) or\
                            isinstance( eval( src( v ) ), types.BuiltinMethodType ) :
                            print( "Object "+repr( v )+" collides with a built-in method of function! "+\
                                   "To avoid collisions, please modify the object, or the object from "+\
                                   "which it is derived to avoid collisions." )
                            try:
                                print( "src([obj]):   "+src( v ) )
                                print( "[obj].name(): "+v.name() )
                            except AttributeError:
                                # If v does not have a __src__ or name method, then expect NameError.
                                pass

                    raise
        elif isinstance(expr,dsl_unop):
            expr = dsl_unop( expr.op(), rough_parse( expr.arg() ) )
        else:
            pass
        return expr
    def tidy_zeroes(expr):
        if isinstance(expr,dsl_binop):
            expr = dsl_binop( expr.op(), tidy_zeroes( expr.left() ),\
                                  tidy_zeroes( expr.right() ) )
            if isinstance(expr.left(),dsl_zero) and \
               not isinstance(expr.right(),dsl_zero):
                    if expr.op() == op_add:
                        expr = expr.right()
                    elif expr.op() == op_sub:
                        expr = -expr.right()

            elif isinstance(expr.right(),dsl_zero) and \
               not isinstance(expr.left(),dsl_zero):
                    if expr.op() == op_add:
                        expr = expr.left()
                    elif expr.op() == op_sub:
                        expr = expr.left()
        elif isinstance(expr,dsl_unop):
            expr = dsl_unop( expr.op(), tidy_zeroes( expr.arg() ) )
        return expr
    rough_expr  = rough_parse( expr )
    simplified_expr = tidy_zeroes( rough_expr )
    return simplified_expr

def expr_list_objects(expr,object_type):
    """Descend the tree of dsl_binop objects to find unique objects
    and return a list of these.

    Objects of object_type must be hashable and provide an equality operation 
    (__eq__ or __cmp__). dsl_* objects use the default __hash__ and
    __cmp__ methods, which use the object identity (address), so
    can be used here.
    
    N.B. A list of objects of any class inherited by dsl_binop, or dsl_unop
    will exclude any dsl_binop and dsl_unop objects in expr."""
    assert not object_type is dsl_binop, "cannot list dsl_binop objects"
    assert not object_type is dsl_unop, "cannot list dsl_unop objects"
    object_set = set()
    def find_object(expr,object_type):
        nonlocal object_set
        if isinstance(expr,dsl_binop):
            expr = dsl_binop( expr.op(),\
                    find_object( expr.left(), object_type ),\
                    find_object( expr.right(), object_type ) )
        elif isinstance(expr,dsl_unop):
            expr = dsl_unop( expr.op(), find_object( expr.arg(), object_type ) )
        elif isinstance(expr,object_type):
             # Set only allows one of each object (objects must be hashable and comparable)
             object_set.add( expr )
        else:
            pass
    find_object(expr,object_type)
    # Convert set to list for output
    object_list = list( object_set )
    return object_list

# Classes which extend and support the DSL defined in dsl.py
class dsl_rr(classtools.classtools):
    """Class representing recurrence relations."""
    def __init__(self,name,rrtype,changing_index,expr):
        self._name = name
        self._rrtype = rrtype
        self._changing_index = changing_index
        self._expr = expr

        # [ Current implementation restriction ]
        assert self._rrtype in ['vrr','hrr'], 'rrtype must be "vrr" or "hrr"'

    def __str__(self):
        """Returns a string of the RR expression."""
        return str( self._expr )

#    def debug_print(self):
#        """Print object attributes for debugging."""
#        print( self.__repr__() )
#        print( 'Name: ' + str( self._name ) )
#        print( 'RR type: ' + str( self._rrtype ) )
#        print( 'Unrolled index: ' + str( self._unrolled_index.name() ) )
#        print( 'Changing index: ' + str( self._changing_index.name() ) )
#        print( 'RR expression:  ' + str( self._expr ) )

    def name(self):
        return self._name

    def rrtype(self):
        return self._rrtype

    def changing_index(self):
        return self._changing_index

    def expr(self):
        return self._expr

class dsl_integral(dsl_base):
    def __init__(self, *args, name=None, desc=None, translational_invariance = None,\
                 indextypes=None, auxtypes=None, argtypes=None, copy_from = None):
        """
        *args is a set of arguments which must be parsed by parse_arg_list to 
        populate the correct attributes of the dsl_integral object.
        Rules for argument input:
        [ to be added ]

        If copy_from is set to another dsl_integral object, then keyword variables
        from that instance of dsl_integral will be copied into the new object.
        Non-keyword variables (passed via *args) will not be copied, and must be
        set manually via *args.
        Variables not set through __init__() will also be copied if copy_from is set.
        """
        # Check copy_from is a dsl_integral object
        if copy_from != None:
            assert isinstance(copy_from,dsl_integral), 'copy_from must be a '+\
                    'dsl_integral object'

        ### Initialize non-keyword arguments ###
        index_list = []
        index_dict = {}
        aux_list   = []
        aux_dict   = {}
        arg_list   = []
        arg_dict   = {}
        if copy_from != None:
            for l, d, cfl, cfd in \
                    [ (index_list, index_dict, copy_from._index_list, copy_from._index_dict ),\
                    (aux_list, aux_dict, copy_from._aux_list, copy_from._aux_dict) ]:
                for x in cfl:
                    # append to copy of list
                    l.append( x )
                    item = cfd[ x.name() ] 
                    # create new dsl_binop in dict
                    d[ x.name() ] = dsl_binop( item.op(), item.left(), item.right() )
            arg_list = copy_from._arg_list[:]
            for a in arg_list:
                arg_dict[a.name()] = copy_from._arg_dict[a.name()]
        self._index_list = index_list
        self._index_dict = index_dict
        self._aux_list   = aux_list  
        self._aux_dict   = aux_dict  
        self._arg_list   = arg_list  
        self._arg_dict   = arg_dict  

        ### Initialize from keyword arguments ###
        if copy_from != None:
            name = copy_from._name
            desc = copy_from._desc
            translational_invariance = copy_from._translational_invariance
            indextypes = copy_from._indextypes
            auxtypes   = copy_from._auxtypes
            argtypes   = copy_from._argtypes
        self._name    = name
        self._desc    = desc
        # translational_invariance has 3 settings: True, False and None.
        #   True:   user asserts that integral is translationally invariant
        #   False:  user asserts integral is not translationally invariant
        #   None:   no user assertion
        # Currently, a setting of None tells Intception to treat the integral
        # as if it is not translationally invariant (do not apply simplifications
        # that apply where translational invariance is known).
        assert isinstance( translational_invariance, bool ) or\
                translational_invariance == None,\
                "translational_variance must be None or bool"
        self._translational_invariance = translational_invariance
        # If index, auxiliary and argument types specified, set 
        # object attributes, else set automatically when
        # parsing args
        self._indextypes = []
        self._auxtypes  = []
        self._argtypes = []
        arg_types = [ indextypes, auxtypes, argtypes ]
        attr_types = [ self._indextypes, self._auxtypes, self._argtypes ]
        for attr in attr_types:
            attr = []
        for arg, attr in zip( arg_types, attr_types ):
            if arg != None:
                if isinstance(arg,list):
                    # List of types, simply assign attr to arg
                    attr = arg
                else:
                    # Assume a single type, append to attr list
                    attr.append(arg)
        self.set_allowed_types()

        ### Initialize attributes not set by __init__() ###
        base = None
        rr_list = []
        if copy_from != None:
            # Copy attributes from copy_from
            base       = copy_from._base
            rr_list    = copy_from._rr_list[:]
        self._base    = base
        self._rr_list = rr_list

        # Parse argument list and set object attributes accordingly
        self.parse_args(args, copy_from)

        assert self._name != None and self._name != '',\
                "dsl_integral must be given unique descriptive name"
        # Auto-generate unique name if copy_from == None
        if copy_from == None:
            self._name = self.generate_unique_name()

#    def debug_print(self):
#        """Print detailed object information for debugging."""
#        print( self.__repr__() )
#        print( 'Name: ' + str( self._name ) )
#        print( 'Description: ' + str( self._desc ) )
#        for n, t, l, d in zip( ['Indexes', 'Auxiliaries', 'Arguments'],\
#                               [ self._indextypes, self._auxtypes, self._argtypes ],\
#                               [ self._index_list, self._aux_list, self._arg_list ],\
#                               [ self._index_dict, self._aux_dict, self._arg_dict ] ) :
#            print( str(n)+':' )
#            print( ' allowed types: ', [ x.__name__ for x in t ] )
#            print( ' list:          ', [ x.name() for x in l ] )
#            print( ' dict:          ', [ str( x ) for x in list( d.values() ) ] )
#        print('Base expr:' + str( self._base ) )
#        print('RRs:')
#        for rr in self._rr_list:
#            print( rr.name() + str( rr.expr() ) )
        
    def set_allowed_types(self):
        """If the allowed index, auxiliary and argument types are not
        manually specified, then manually set these attributes.

        Implementation notes:
        * Only 1 index type and 1 auxiliary type are allowed currently.
        * The auxiliary type must be dsl_integer_index.
        * The index type must be dsl_cartesian_gaussian.
        """
        allowed_indextypes = [ dsl_cartesian_gaussian ]
        allowed_auxtypes = [ dsl_integer_index ]
        allowed_argtypes = [ dsl_scalar, dsl_position ]
        if self._indextypes == []:
            self._indextypes = allowed_indextypes
        if self._auxtypes == []:
            self._auxtypes = allowed_auxtypes
        if self._argtypes == []:
            self._argtypes = allowed_argtypes

        # [ Current implementation restriction ]
        assert len(self._indextypes) == 1, 'only one index type allowed'
        assert self._indextypes[0] is dsl_cartesian_gaussian,\
                'only dsl_cartesian_gaussian type index allowed'
        assert len(self._auxtypes) == 1, 'only one index type allowed'
        assert self._auxtypes[0] is dsl_integer_index,\
                'only dsl_integer_index type auxiliary allowed'
        assert len(self._argtypes) == 2, 'only two argument types allowed'
        for argtype in [ dsl_scalar, dsl_position ]:
            assert argtype in self._argtypes,\
                    'only dsl_scalar, dsl_position type arguments allowed'

    def parse_args(self,args,copy_from = None):
        """Goes through a tuple of arguments passed to the __init__() method
        and updates the corresponding object attributes.

        If copy_from is set to a dsl_integral object, parse_args additionally
        checks that the index, auxilary and argument objects in args are 
        in copy_from._index_list, copy_from._aux_list and copy_from._arg_list,
        i.e. if the new object is being copied from another, no new indexes,
        auxiliary or arguments may be added.
        """

        def arg_is_dsl_binop(arg,upd=False):
            """If arg is dsl_binop, then the dsl_binop should represent the
            addition or subtraction of an integer literal from a index or
            auxiliary."""
            if arg.op() is op_add or arg.op() is op_sub:
                if isinstance_list(arg.left(),self._indextypes ):
                    if isinstance(arg.right(),int):
                        c = arg.right()
                        a = arg.left()
                        op = arg.op()
                        self.add_index(a,c,op,update=upd)
                        return None
                    else:
                        raise Exception('Only integer literals '\
                                       +'can be added to an index')
                elif isinstance_list(arg.left(),self._auxtypes):
                    if isinstance(arg.right(),int):
                        c = arg.right()
                        a = arg.left()
                        op = arg.op()
                        self.add_auxiliary(a,c,op,update=upd)
                        return None
                    else:
                        raise Exception('Only integer literals '\
                                       +'can be added to an index')
                else:
                    raise Exception('binary operation argument must take the form '+\
                                    '[index_object] + [integer] or '+\
                                    '[auxiliary_object] + [integer]')
            else:
                raise Exception('Can only accept addition or subtraction '\
                               +'for increments')

        def copy_from_check(obj):
            """Raises and AssertionError if obj is not present in the the
            _index_list, _aux_list or _arg_list attributes of copy_from."""
            all_obj_list = copy_from._index_list + copy_from._aux_list +\
                           copy_from._arg_list
            assert obj in all_obj_list,\
                obj.name()+' not present in copy_from. If copy_from is set, '+\
                'obj must have been added to the copy_from dsl_integral instance.'
        if copy_from != None: 
            upd = True
        else:
            upd = False
        i = 1
        for arg in args:
            obj = arg
            if isinstance(arg,dsl_binop):
                arg_is_dsl_binop(arg,upd)
                obj = arg.left()
            elif isinstance_list(arg,self._indextypes):
                # No change in a Cartesian Gaussian
                self.add_index(arg,update=upd)
            elif isinstance_list(arg,self._auxtypes):
                # No change in an auxiliary index
                self.add_auxiliary(arg,update=upd)
            elif isinstance_list(arg,self._argtypes):
                self.add_argument(arg,update=upd)
            else:
                raise Exception('Argument number '+str(i)+': '+str(arg)+\
                                ' not understood.')
            if upd == True: copy_from_check(obj)
            i +=1

        if upd==True:
            # Additionally ensure that the same number of indexes, auxiliaries,
            # and arguments are present in a copy and original
            assert len( copy_from._index_list ) == len( self._index_list ),\
                'cannot change number of indexes in copy'
            assert len( copy_from._aux_list ) == len( self._aux_list ),\
                'cannot change number of auxiliaries in copy'
            assert len( copy_from._arg_list ) == len( self._arg_list ),\
                'cannot change number of arguments in copy'


    def generate_unique_name(self):
        """Returns a unique string based on the set of indices, auxiliaries and
        arguments of a dsl_integral object, which can be used to label a function
        for generating that integral class in C source code."""
        list_of_lists = [ self._index_list, self._aux_list, self._arg_list ]
        if sum( [ len(l) for l in list_of_lists ] ) == 0:
            return None # don't set a name if no indexes, auxiliaries or args
        out = []
        out.append('intception_integral_')
        out.append( self._name+'_' ) 
        for l in list_of_lists:
            for obj in l:
                out.append( str( obj.name() ) )
            if len(l) > 0: out.append('_')
        del out[-1]
        return ''.join(out)

    def __str__(self):
        return dsl_base.__str__(self)

    def src(self): 
        """There is not sufficient information in this abstract class to output a string which
        fully represents an integral in source code. This should be overloaded in a derived class."""
        return self._name

    def name(self):
        """Returns a string that represents the integral class in C source code."""
        return self._name

    def description(self):
        """Returns a string that describes the integral in human-readable format."""
        return self._desc

    def translational_invariance(self):
        """Returns a bool (or None) indidcating whether the integral is translationally
        invariant:
           True:   user asserts that integral is translationally invariant
           False:  user asserts integral is not translationally invariant
           None:   no user assertion
        """
        # Currently, a setting of None tells Intception to treat the integral
        # as if it is not translationally invariant (do not apply simplifications
        # that apply where translational invariance is known).
        # In the future, Intception may attempt to automatically decide whether
        # an integral type is translationally invariant.
        return self._translational_invariance

    def base(self):
        """Returns the DSL expression which defines the zero-angular momentum case for
        the integral class. For integral classes with auxiliary indices, the base case
        may cover a range of auxiliary index values and may include functions like the
        Boys function. This expression (or set of expressions) is evaluated prior to
        the application of recurrence relation, and thus 'seeds' the recurrence relations."""
        return self._base

    def indextypes(self):
        """Returns list of allowed index types."""
        return self._indextypes

    def auxtypes(self):
        """Returns list of allowed auxiliary types."""
        return self._auxtypes

    def argtypes(self):
        """Returns list of allowed arguments types."""
        return self._argtypes

    def index_list(self):
        """Returns list of indexes."""
        return self._index_list

    def aux_list(self):
        """Returns list of auxiliaries."""
        return self._aux_list

    def arg_list(self):
        """Returns list of arguments."""
        return self._arg_list

    def index_binop(self,i):
        """Returns binary operation corresponding to index."""
        return self._index_dict[i.name()]
        
    def aux_binop(self,a):
        """Returns binary operation corresponding to auxiliary."""
        return self._aux_dict[a.name()]

    def rr(self,name,index):
        """Returns the rr with name and index specified."""
        for rr in self._rr_list:
            if rr.name() == name and rr.changing_index() is index:
                return rr
        return None

    def rr_list(self):
        """Returns a list of RRs associated with the dsl_integral object."""
        return self._rr_list

    def set_base(self,expr):
        """
        Sets the expression for the base (typicalluy zero angular momentum) case of the
        integral class.

        This expression may not contain dsl_integral objects (it must be defined only
        in terms of the constants and parameters relating to the integral class, e.g.
        exponents, centres).

        The expression may contain objects which are indexed by auxiliary indexes, e.g.
        the Boys function as a dsl_function object.

        expr:   expression consisting of DSL objects which describes the zero-angular momentum
                base class (zero angular momentum starting point for RRs).
        """
        assert isinstance(expr,dsl_base), "base expr must be a expression consisting of "+\
                "classes that inherit from dsl_base"
        # Check for dsl_integral objects and reject if present
        integral_list = expr_list_objects( expr, dsl_integral )
        assert len( integral_list ) == 0, "base expr cannot contain dsl_integral objects"
        # Set attribute
        self._base = expr

    def _add_obj(self,l_obj,l,d_obj,d,update):
        """Adds a objects to a list, l, and a dictionary, d, with key l_obj.name().
        If update == True, then dictionary d is updated, but list l is unchanged."""
        assert isinstance(l_obj.name(),str), str(l_obj)+' name() method does not return '+\
                'a string'
        if update == False:
            assert l_obj not in l, 'duplicate objects not allowed when update == False'
            l.append(l_obj)
        else:
            assert l_obj in l, 'update == True, but '+str(l_obj)+' not in list'
        d[ l_obj.name() ] = d_obj


    def add_index(self,index,change=0,op=op_add,update=False):
        """Adds a DSL object (e.g. dsl_cartesian_gaussian) to a list of 
        integral indices and a dictionary of indices.These are non-auxiliary indices, 
        i.e. they define the desired end-result class of integrals.

        If update == True, do not add a new index, but modify and existing dsl_binop
        in self._index_dict with updated change and op.

        Implementation notes: 
        * Currently this set of indices must be all of the same type, e.g.
          all dsl_cartesian_gaussian objects. In future this may be expanded to allow for
          integrals with mixed index types."""
        assert isinstance_list(index,self._indextypes),\
                'index must be one of the allowed types in self._indextypes'
        self._add_obj(index,self._index_list,\
                      dsl_binop(op,index,change),self._index_dict,\
                      update)
#        self._index_list.append(index)
#        self._index_dict[ index.name() ] = dsl_binop(op,index,change)

    def add_auxiliary(self, aux,change=0,op=op_add,update=False):
        """Adds a DSL object (e.g. dsl_integer_index) to a list of 
        integral auxiliary indices and a dictionary of auxiliary indices.
        These are indices that define integral intermediates, but are not required
        to define the final end-result class of integrals.

        If update == True, do not add a new index, but modify and existing dsl_binop
        in self._index_dict with updated change and op.
        
        Implementation notes: 
        * Currently this set of indices must be all of the same type, e.g.
          all dsl_integer_index objects. In future this may be expanded to allow for
          integrals with mixed auxiliary index types.
        * Only a single dsl_integer_index object is currently supported, and this can
          only be incremented in RR expressions (no decrements).
        * It is assumed that the auxiliary indices are equal to zero in the 
          final end-result class of integrals.
        """
        assert isinstance_list(aux,self._auxtypes),\
                'aux must be one of the allowed types in self._auxtypes'
        self._add_obj(aux,self._aux_list,\
                      dsl_binop(op,aux,change),self._aux_dict,\
                      update)
#        self._aux_list.append(aux)
#        self._aux_dict[ aux.name() ] = dsl_binop(op,aux,change)

        # [ Current implementation restriction ]
        assert len( self._aux_list ) <= 1, 'only one auxiliary allowed'
        assert op == op_add, 'only increments allowed in auxiliary index'
    
    def add_argument(self, arg, update=False):
        """Adds a DSL object (e.g. dsl_scalar) to a list and dictionary of 
        arguments required to define the integral class in addition to the 
        indices and auxiliaries. For example, a program for evaluating 
        nuclear attraction integrals ( a | 1/rC | b ) requires an argument
        specifying the nuclear centre 'C', as well as arguments associated with the
        two indices 'a' and 'b'. 
        A DSL object added to the list/dictionary of arguments will be appended to 
        the arguments for the subroutine used to generate the integral class and can 
        therefore be passed in by the calling unit.
        """
        assert isinstance_list(arg,self._argtypes),\
                'arg must be one of the allowed types in self._argtypes'
        self._add_obj(arg,self._arg_list,\
                      arg,self._arg_dict,\
                      update)
#        self._arg_list.append(arg)
#        self._arg_dict[ arg.name() ] = arg

    def _add_rr(self,rr):
        """Adds a dsl_rr object representing a recurrence relation (RR) to a list of
        RR operations. This generally should not be called by a user, but by one of
        add_vrr() or add_hrr(), or via __init__().
        """
        self._rr_list.append( rr )
    
    def add_vrr(self,name,index,expr): 
        """
        Takes a set of arguments nd creates a dsl_rr object appropriate to the 
        RR specified. This is passed to _add_rr which then processes the object and adds
        it to a list of RR operations associated with the integral type.
        

        Implementation notes:
        * VRRs are currently defined as (and limited to), RR expressions in which
          intermediate integrals have:
          - decrements in non-auxiliary (e.g. dsl_cartesian_gaussian) indices,
          - increments in auxiliary indices (e.g. dsl_integer_index), if present.
        * Permutation of indexes in VRRs is disabled
        """
        ### 22/01/15 ###
        # Permutation of indexes for VRR is currently disabled, since this causes
        # the creation of multiple variables with expressions that evaluate to the same
        # value. The equality of expressions is not easy to determine during code generation
        # so this feature has been disabled.
        #If the set of arguments corresponds to an existing VRR, but with a different
        #changing index (and a new expr has not been specified) then the indices in the 
        #existing RR expr are permuted to create a new dsl_rr object with a permuted expr.
        #rr_already_defined = False
        #for rr in self._rr_list:
        #    if name == rr.name():
        #        # rr already defined
        #        assert expr == None, 'RR '+name+' already defined. '+\
        #                             'New expr cannot be added.'
        #        assert rr.rrtype() == 'vrr', 'Only VRR exprs can be permuted.'
        #        rr_already_defined = True
        #        expr_to_permute = rr.expr()
        #        index_to_permute = rr.changing_index()
        #        break
        #if rr_already_defined == True:
        #    # Permute RR expr from already-added RR
        #    rr_expr = expr_permute( expr_to_permute, index, index_to_permute )
        #else:
        ##     assert len( self._rr_list ) == 0,\
        ##             "Currently only one unique VRR expression is allowed."
        #    rr_expr = expr
        ###
        assert index in self.index_list(),\
                "index being incremented must be in integral definition"
        rr_expr = expr
        # Add using method for arbitrary dsl_rr object
        self._add_rr( dsl_rr( name,'vrr',index,rr_expr) )

        # [ Current implementation restrictions ]
        # * Only decrements allowed in index objects
        # * Only increments allowed in auxiliary index objects
        # * Only changes of 0 and 1 allowed in auxiliary index objects
        # * VRRs may not follow HRRs
        for rr in self._rr_list:
            assert rr.rrtype() != 'hrr', "A VRR may not follow a HRR."
        integrals = expr_list_objects(rr_expr,dsl_integral)
        for integral in integrals:
            for i_binop in integral._index_dict.values():
                if i_binop.op() == op_add:
                    assert i_binop.right() == 0, \
                            'VRRs cannot contain increments in index objects.'
                assert i_binop.left() in self._index_list,\
                            'index '+i_binop.left().name()+' not present in integral definition'
            for a_binop in integral._aux_dict.values():
                if a_binop.op() == op_sub:
                    assert a_binop.right() == 0, \
                            'VRRs cannot contain decrements in auxiliary objects.'
                elif a_binop.op() == op_add:
                    assert a_binop.right() <= 1,\
                            'VRRs can only contain increments of up to 1 unit in auxiliary indexes'
                assert a_binop.left() in self._aux_list,\
                            'auxiliary '+a_binop.left().name()+' not present in integral definition'

    def add_trr(self,name,move_to):
        """
        Create a dsl_rr object corresponding to the transfer equation, or translational
        recurrence relation (TRR).

        Currently not implemented.
        """
        raise Exception("TRR (transfer equation) not implemented.")

    def add_hrr(self,name,move_from,move_to):
        """Takes a set of arguments and creates a dsl_rr object appropriate to the 
        RR specified. This is passed to _add_rr which then processes the object and adds
        it to a list of RR operations associated with the integral type.
        Uses the Head-Gordon-Pople type HRR expression for moving angular momentum between 
        integral indexes over a single electron coordinate:
             ( a + 1i | b ... n ) = ( a | b + 1i ... n ) - ABi ( a | b ... n )
        This simplifies RR input, since this expression is common to many integral types.

        Angular momentum is moved from the move_from index, to the move_to index.

        Implementation notes:
        * HRRs are currently defined as (and limited to), RR expressions in which 
          intermediate integrals have:
          - increments in one or more non-auxiliary (e.g. dsl_cartesian_gaussian) indices,
          - no change in auxiliary indices.
        * Only one type of HRR is currently supported (the HRR described above).
        * In order for a HRR to be defined, a VRR must be defined for move_from.
        """
        Rfromto = dsl_position( expr = move_from.cen() - move_to.cen() )
        # User defined variable for Rfromto, only add
        rr_expr = self.int( move_from+1,move_to-1 ) +\
                    Rfromto * self.int(move_from,move_to-1 ) 
        # Add using method for arbitrary dsl_rr object
        self._add_rr( dsl_rr( name,'hrr',move_to,rr_expr) )

        # [ Current implementation restrictions ]
        move_from_already_incremented = False
        for rr in self._rr_list:
            if rr.rrtype() == 'vrr' and rr.changing_index() is move_from:
                move_from_already_incremented = True
            elif rr.rrtype() == 'hrr' and rr.changing_index() is move_from:
                move_from_already_incremented = True
        assert move_from_already_incremented == True,\
                "move_from must be incremented by a preceding RR operation"


        integrals = expr_list_objects(rr_expr,dsl_integral)
        index_incremented = False
        for integral in integrals:
            for i_binop in integral._index_dict.values():
                if i_binop.op() == op_add:
                    if i_binop.right() > 0: 
                        index_incremented = True
            for a_binop in integral._aux_dict.values():
                assert a_binop.right() == 0, \
                            'HRRs cannot contain have changes in auxiliary objects.'
        assert index_incremented == True, 'HRRs must have at least one index incremented.'

    def int(self,*args):
        """Returns a new instance of dsl_integral, checking that the Cartesian
        Gaussians and indices given are the same as in the current instance, 
        though the changes in those Cartesian Gaussians and indices may change.
        Copy references to RRs and base class to new instance, and copy a 
        reference to the dsl_direction object which determines the Cartesian
        component for RRs."""
        # Checking that args are present in self and correct initialization of
        # variables in the new dsl_integral instance is done by __init__() and
        # parse_args.
        return dsl_integral( *args, copy_from = self )

