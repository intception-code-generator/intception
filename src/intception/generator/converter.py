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
# General intception modules

from intception.dsl import *
from intception.generator.src_dsl import *
from intception.generator.arrays import *
from intception.dsl_extensions import expr_deepish_copy, expr_find_and_replace, expr_simplify


class converter(classtools.classtools):
    """
    Provides an expr_to_src_convert method, for converting dsl_* expressions to
    src_dsl_* expressions.
    
    Carries with it information required for the conversion from dsl_* to 
    src_dsl_*.
    """
    def __init__(self,wintegral,direction,contraction_info,gen):
        """
        wintegral:  integral_wrapper object
        direction:  dsl_direction object
        contraction_info: contraction_info object
        gen:        generator object
        """
        self._wintegral = wintegral
        self._direction = direction
        self._contraction_info = contraction_info
        self._generator = gen

        self._contract_convert = True 
        # allows temporary suppression of contraction conversion

    def wintegral(self):
        return self._wintegral

    def direction(self):
        return self._direction

    def contraction_info(self):
        return self._contraction_info

    def generator(self):
        return self._generator

    def contract_convert(self):
        return self._contract_convert

    def set_contract_convert(self, contract_convert ):
        """
        Sets value of self._contract_convert (must be bool).
        
        If True (default), then expr_full_convert will call 
        expr_to_contract_convert during the conversion process.
        If False, expr_full_convert only calls expr_to_src_convert.
        This is useful for temporary supression of conversion of
        dsl_* objects to their contracted equivalents.
        """
        assert isinstance( contract_convert, bool ),\
                "contract_convert must be of type bool"
        self._contract_convert = contract_convert


    def expr_full_convert(self,expr):
        """
        Returns a src_dsl_* expression, with any subsitutions required for
        contracted integrals processed by running expr_to_contract_convert and
        expr_to_src_convert on expr.

        expr:   dsl_* object (can be dsl_binop/unop tree)

        If self._contract_convert == False, then expr_to_contract_convert
        is not called during conversion process.
        """
        # If contracted integrals required, run expr_to_contract_convert
        if self.generator().generator_options().contracted() == True and \
                self.contraction_info() != None and self.contract_convert() == True:
            # if self.contraction_info == None, assume contraction not required
            expr1 = self.expr_to_contract_convert( expr )
        else:
            expr1 = expr
        # Always run expr_to_src_convert
        expr2 = self.expr_to_src_convert( expr1 )
        return expr2

    def expr_to_contract_convert(self,expr):
        """
        Takes a DSL expr containing dsl_* objects and substitutes any
        objects which match keys in contraction_info.contraction_convert_dict
        for the corresponding values in the dictionary.
        """
        assert self.contraction_info() != None,\
                "a contraction_info object must be attached to converter for "+\
                "expr_to_contract_convert to function"
        contract_expr = expr_deepish_copy( expr )
        if not isinstance( contract_expr, dsl_binop ) and\
                not isinstance( contract_expr, dsl_unop ):
            # If expr is not a dsl_* expression tree, simply check dictionary for
            # matching object
            try:
                return self.contraction_info().contraction_convert_dict[ contract_expr ]
            except KeyError:
                # If no match found, return the original object unchanged
                return expr
        else:
            for k, v in self.contraction_info().contraction_convert_dict.items():
                contract_expr = expr_find_and_replace( contract_expr, k, v )
            return contract_expr

    def expr_to_src_convert(self,expr):
        """
        Takes a DSL expr containing dsl_* objects and replaces the 
        dsl_* objects with corresponding src_dsl_* wrappers.

        Also removes any terms in the src_dsl_* expression which unambiguously
        evaluate to zero.
        """
        if not isinstance(expr,dsl_binop) and\
            not isinstance(expr,dsl_unop):
            # If expr is not an expression, but just a single object, call faster
            # routine
            return self.dsl_to_src_dsl_convert(expr)
        wintegral  = self._wintegral
        direction  = self._direction
        #gen        = self._generator
        src_expr = expr_deepish_copy( expr )
        # list dsl_integral objects in expr
        dintegral_list = expr_list_objects( src_expr, dsl_integral )
        # list dsl_index objects in expr
        dindex_list    = expr_list_objects( src_expr, dsl_index )
        # list dsl_variable objects in expr
        dvariable_list = expr_list_objects( src_expr, dsl_variable )
        # list dsl_value objects in expr
        dvalue_list    = expr_list_objects( src_expr, dsl_value )
        # list dsl_function objects in expr
        dfunction_list = expr_list_objects( src_expr, dsl_function )
        # list dsl_pointer objects in expr
        dpointer_list  = expr_list_objects( src_expr, dsl_pointer )
        # additionally, check for general_array objects and list these
        general_array_list = expr_list_objects( src_expr, general_array )
        # swap dsl_integral objects for src_dsl_integral objects
        for dintegral in dintegral_list:
            if wintegral != None:
                # if wintegral is defined...
                # all dsl_integral objects should be the same class as the integral_wrapper
                # object with this method
                assert dintegral.name() == wintegral.integral().name(),\
                    "all dsl_integral objects in expr should be of same class as "+\
                    "the integral class being generated."
            #sintegral = src_dsl_integral( \
            #                dintegral, wintegral, direction, gen )
            sintegral = src_dsl_integral( dintegral, self )
            src_expr = expr_find_and_replace( src_expr, dintegral, sintegral )
        # swap dsl_* objects for src_dsl_* objects
        # src_dsl_{index,variable,value} objects have the same interface
#        for dobj_list, src_obj_type in [\
#                ( dindex_list, src_dsl_index ),\
#                ( dvariable_list, src_dsl_variable ),\
#                ( dvalue_list, src_dsl_value ),\
#                ( dfunction_list, src_dsl_function ), \
#                ( dpointer_list, src_dsl_pointer ) ]:
        #    for dobj in dobj_list:
        #        sobj = src_obj_type( dobj, self )
        #        src_expr = expr_find_and_replace( src_expr, dobj, sobj )
        dobj_list = dindex_list + dvariable_list + dfunction_list + dpointer_list + \
                    general_array_list
        for dobj in dobj_list:
            sobj = self.dsl_to_src_dsl_convert( dobj )
            src_expr = expr_find_and_replace( src_expr, dobj, sobj )
        # swap general_array objects for src_general_array objects
#        for ga in general_array_list:
#            src_ga = src_general_array( ga, self )
#            src_expr = expr_find_and_replace( src_expr, ga, src_ga )
        # Remove terms that unambiguously evaluate to zero
        src_expr = expr_simplify( src_expr )
        return src_expr

    def dsl_to_src_dsl_convert(self,dobj):
        """
        Takes a single dsl_* objects returns the  corresponding src_dsl_* wrapper.
        Prefer this to expr_to_src_convert where possible, since it should be faster.
        If dobj is not one of the convertable objects, dobj is returned unmodified.

        dobj: dsl_* object to be converted.

        Note: dobj cannot be a dsl_binop or dsl_unop, only an object with a corresponding
              src_dsl_* class.
        """
        assert not isinstance(dobj,dsl_binop) and\
                not isinstance(dobj,dsl_unop),\
                "dsl_to_src_dsl_convert cannot process dsl_binop/dsl_unop expression trees -- "+\
                "use expr_to_src_convert instead."
        # dsl_classname to dsl_* and src_dsl_* dictionary
        d = { dsl_integral : src_dsl_integral,\
              dsl_index    : src_dsl_index,\
              dsl_variable : src_dsl_variable,\
              dsl_value    : src_dsl_value,\
              dsl_function : src_dsl_function,\
              dsl_pointer  : src_dsl_pointer,\
              general_array: src_general_array }
        src_class = None
        for dclass in d.keys():
            if isinstance(dobj,dclass):
                src_class = d[dclass]
                break
        if src_class != None:
            # If a corresponding src_dsl_* class is found, return instance of this
            return src_class( dobj, self )
        else:
            # If no conversion is possible, just return original object unchanged
            return dobj
        
        

