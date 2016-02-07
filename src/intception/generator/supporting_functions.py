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
from intception.dsl import dsl_variable
from intception.dsl_extensions import expr_list_objects

# Subroutines supporting source code generation
def expr_recursive_search( var ):
    """
    Recursively search dsl_variable object expr() values for additional
    dsl_variable objects and list dsl_variable objects found.

    The list returned is of unique dsl_variable objects, and since it
    is internally represented by a set, the dsl_variable class should
    provide __hash__ and comparison methods (currently uses the Python
    default, which is based on object identity).
    """
    var_set = set()
    def recursive_search( var ):
        nonlocal var_set
        if var.expr() != None:
            l = expr_list_objects( var.expr(), dsl_variable )
            for v in l:
                var_set.add(v)
                recursive_search( v )
        elif var.expr() == None and var.component_of() != None:
            var_set.add( var.component_of() ) 
            recursive_search( var.component_of() )
    recursive_search( var )
    var_list = list( var_set )
    return var_list

def determine_assignment_order(variable_list):
    """Checks the expr values of dsl_scalars and dsl_positions with autoassign == True for
    interdependencies, and reorders the self._var_order and self._position_order lists to
    make sure that during assignment, dependencies are satisfied.
    A dependency graph is generated, followed by a topological sort.
    
    Additionally, if items in variable_list are components of other dsl_variable objects
    (i.e. dsl_variable.component_of() != None), then the components are removed from the
    list and replaced with the dsl_variable they are a component of. If variable_list
    contained dsl_scalar x[0], then the dsl_position x would replace it in the output, since
    this is what would be assigned in the source code."""
    assign_list = []
    class rnode:
        def __init__(self,value,parents = [],children = []):
            self._value = value
            self._parents = []
            self._children = []
    

    nodes = {}
#    for v in list( self._vars.values() ) + list( self._posns.values() ):
    for v in variable_list:
#        if v.autoassign() == True:
#            if v.expr() == None:
#                raise Exception('autoassign set to True, but no expr set for '+v.name() )
#            nodes[ v.name() ] =  rnode( v )
        if v.component_of() != None:
            # Replace vector components with the vector itself, e.g. dsl_scalar objects
            # with component_of() != None will be replaced with the dsl_position to which
            # component_of() is set.
            assert v.expr() == None, "when component_of() == True, expr() must be None"
            nodes[ v.component_of().name() ] = rnode( v.component_of() )
        else:
            nodes[ v.name() ] =  rnode( v )
            
    
    # Build dependency graph
    for v in variable_list:
        if v.expr() != None: 
            #var_list = expr_list_objects( v.expr(), dsl_variable )
            var_set = set( expr_recursive_search( v )  )
            # check for vector components and add corresponding positions posn_list if present
            remove_set = set()
            for var in var_set:
                if var.component_of() != None:
                    assert isinstance(var.component_of(),dsl_variable),\
                            "a dsl_variable can only be a component of a dsl_variable"
                    assert var.component_of() in var_set,\
                            "component found, but component_of() object not in var_list"
                    remove_set.add( var )
            for var in remove_set:
                var_set.discard( var )
            # DEBUG
            #print( v.name(), [ x.name() for x in var_list ] )
                    
            for var in var_set:
                if var in variable_list:
                    nodes[ v.name() ]._parents.append( var ) 
                    nodes[ var.name() ]._children.append( v )

#   print('Dependency graph')
#   for node in nodes.values():
#      print( node._value.name()+':')
#      print( 'Parents:', [ n.name() for n in node._parents ] )
#      print( 'Children:', [ n.name() for n in node._children ] )

    # Topological ordering algorithm based on 
    # http://architects.dzone.com/articles/algorithm-week-topological
    # NB. This will fail if the dependency graph is not a directed acyclic
    # graph. However, Python should not allow cycles (circular dependencies)
    # since it will not allow you to use a variable that has not yet been
    # declared.
    orphans = []
    # Collect all variables and positions which have no parents
    for n in nodes.values():
        if len( n._parents ) == 0 :
            orphans.append( n )
    stack = orphans
    ordered_nodes = []
#    print( 'stack:  ',[ n._value.name() for n in stack ] )
#    print( 'ordered:',[ n._value.name() for n in ordered_nodes ] )
    while len( stack ) > 0:
        node = stack.pop()
        ordered_nodes.append( node )
        for child in node._children:
            child_node = nodes[child.name()]
            child_node._parents.remove( node._value )
            if len( child_node._parents ) == 0:
                stack.append( child_node )
        node._children = []
#        print( 'stack:  ',[ n._value.name() for n in stack ] )
#        print( 'ordered:',[ n._value.name() for n in ordered_nodes ] )

    for node in nodes.values():
        # Check that the dependency graph has no branches, even though this should not
        # occur due to Python preventing the use of undeclared variables.
        # If the sort fails, then it will end before all the branches have been removed.
        # As a result, some nodes will still have parents, indicating a cyclic graph.
        if len( node._parents ) != 0:
            raise Exception('Topological sort failed! Possible circular dependencies.')

    for node in ordered_nodes:
        assign_list.append( node._value )
#    print( [ item.name() for item in assign_list ] )

    return assign_list

