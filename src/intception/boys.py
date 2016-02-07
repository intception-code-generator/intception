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
import os
import sys
import math
from intception.dsl import *
from intception.printer import printer
from intception.read_json import read_json

class boys_function:
    def __init__(self,name='intception_boys',f_label='f',nmin_label='nmin',n_label='nmax',x_label='x',\
                 smallx_data_file = None,mediumx_data_file = None):
        # Set up some supporting information
        if smallx_data_file == None:
            raise Exception('Tabulated data for small x required.')
        else:
            smallx_data = read_json(smallx_data_file)
            self._smallx_xmin = dsl_value('small_xmin',smallx_data['xmin'],vartype='double')
            self._smallx_xmax = dsl_value('small_xmax',smallx_data['xmax'],vartype='double')
            # xmin, xmax are stored as strings containing the desire double precision numbers
            self._smallx_xstep = dsl_value('small_xstep',smallx_data['xstep'],vartype='int')
            self._smallx_xpoints = dsl_value('small_xpoints',smallx_data['xpoints'],vartype='int')
            self._smallx_nvals = dsl_value('small_nvals',smallx_data['nvals'],vartype='int')
            self._smallx_ntaylor = dsl_value('small_ntaylor',smallx_data['ntaylor'],vartype='int')
            # xstep, xpoints, nvals, ntaylor are stored as integers, so can be manipulated in Python
            # as numbers
            self._smallx_boys_data = smallx_data['boys data']
            self._smallx_pointer = dsl_pointer('boys_data_smallx',vartype='double', constant = True)
        if mediumx_data_file == None:
            raise Exception('Tabulated data for medium x required.')
        else:
            mediumx_data = read_json(mediumx_data_file)
            self._mediumx_xmin = dsl_value('medium_xmin',mediumx_data['xmin'],vartype='double')
            self._mediumx_xmax = dsl_value('medium_xmax',mediumx_data['xmax'],vartype='double')
            # xmin, xmax are stored as strings containing the desire double precision numbers
            self._mediumx_xstep = dsl_value('medium_xstep',mediumx_data['xstep'],vartype='int')
            self._mediumx_xpoints = dsl_value('medium_xpoints',mediumx_data['xpoints'],vartype='int')
            self._mediumx_nvals = dsl_value('medium_nvals',mediumx_data['nvals'],vartype='int')
            self._mediumx_ntaylor = dsl_value('medium_ntaylor',mediumx_data['ntaylor'],vartype='int')
            # xstep, xpoints, nvals, ntaylor are stored as integers, so can be manipulated in Python
            # as numbers
            self._mediumx_boys_data = mediumx_data['boys data']
            self._mediumx_pointer = dsl_pointer('boys_data_mediumx',vartype='double', constant = True)
        self._largex_xmin = dsl_value('xmin',self._mediumx_xmax)
        # Build one over factorial array (oof) 
        self._oof_pointer = dsl_pointer('oof',vartype='double', constant = True)
        self._oof_data = []
        self._oof_length = max( self._smallx_ntaylor.value(), self._mediumx_ntaylor.value() )
        self._oof_data.append( '1.0' ) # 1/(0!)
        self._oof_data.append( '1.0' ) # 1/(1!)
        for i in range( self._oof_length)[2:]:
            self._oof_data.append( '1.0/'+str( float(math.factorial( i ) ) ) )
        # Set up dsl function object
        f = dsl_pointer(f_label,vartype='double') # output
        ifskip = dsl_scalar('i'+f_label+'skip',vartype='int')
        nmin = dsl_scalar(nmin_label,vartype='int')
        n = dsl_scalar(n_label,vartype='int')
        x = dsl_scalar(x_label,vartype='double') 
        args = [ f, ifskip, nmin, n, x ]
        xn = dsl_scalar('xn',vartype='double') 
        dx = dsl_scalar('dx',vartype='double') 
        ii = dsl_scalar('ii',vartype='int')
        i  = dsl_scalar('i',vartype='int')
        expx = dsl_scalar('expx',vartype='double',expr=dsl_unop(op_exp, -x ) ) 
        expx.set_autoassign(False)
        oox  = dsl_scalar('oox',vartype='double',expr= 1 / x )
        oox.set_autoassign(False)
        ooxstep  = dsl_scalar('ooxstep',vartype='double' )
        ooxstep.set_autoassign(False)
        local_vars = [ xn, dx, ii, i, expx, oox, ooxstep ]
        # Dictionary for easy access of individual local variables
        self._boys_local_variable_dict = {}
        for local_var in local_vars:
            self._boys_local_variable_dict[ local_var.name() ] = local_var
        self._boys_dsl_function = dsl_function(name, args, vartype='void' )
        self._boys_xarg    = x
        self._boys_nminarg = nmin
        self._boys_narg    = n
        self._boys_output_pointer = f
        self._boys_ifskiparg = ifskip

    def boys_local_variable_dict(self):
        # Python dictionary
        return self._boys_local_variable_dict

    def boys_xarg(self):
        # dsl_scalar
        return self._boys_xarg

    def boys_nminarg(self):
        #dsl_scalar
        return self._boys_nminarg

    def boys_narg(self):
        # dsl_scalar
        return self._boys_narg

    def boys_ifskiparg(self):
        # dsl_scalar
        return self._boys_ifskiparg

    def boys_output_pointer(self):
        # dsl_pointer
        return self._boys_output_pointer

    
    def oof_pointer(self):
        # dsl_pointer
        return self._oof_pointer

    def oof_data(self):
        # Python list
        return self._oof_data

    def oof_length(self):
        # Python integer
        return self._oof_length
        
    def smallx_pointer(self):
        # dsl_pointer
        return self._smallx_pointer

    def smallx_boys_data(self):
        # Python list (nested)
        return self._smallx_boys_data

    def smallx_xmin(self):
        # dsl_value
        return self._smallx_xmin

    def smallx_xmax(self):
        # dsl_value
        return self._smallx_xmax
    
    def smallx_xstep(self):
        # dsl_value
        return self._smallx_xstep

    def smallx_xpoints(self):
        # dsl_value
        return self._smallx_xpoints

    def smallx_nvals(self):
        # dsl_value
        return self._smallx_nvals

    def smallx_ntaylor(self):
        # dsl_value
        return self._smallx_ntaylor

    def mediumx_pointer(self):
        # dsl_pointer
        return self._mediumx_pointer

    def mediumx_boys_data(self):
        # Python list (nested)
        return self._mediumx_boys_data

    def mediumx_xmin(self):
        # dsl_value
        return self._mediumx_xmin

    def mediumx_xmax(self):
        # dsl_value
        return self._mediumx_xmax
    
    def mediumx_xstep(self):
        # dsl_value
        return self._mediumx_xstep

    def mediumx_xpoints(self):
        # dsl_value
        return self._mediumx_xpoints

    def mediumx_nvals(self):
        # dsl_value
        return self._mediumx_nvals

    def mediumx_ntaylor(self):
        # dsl_value
        return self._mediumx_ntaylor

    def boys_dsl_function(self):
        # dsl_function
        return self._boys_dsl_function

# boys_function object, which can be used in input script
datadir = os.path.dirname( os.path.realpath(__file__) ) + "/boys_data/"
boys_default = boys_function(smallx_data_file=datadir+'boys_data_smallx.json',\
                     mediumx_data_file=datadir+'boys_data_mediumx.json')

class boys_f(dsl_base):
    """Class which enables the Boys function to be manipulated like other
    DSL objects. By default, uses the above instance of boys_function"""
    # The default boys_function_obj is a class variable, rather than an object variable
    # so can be accessed without creating an instance of the class.
    # This is the default Boys function implementation used by Intception.
    boys_function_obj = boys_default
    def __init__(self,m,x,varname='boys_f'):
        self._m = m # should be dsl_index
        self._x = x # should be dsl_scalar, or dsl_binop
        self._varname = varname

    def __str__(self):
        out = []
        out.append(self._varname)
        out.append('[ ')
        out.append( str(self._m) )
        out.append(' ]')
        return ''.join(out)

    def m(self):
        return self._m

    def x(self):
        return self._x

# Test Boys function output code
if __name__ == "__main__":
     p = printer(endl=';',output_file=sys.stdout)
     boys = boys_function(smallx_data_file='boys_data_smallx.json',\
                          mediumx_data_file='boys_data_mediumx.json')
     p.out( '#include<stdio.h>',endl='' )
     p.out( '#include<math.h>',endl='' )
     p.out( '#include"boys.h"',endl='' )
     boys.datadump(p)
     boys.funcdump(p)
