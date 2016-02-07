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
from intception import classtools
from intception.dsl import dsl_base
from intception.generator.arrays import general_array
from intception.printer import printer, src

class cprinter(classtools.classtools):
    """
    Output of source code from objects in the DSL language in C99 format.

    Automatically converts dsl_* objects to the src_dsl* representation and 
    provides additional language-specific methods of source code output.
    """
    def __init__(self,converter,output_file=sys.stdout,preamble_comment_str=None):
        """
        converter:      converter object providing expr_full_convert method
        output_file:    file object or file-like object to be used for text output
                        (must be open for writing)
        preamble_comment_str:
                        if not None, then immediately output a comment line containing
                        this string to the top of the output_file

        Note: the expr_full method should have the definition

            expr_full_convert(expr)

        where expr is a dsl_* expression (comprised of dsl_binop and dsl_unop objects), 
        or a dsl_* object, and which returns corresponding expr or object with all
        dsl_* objects converted to src_dsl_* objects.
        If expr is not a dsl_* expression or object (e.g. it is a str, or int), then
        it should pass through expr_full without modification.
        """
        self.converter   = converter
        self.printer_obj = printer(endl=';', output_file=output_file)
        if preamble_comment_str != None:
            assert isinstance( preamble_comment_str, str ), \
                    "preamble_comment_str must be instance of str" 
            self.comment( preamble_comment_str )

    def out(self,arg,endl=None,output_file=None):
        """
        Converts arg to src_dsl_* representation (if appropriate) and then
        calls self.printer_obj.out() to output text to a file or file-like object.
        """
        if isinstance( arg, dsl_base ) or isinstance( arg, general_array ):
            src_arg = self.converter.expr_full_convert( arg )
        else:
            src_arg = arg
        self.printer_obj.out(src_arg,endl,output_file)

    def indent(self):
        """
        Calls self.printer_obj.indent() to increment indentation level.
        """
        self.printer_obj.indent()

    def outdent(self):
        """
        Calls self.printer_obj.outdent() to decrement indentation level.
        """
        self.printer_obj.outdent()

    def blankline(self,output_file=None):
        """
        Calls self.printer_obj.blankline() to output a blank line to output_file.
        """
        self.printer_obj.blankline()

    def comment(self,comment_str):
        """
        Using self.printer_obj, output a C // comment line.
        """
        self.printer_obj.out( '// '+str(comment_str), endl='' )

    def include(self,header_filename,quotes=False):
        """
        Using self.printer_obj, output a C #include line.

        If quotes == False, then the output uses angle brackets (search directories in
        include path list first, typically):
            #include<[header_filename]>
        else, quotes are used (search directory containing source file first, typically)a:
            #include"header_filename"
        """
        if quotes == True: 
            self.printer_obj.out( '#include"'+str(header_filename)+'"', endl='' )
        else:
            self.printer_obj.out( '#include<'+str(header_filename)+'>', endl='' )

    def funcdef(self,dobj, const_dict = {}):
        """
        Coverts dsl_* object dobj to corresponding src_dsl_* object and then
        calls src_dsl_*.prototype(), adding a curly bracket and indent to
        open a function definition block.
        """
        src_obj = self.converter.expr_full_convert( dobj )
        assert hasattr(src_obj,'prototype'), "src_dsl_* object must have prototype() method"
        self.printer_obj.out( src_obj.prototype(const_dict) + '{ ', endl='' )
        self.printer_obj.indent()

    def endfuncdef(self):
        """
        Closes a function definition block.
        """
        self.endblock()

    def funccall(self,dfunc,darg_list=None):
        """
        Converts dsl_* object dfunc to corresponding src_dsl_* object, converts
        all darg_list dsl_* objects to src_dsl_* objects and then calls
        src_dsl_*.call(src_arg_list). The result of this call is output to the 
        file/file-like object through self.printer_obj.

        The src_dsl_*.call() method should output a valid function call for dfunc and
        be able to accept the argument src_arg_list, so that an alternative set of
        function arguments can be provided. If src_arg_list == None, then the default
        set of arguments is used (if any).
        """
        src_func = self.converter.expr_full_convert( dfunc )
        assert hasattr(src_func,'call'), "src_dsl_* object must have call() method"
        if darg_list == None:
            darg_list = dfunc.args()
        src_arg_list = []
        for dobj in darg_list:
            if isinstance(dobj,dsl_base) or isinstance(dobj,general_array): 
                src_arg = self.converter.expr_full_convert( dobj )
            else:
                src_arg = dobj
            src_arg_list.append( src_arg )
        self.printer_obj.out( src_func.call(src_arg_list) )

    def funcprototype(self,dfunc,const_dict = {}):
        """
        Converts dsl_* object dfunc to corresponding src_dsl_* object, converts
        all darg_list dsl_* objects to src_dsl_* objects and then calls
        src_dsl_*.prototype(). The result of this call is output to the 
        file/file-like object through self.printer_obj.
        """
        src_func = self.converter.expr_full_convert( dfunc )
        self.printer_obj.out( src_func.prototype(const_dict) )

    def returnstmt(self,arg=None):
        """
        Converts arg from dsl_* to src_dsl_* object,if appropriate and outputs a
        return statement. 

        If arg == None, then the return statement has no argument.
        """
        if arg == None:
            self.printer_obj.out('return')
        else:
            if isinstance(arg,dsl_base) or isinstance(arg,general_array):
                src_arg = self.converter.expr_full_convert( arg )
            else:
                src_arg = arg
            self.printer_obj.out('return '+src( src_arg ) )

    def declaration(self,dobj):
        """
        Coverts dsl_* object dobj to corresponding src_dsl_* object and then
        calls src_dsl_*.declaration()
        """
        src_obj = self.converter.expr_full_convert( dobj )
        assert hasattr(src_obj,'declaration'), "src_dsl_* object must have declaration() method"
        self.printer_obj.out( src_obj.declaration() )
        
    def forloop(self,start,condition,iteration):
        """
        Converts dsl_* objects, start, condition, iteration to src_dsl_* object and
        outputs the top part of a C for loop based on these, including a single indent.
        """
        src_start     = self.converter.expr_full_convert( start )
        src_condition = self.converter.expr_full_convert( condition )
        src_iteration = self.converter.expr_full_convert( iteration )
        self.printer_obj.out('for( '+src( src_start )+' ; '+\
                             src( src_condition )+' ; '+\
                             src( src_iteration )+' ) {', endl='' )
        self.indent()

    def endforloop(self):
        """
        Closes a C for loop.
        """
        self.endblock()

    def ifblock(self,condition):
        """
        Converts dsl_* object condition to src_dsl_* object and
        outputs the top part of a C if block including a single indent.
        """
        src_condition = self.converter.expr_full_convert( condition )
        self.printer_obj.out('if( '+src( src_condition )+' ) {',endl='')
        self.indent()

    def elifblock(self,condition):
        """
        Converts dsl_* object condition to src_dsl_* object and
        outputs the following:
            - close of preceding if or elif block, with single outdent,
            - top part of an elif block including a single indent.
        """
        src_condition = self.converter.expr_full_convert( condition )
        self.endblock()
        self.printer_obj.out('else if( '+src( src_condition )+' ) {',endl='')
        self.indent()

    def elseblock(self):
        """
        Outputs the following:
            - close of preceding if or elif block, with single outdent,
            - top part of an else block including a single indent.
        """
        self.endblock()
        self.printer_obj.out('else {',endl='')
        self.indent()

    def endifblock(self):
        """
        Closes a C if (..elif..else..) block.
        """
        self.endblock()

    def endblock(self,endl=None):
        """
        Close a block with a curly bracket and unindent one tab.
        """
        self.outdent()
        self.printer_obj.out('}',endl='')
        


        






