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
# GNU Makefile for compilation of libintception to a static library
# James C. Womack, 2014-2016
# 
# Place the output of Intception in the directory ./src and execute
# 	make
# to build a static library libintception.a.
SHELL = /bin/sh

UNAME := $(shell uname)
CC := gcc
CSTD := --std=c99
COPT := -O3
WARN := -Wall -Wno-maybe-uninitialized
SRCDIR := ./output
SRC := $(wildcard ${SRCDIR}/*.c)
OBJ := $(patsubst %.c,%.o,${SRC})
DEP := $(patsubst %.c,%.d,${SRC})
NAME := intception
LIBNAME := lib${NAME}.a
# BLASLIB and BLASLINK must be set on command line, if required, e.g. for OS X, with CBLAS
# BLASLIB="-L/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current"
# BLASINC="-I/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers/"
# BLASLINK="-lcblas"
#BLASLIB := 
#BLASINC := 
#BLASLINK := 

# Compile shared library and use real name convention
# http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html
${LIBNAME} : ${OBJ}
	ar rcs ${LIBNAME} ${OBJ}

# The -MMD flag instructs GCC to output a %.d file in addition to the
# %.o file, containing dependency rules for that object file
${SRCDIR}/%.o : ${SRCDIR}/%.c
ifdef BLASLIB
	${CC} ${CSTD} ${COPT} ${BLASLIB} ${BLASINC} ${WARN} -MMD -c $< ${BLASLINK} -o $@ 
else
	${CC} ${CSTD} ${COPT} ${WARN} -MMD -c $< -o $@ 
endif

.PHONY : clean

clean :
	rm -f ${OBJ} ${DEP} ${LIBNAME}

# Include (if present) %.d dependency rules
-include ${DEP}

