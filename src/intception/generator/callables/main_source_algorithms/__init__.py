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
# Define set of algorithms that are imported when 
#    import intception.generator.callables.main_source_algorithms
# is use
from . import primitive_vrr_only_no_auxiliary
from . import primitive_vrr_only_with_auxiliary
from . import primitive_vrr_hrr_no_auxiliary
from . import primitive_vrr_hrr_with_auxiliary
from . import contracted_vrr_only_no_auxiliary
from . import contracted_vrr_only_with_auxiliary
from . import contracted_vrr_hrr_no_auxiliary
from . import contracted_vrr_hrr_with_auxiliary

# Create a dictionary of algo_descriptors that maps descriptors to the correct module
algo_desc_dict = {\
        primitive_vrr_only_no_auxiliary.algo_desc : primitive_vrr_only_no_auxiliary, \
        primitive_vrr_only_with_auxiliary.algo_desc : primitive_vrr_only_with_auxiliary, \
        primitive_vrr_hrr_no_auxiliary.algo_desc : primitive_vrr_hrr_no_auxiliary, \
        primitive_vrr_hrr_with_auxiliary.algo_desc : primitive_vrr_hrr_with_auxiliary, \
        contracted_vrr_only_no_auxiliary.algo_desc : contracted_vrr_only_no_auxiliary, \
        contracted_vrr_only_with_auxiliary.algo_desc : contracted_vrr_only_with_auxiliary, \
        contracted_vrr_hrr_no_auxiliary.algo_desc : contracted_vrr_hrr_no_auxiliary, \
        contracted_vrr_hrr_with_auxiliary.algo_desc : contracted_vrr_hrr_with_auxiliary \
        }

