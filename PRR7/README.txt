Codes for analyzing the data of PRR7 in the Arabidopsis circadian system

"get_preprocess.py"

: This python file reads the data files for the curves in pixels and transforms them into the normalized curves as well as their splined curves. For the details, see the comments in the python file.

"get_analysis.py"

: This python file reads the data files for splined x(t) and g_m(t) to estimate r(t) using the experimentally known values of r(t) at a few times. Then it calculates the local cost and global cost. For the details, see the comments in the python file.
