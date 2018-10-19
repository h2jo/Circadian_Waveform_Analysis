Codes for analyzing the data of PER2 in the mouse circadian system

"get_preprocess.py"

: This python file reads the data files for the curves and do the moving average with a given time window. For the details, see the comments in the python file.

"get_analysis.py"

: This python file reads the data files for the moving-averaged x(t) to estimate r(t) using the experimentally known values of r(t) at a few times. Then it calculates the local cost and global cost. For the details, see the comments in the python file.
