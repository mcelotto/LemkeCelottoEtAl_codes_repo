% Add all folders required by PID module to MATLAB path.
addpath(genpath('$PID_MATLAB_ROOT'));

% add paths to Python search paths
insert(py.sys.path, int64(0), '$II_PYTHON_INTERFACE_ROOT');
insert(py.sys.path, int64(0), '$PID_BROJA_2PID_ROOT');

% load Matlab/Python interface modules
pymod = py.importlib.import_module('ComputePID');
