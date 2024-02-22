% Add all folders required by PID module to MATLAB path.
addpath(genpath('C:\Users\slemke\Documents\InfoFlowProject\nit\src\PID'));

% add paths to Python search paths
insert(py.sys.path, int64(0), 'C:\Users\slemke\Documents\InfoFlowProject\nit\src\PID\PIDPyInterface');
insert(py.sys.path, int64(0), 'C:\\Users\\slemke\\Documents\\InfoFlowProject\\nit\\extern\\BROJA_2PID');

% load Matlab/Python interface modules
pymod = py.importlib.import_module('ComputePID');
