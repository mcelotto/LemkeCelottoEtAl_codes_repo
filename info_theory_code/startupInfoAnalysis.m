% setup PID
% path to BROJA2_PID package
BROJA_2PID_path = '/$BROJA_2PID_PATH';
% path to PIDPyInterface
PIDPythonInterface_path = '$PIDPyInterface_PATH';
insert(py.sys.path, int64(0), PIDPythonInterface_path);
insert(py.sys.path, int64(0), BROJA_2PID_path);

% load Matlab/Python interface modules
pymod = py.importlib.import_module('ComputePID');