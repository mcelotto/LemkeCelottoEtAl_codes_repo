% setup PID
% path to BROJA2_PID package
% BROJA_2PID_path = '/$BROJA_2PID_PATH';
BROJA_2PID_path = 'C:\Users\mcelotto\Desktop\Neural Computation\PhD\Lemke\Reach_to_grasp_project\Scripts\code_from_Stefan_220224\BROJA_2PID-master';
% path to PIDPyInterface
%PIDPythonInterface_path = '$PIDPyInterface_PATH';
PIDPythonInterface_path = 'C:\Users\mcelotto\Desktop\Neural Computation\PhD\Sur_project\scripts\github_sync\LemkeCelottoEtAl_codes_repo\info_theory_code\PIDPyInterface';
insert(py.sys.path, int64(0), PIDPythonInterface_path);
insert(py.sys.path, int64(0), BROJA_2PID_path);

% load Matlab/Python interface modules
pymod = py.importlib.import_module('ComputePID');