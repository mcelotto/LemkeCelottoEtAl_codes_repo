disp('Building GLMnet mex files:')

GLMnetDir = 'C:\Users\slemke\Documents\InfoFlowProject\nit\extern\glmnet_matlab';
glmnet_file1 = fullfile(GLMnetDir, 'glmnetMex.F');
glmnet_file2 = fullfile(GLMnetDir, 'GLMnet.f');
arch = computer('arch');
if 0
	compileOpts = {'-g'};
else
	compileOpts = {'-O'};
end

try
    if strcmp(arch,'maci64') 
        mex('FFLAGS= -fdefault-real-8 -ffixed-form -fPIC -w', compileOpts{:}, glmnet_file1, glmnet_file2)
    elseif strcmp(arch,'glnxa64')
        mex('FFLAGS= -fdefault-real-8 -ffixed-form -fPIC -w', compileOpts{:}, glmnet_file1, glmnet_file2)
    elseif strcmp(arch,'win64')
        mex('FFLAGS= -fPIC -w', '-largeArrayDims', compileOpts{:}, glmnet_file1, glmnet_file2)
    elseif strcmp(arch,'win32')
        mex('FFLAGS= -fPIC -w', compileOpts{:}, glmnet_file1, glmnet_file2)
    end
    testflag = 1;
    disp('Build OK');
catch
    testflag = 0;
    warning('MATLAB did not detect any fortran compiler. GLMnet cannot be used with this installation')
end