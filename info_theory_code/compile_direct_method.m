% Set MATLAB to current folder (containing this script)

    cd(fileparts(matlab.desktop.editor.getActiveFilename))

% add path to code for information theoretic analysis

    path_to_info_analysis_code = fullfile(pwd,'Code');
    addpath(genpath(path_to_info_analysis_code));

% compile code
    mex('-g', '-largeArrayDims', '-outdir', path_to_info_analysis_code,...
        fullfile(pwd,'info_theory_code','MI','direct_method.c'));