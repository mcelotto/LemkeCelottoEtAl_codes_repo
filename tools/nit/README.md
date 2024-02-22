# ]:[ niːt neuroscience information toolbox
## Software requirements
#### Linux:
- MATLAB supported C/C++ compiler (tested with gcc version > 6.3.1)
- MATLAB supported Fortran compiler (only necessary for the GLM module, tested with gfortran version > 6.3.1)
- MATLAB (tested with versions > r2019b but earlier version should be compatible)
  - Statistics and Machine Learning toolbox
  - Optimization toolbox
- Python > 3.6
- MATLAB Engine API for Python (see installation instructions install [here](https://it.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html))
#### MacOS:
- MATLAB supported C/C++compiler
- MATLAB supported Fortran compiler (only necessary for the GLM module)
- MATLAB (tested with versions > r2019b but earlier version should be compatible)
  - Statistics and Machine Learning toolbox
  - Optimization toolbox
- Python > 3.6
- MATLAB Engine API for Python (see installation instructions [here](https://it.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html))
#### Windows:
- MATLAB supported C/C++ compiler (tested with MinGW-w64, see installation instructions [here](https://it.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler))
- MATLAB supported Fortran compiler (only necessary for the GLM module)
- Microsoft Visual C++
- Microsoft Visual C++ build tools
- MATLAB (tested with versions > r2019b but earlier version should be compatible)
  - Statistics and Machine Learning toolbox
  - Optimization toolbox
- Python > 3.6 (tested using Anaconda prompt to execute BuildAndTest.py, without the need for WSL)
- MATLAB Engine API for Python (see installation instructions [here](https://it.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html))
## Build and test instructions
Before starting the installation make sure that all software requirements described [above](software-requirements) are satisfied. After this please proceed with the following steps. (**A note for Windows users:** *it is recommended to execute all command line commands below from an Anaconda prompt rather than from Windows command line*)

1. Download the software and copy it into the install directory of your choice,

2. Open a command line terminal (*or an Anaconda prompt for Windows users*) 

3. Execute the build and test script:

   ​	`python BuildAndTest.py -d full/path/to/nit_root_dir`

   where `full/path/to/nit_root_dir` is the full path of the root directory of the software. 

The script compiles the code that needs compiling, adds the required source files to MATLAB search path by adding a `StartupNIT.m` file in the MATLAB user home folder and by editing the default user `Startup.m` file and performs the required install tests. If all software requirements are satisfied the process should seamlessly finish without errors.

Modules which are not required can be excluded from install by specifying, through the `-m` optional argument (see table below), the list of modules to be installed. 

A full list of arguments accepted by the build and test script is given below:

| Argument     | Description                                                  | `Value(s)`                                                   | Default | Mandatory |
| ------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------- | --------- |
| -d           | specifies the full path of the root directory of the software | `full/path/to/nit_root_dir`                                  | N/A     | yes       |
| --debug      | allows to compile the MEX files in debug mode                | N/A                                                          | false   | no        |
| -m --modules | list of modules to be included in the install                | `all`<br>`svm`<br>`mi_binmethods`<br>`mi_conmethods`<br>`glm`<br>`uml`<br>`pid`<br>`tools` | `all`   | no        |

If no `-m` option is given, by default the software installs/enables all modules.
## How to use NIT
A collection of demo scripts for each module can be found in the `src/Tutorials` directory.

**A(nother) note for windows users:** the *PID* module is based on Python. When building or using NIT on a Windows machine with the Intersection Information module active (i.e. if you ran `BuildAndTest.py` without any `-m` option or if you explicitly selected the install of the PID module through `-m pid`), it is recommended to launch MATLAB directly from an Anaconda prompt by opening the prompt and typing `matlab`.

## License

