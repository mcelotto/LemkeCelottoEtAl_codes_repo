This code package quantifies the information-theoretic quantities we use in the associated study: 

	Lemke SM, Celotto M, Maffulli R, Ganguly K, Panzeri S. Information flow between motor cortex and striatum reverses during skill learning.

Below are instructions to set up the sofware and paths to run the ComputeInformationQuantaties script, folowed by a brief description of the analysis and the information-theoretic quantities that we compute.

# Prerequisites and Set Up

## Prerequisites

1. MATLAB supported C/C++ compiler: see www.mathworks.com/support/requirements/supported-compilers.html for installation
2. Python (Tested with Versions > 3.6) [^+]
3. MATLAB Engine API for Python: see installation instructions at www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html [^+]
4. BROJA_2PID package: can be downloaded from: https://github.com/Abzinger/BROJA_2PID. [^+]

## Setup instructions
1. Edit `ComputeInfoQuantatiesPackage/Code/startupInfoAnalysis.m`
    - edit `$BROJA_2PID_path` to full path where BROJA_2PID package has been downloaded
    - edit `$PIDPythonInterface_path` to full path where PIDPyInterface folder is located (found in ComputeInfoQuantatiesPackage/Code folder)
2. Add execution of file `./Code/startupInfoAnalysis.m` to the MATLAB startup file. [^+]
    - Open startup.m file in MATLAB using the following Matlab command: `edit(fullfile(userpath,'startup.m'))`: see https://it.mathworks.com/help/matlab/ref/startup.html for information on startup.m file
    - Add the following line to `startup.m`: `run 'PATH\Code\startupInfoAnalysis.m'`
      - where PATH is the full path to the downloaded `ComputeInfoQuantitiesPackage`
      - remember to use the correct forward or backslashes ("\" for windows, "/" for linux/mac)
    - Restart MATLAB
3. Run "Setup" section in `ComputeInformationQuantaties.m` script

Now you can succesfully run the main script!

[^+]: These steps are only required to run Partial Information Decomposition section. The other sections can be run without these prerequisites and setup steps.

# Description of Analyses
In the main script, we provide an example of the information-theoretic analyses we use in the associated study.
To keep computations simple, we analyze a single pair of LFP channels recorded simultaneously from M1 and DLS during a day of reach-to-grasp training. We include one reach feature, the maximum reaching velocity  during each trial. In the study we compute a null hypothesis for each  information quantity via trial-shuffling which we have not included for simplicity.
The script will compute these information quantities:
- Mutual Information (MI) between a neural signal and the included  reach feature, computed across trials for each neural signal separately 
  at each time bin.
- Shared or redundant information that a pair of neural signals carry  about the included reach feature. Shared information is computed using 
  the partial information decomposition (PID) framework. For each time  point, shared information is computed for a range of temporal delays  between the neural signals, to explore the possibility of one signal  carrying movement information that is present in the other signal at a  later time, and vice versa.
- The Directed Information (DI) between a pair of neural signals. DI is computed at each time point, across temporal delays, to capture the the total amount of information transferred between two neural signals.
- The Feature-specific Information Transfer (FIT) between a pair of neural signals. FIT is computed at each time point, across temporal delays,  capturing the movement-specific information transferred between two neural signals.
