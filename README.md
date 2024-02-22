This repository contains the code for the paper: 

	Lemke SM, Celotto M, Maffulli R, Ganguly K, Panzeri S. Information flow between motor cortex and striatum reverses during skill learning.

Below are instructions to set up the sofware and paths to run the routines in "info_theory_code" script, folowed by a brief description of the analysis and the information-theoretic quantities that we compute.

# Prerequisites and Set Up

## Prerequisites

1. MATLAB (Tested with Versions > 2023a)
2. MATLAB supported C/C++ compiler: see www.mathworks.com/support/requirements/supported-compilers.html for installation
3. Python (Tested with Versions > 3.10) [^+]
4. BROJA_2PID package: can be downloaded from: https://github.com/Abzinger/BROJA_2PID. [^+]
5. ecos package: can be installed as 'pip install ecos' [^+]

## Setup instructions
1. Edit `ComputeInfoQuantatiesPackage/Code/startupInfoAnalysis.m`
    - edit `$BROJA_2PID_path` to full path where BROJA_2PID package has been downloaded
    - edit `$PIDPythonInterface_path` to full path where PIDPyInterface folder is located (found in ComputeInfoQuantatiesPackage/Code folder)
2. Add execution of file `./info_theory_code/startupInfoAnalysis.m` to the MATLAB startup file. [^+]
    - Open startup.m file in MATLAB using the following Matlab command: `edit(fullfile(userpath,'startup.m'))`: see https://it.mathworks.com/help/matlab/ref/startup.html for information on startup.m file
    - Add the following line to `startup.m`: `run 'PATH\Code\startupInfoAnalysis.m'`
      - where PATH is the full path to the downloaded `ComputeInfoQuantitiesPackage`
      - remember to use the correct forward or backslashes ("\" for windows, "/" for linux/mac)
    - Restart MATLAB
3. Run "compile_direct_method" script

Now you can succesfully run information-theoretic analyses!

[^+]: These steps are only required to run Partial Information Decomposition section. The other sections can be run without these prerequisites and setup steps.

# Description of Analyses
The main scripts to run and plot all analyses in the paper is 1_masterScript\masterScript

The masterScript will:
- Compute and binarize reaching features and neural signals (LFP and spiking activity)
- Compute Mutual Information (MI) between a neural signal and reach feature, computed across trials for each neural signal separately 
  at each time bin.
- Compute shared information that a pair of neural signals carry  about reach features. Shared information is computed using 
  the partial information decomposition (PID) framework. For each time  point, shared information is computed for a range of temporal delays  between the neural signals, to explore the possibility of one signal  carrying movement information that is present in the other signal at a  later time, and vice versa.
- The Transfer Entropy Information (TE) between a pair of neural signals. TE is computed at each time point, across temporal delays, to capture the the total amount of information transferred between two neural signals.
- The Feature-specific Information Transfer (FIT) between a pair of neural signals. FIT is computed at each time point, across temporal delays,  capturing the movement-specific information transferred between two neural signals.

For any issue in setting up information theoretic software, or in running the analysis scripts, please contact marco.celotto[at]iit.it