This toolbox allows computing, in MATLAB, the intersection information as defined in the framework of Partial Information Decomposition. 
The toolbox uses the 'BROJA_2PID' Python package (https://github.com/Abzinger/BROJA_2PID) to perform the PID decomposition as defined in Bertshinger et al. 'Quantifying unique information' (2014).

The function 'II_BROJA2.m' computes the intersection information for a set of N neural signals recorded over T time steps for Ntr trials.

The function 'II_bias_corr.m'computes the limited sampling bias-corrected intersection information for a set of N neural signals recorded over T time steps for Ntr trials. 
The bias correction method can be either 'linear' or 'quadratic' extrapolation. Note that this function returns even the II_biased array, which is the same output of 'II.m'.

-------------------------------------------------------------------

To make these scripts work properly, follow the steps below:

Install Python 3.6 in '/usr/bin/':
>> sudo apt update
>> sudo apt-get install python3.6

Install pip:
>> sudo apt install python3-pip

Install the Embedded Conic Solver (ECOS):
>> pip3 install ecos

In 'BROJA2_for_matlab' open 'auxiliary_functions' and modify:

- In 'python.m' line 84 the pythonInst variable into the PATH where the python3.6 executable is located
(If you are using Windows you have to modify 'python3.6' --> 'python.exe' at line 80 too)

Then, in the main 'BROJA2_for_matlab' folder:

- In 'II_BROJA2.m' modify line 37 with the PATH where you put the 'BROJA2_for_matlab' folder
- In 'II_BROJA2_bias_corr.m' modify line 26 with the PATH where you put the 'BROJA2_for_matlab' folder

To test if everything works appropriately go to 'test' and run 'test.m'

Now you can call the 'II_BROJA2.m' and 'II_BROJA2_bias_corr.m' script giving as input the stimulus, responses and choice arrays for an experimental session.

