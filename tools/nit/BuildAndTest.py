#
#  This source code is part of:
#  NIT - Neuroscience Information Toolbox
#  Copyright (C) 2020  Roberto Maffulli, Miguel Angel Casal Santiago
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import os
import io
import sys
import shutil
import logging
import argparse
import subprocess
import matlab.engine
from io import StringIO


def config_applogger(thelogger, modname):
    """Configure the logger for each specific module"""
    modname = {'app_name': modname}
    modlogger = logging.LoggerAdapter(thelogger, modname)
    return modlogger


def find_and_replace_in(templatefile, newfile, findstr, replacestr):
    replace = dict()
    # if findstr and replacestr are single strings and not lists of
    # strings turn them into lists of just one string for the next
    # for loop to work properly
    if type(findstr) is not list:
        findstr = [findstr]
    if type(replacestr) is not list:
        replacestr = [replacestr]
    for fs, rs in zip(findstr, replacestr):
        replace[fs] = rs
    with open(templatefile, 'r') as infile, open(newfile, 'w') as outfile:
        for line in infile:
            for src, target in replace.items():
                line = line.replace(src, target)
            outfile.write(line)


def add_to_startup(startupfile):
    # appends the configured startupfile to userpath/startup.m
    eng = matlab.engine.start_matlab("-nodesktop")
    userpath = eng.eval('userpath')
    eng.quit()
    stfile = open(startupfile, "r")
    stfilelines = stfile.read()
    stfile.close()
    modlogger.info("Appending " + os.path.split(startupfile)[-1] + " to: " + os.path.join(userpath, 'startupNIT.m'))
    stfileout = open(os.path.join(userpath, 'startupNIT.m'), "a")
    stfileout.write(stfilelines)
    stfileout.close()


def svm_bt(configdir, libsvmdir, mldir):
    """Configure and build the SVM module"""
    # configure local logger
    modlogger = config_applogger(logger, 'SVM BT')

    # generates config files based on local directories
    modlogger.info("Generating LIBSVM config files")
    makefile = os.path.join(libsvmdir, 'matlab/Makefile')
    startupfile = os.path.join(configdir, 'startup_svmlib.m')
    startupfiletemplate = os.path.join(configdir, 'startup_svmlib_TEMPLATE.m')
    if os.name == 'nt':
        libsvmmexdir = os.path.join(libsvmdir, 'windows')
    else:
        libsvmmexdir = os.path.join(libsvmdir, 'matlab')
    find_and_replace_in(startupfiletemplate, startupfile, ['$LIBSVM_MATLAB_ROOT', '$SML_SRC_ROOT'], [str(libsvmmexdir), str(mldir)])

    # append the configured startupfile to userpath/startup.m
    add_to_startup(startupfile)

    # compile LibSVM MEX library and test
    modlogger = config_applogger(logger, 'SVM BT - MATLAB')
    matout = io.StringIO()
    materr = io.StringIO()
    eng = matlab.engine.start_matlab("-nodesktop")
    eng.cd(libsvmmexdir, stdout=matout, stderr=materr)
    if os.name != 'nt':
        try:
            makefiletemplate = os.path.join(libsvmdir, 'matlab/Makefile_TEMPLATE')
            matlabexe = subprocess.run(['which', 'matlab'], stdout=subprocess.PIPE).stdout
            find_and_replace_in(makefiletemplate, makefile, '$MATLAB_EXE', str(matlabexe))
            eng.make(nargout=0, stdout=matout, stderr=materr)
            for theline in matout.getvalue().splitlines():
                modlogger.info(theline)
        except Exception as err:
            eng.quit()
            for theline in materr.getvalue().splitlines():
                modlogger.error(theline)
            quit()
    try:
        modlogger.info("Testing library:")
        matout.truncate(0)
        matout.seek(0)
        eng.svm_test(nargout=0, stdout=matout, stderr=materr)
        for theline in matout.getvalue().splitlines():
            modlogger.info(theline)
    except Exception as err:
            eng.quit()
            for theline in materr.getvalue().splitlines():
                modlogger.error(theline)
            quit()
    eng.quit()


def mi_binmethods_bt(configdir, srcdir, debugflag):
    """Configure and build the MI Binned Methods module"""
    # configure local logger
    modlogger = config_applogger(logger, 'MI Binned Methods BT')

    # generates config files based on local directories
    modlogger.info("Generating MI Binned Methods config files")
    makefiletemplate = os.path.join(configdir, 'make_MI_binmethods_TEMPLATE.m')
    startupfiletemplate = os.path.join(configdir, 'startup_MI_binmethods_TEMPLATE.m')
    makefile = os.path.join(configdir, 'make_MI_binmethods.m')
    startupfile = os.path.join(configdir, 'startup_MI_binmethods.m')
    find_and_replace_in(makefiletemplate, makefile,
                        ['$MI_binmethods_ROOT', '$MEX_debug'],
                        [str(srcdir), str(int(debugflag))])
    find_and_replace_in(startupfiletemplate, startupfile, '$MI_binmethods_ROOT', str(srcdir))
    add_to_startup(startupfile)

    # compile MI - Binned Methods MEX library and test
    modlogger = config_applogger(logger, 'MI Binned Methods BT - MATLAB')
    matout = io.StringIO()
    materr = io.StringIO()
    try:
        eng = matlab.engine.start_matlab("-nodesktop")
        eng.cd(configdir, stdout=matout, stderr=materr)
        eng.make_MI_binmethods(nargout=0, stdout=matout, stderr=materr)
        for theline in matout.getvalue().splitlines():
            modlogger.info(theline)
        modlogger.info("Testing MI Binned Methods library:")
        matout.truncate(0)
        matout.seek(0)
        eng.MI_binmethods_test(nargout=0, stdout=matout, stderr=materr)
        eng.quit()
        for theline in matout.getvalue().splitlines():
            modlogger.info(theline)
    except Exception:
        eng.quit()
        for theline in materr.getvalue().splitlines():
            modlogger.error(theline)
        quit()


def mi_conmethods_bt(configdir, srcdir, debugflag):
    """Configure and build the MI Continuous Methods module"""
    # configure local logger
    modlogger = config_applogger(logger, 'MI Continuous Methods BT')

    # generates config files based on local directories
    modlogger.info("Generating MI Continuous Methods config files")
    startupfiletemplate = os.path.join(configdir, 'startup_MI_conmethods_TEMPLATE.m')
    startupfile = os.path.join(configdir, 'startup_MI_conmethods.m')
    find_and_replace_in(startupfiletemplate, startupfile, '$MI_conmethods_ROOT', str(srcdir))

    # generates make files for NPC based on local directories
    makefiletemplate = os.path.join(configdir, 'make_MI_conmethods_NPC_TEMPLATE.m')
    makefile = os.path.join(configdir, 'make_MI_conmethods_NPC.m')
    find_and_replace_in(makefiletemplate, makefile,
                        ['$MI_conmethods_ROOT', '$MEX_debug'],
                        [str(srcdir), str(int(debugflag))])
    add_to_startup(startupfile)

    # test MI - Continuous Methods (MVC)
    modlogger = config_applogger(logger, 'Testing MI Continuous Methods (MVC) - MATLAB')
    matout = io.StringIO()
    materr = io.StringIO()
    try:
        eng = matlab.engine.start_matlab("-nodesktop")
        eng.cd(configdir, stdout=matout, stderr=materr)
        eng.MI_conmethods_MVC_test(nargout=0, stdout=matout, stderr=materr)
        eng.quit()
        for theline in matout.getvalue().splitlines():
            modlogger.info(theline)
    except Exception as err:
        eng.quit()
        for theline in materr.getvalue().splitlines():
            modlogger.error(theline)
        quit()

    # compile MI - Continuous Methods (NPC) MEX library and test
    modlogger = config_applogger(logger, 'MI Continuous Methods (NPC) BT - MATLAB')
    matout = io.StringIO()
    materr = io.StringIO()
    try:
        eng = matlab.engine.start_matlab("-nodesktop")
        eng.cd(configdir, stdout=matout, stderr=materr)
        eng.make_MI_conmethods_NPC(nargout=0, stdout=matout, stderr=materr)
        for theline in matout.getvalue().splitlines():
            modlogger.info(theline)
        modlogger.info("Testing MI Continuous Methods (NPC) library:")
        matout.truncate(0)
        matout.seek(0)
        eng.MI_conmethods_NPC_test(nargout=0, stdout=matout, stderr=materr)
        eng.quit()
        for theline in matout.getvalue().splitlines():
            modlogger.info(theline)
    except Exception:
        eng.quit()
        for theline in materr.getvalue().splitlines():
            modlogger.error(theline)
        quit()


def glm_bt(configdir, glmdir, mldir, debugflag):
    """Configure and build the GLM module"""
    # configure local logger
    modlogger = config_applogger(logger, 'GLM BT')

    # generates config files based on local directories
    modlogger.info("Generating GLMnet config files")
    startupfile = os.path.join(configdir, 'startup_glmnet.m')
    startupfiletemplate = os.path.join(configdir, 'startup_glmnet_TEMPLATE.m')
    find_and_replace_in(startupfiletemplate, startupfile, ['$GLMNET_MATLAB_ROOT', '$SML_SRC_ROOT'], [str(glmdir), str(mldir)])
    add_to_startup(startupfile)

    modlogger = config_applogger(logger, 'GLMnet BT - MATLAB')
    matout = io.StringIO()
    materr = io.StringIO()
    # generates make files for GLMnet based on local directories
    makefiletemplate = os.path.join(configdir, 'make_GLMnet_TEMPLATE.m')
    makefile = os.path.join(configdir, 'make_GLMnet.m')
    find_and_replace_in(makefiletemplate, makefile,
                        ['$GLMnet_ROOT', '$MEX_debug'],
                        [str(glmdir), str(int(debugflag))])

    # compile GLMnet MEX library  and test
    try:
        eng = matlab.engine.start_matlab("-nodesktop")
        eng.cd(glmdir, stdout=matout, stderr=materr)
        eng.cd(glmdir, stdout=matout, stderr=materr)
        eng.make_GLMnet(nargout=0, stdout=matout, stderr=materr)
        testflag = eng.workspace['testflag']
        for theline in matout.getvalue().splitlines():
            modlogger.info(theline)
        matout.truncate(0)
        matout.seek(0)
        if testflag:
            eng.startup
            eng.GLMnet_test(nargout=0, stdout=matout, stderr=materr)
            for theline in matout.getvalue().splitlines():
                modlogger.info(theline)
            eng.quit()
    except Exception as err:
        eng.quit()
        for theline in materr.getvalue().splitlines():
            modlogger.error(theline)
        quit()


def pid_bt(configdir, srcdir, brojadir):
    """Configure and build the Intersection Information module"""
    # configure local logger
    modlogger = config_applogger(logger, 'PID:II&FIT BT')

    # install dependencies
    modlogger.info("Installing python dependencies")
    modlogger = config_applogger(logger, 'PID:II&FIT BT - PIP')
    ii_pip = subprocess.Popen(["pip", "install", "ecos"],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    output, error = ii_pip.communicate()
    if error:
        for theline in StringIO(error):
            modlogger.error(theline.rstrip())
        quit()
    for theline in StringIO(output):
        modlogger.info(theline.rstrip())

    # generates config files based on local directories
    modlogger = config_applogger(logger, 'PID:II&FIT BT')
    modlogger.info("Generating PID:II&FIT config files")
    startupfile = os.path.join(configdir, 'startup_PID.m')
    startupfiletemplate = os.path.join(configdir, 'startup_PID_TEMPLATE.m')
    find_and_replace_in(startupfiletemplate, startupfile,
                        ['$PID_MATLAB_ROOT', '$II_PYTHON_INTERFACE_ROOT', '$PID_BROJA_2PID_ROOT'],
                        [srcdir, os.path.join(srcdir, 'PIDPyInterface'), brojadir])
    add_to_startup(startupfile)

    # adds correct path in src files based on local directories
    pidpyinterfacedir = os.path.join(srcdir, 'PIDPyInterface')
    pyinterfacefiletemplate = os.path.join(pidpyinterfacedir, 'ComputePID_TEMPLATE.py')
    pyinterfacefile = os.path.join(pidpyinterfacedir, 'ComputePID.py')
    find_and_replace_in(pyinterfacefiletemplate, pyinterfacefile, '$BROJA_2PID_PATH', brojadir)


    # test II module
    modlogger = config_applogger(logger, 'PID: II BT - MATLAB')
    matout = io.StringIO()
    materr = io.StringIO()
    try:
        eng = matlab.engine.start_matlab("-nodesktop")
        eng.cd(configdir, stdout=matout, stderr=materr)
        eng.PID_II_test(nargout=0, stdout=matout, stderr=materr)
        for theline in matout.getvalue().splitlines():
            modlogger.info(theline)
        eng.quit()
    except Exception as err:
        eng.quit()
        for theline in materr.getvalue().splitlines():
            modlogger.error(theline)
        quit()

    # test PID module
    modlogger = config_applogger(logger, 'PID: PID BT - MATLAB')
    matout = io.StringIO()
    materr = io.StringIO()
    try:
        eng = matlab.engine.start_matlab("-nodesktop")
        eng.cd(configdir, stdout=matout, stderr=materr)
        eng.PID_PID_test(nargout=0, stdout=matout, stderr=materr)
        for theline in matout.getvalue().splitlines():
            modlogger.info(theline)
        eng.quit()
    except Exception as err:
        eng.quit()
        for theline in materr.getvalue().splitlines():
            modlogger.error(theline)
        quit()

    # # test FIT module
    # modlogger = config_applogger(logger, 'PID: FIT BT - MATLAB')
    # matout = io.StringIO()
    # materr = io.StringIO()
    # try:
    #     eng = matlab.engine.start_matlab("-nodesktop")
    #     eng.cd(configdir, stdout=matout, stderr=materr)
    #     eng.PID_FIT_test(nargout=0, stdout=matout, stderr=materr)
    #     for theline in matout.getvalue().splitlines():
    #         modlogger.info(theline)
    #     eng.quit()
    # except Exception as err:
    #     eng.quit()
    #     for theline in materr.getvalue().splitlines():
    #         modlogger.error(theline)
    #     quit()


def tools_bt(configdir, srcdir):
    """Builds and tests the tools module"""
    # configure local logger
    modlogger = config_applogger(logger, 'Tools BT')

    # generates config files based on local directories
    modlogger.info("Generating Tools config files")
    startupfiletemplate = os.path.join(configdir, 'startup_tools_TEMPLATE.m')
    startupfile = os.path.join(configdir, 'startup_tools.m')
    find_and_replace_in(startupfiletemplate, startupfile, '$TOOLS_SRC_ROOT', str(srcdir))
    add_to_startup(startupfile)


def uml_bt(configdir, srcdir):
    """Builds and tests the PCA module"""
    # configure local logger
    modlogger = config_applogger(logger, 'UML BT')

    # generates config files based on local directories
    modlogger.info("Generating UML config files")
    startupfiletemplate = os.path.join(configdir, 'startup_uml_TEMPLATE.m')
    startupfile = os.path.join(configdir, 'startup_uml.m')
    find_and_replace_in(startupfiletemplate, startupfile, '$UML_SRC_ROOT', str(srcdir))
    add_to_startup(startupfile)


def main():

    # configure local logger
    modlogger = config_applogger(logger, "BuildAndTest")
    
    # configuration
    parser = argparse.ArgumentParser(description='Build script for the Neural Information Toolbox software')
    parser.add_argument('-d', '--dir', help='path to the root directory of the toolbox (mandatory argument)', type=str,
                        required=True)
    parser.add_argument('--debug', help='flag to compile MATLAB mex files in debug mode (optional argument)', type=bool,
                        required=False, default=False)
    parser.add_argument('-m', '--modules', nargs='+',
                        help='list of modules to be included in the toolbox installation (optional argument)', type=str,
                        required=False, default=["all"])
    args = parser.parse_args()
    modlogger.info("Root directory: %s" % args.dir)

    # defines directories for system set-up
    modlogger.info("Define local directories:")
    rootdir = os.path.abspath(args.dir)
    if not os.path.isdir(rootdir):
        raise Exception('Specified path for the root directory -d/--dir is non-existing. Aborting.')

    configdir = os.path.join(rootdir, 'config')
    externdir = os.path.join(rootdir, 'extern')
    srcdir = os.path.join(rootdir, 'src')
    modlogger.info(" - Configuration directory: " + configdir)
    modlogger.info(" - Extern directory:        " + externdir)
    modlogger.info(" - Source directory:        " + srcdir)

    
    # define modules to be compiled/set-up
    allowed_modules = ["all", "svm", "mi_binmethods", "mi_conmethods", "glm", "pid", "tools", "uml"]
    for mod in args.modules:
        if mod not in allowed_modules:
            raise Exception("Module argument " + mod + " not recognized. Allowed module keywords are: '"
                            + "' '".join(map(str, allowed_modules)) + "'. Aborting.")
    buildall = buildsvm = buildmi_binned = buildmi_cont = buildglm = buildpid = buildtools = builduml = False
    if "all" in args.modules:
        buildall = True
    else:
        if "svm" in args.modules:
            buildsvm = True
        if "mi_binmethods" in args.modules:
            buildmi_binned = True
            if not buildtools:
                modlogger.info("Activating Tools module as it is required by MI binned methods module.")
                buildtools = True
        if "mi_conmethods" in args.modules:
            buildmi_cont = True
        if "glm" in args.modules:
            buildglm = True
        if "pid" in args.modules:
            buildpid = True
            if not buildtools:
                modlogger.info("Activating Tools module as it is required by II module.")
                buildtools = True
        if "tools" in args.modules:
            buildtools = True
        if "uml" in args.modules:
            builduml = True
    
    # assert python version
    if sys.version_info < (3, 6, 0):
        raise Exception("Minimum Python version required is 3.6")

    # start adding requirements to matlab startup
    eng = matlab.engine.start_matlab("-nodesktop")
    userpath = eng.eval('userpath')
    eng.quit()
    files = os.listdir(userpath)

    if 'startupNIT.m' in files:  # check if startupNIT exists
        os.remove(os.path.join(userpath, 'startupNIT.m'))  # remove startupNIT
    stfileout = open(os.path.join(userpath, 'startupNIT.m'), "a")
    stfileout.write("%% execute requirements of the nit toolbox \n"
        "addpath(genpath(\"" + configdir + "\"))\n")
    stfileout.close()
    stfileout = open(os.path.join(userpath, 'startup.m'), "a")  # add startupNIT to matlab startup
    stfileoutread = open(os.path.join(userpath, 'startup.m'), "r")
    stfileoutlines = stfileoutread.readlines()
    if not "startupNIT \n" in stfileoutlines:
        stfileout.write("\n%% execute requirements of the nit toolbox \n"
            "startupNIT \n")
    stfileout.close()

    # build and test tools module
    if buildall or buildtools:
        toolsdir = os.path.join(srcdir, 'tools')
        tools_bt(configdir, toolsdir)

    # build and test uml module
    if buildall or builduml:
        umldir = os.path.join(srcdir, 'UML')
        uml_bt(configdir, umldir)

    # build and test MI - BinnedMethods
    if buildall or buildmi_binned:
        srcdirmi = os.path.join(srcdir, 'MI', 'BinnedMethods')
        mi_binmethods_bt(configdir, srcdirmi, args.debug)

    # build and test MI - ContinuousMethods
    if buildall or buildmi_cont:
        srcdirmi = os.path.join(srcdir, 'MI', 'ContinuousMethods')
        mi_conmethods_bt(configdir, srcdirmi, args.debug)

    # build and test svm module
    if buildall or buildsvm:
        libsvmdir = os.path.join(externdir, 'libsvm-3.24')
        mldir = os.path.join(srcdir,'SML')
        svm_bt(configdir, libsvmdir, mldir)

    # build and test glm module
    if buildall or buildglm:
        glmdir = os.path.join(externdir, 'glmnet_matlab')
        mldir = os.path.join(srcdir,'SML')
        glm_bt(configdir, glmdir, mldir, args.debug)

    # build and test PID module
    if buildall or buildpid:
        srcdirpid = os.path.join(srcdir, 'PID')
        brojadir = os.path.join(externdir, 'BROJA_2PID')
        # escape the backslash of the dir delimiter in windows as python complains (rightfully)
        if os.name == 'nt':
            brojadir = brojadir.replace('\\', '\\\\')
        pid_bt(configdir, srcdirpid, brojadir)

    # ending adding requirements to the startup file
    eng = matlab.engine.start_matlab("-nodesktop")
    userpath = eng.eval('userpath')
    eng.quit()
    stfileout = open(os.path.join(userpath, 'startupNIT.m'), "a")
    stfileout.write("rmpath(genpath(\""+ configdir + "\"))\n"
        "%% finish adding nit requirements\n")
    stfileout.close()
    stfileout = open(os.path.join(userpath, 'startup.m'), "a")
    stfileoutread = open(os.path.join(userpath, 'startup.m'), "r")
    stfileoutlines = stfileoutread.readlines()
    if not "%% finish adding nit requirements\n" in stfileoutlines:
        stfileout.write("%% finish adding nit requirements\n")
    stfileout.close()


if __name__ == "__main__":
    
    # configure logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(app_name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    modlogger = config_applogger(logger, 'BuildAndTest')
    try:
        main()
    except Exception as e:
        modlogger.exception(e)


