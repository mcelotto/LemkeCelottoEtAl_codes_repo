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
import re
import sys
import git
import shutil
import logging
import argparse
from datetime import datetime

# name of temporary directory for documentation
docsdirname = 'TMP_Documentation'
# set of delimiters that this script is going to use to extract comments from source files
comments_delim = ['%%%']
# list of source files extensions to be considered when creating the documentation
src_files_extensions = ['.m']


def create_docs_from_src(srcfilepath, docsfilepath, identifiers):
    """ The function takes all lines starting with any of identifiers from srcfile and copies them into a docsfile """
    createfile = False  # create doc file only if the corresponding source has comments lines to be inserted
    srcfile = open(srcfilepath, "r")
    for line in srcfile.readlines():
        for ii in identifiers:
            if line.startswith(ii):
                if not createfile:
                    os.makedirs(os.path.dirname(docsfilepath), exist_ok=True)
                    docsfile = open(docsfilepath, "w")
                    createfile = True
                targetregex = '^' + ii
                rgx = re.compile(targetregex)
                docline = rgx.sub('', line)
                docsfile.write(docline)
    srcfile.close()
    if createfile:
        docsfile.close()
    return createfile


def build_sidebar(root_docs_path):
    sidebarfile = open(os.path.join(root_docs_path, "_Sidebar.md"), "w")
    rgx = re.compile('.md$')
    for root, dirs, files in os.walk(root_docs_path):
        level = root.replace(root_docs_path, '').count(os.sep) - 1
        indent = "\t" * level
        subindent = "\t" * (level + 1)
        dirname = os.path.basename(root)
        # handle directories (skipping the home directory)
        if dirname != docsdirname:
            # dirs that do not have necessarily a .md file associated
            if not os.path.exists(os.path.join(root, os.path.basename(root) + '.md')):
                sidebarfile.write('{}* {}\n'.format(indent, os.path.basename(root)))
            # dirs that contain a .md file with their same name inside them
            else:
                tocentry = os.path.basename(root)
                filepath = os.path.join(os.path.relpath(root, root_docs_path), os.path.basename(root))
                sidebarfile.write('{}* [{}]({})\n'.format(indent, tocentry, filepath))
        for f in files:
            # build TOC for existing markdown files but exclude the _Sidebar.md itself and files that have the same name
            # as the directory they are in (already included in the parent indent level)
            filename = os.path.basename(f)
            tocentry = rgx.sub('', filename)#.capitalize()
            filepath = os.path.join(os.path.relpath(root, root_docs_path), tocentry)
            if filename != "_Sidebar.md":
                if tocentry != dirname:
                    sidebarfile.write('{}* [{}]({})\n'.format(subindent, tocentry, filepath))
    sidebarfile.close()


def main():
    # configuration
    parser = argparse.ArgumentParser(description='Builds the documentation for the Neural Information Toolbox software '
                                                 'by parsing the project directories and pushing the resulting docs '
                                                 'to the wiki repository for the gitlab project.\n'
                                                 'The script builds a documentation tree having the same structure as '
                                                 'the directory tree in the project source. The *.md files that are '
                                                 'found in the directories are copied straight to the documentation '
                                                 'repository, while the comments preceded by \'%%%\' in *.m files '
                                                 'are stripped and used to generate the documentation of the '
                                                 'corresponding Matlab functions.\n'
                                                 'Files matching regular expressions contained in .docsignore files '
                                                 'that can be present in each directory, are excluded from the '
                                                 'documentation.')
    parser.add_argument('-r', '--repo', help='git SSH url for the project\'s Wiki (mandatory argument)',
                        type=str, required=True)
    args = parser.parse_args()
    repurl = args.repo
    # check that remote repository exists and that user has access rights
    git.cmd.Git().ls_remote(repurl)
    logger.info("git SSH url for the project\'s Wiki: %s" % repurl)

    # remove docs dir if already present and create a blank documentation dir
    docsdir = os.path.join(os.getcwd(), docsdirname)
    if os.path.exists(docsdir):
        shutil.rmtree(docsdir)
    os.mkdir(docsdir)

    # loop over directories to find markdown and source files
    nitrootdir = os.getcwd()
    srcdir = os.path.join(nitrootdir, 'src')
    dirs = [srcdir]
    markdown_files = []
    ignorelist = []
    while dirs:
        dirname = dirs.pop()
        # read .docsignore file if present in current dir and store entries
        docsignore = os.path.join(dirname,'.docsignore')
        if os.path.isfile(docsignore):
            ignorelist = [line.rstrip('\n') for line in open(docsignore)]
            # escape * character in regexps (fixes a bur in python regexp
            # see: https://stackoverflow.com/questions/3675144/regex-error-nothing-to-repeat)
            ignorelist = [il.replace('*', '\\\*') for il in ignorelist]
            ignoreregexps = "|".join(ignorelist)
        # scan directories recursively
        with os.scandir(dirname) as files:
            for entry in files:  # loop through the folder
                # markdown files are copied straight to the documentation
                if entry.is_dir() and entry.name != docsdirname:
                    dirs.append(entry.path)
                # write documentation if ignorelist exists and file is not in ignore list or
                # if file has extension .md if ignorelist is not existing
                elif (entry.name.endswith('.md') and ignorelist and not re.search(ignoreregexps, entry.name)) or \
                        (entry.name.endswith('.md') and not ignorelist):
                    relativepath = os.path.relpath(entry.path, srcdir)
                    absolutedocpath = os.path.join(docsdir, relativepath)
                    os.makedirs(os.path.dirname(absolutedocpath), exist_ok=True)
                    shutil.copyfile(entry.path, absolutedocpath)
                    markdown_files.append(absolutedocpath)
                # source files are processed to strip the relevant lines to be added to the documentation
                else:
                    for ex in src_files_extensions:
                        # write documentation if ignorelist exists and file is not in ignore list or
                        # if file has extension ex if ignorelist is not existing
                        if entry.name.endswith(ex) and ignorelist and not re.search(ignoreregexps, entry.name) or \
                                (entry.name.endswith(ex) and not ignorelist):
                            # code_files.append(entry.name)
                            rgx = re.compile(ex + '$')
                            relativepath = os.path.relpath(entry.path, srcdir)
                            docsfilename = rgx.sub('', relativepath) + '.md'
                            absolutedocpath = os.path.join(docsdir, docsfilename)
                            if create_docs_from_src(entry.path, absolutedocpath, comments_delim):
                                markdown_files.append(absolutedocpath)
        # empty ignorelist for next directory
        if not ignorelist:
            ignorelist = []

    # copy separately the markdown files that are in the nitrootdir directory
    for file in os.listdir(nitrootdir):
        filename = os.fsdecode(file)
        if filename.endswith(".md"):
            # treat separately the readme file which becomes the Wiki homepage
            if filename == 'README.md':
                shutil.copyfile(os.path.join(nitrootdir, filename), os.path.join(docsdir, 'home.md'))
            else:
                shutil.copyfile(os.path.join(nitrootdir, filename), os.path.join(docsdir, filename))

    # create _Sidebar.md file for TOC
    build_sidebar(docsdir)
    # push changes to the Wiki repository
    os.chdir(docsdirname)
    wiki_rep = git.Repo.init()
    wiki_rep.git.add(A=True)
    when = datetime.now().strftime("%B %d, %Y - %H:%M:%S")
    wiki_rep.git.commit('-m', 'Built the documentation on ' + when)
    remote = wiki_rep.create_remote('NIT-Wiki', url=repurl)
    remote.fetch()
    #remote.push(u=True,force=True)
    wiki_rep.git.push('-u', 'NIT-Wiki', 'master','--force')
    # cleanup of temporary documentation directory
    shutil.rmtree(docsdir)


if __name__ == "__main__":
    # configure logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - BuildDocumentation - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    try:
        main()
    except Exception as e:
        logger.error(e)
        logger.exception(e)
