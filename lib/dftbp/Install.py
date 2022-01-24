#!/usr/bin/env python

"""
Install.py tool to download, unpack, build, and link to the DFTBP library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, subprocess, shutil, tarfile
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, checkmd5sum

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

version = '1.2.1'
suffix = 'gfortran'
#gfortran was not tested with DFTBP 


# known checksums for different DFTBP versions. used to validate the download.
checksums = { \
        '1.1.0' : '533635721ee222d0ed2925a18fb5b294', \
        '1.2.0' : '68bf0db879da5e068a71281020239ae7', \
        '1.2.1' : '85ac414fdada2d04619c8f936344df14', \
        }

# help message

HELP = """
Syntax from src dir: make lib-dftbp args="-b"
                 or: make lib-dftbp args="-p /usr/local/dftbp"
                 or: make lib-dftbp args="-m gfortran"
                 or: make lib-dftbp args="-b -v 1.2.1"

Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/dftbp
                 or: python Install.py -m gfortran
                 or: python Install.py -v 1.2.1 -b

Example:

make lib-dftbp args="-b -m gfortran"   # download/build in lib/dftbp
make lib-dftbp args="-p $HOME/dftbp"   # use existing DFTBP installation
"""


pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the DFTBP library")
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing DFTBP installation")
parser.add_argument("-m", "--machine", choices=['gfortran', 'ifort', 'linalg', 'serial', 'mpi'],
                    help="suffix of a Makefile.lammps.* file used for linking LAMMPS with this library")
parser.add_argument("-v", "--version", default=version,
                    help="set version of DFTBP to download and build (default: %s)" % version)

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path:
  parser.print_help()
  sys.exit(HELP)

homepath = fullpath(".")

buildflag = args.build
pathflag = args.path is not None
version = args.version
suffixflag = args.machine is not None
suffix = args.machine

if pathflag:
  dftbpdir = args.path
  if not os.path.isdir(dftbpdir):
    sys.exit("DFTBP path %s does not exist" % dftbpdir)
  dftbpdir = fullpath(dftbpdir)

homedir = "DFTBP-%s" % version

if buildflag:
  url = "https://github.com/gsjung0419/DFTBP/archive/v%s.tar.gz" % version
  dftbppath = fullpath(homepath)
  dftbpdir = os.path.join(dftbppath, homedir)

# download and unpack DFTBP tarball

if buildflag:
  print("Downloading DFTBP ...")
  geturl(url, "DFTBP.tar.gz")

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version], 'DFTBP.tar.gz'):
      sys.exit("Checksum for DFTBP library does not match")

  print("Unpacking DFTBP ...")
  if os.path.exists(dftbpdir):
    shutil.rmtree(dftbpdir)
  if tarfile.is_tarfile('DFTBP.tar.gz'):
    tgz = tarfile.open('DFTBP.tar.gz')
    tgz.extractall()
    os.remove('DFTBP.tar.gz')
  else:
    sys.exit("File DFTBP.tar.gz is not a supported archive")

  # build DFTBP
  print("Building DFTBP ...")
  cmd = 'cd "%s"; make' % dftbpdir
  try:
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    sys.exit("Make failed with:\n %s" % e.output.decode('UTF-8'))

# create 3 links in lib/dftbp to DFTBP dirs

print("Creating links to DFTBP files")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
if os.path.isfile("filelink.o") or os.path.islink("filelink.o"):
  os.remove("filelink.o")
os.symlink(os.path.join(dftbpdir, 'src'), 'includelink')
os.symlink(dftbpdir, 'liblink')
os.symlink(os.path.join(dftbpdir, 'src', 'dftbp_c_bind.o'), 'filelink.o')

# copy Makefile.lammps.suffix to Makefile.lammps

if suffixflag or not os.path.exists("Makefile.lammps"):
  if suffix is None:
    suffix = 'gfortran'
  print("Creating Makefile.lammps")
  if os.path.exists("Makefile.lammps.%s" % suffix):
    shutil.copyfile("Makefile.lammps.%s" % suffix, 'Makefile.lammps')
