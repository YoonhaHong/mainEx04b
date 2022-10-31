# PYTHIA configuration file.
# Generated on 2020. 11. 07. (토) 23:33:35 KST with the user supplied options:
# --prefix=/opt/install/pythia8303
# --with-hepmc2=/opt/install/HepMC-2.06.11
# --with-fastjet3=/opt/install/fastjet-3.3.4
# --with-root=/opt/install/root-6.18.04
# --with-hepmc3=/opt/install/HepMC3-3.2.2
# --with-lhapdf6=/opt/install/LHAPDF-6.3.0

# Install directory prefixes.
PREFIX_BIN=/opt/install/pythia8303/bin
PREFIX_INCLUDE=/opt/install/pythia8303/include
PREFIX_LIB=/opt/install/pythia8303/lib
PREFIX_SHARE=/opt/install/pythia8303/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/opt/install/fastjet-3.3.4/bin/
FASTJET3_INCLUDE=-I/opt/install/fastjet-3.3.4/include
FASTJET3_LIB=-L/opt/install/fastjet-3.3.4/lib -Wl,-rpath,/opt/install/fastjet-3.3.4/lib -lfastjet

HEPMC2_USE=true
HEPMC2_CONFIG=
HEPMC2_BIN=/opt/install/HepMC-2.06.11/
HEPMC2_INCLUDE=-I/opt/install/HepMC-2.06.11/include
HEPMC2_LIB=-L/opt/install/HepMC-2.06.11/lib -Wl,-rpath,/opt/install/HepMC-2.06.11/lib -lHepMC

HEPMC3_USE=true
HEPMC3_CONFIG=HepMC3-config
HEPMC3_BIN=/opt/install/HepMC3-3.2.2/bin/
HEPMC3_INCLUDE=-I/opt/install/HepMC3-3.2.2/include
HEPMC3_LIB=-L/opt/install/HepMC3-3.2.2/lib -Wl,-rpath,/opt/install/HepMC3-3.2.2/lib -lHepMC3

LHAPDF5_USE=false
LHAPDF5_CONFIG=HepMC3-config
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/opt/install/LHAPDF-6.3.0/bin/
LHAPDF6_INCLUDE=-I/opt/install/LHAPDF-6.3.0/include
LHAPDF6_LIB=-L/opt/install/LHAPDF-6.3.0/lib -Wl,-rpath,/opt/install/LHAPDF-6.3.0/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=lhapdf-config
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=lhapdf-config
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/opt/install/root-6.18.04-20210130/bin/
ROOT_INCLUDE=-I/opt/install/root-6.18.04-20210130/include
ROOT_LIB=-L/opt/install/root-6.18.04-20210130/lib -Wl,-rpath,/opt/install/root-6.18.04-20210130/lib -lCore

YODA_USE=false
YODA_CONFIG=root-config
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

GZIP_USE=false
GZIP_CONFIG=root-config
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=root-config
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=root-config
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=root-config
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=
