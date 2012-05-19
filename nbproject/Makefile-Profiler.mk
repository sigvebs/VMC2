#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=mpicxx
CXX=mpicxx
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Profiler
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/QD/QDJastrow.o \
	${OBJECTDIR}/includes/zigrandom.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/QD/QDOrbital.o \
	${OBJECTDIR}/OneBodyDensity.o \
	${OBJECTDIR}/Slater.o \
	${OBJECTDIR}/Hamiltonian.o \
	${OBJECTDIR}/includes/lib.o \
	${OBJECTDIR}/SGD.o \
	${OBJECTDIR}/DMC.o \
	${OBJECTDIR}/includes/zignor.o \
	${OBJECTDIR}/Blocking.o \
	${OBJECTDIR}/QD/QDHamiltonian.o \
	${OBJECTDIR}/Hermite.o \
	${OBJECTDIR}/Jastrow.o \
	${OBJECTDIR}/ComputeOneWF.o \
	${OBJECTDIR}/Orbital.o \
	${OBJECTDIR}/includes/ziggurat.o \
	${OBJECTDIR}/includes/ini.o \
	${OBJECTDIR}/WaveFunction.o \
	${OBJECTDIR}/ComputeGrid.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-pg -O3 -DMPICH_IGNORE_CXX_SEEK -llapack -lblas -larmadillo
CXXFLAGS=-pg -O3 -DMPICH_IGNORE_CXX_SEEK -llapack -lblas -larmadillo

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-larmadillo

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${TESTDIR}/TestFiles/f2

${TESTDIR}/TestFiles/f2: ${OBJECTFILES}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc} -o ${TESTDIR}/TestFiles/f2 ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/QD/QDJastrow.o: QD/QDJastrow.cpp 
	${MKDIR} -p ${OBJECTDIR}/QD
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/QD/QDJastrow.o QD/QDJastrow.cpp

${OBJECTDIR}/includes/zigrandom.o: includes/zigrandom.c 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.c) -g -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/zigrandom.o includes/zigrandom.c

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/QD/QDOrbital.o: QD/QDOrbital.cpp 
	${MKDIR} -p ${OBJECTDIR}/QD
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/QD/QDOrbital.o QD/QDOrbital.cpp

${OBJECTDIR}/OneBodyDensity.o: OneBodyDensity.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/OneBodyDensity.o OneBodyDensity.cpp

${OBJECTDIR}/Slater.o: Slater.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Slater.o Slater.cpp

${OBJECTDIR}/Hamiltonian.o: Hamiltonian.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Hamiltonian.o Hamiltonian.cpp

${OBJECTDIR}/includes/lib.o: includes/lib.cpp 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/lib.o includes/lib.cpp

${OBJECTDIR}/SGD.o: SGD.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/SGD.o SGD.cpp

${OBJECTDIR}/DMC.o: DMC.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/DMC.o DMC.cpp

${OBJECTDIR}/includes/zignor.o: includes/zignor.c 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.c) -g -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/zignor.o includes/zignor.c

${OBJECTDIR}/Blocking.o: Blocking.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Blocking.o Blocking.cpp

${OBJECTDIR}/QD/QDHamiltonian.o: QD/QDHamiltonian.cpp 
	${MKDIR} -p ${OBJECTDIR}/QD
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/QD/QDHamiltonian.o QD/QDHamiltonian.cpp

${OBJECTDIR}/Hermite.o: Hermite.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Hermite.o Hermite.cpp

${OBJECTDIR}/Jastrow.o: Jastrow.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Jastrow.o Jastrow.cpp

${OBJECTDIR}/ComputeOneWF.o: ComputeOneWF.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/ComputeOneWF.o ComputeOneWF.cpp

${OBJECTDIR}/Orbital.o: Orbital.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Orbital.o Orbital.cpp

${OBJECTDIR}/includes/ziggurat.o: includes/ziggurat.cpp 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/ziggurat.o includes/ziggurat.cpp

${OBJECTDIR}/includes/ini.o: includes/ini.cpp 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/ini.o includes/ini.cpp

${OBJECTDIR}/WaveFunction.o: WaveFunction.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/WaveFunction.o WaveFunction.cpp

${OBJECTDIR}/ComputeGrid.o: ComputeGrid.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/ComputeGrid.o ComputeGrid.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${TESTDIR}/TestFiles/f2

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
