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
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/QD/QDJastrow.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/QD/QDOrbital.o \
	${OBJECTDIR}/includes/zigrandom.o \
	${OBJECTDIR}/Slater.o \
	${OBJECTDIR}/Hamiltonian.o \
	${OBJECTDIR}/includes/lib.o \
	${OBJECTDIR}/includes/zignor.o \
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
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/vmc2

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/vmc2: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/vmc2 ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/QD/QDJastrow.o: QD/QDJastrow.cpp 
	${MKDIR} -p ${OBJECTDIR}/QD
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/QD/QDJastrow.o QD/QDJastrow.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/QD/QDOrbital.o: QD/QDOrbital.cpp 
	${MKDIR} -p ${OBJECTDIR}/QD
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/QD/QDOrbital.o QD/QDOrbital.cpp

${OBJECTDIR}/includes/zigrandom.o: includes/zigrandom.c 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.c) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/zigrandom.o includes/zigrandom.c

${OBJECTDIR}/Slater.o: Slater.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/Slater.o Slater.cpp

${OBJECTDIR}/Hamiltonian.o: Hamiltonian.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/Hamiltonian.o Hamiltonian.cpp

${OBJECTDIR}/includes/lib.o: includes/lib.cpp 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/lib.o includes/lib.cpp

${OBJECTDIR}/includes/zignor.o: includes/zignor.c 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.c) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/zignor.o includes/zignor.c

${OBJECTDIR}/QD/QDHamiltonian.o: QD/QDHamiltonian.cpp 
	${MKDIR} -p ${OBJECTDIR}/QD
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/QD/QDHamiltonian.o QD/QDHamiltonian.cpp

${OBJECTDIR}/Hermite.o: Hermite.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/Hermite.o Hermite.cpp

${OBJECTDIR}/Jastrow.o: Jastrow.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/Jastrow.o Jastrow.cpp

${OBJECTDIR}/ComputeOneWF.o: ComputeOneWF.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/ComputeOneWF.o ComputeOneWF.cpp

${OBJECTDIR}/Orbital.o: Orbital.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/Orbital.o Orbital.cpp

${OBJECTDIR}/includes/ziggurat.o: includes/ziggurat.cpp 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/ziggurat.o includes/ziggurat.cpp

${OBJECTDIR}/includes/ini.o: includes/ini.cpp 
	${MKDIR} -p ${OBJECTDIR}/includes
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/includes/ini.o includes/ini.cpp

${OBJECTDIR}/WaveFunction.o: WaveFunction.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/WaveFunction.o WaveFunction.cpp

${OBJECTDIR}/ComputeGrid.o: ComputeGrid.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/ComputeGrid.o ComputeGrid.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/vmc2

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
