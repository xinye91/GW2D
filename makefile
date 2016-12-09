FCM = mpiifort


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
all: petsc2D.o linkVdataFromGD.o flow.o gd.o displayMod.o vdata.o ncTable.o proc_mod.o pro.f90 chkopts
	-${FLINKER} -o gdNCp.out petsc2D.o linkVdataFromGD.o flow.o gd.o displayMod.o vdata.o ncTable.o proc_mod.o pro.f90 -I/${NFDIR}/include -L/${NFDIR}/lib -lnetcdff -fpp ${PETSC_KSP_LIB}
petsc2D.o : petsc2Dsolver.f90 proc_mod.o
	-${FLINKER} -o petsc2D.o -c petsc2Dsolver.f90 proc_mod.o -fpp ${PETSC_FC_INCLUDES} -L${PETSC_LIB_DIR} ${PETSC_KSP_LIB_BASIC}
linkVdataFromGD.o : linkVdataFromGd.f90 gd.o vdata.o
	${FCM} -o linkVdataFromGD.o -c linkVdataFromGd.f90 -fpp
flow.o : flow.f90 displayMod.o vdata.o proc_mod.o
	${FCM} -o flow.o -c flow.f90 displayMod.o vdata.o proc_mod.o -fpp
vdata.o : vdata.f90 displayMod.o gd.o
	${FCM} -o vdata.o -c vdata.f90 gd.o -fpp
proc_mod.o : proc_mod.f90 gd.o
	${FCM} -o proc_mod.o -c proc_mod.f90 gd.o -fpp 
ncTable.o : ncTableHierarchy.f90 gd.o
	${FCM} -o ncTable.o -c ncTableHierarchy.f90 gd.o -I/${NFDIR}/include -L/${NFDIR}/lib -lnetcdff -fpp 
gd.o : gd.f90
	${FCM} -o gd.o -c gd.f90 -fpp
displayMod.o : ch_assert.f90
	${FCM} -o displayMod.o -c ch_assert.f90 -fpp


#include ${PETSC_DIR}/lib/petsc/conf/test
