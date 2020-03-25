import os

os.system('rm -f *.o')
os.system('rm -f *.out')

CPC='mpiicpc'
CPC_FLAGS='-I/root/PetscGit/include -I/root/intel/impi/2018.3.222/include64 -I/root/intel/mkl/include -O3 -static-intel -qopenmp -DAdd_ -Wno-unknown-pragmas'
LIB_FLAGS='${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_cdft_core.a ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a /root/intel/compilers_and_libraries_2018/linux/lib/intel64_lin/libirc.a /root/intel/compilers_and_libraries_2018/linux/lib/intel64_lin/libiomp5.a /root/PetscGit/lib/libmumps_common.a /root/PetscGit/lib/libpetsc.a /root/PetscGit/lib/libptscotch.a /root/PetscGit/lib/libscotcherrexit.a /root/PetscGit/lib/libptscotcherr.a /root/PetscGit/lib/libscotcherr.a /root/PetscGit/lib/libscotchmetis.a /root/PetscGit/lib/libdmumps.a /root/PetscGit/lib/libptscotchparmetis.a /root/PetscGit/lib/libdmumps.a /root/PetscGit/lib/libpord.a /root/PetscGit/lib/libmetis.a /root/PetscGit/lib/libsuperlu_dist.a  /root/PetscGit/lib/libparmetis.a /root/PetscGit/lib/libscalapack.a /root/PetscGit/lib/libesmumps.a /root/intel/impi/2018.3.222/lib64/libmpi.a /root/PetscGit/lib/libscotch.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -lX11 -lz -O3 -static-intel -qopenmp  -lifcore -limf  -wd10237'



Source_Files=['Globals.cpp', 'Petsc_Functions.cpp', 'Data_Input.cpp', 'Minor_NR.cpp', 'PrimeSim.cpp', 'PrimeSim2.cpp', 'Init_Calcs.cpp', 'MixedFunctions.cpp', 'Do_Cycle.cpp', 'MKL_Solver.cpp', 'ViewAssist.cpp', 'Do_Time_Step.cpp', 'Mass_Tran.cpp', 'Numeric.cpp', 'MIfstream.cpp', 'PetscAUX.cpp']

Outname=''
for FName in Source_Files:
	Compiler_Run=CPC+' -c '+FName+' '+CPC_FLAGS
	#print(Compiler_Run)
	os.system(Compiler_Run)	
	j=0
	FLen=len(FName)-3
		
	for ch in FName:
		if j<FLen:
			Outname=Outname+ch
		j=j+1
	Outname=Outname+'o '

Linker_Run=CPC+' -o osiris.out'+' '+Outname+' '+LIB_FLAGS
#print(Linker_Run)
os.system(Linker_Run)
