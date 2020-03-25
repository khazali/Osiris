#include "Functions.h"
#include <petscsnes.h>

#include <petscsnes.h>
#include <petscsys.h>
#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <HYPRE_config.h>
//#include <HYPRE_utilities.h>

PetscErrorCode SetUpPetscEnv(int *argc, char ***argv, int ProblemSize, int BlockSize) {
	int i;

	PetscOptionsSetValue(NULL, "-sub_pc_type", "lu");
	PetscOptionsSetValue(NULL, "-sub_pc_factor_mat_solver_type", "mkl_pardiso");
	//PetscOptionsSetValue(NULL, "-ksp_gmres_restart", "60");
	//PetscOptionsSetValue(NULL, "-pc_type", "lu");
	//PetscOptionsSetValue(NULL, "-sub_pc_factor_mat_solver_type", "mkl_pardiso");
	/*PetscOptionsSetValue(NULL, "-mat_mumps_icntl_13", NULL);
	PetscOptionsSetValue(NULL, "-mat_mumps_icntl_23", "1000");
	PetscOptionsSetValue(NULL, "-mat_mumps_icntl_24", "1");
	PetscOptionsSetValue(NULL, "-mat_mumps_icntl_25", NULL);
	PetscOptionsSetValue(NULL, "-mat_mumps_icntl_28", "2");
	PetscOptionsSetValue(NULL, "-mat_mumps_icntl_29", "2");
	PetscOptionsSetValue(NULL, "-mat_mumps_cntl_3", "1e-20");
	PetscOptionsSetValue(NULL, "-mat_mumps_cntl_5", "1e-19");*/

	/*PetscOptionsSetValue(NULL, "-snes_linesearch_minlambda", "1e-19");
	PetscOptionsSetValue(NULL, "-snes_linesearch_rtol", "1e-16");
	PetscOptionsSetValue(NULL, "-snes_linesearch_atol", "1e-25");
	PetscOptionsSetValue(NULL, "-snes_linesearch_ltol", "1e-16");
	PetscOptionsSetValue(NULL, "-snes_linesearch_max_it", "20");
	PetscOptionsSetValue(NULL, "-snes_linesearch_keeplambda", "true");
	PetscOptionsSetValue(NULL, "-snes_linesearch_order", "3");*/

	for (i = 0; i < NOOptionsDB; i++) {
		if (strcmp(&(OptionsDBSecond[i][0]),"NULL")) PetscOptionsSetValue(NULL, &(OptionsDBFirst[i][0]), &(OptionsDBSecond[i][0]));
		else PetscOptionsSetValue(NULL, &(OptionsDBFirst[i][0]), NULL);
	}
	

	ierr = PetscInitialize(argc, argv, (char*)0, NULL);
	if (ierr) exit(ierr);

	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &MPIsize);
	CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &MPIRank);

	//if ((VectorMap = (int *)malloc(MPIsize * sizeof(int))) == NULL) TerM("Can not allocate memory for VectorMap");
	//if ((DisplacementMap = (int *)malloc(MPIsize * sizeof(int))) == NULL) TerM("Can not allocate memory for DisplacementMap");
	return ierr;

}

PetscErrorCode FormFunction(SNES snes, Vec Petsc_x, Vec Petsc_f, void *aptr) {
	int VecStart, VecEnd;
	PetscErrorCode ierr;
	register int i;
	AppCtx *actx= ((AppCtx*)(aptr));
	FType *Allxs;
	
	MethQM = 0;
	//MethQF = 0;
	BMethQF = 0;
	BMethQM = 0;

	actx->SimulationTime = SSimTime;
	ierr = VecGetOwnershipRange(Petsc_x, &VecStart, &VecEnd);
	CHKERRQ(ierr);
	if (MPIsize > 1) {		
		UpdateAll(Petsc_x, false);
		//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
		for (i = VecStart; i < VecEnd; i++) {		
			actx->GlobalIndex = i;
			actx->LocalIndex = i - VecStart;
			IteratorFunction(Petsc_x, Petsc_f, actx);
		}
		ierr = PetscBarrier(NULL); CHKERRQ(ierr);
	}
	else {
		UpdateAll(Petsc_x, false);
		for (i = 0; i < MatrixSize; i++) {
			actx->GlobalIndex = i;
			actx->LocalIndex = i;
			IteratorFunction(Petsc_x, Petsc_f, actx);
		}
	}
	
	
	ierr = VecGetArrayRead(Petsc_x, (const PetscScalar **)(&Allxs));
	CHKERRQ(ierr);
	ierr = VecRestoreArrayRead(Petsc_x, (const PetscScalar **)(&Allxs));
	CHKERRQ(ierr);
	ierr = VecGetArrayRead(Petsc_f, (const PetscScalar **)(&Allxs));
	CHKERRQ(ierr);
	ierr = VecRestoreArrayRead(Petsc_f, (const PetscScalar **)(&Allxs));
	CHKERRQ(ierr);

	return 0;
}

PetscErrorCode FormJacobian(SNES snes, Vec Petsc_x, Mat jac, Mat Petsc_B, void *aptr) {
	int VecStart, VecEnd;
	PetscErrorCode ierr;
	register int i;
	AppCtx *actx = ((AppCtx*)(aptr));
	PetscInt j;

	//double *gettest = new double[MatrixSize];
	//int *getitest = new int[MatrixSize];
	//int roww;

	//for (i = 0; i < MatrixSize; i++) getitest[i] = i;

	actx->SimulationTime = SSimTime;
	ierr = VecGetOwnershipRange(Petsc_x, &VecStart, &VecEnd);
	CHKERRQ(ierr);
	j = VecEnd - VecStart;

	if (MPIsize > 1) {		
		
		UpdateAll(Petsc_x, false);
		//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
		MatZeroEntries(Petsc_B);
		for (i = VecStart; i < VecEnd; i++) {
			actx->GlobalIndex = i;
			IteratorJacobian(Petsc_x, Petsc_B, actx);
		}
		//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
	}
	else {
		MatZeroEntries(Petsc_B);
		UpdateAll(Petsc_x, false);
		//for (i = VecStart; i < VecEnd; i++) {
		for (i = 0; i < MatrixSize; i++) {
			actx->GlobalIndex = i;
			IteratorJacobian(Petsc_x, Petsc_B, actx);
		}
	}

	
	jac = Petsc_B;
	ierr = MatAssemblyBegin(Petsc_B, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Petsc_B, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
	ierr = MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
	ierr = MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);


	//roww = 1;
	//MatGetValues(jac, 1, &roww, MatrixSize, getitest, gettest);
	
	ierr = MatSetVariableBlockSizes(jac, NewBlockSizes, (PetscInt*)(&(SadeqSize[SadeqIndex])));
	CHKERRQ(ierr);
	ierr = MatSetVariableBlockSizes(Petsc_B, NewBlockSizes, (PetscInt*)(&(SadeqSize[SadeqIndex])));
	CHKERRQ(ierr);
	/*ierr = MatSetMRLine(jac, NewBlockSizes, (PetscInt*)(&(PressureRow[SadeqIndex])));
	CHKERRQ(ierr);
	ierr = MatSetMRLine(Petsc_B, NewBlockSizes, (PetscInt*)(&(PressureRow[SadeqIndex])));
	CHKERRQ(ierr);*/
	

	//delete[] getitest;
	//delete[] getitest;
	return 0;
}


PetscInt UpdateAll(Vec XVec, bool FullTimeStep) {
	int Ix, Iy, Iz, i, BLNo;
	FType normSuml, normSumg;
	PetscScalar *XNew;
	Vec LovalXVec;
	VecScatter ctx;

	
	

	if (MPIsize > 1) {
		//ierr = PetscBarrier((PetscObject)(XVec)); CHKERRQ(ierr);
		//ierr= MPI_Allgatherv(XNew2, VectorMap[MPIRank], MPI_DOUBLE, SolvedXs, VectorMap, DisplacementMap, MPI_DOUBLE, PETSC_COMM_WORLD);		
		//XNew = SolvedXs;
		//ierr = VecCreateSeq(PETSC_COMM_SELF, MatrixSize, &LovalXVec);
		//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
		ierr = VecScatterCreateToAll(XVec, &ctx, &LovalXVec); CHKERRQ(ierr);
		ierr = VecScatterBegin(ctx, XVec, LovalXVec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
		ierr = VecScatterEnd(ctx, XVec, LovalXVec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
		ierr = VecGetArrayRead(LovalXVec, (const PetscScalar **)(&XNew));
		CHKERRQ(ierr);
		//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
	}
	else {
		ierr = VecGetArrayRead(XVec, (const PetscScalar **)(&XNew));
		CHKERRQ(ierr);
		//XNew = XNew2;
	}

	BLNo = 0;
	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				//if (porosity[Ix][Iy][Iz] > POR_TRESHOLD) {
					if (phaseStat[Ix][Iy][Iz] == 1) {
						normSuml = 0;
						for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][0] = XNew[BLNo];
							if (comp[Ix][Iy][Iz][i][0] < EPS_COMP) comp[Ix][Iy][Iz][i][0] = 0;
							else if (comp[Ix][Iy][Iz][i][0] > 1) comp[Ix][Iy][Iz][i][0] = 1;

							normSuml += comp[Ix][Iy][Iz][i][0];
							BLNo++;
						}
						if (normSuml) for (i = 0; i < Nc; i++) {
							 comp[Ix][Iy][Iz][i][0] /= normSuml;
						}
						else for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][0] = 1.0 / Nc;
						}

						if ((fabs(PCoeff*XNew[BLNo])) > (P[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) P[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF * dSgn(XNew[BLNo] * PCoeff - P[Ix + 1][Iy + 1][Iz + 1][1]);
						else P[Ix + 1][Iy + 1][Iz + 1][1] = XNew[BLNo] * PCoeff;
						if (P[Ix + 1][Iy + 1][Iz + 1][1] < MINPPP) P[Ix + 1][Iy + 1][Iz + 1][1] = MINPPP;
						BLNo++;

						sat[Ix][Iy][Iz][0] = XNew[BLNo];
						BLNo++;

						sat[Ix][Iy][Iz][1] = XNew[BLNo];
						BLNo++;

						if (sat[Ix][Iy][Iz][0] < 0) sat[Ix][Iy][Iz][0] = 0;
						if (sat[Ix][Iy][Iz][0] > 1) sat[Ix][Iy][Iz][0] = 1;
						if (sat[Ix][Iy][Iz][1] > 1) sat[Ix][Iy][Iz][1] = 1;
						if (sat[Ix][Iy][Iz][1] < 0) sat[Ix][Iy][Iz][1] = 0;

						normSuml = sat[Ix][Iy][Iz][0] + sat[Ix][Iy][Iz][1];
						if (normSuml) sat[Ix][Iy][Iz][0] /= normSuml;
						if (normSuml) sat[Ix][Iy][Iz][1] /= normSuml;
					}

					else if (phaseStat[Ix][Iy][Iz] == -1) {
						normSumg = 0;
						for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][1] = XNew[BLNo];
							if (comp[Ix][Iy][Iz][i][1] < EPS_COMP) comp[Ix][Iy][Iz][i][1] = 0;
							else if (comp[Ix][Iy][Iz][i][1] > 1) comp[Ix][Iy][Iz][i][1] = 1;
							normSumg += comp[Ix][Iy][Iz][i][1];
							BLNo++;
						}
						if (normSumg) for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][1] /= normSumg;
						}
						else for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][1] = 1.0 / Nc;
						}

						if ((fabs(PCoeff*XNew[BLNo])) > (P[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) P[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF * dSgn(XNew[BLNo] * PCoeff - P[Ix + 1][Iy + 1][Iz + 1][1]);
						else P[Ix + 1][Iy + 1][Iz + 1][1] = XNew[BLNo] * PCoeff;
						if (P[Ix + 1][Iy + 1][Iz + 1][1] < MINPPP) P[Ix + 1][Iy + 1][Iz + 1][1] = MINPPP;
						BLNo++;

						sat[Ix][Iy][Iz][0] = XNew[BLNo];
						BLNo++;

						sat[Ix][Iy][Iz][2] = XNew[BLNo];
						BLNo++;

						if (sat[Ix][Iy][Iz][0] < 0) sat[Ix][Iy][Iz][0] = 0;
						if (sat[Ix][Iy][Iz][0] > 1) sat[Ix][Iy][Iz][0] = 1;
						if (sat[Ix][Iy][Iz][2] > 1) sat[Ix][Iy][Iz][2] = 1;
						if (sat[Ix][Iy][Iz][2] < 0) sat[Ix][Iy][Iz][2] = 0;


						normSuml = sat[Ix][Iy][Iz][0] + sat[Ix][Iy][Iz][2];
						if (normSuml) sat[Ix][Iy][Iz][0] /= normSuml;
						if (normSuml) sat[Ix][Iy][Iz][2] /= normSuml;
					}
					else if (phaseStat[Ix][Iy][Iz] == 0) {
						normSuml = 0;
						for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][0] = XNew[BLNo];
							if (comp[Ix][Iy][Iz][i][0] < EPS_COMP) comp[Ix][Iy][Iz][i][0] = 0;
							else if (comp[Ix][Iy][Iz][i][0] > 1) comp[Ix][Iy][Iz][i][0] = 1;
							normSuml += comp[Ix][Iy][Iz][i][0];
							BLNo++;
						}
						normSumg = 0;
						for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][1] = XNew[BLNo];
							if (comp[Ix][Iy][Iz][i][1] < EPS_COMP) comp[Ix][Iy][Iz][i][1] = 0;
							else if (comp[Ix][Iy][Iz][i][1] > 1) comp[Ix][Iy][Iz][i][1] = 1;
							normSumg += comp[Ix][Iy][Iz][i][1];
							BLNo++;
						}
						if (normSuml) for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][0] /= normSuml;
						}
						else for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][0] = 1.0 / Nc;
						}
						if (normSumg) for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][1] /= normSumg;
						}
						else for (i = 0; i < Nc; i++) {
							comp[Ix][Iy][Iz][i][1] = 1.0 / Nc;
						}

						if ((fabs(PCoeff*XNew[BLNo])) > (P[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) P[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF * dSgn(XNew[BLNo] * PCoeff - P[Ix + 1][Iy + 1][Iz + 1][1]);
						else P[Ix + 1][Iy + 1][Iz + 1][1] = XNew[BLNo] * PCoeff;
						if (P[Ix + 1][Iy + 1][Iz + 1][1] < MINPPP) P[Ix + 1][Iy + 1][Iz + 1][1] = MINPPP;
						BLNo++;

						sat[Ix][Iy][Iz][0] = XNew[BLNo];
						BLNo++;

						sat[Ix][Iy][Iz][1] = XNew[BLNo];
						BLNo++;

						sat[Ix][Iy][Iz][2] = XNew[BLNo];
						BLNo++;

						if (sat[Ix][Iy][Iz][0] < 0) sat[Ix][Iy][Iz][0] = 0;
						if (sat[Ix][Iy][Iz][0] > 1) sat[Ix][Iy][Iz][0] = 1;
						if (sat[Ix][Iy][Iz][1] > 1) sat[Ix][Iy][Iz][1] = 1;
						if (sat[Ix][Iy][Iz][2] > 1) sat[Ix][Iy][Iz][2] = 1;
						if (sat[Ix][Iy][Iz][1] < 0) sat[Ix][Iy][Iz][1] = 0;
						if (sat[Ix][Iy][Iz][2] < 0) sat[Ix][Iy][Iz][2] = 0;

					}

					normSuml = sat[Ix][Iy][Iz][0] + sat[Ix][Iy][Iz][1] + sat[Ix][Iy][Iz][2];
					if (normSuml) sat[Ix][Iy][Iz][0] /= normSuml;
					if (normSuml) sat[Ix][Iy][Iz][1] /= normSuml;
					if (normSuml) sat[Ix][Iy][Iz][2] /= normSuml;
				//}
			}


	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				//if (porosity[Ix][Iy][Iz] > POR_TRESHOLD) {
					if (fphaseStat[Ix][Iy][Iz] == 1) {
						normSuml = 0;
						for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][0] = XNew[BLNo];
							if (fcomp[Ix][Iy][Iz][i][0] < EPS_COMP) fcomp[Ix][Iy][Iz][i][0] = 0;
							else if (fcomp[Ix][Iy][Iz][i][0] > 1) fcomp[Ix][Iy][Iz][i][0] = 1;

							normSuml += fcomp[Ix][Iy][Iz][i][0];
							BLNo++;
						}
						if (normSuml) for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][0] /= normSuml;
						}
						else for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][0] = 1.0 / Nc;
						}
						

						if ((fabs(PCoeff*XNew[BLNo])) > (fP[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) fP[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF * dSgn(XNew[BLNo] * PCoeff - fP[Ix + 1][Iy + 1][Iz + 1][1]);
						else fP[Ix + 1][Iy + 1][Iz + 1][1] = XNew[BLNo] * PCoeff;
						if (fP[Ix + 1][Iy + 1][Iz + 1][1] < MINPPP) fP[Ix + 1][Iy + 1][Iz + 1][1] = MINPPP;
						BLNo++;

						fsat[Ix][Iy][Iz][0] = XNew[BLNo];
						BLNo++;

						fsat[Ix][Iy][Iz][1] = XNew[BLNo];
						BLNo++;

						if (fsat[Ix][Iy][Iz][0] < 0) fsat[Ix][Iy][Iz][0] = 0;
						if (fsat[Ix][Iy][Iz][0] > 1) fsat[Ix][Iy][Iz][0] = 1;
						if (fsat[Ix][Iy][Iz][1] > 1) fsat[Ix][Iy][Iz][1] = 1;
						if (fsat[Ix][Iy][Iz][1] < 0) fsat[Ix][Iy][Iz][1] = 0;

						normSuml = fsat[Ix][Iy][Iz][0] + fsat[Ix][Iy][Iz][1];
						if (normSuml) fsat[Ix][Iy][Iz][0] /= normSuml;
						if (normSuml) fsat[Ix][Iy][Iz][1] /= normSuml;
					}

					else if (fphaseStat[Ix][Iy][Iz] == -1) {
						normSumg = 0;
						for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][1] = XNew[BLNo];
							if (fcomp[Ix][Iy][Iz][i][1] < EPS_COMP) fcomp[Ix][Iy][Iz][i][1] = 0;
							else if (fcomp[Ix][Iy][Iz][i][1] > 1) fcomp[Ix][Iy][Iz][i][1] = 1;
							normSumg += fcomp[Ix][Iy][Iz][i][1];
							BLNo++;
						}
						if (normSumg) for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][1] /= normSumg;
						}
						else for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][1] = 1.0 / Nc;
						}

						if ((fabs(PCoeff*XNew[BLNo])) > (fP[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) fP[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF * dSgn(XNew[BLNo] * PCoeff - fP[Ix + 1][Iy + 1][Iz + 1][1]);
						else fP[Ix + 1][Iy + 1][Iz + 1][1] = XNew[BLNo] * PCoeff;
						if (fP[Ix + 1][Iy + 1][Iz + 1][1] < MINPPP) fP[Ix + 1][Iy + 1][Iz + 1][1] = MINPPP;
						BLNo++;

						fsat[Ix][Iy][Iz][0] = XNew[BLNo];
						BLNo++;

						fsat[Ix][Iy][Iz][2] = XNew[BLNo];
						BLNo++;

						if (fsat[Ix][Iy][Iz][0] < 0) fsat[Ix][Iy][Iz][0] = 0;
						if (fsat[Ix][Iy][Iz][0] > 1) fsat[Ix][Iy][Iz][0] = 1;
						if (fsat[Ix][Iy][Iz][2] > 1) fsat[Ix][Iy][Iz][2] = 1;
						if (fsat[Ix][Iy][Iz][2] < 0) fsat[Ix][Iy][Iz][2] = 0;

						normSuml = fsat[Ix][Iy][Iz][0] + fsat[Ix][Iy][Iz][2];
						if (normSuml) fsat[Ix][Iy][Iz][0] /= normSuml;
						if (normSuml) fsat[Ix][Iy][Iz][2] /= normSuml;
					}
					else if (fphaseStat[Ix][Iy][Iz] == 0) {
						normSuml = 0;
						for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][0] = XNew[BLNo];
							if (fcomp[Ix][Iy][Iz][i][0] < EPS_COMP) fcomp[Ix][Iy][Iz][i][0] = 0;
							else if (fcomp[Ix][Iy][Iz][i][0] > 1) fcomp[Ix][Iy][Iz][i][0] = 1;
							normSuml += fcomp[Ix][Iy][Iz][i][0];
							BLNo++;
						}
						normSumg = 0;
						for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][1] = XNew[BLNo];
							if (fcomp[Ix][Iy][Iz][i][1] < EPS_COMP) fcomp[Ix][Iy][Iz][i][1] = 0;
							else if (fcomp[Ix][Iy][Iz][i][1] > 1) fcomp[Ix][Iy][Iz][i][1] = 1;
							normSumg += fcomp[Ix][Iy][Iz][i][1];
							BLNo++;
						}
						if (normSuml) for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][0] /= normSuml;
						}
						else for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][0] = 1.0 / Nc;
						}
						if (normSumg) for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][1] /= normSumg;
						}
						else for (i = 0; i < Nc; i++) {
							fcomp[Ix][Iy][Iz][i][1] = 1.0 / Nc;
						}
						if ((fabs(PCoeff*XNew[BLNo])) > (fP[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) fP[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF * dSgn(XNew[BLNo] * PCoeff - fP[Ix + 1][Iy + 1][Iz + 1][1]);
						else fP[Ix + 1][Iy + 1][Iz + 1][1] = XNew[BLNo] * PCoeff;
						if (fP[Ix + 1][Iy + 1][Iz + 1][1] < MINPPP) fP[Ix + 1][Iy + 1][Iz + 1][1] = MINPPP;
						BLNo++;

						fsat[Ix][Iy][Iz][0] = XNew[BLNo];
						BLNo++;

						fsat[Ix][Iy][Iz][1] = XNew[BLNo];
						BLNo++;

						fsat[Ix][Iy][Iz][2] = XNew[BLNo];
						BLNo++;

						if (fsat[Ix][Iy][Iz][0] < 0) fsat[Ix][Iy][Iz][0] = 0;
						if (fsat[Ix][Iy][Iz][0] > 1) fsat[Ix][Iy][Iz][0] = 1;
						if (fsat[Ix][Iy][Iz][1] > 1) fsat[Ix][Iy][Iz][1] = 1;
						if (fsat[Ix][Iy][Iz][2] > 1) fsat[Ix][Iy][Iz][2] = 1;
						if (fsat[Ix][Iy][Iz][1] < 0) fsat[Ix][Iy][Iz][1] = 0;
						if (fsat[Ix][Iy][Iz][2] < 0) fsat[Ix][Iy][Iz][2] = 0;

						normSuml = fsat[Ix][Iy][Iz][0] + fsat[Ix][Iy][Iz][1] + fsat[Ix][Iy][Iz][2];
						if (normSuml) fsat[Ix][Iy][Iz][0] /= normSuml;
						if (normSuml) fsat[Ix][Iy][Iz][1] /= normSuml;
						if (normSuml) fsat[Ix][Iy][Iz][2] /= normSuml;
					}
				//}
			}

			
			if (MPIsize > 1) {
				ierr = VecRestoreArrayRead(LovalXVec, (const PetscScalar **)(&XNew));
				CHKERRQ(ierr);				
				ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
				ierr = VecDestroy(&LovalXVec); CHKERRQ(ierr);
			}
			else {
				ierr = VecRestoreArrayRead(XVec, (const PetscScalar **)(&XNew));
				CHKERRQ(ierr);
			}
			//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
			UpdateDependencies(FullTimeStep);

			return 0;
}

void UpdateDependencies(bool FullTimeStep) {

	if (FullTimeStep) {
		PhaseCtrl();
		fPhaseCtrl();
		CalcPCoeff();
		PJac();
	}
	/*
	MixedCalcIFT();
	MixedCalcSatFs();
	MixedCalcP();
	MixedCalcVisco();
	MixedPreMul();
	MixedCalcDiffusion();
	MixedCalcGrav();
	MixedCalcTranses();
	*/

	if (MPIsize > 1) {
		MixedCalcIFT();
		MixedCalcSatFs();
		MixedCalcP();
		MixedCalcVisco();
		MixedPreMul();
		MixedCalcDiffusion();
		MixedCalcGrav();		
		MixedCalcTranses();		
	}
	else {
		CalcIFT();
		CalcSatFs();
		Calc_Phase_P();
		Calc_Visco();
		Pre_Mul();
		fCalcIFT();
		fCalcSatFs();
		fCalc_Phase_P();
		fCalc_Visco();
		fPre_Mul();
		CalcDiffusion();
		//CalcPCoeff();
		CalcSatFs();
		CalcPhaseGrav();
		Calc_Transes();
		fCalcDiffusion();
		fCalcSatFs();
		fCalcPhaseGrav();
		fCalc_Transes();
	}
	

	Qoprod = 0;
	Qginj = 0;
	WellCompSadeq[0] = 0;
	WellCompSadeq[1] = 0;
	WellCompSadeq[2] = 0;
	WellCompSadeq2[0] = 0;
	WellCompSadeq2[1] = 0;
	WellCompSadeq2[2] = 0;
	
	ZeroWellQ();
	
	gasInjStat = 0;
}

PetscErrorCode SolveNonLinear(void) {
	int VecStart, VecEnd, j;
	PetscInt IterNoX, IterNoL;
	KSPConvergedReason KReason;
	SNESLineSearch    linesearch;
	SNES snes;
	KSP Petsc_ksp;
	PC Petsc_pc;
	Vec Petsc_X, Petsc_R, solup;
	Mat Petsc_J;
	register int i;
	//static bool FirstTime = true;
	PetscScalar vnorm, *XNew2;

	
	ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
	CHKERRQ(ierr);
	//ierr = PetscBarrier((PetscObject)(snes)); CHKERRQ(ierr);

	/*if (FirstTime) {
		UpdateDependencies(false);
		FirstTime = false;
	}*/
	j = CreateAllocationTable();
	//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, j, MatrixSize, &Petsc_X);

	CHKERRQ(ierr);
	ierr = VecDuplicate(Petsc_X, &Petsc_R);
	CHKERRQ(ierr);

	ierr = MatCreateAIJ(PETSC_COMM_WORLD, j, j, MatrixSize, MatrixSize, 10*(2 * Nc + 4), NULL, 10*(2 * Nc + 4), NULL, &Petsc_J);
	//ierr = MatCreateBAIJ(PETSC_COMM_WORLD, j, j, MatrixSize, MatrixSize, (3 * Nc + 4), NULL, 7 * (3 * Nc + 4), NULL, &Petsc_J);
	
	CHKERRQ(ierr);
	ierr = MatSetOption(Petsc_J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	CHKERRQ(ierr);

	VecAssemblyBegin(Petsc_X);
	VecAssemblyEnd(Petsc_X);
	VecAssemblyBegin(Petsc_R);
	VecAssemblyEnd(Petsc_R);
	//InitGuess2(snes, Petsc_X, NULL);
	//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(Petsc_X, &VecStart, &VecEnd);
	CHKERRQ(ierr);
	//j = VecEnd - VecStart;
	
	//MPI_Allgather(&j, 1, MPI_INT, VectorMap, 1, MPI_INT, PETSC_COMM_WORLD);
	//DisplacementMap[0] = 0;
	//for (i = 0; i < MPIsize - 1; i++) DisplacementMap[i + 1] = DisplacementMap[i] + VectorMap[i];
	MakeList(VecStart, VecEnd);

	SNESSetFunction(snes, Petsc_R, FormFunction, (void *)(&userData));
	SNESSetJacobian(snes, Petsc_J, Petsc_J, FormJacobian, (void *)(&userData));
	//SNESSetJacobian(snes, NULL, NULL, SNESComputeJacobianDefault, NULL);
	SNESSetComputeInitialGuess(snes, InitGuess2, NULL);
	
	CSRrow[totalRow] = totalJac;


	
	ierr = SNESGetKSP(snes, &Petsc_ksp);
	CHKERRQ(ierr);
	ierr = KSPGetPC(Petsc_ksp, &Petsc_pc);
	CHKERRQ(ierr);
	ierr = KSPSetTolerances(Petsc_ksp, KSPRelTol, KSPAbsTol, KSPDivTol, MaxIters);
	CHKERRQ(ierr);
	ierr = SNESSetTolerances(snes, 1e-25, 1e-17, 1e-25, 2000, 2000);
	CHKERRQ(ierr);
	ierr = SNESSetType(snes, SNESKSPONLY);
	CHKERRQ(ierr);
	/*ierr = SNESGetLineSearch(snes, &linesearch);
	CHKERRQ(ierr);
	ierr = SNESLineSearchSetType(linesearch, SNESLINESEARCHBT);
	CHKERRQ(ierr);*/
	//if (Dt > 10 * MINDT) {
		ierr = KSPSetType(Petsc_ksp, SolverType);
		CHKERRQ(ierr);
		ierr = PCSetType(Petsc_pc, PCVPBJACOBI);
		CHKERRQ(ierr);
	/*}
	else {
		ierr = KSPSetType(Petsc_ksp, KSPBCGS);
		CHKERRQ(ierr);
		ierr = PCSetType(Petsc_pc, PCMINIMALRESIDUAL);
		CHKERRQ(ierr);
		ierr = PCMinimalResidualSetLines(Petsc_pc, 2*Nx*Ny*Nz, SadeqIndex, (PetscInt*)(PressureRow));
		CHKERRQ(ierr);
		ierr = PCMinimalResidualSetInnerIterations(Petsc_pc, 2);
		CHKERRQ(ierr);
		ierr = PCMinimalResidualSetNNZ(Petsc_pc, MatrixSize);
		CHKERRQ(ierr);
	}*/
	
	
	
	
	/*ierr = SNESLineSearchSetFromOptions(linesearch);
	CHKERRQ(ierr);*/
	ierr = SNESSetFromOptions(snes);
	CHKERRQ(ierr);
		
	ierr = SNESSolve(snes, NULL, Petsc_X);
	CHKERRQ(ierr);
	//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
	ierr = SNESGetSolutionUpdate(snes, &solup);
	CHKERRQ(ierr);
	ierr = VecNorm(solup, NORM_1, &vnorm);
	vnorm /= MatrixSize;
	/*ierr = VecGetArrayRead(solup, (const PetscScalar **)(&XNew2));
	CHKERRQ(ierr);
	ierr = VecRestoreArrayRead(solup, (const PetscScalar **)(&XNew2));
	CHKERRQ(ierr);*/

	//SNESComputeJacobian(snes, Petsc_X, Petsc_J, Petsc_J); CHKERRQ(ierr);

	//ierr = KSPView(Petsc_ksp, PETSC_VIEWER_STDOUT_WORLD);
	
	ierr = KSPGetConvergedReason(Petsc_ksp, &KReason);
	CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(Petsc_ksp, &IterNoL);
	CHKERRQ(ierr);
	ierr = SNESGetConvergedReason(snes, &snesreason);
	CHKERRQ(ierr);
	ierr = SNESGetIterationNumber(snes, &IterNoX);
	CHKERRQ(ierr);
	//PCView(Petsc_pc, PETSC_VIEWER_STDOUT_WORLD);
	//KSPView(Petsc_ksp, PETSC_VIEWER_STDOUT_WORLD);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\nSNES iterations = %D and SNES Convergence Reason=%D and Linear Solver Iterations=%d and Linear Solver Convergence Reason=%D\n", IterNoX, snesreason, IterNoL, KReason);
	CHKERRQ(ierr);
	//ierr = SNESView(snes, PETSC_VIEWER_STDOUT_WORLD);
	//CHKERRQ(ierr);
	
	

	//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
	if ((snesreason <= 0) || (KReason <= 0) || (vnorm > TOLMAIN) || ((vnorm - vnorm) != 0)) {
		ierr = VecDestroy(&Petsc_X);
		CHKERRQ(ierr);
		ierr = VecDestroy(&Petsc_R);
		CHKERRQ(ierr);
		//ierr = VecDestroy(&solup);
		//CHKERRQ(ierr);
		ierr = MatDestroy(&Petsc_J);
		CHKERRQ(ierr);
		ierr = SNESDestroy(&snes);
		CHKERRQ(ierr);
		return 0;
	}
	else {
		ierr = PetscBarrier(NULL); CHKERRQ(ierr);
		UpdateAll(Petsc_X, true);
		ierr = VecDestroy(&Petsc_X);
		CHKERRQ(ierr);
		ierr = VecDestroy(&Petsc_R);
		CHKERRQ(ierr);
		//ierr = VecDestroy(&solup);
		//CHKERRQ(ierr);
		ierr = MatDestroy(&Petsc_J);
		CHKERRQ(ierr);
		ierr = SNESDestroy(&snes);
		CHKERRQ(ierr);
		return 1;
	}
}

PetscErrorCode PetscFinishAll(void) {	
	//free(VectorMap);
	//free(DisplacementMap);
	ierr = PetscFinalize();

	return ierr;
}

PetscErrorCode InitGuess(SNES snes, Vec InitialGuess, void *ptr) {
	register int Ix, Iy, Iz, i, BLNo;
	FType *XNew;
	


	ierr = VecGetArray((InitialGuess), &XNew);
	CHKERRQ(ierr);
	

	BLNo = 0;
	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				if (phaseStat[Ix][Iy][Iz] == 1) {
					for (i = 0; i < Nc; i++) {
						XNew[BLNo] = comp[Ix][Iy][Iz][i][0];						
						BLNo++;
					}

					XNew[BLNo] = P[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
					BLNo++;

					XNew[BLNo] = sat[Ix][Iy][Iz][0];
					BLNo++;

					XNew[BLNo] = sat[Ix][Iy][Iz][1];
					BLNo++;
				}

				else if (phaseStat[Ix][Iy][Iz] == -1) {
					
					for (i = 0; i < Nc; i++) {
						XNew[BLNo] = comp[Ix][Iy][Iz][i][1];						
						BLNo++;
					}
					
					XNew[BLNo] = P[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
					BLNo++;

					XNew[BLNo] = sat[Ix][Iy][Iz][0];
					BLNo++;

					XNew[BLNo] = sat[Ix][Iy][Iz][2];
					BLNo++;
				}
				else if (phaseStat[Ix][Iy][Iz] == 0) {
					for (i = 0; i < Nc; i++) {
						XNew[BLNo] = comp[Ix][Iy][Iz][i][0];
						BLNo++;
					}
					
					for (i = 0; i < Nc; i++) {
						XNew[BLNo] = comp[Ix][Iy][Iz][i][1];						
						BLNo++;
					}
					
					XNew[BLNo] = P[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
					BLNo++;

					XNew[BLNo] = sat[Ix][Iy][Iz][0];
					BLNo++;

					XNew[BLNo] = sat[Ix][Iy][Iz][1];
					BLNo++;

					XNew[BLNo] = sat[Ix][Iy][Iz][2];
					BLNo++;
				}				
			}


	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				if (fphaseStat[Ix][Iy][Iz] == 1) {
					for (i = 0; i < Nc; i++) {
						XNew[BLNo] = fcomp[Ix][Iy][Iz][i][0];
						BLNo++;
					}

					XNew[BLNo] = fP[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
					BLNo++;

					XNew[BLNo] = fsat[Ix][Iy][Iz][0];
					BLNo++;

					XNew[BLNo] = fsat[Ix][Iy][Iz][1];
					BLNo++;
				}

				else if (fphaseStat[Ix][Iy][Iz] == -1) {

					for (i = 0; i < Nc; i++) {
						XNew[BLNo] = fcomp[Ix][Iy][Iz][i][1];
						BLNo++;
					}

					XNew[BLNo] = fP[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
					BLNo++;

					XNew[BLNo] = fsat[Ix][Iy][Iz][0];
					BLNo++;

					XNew[BLNo] = fsat[Ix][Iy][Iz][2];
					BLNo++;
				}
				else if (fphaseStat[Ix][Iy][Iz] == 0) {
					for (i = 0; i < Nc; i++) {
						XNew[BLNo] = fcomp[Ix][Iy][Iz][i][0];
						BLNo++;
					}

					for (i = 0; i < Nc; i++) {
						XNew[BLNo] = fcomp[Ix][Iy][Iz][i][1];
						BLNo++;
					}

					XNew[BLNo] = fP[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
					BLNo++;

					XNew[BLNo] = fsat[Ix][Iy][Iz][0];
					BLNo++;

					XNew[BLNo] = fsat[Ix][Iy][Iz][1];
					BLNo++;

					XNew[BLNo] = fsat[Ix][Iy][Iz][2];
					BLNo++;
				}
			}

	ierr = VecRestoreArray((InitialGuess), &XNew);

	return 0;
}

PetscErrorCode InitGuess2(SNES snes, Vec InitialGuess, void *ptr) {
	register int Ix, Iy, Iz, i, EqNo, qxtemp, iix, BNo, iIndex;
	FType *XNew;
	PetscInt AStart, AEnd;
	bool IsBFrac;


	//ierr = PetscBarrier(NULL); CHKERRQ(ierr);
	ierr = VecGetArray((InitialGuess), &XNew);
	CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(InitialGuess, &AStart, &AEnd);
	CHKERRQ(ierr);
	//ALength = AEnd - AStart;

	for (iix = AStart; iix < AEnd; iix++) {
		iIndex = iix - AStart;

		for (i = 1; i < (2 * Nx*Ny*Nz); i++) if (iix < preConRow[i]) break;
		BNo = i - 1;
		EqNo = iix - preConRow[BNo];

		if (BNo >= (Nx*Ny*Nz)) {
			BNo -= Nx * Ny*Nz;
			IsBFrac = true;
		}
		else IsBFrac = false;

		Iz = BNo / (Nx*Ny);
		qxtemp = BNo % (Nx*Ny);
		Iy = qxtemp / Nx;
		Ix = qxtemp % Nx;

		if (IsBFrac) {
			if (fphaseStat[Ix][Iy][Iz] == 1) {
				if (EqNo < Nc) XNew[iIndex] = fcomp[Ix][Iy][Iz][EqNo][0];
				else if (EqNo == Nc) XNew[iIndex] = fP[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
				else if (EqNo == (Nc + 1)) XNew[iIndex] = fsat[Ix][Iy][Iz][0];
				else XNew[iIndex] = fsat[Ix][Iy][Iz][1];
			}

			else if (fphaseStat[Ix][Iy][Iz] == -1) {
				if (EqNo < Nc) XNew[iIndex] = fcomp[Ix][Iy][Iz][EqNo][1];
				else if (EqNo == Nc) XNew[iIndex] = fP[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
				else if (EqNo == (Nc + 1)) XNew[iIndex] = fsat[Ix][Iy][Iz][0];
				else XNew[iIndex] = fsat[Ix][Iy][Iz][2];
			}
			else if (fphaseStat[Ix][Iy][Iz] == 0) {
				if (EqNo < Nc) XNew[iIndex] = fcomp[Ix][Iy][Iz][EqNo][0];
				else if (EqNo < 2 * Nc) XNew[iIndex] = fcomp[Ix][Iy][Iz][EqNo - Nc][1];
				else if (EqNo == 2 * Nc) XNew[iIndex] = fP[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
				else if (EqNo == (2 * Nc + 1)) XNew[iIndex] = fsat[Ix][Iy][Iz][0];
				else if (EqNo == (2 * Nc + 2)) XNew[iIndex] = fsat[Ix][Iy][Iz][1];
				else XNew[iIndex] = fsat[Ix][Iy][Iz][2];
			}
		}
		else {
			if (phaseStat[Ix][Iy][Iz] == 1) {
				if (EqNo < Nc) XNew[iIndex] = comp[Ix][Iy][Iz][EqNo][0];
				else if (EqNo == Nc) XNew[iIndex] = P[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
				else if (EqNo == (Nc + 1)) XNew[iIndex] = sat[Ix][Iy][Iz][0];
				else XNew[iIndex] = sat[Ix][Iy][Iz][1];
			}

			else if (phaseStat[Ix][Iy][Iz] == -1) {
				if (EqNo < Nc) XNew[iIndex] = comp[Ix][Iy][Iz][EqNo][1];
				else if (EqNo == Nc) XNew[iIndex] = P[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
				else if (EqNo == (Nc + 1)) XNew[iIndex] = sat[Ix][Iy][Iz][0];
				else XNew[iIndex] = sat[Ix][Iy][Iz][2];			
			}
			else if (phaseStat[Ix][Iy][Iz] == 0) {
				if (EqNo < Nc) XNew[iIndex] = comp[Ix][Iy][Iz][EqNo][0];
				else if (EqNo < 2 * Nc) XNew[iIndex] = comp[Ix][Iy][Iz][EqNo - Nc][1];
				else if (EqNo == 2 * Nc) XNew[iIndex] = P[Ix + 1][Iy + 1][Iz + 1][1] / PCoeff;
				else if (EqNo == (2 * Nc + 1)) XNew[iIndex] = sat[Ix][Iy][Iz][0];
				else if (EqNo == (2 * Nc + 2)) XNew[iIndex] = sat[Ix][Iy][Iz][1];
				else XNew[iIndex] = sat[Ix][Iy][Iz][2];
			}
		}
		
	}

	ierr = VecRestoreArray((InitialGuess), &XNew);
	CHKERRQ(ierr);
	ierr = PetscBarrier(NULL); CHKERRQ(ierr);

	return 0;
}

void MakeList(int Begining, int Finishing) {
	int Ix, Iy, Iz, itemp, C2UMSize, IsFracInt;
	register int i, j, k, p, r, d, q;
	

	for (i = 1; i < (2 * Nx*Ny*Nz); i++) if (Begining < preConRow[i]) break;
	Begining = i - 1;

	for (; i < (2 * Nx*Ny*Nz); i++) if ((Finishing - 1) < preConRow[i]) break;
	Finishing = i;

	SchindlerListLength = 0;
	for (i = 0; i < wellNO; i++) {
		for (j = welli[i][3]; j <= welli[i][4]; j++) {
			ShiftArray[SchindlerListLength] = j * (Nx*Ny) + welli[i][2] * Nx + welli[i][1];
			SchindlerListLength++;
			ShiftArray[SchindlerListLength] = ShiftArray[SchindlerListLength - 1] + Nx * Ny*Nz;
			SchindlerListLength++;
		}
	}

	for (i = Begining; i < Finishing; i++) {
		//if (porosity[Ix][Iy][Iz] > POR_TRESHOLD) {
			if (i < Nx*Ny*Nz) {
				p = i;
				ShiftArray[SchindlerListLength] = p + Nx * Ny*Nz;
				IsFracInt = 0;
			}
			else {
				p = i - Nx * Ny*Nz;
				ShiftArray[SchindlerListLength] = p;
				IsFracInt = 1;
			}
			SchindlerListLength++;



			Iz = p / (Nx*Ny);
			q = p % (Nx*Ny);
			Iy = q / Nx;
			Ix = q % Nx;

			ShiftArray[SchindlerListLength] = i;
			SchindlerListLength++;
			if (Ix > 0) {
				ShiftArray[SchindlerListLength] = Iz * (Nx*Ny) + Iy * Nx + (Ix - 1) + IsFracInt * Nx*Ny*Nz;
				SchindlerListLength++;
			}
			if (Ix < (Nx - 1)) {
				ShiftArray[SchindlerListLength] = Iz * (Nx*Ny) + Iy * Nx + Ix + 1 + IsFracInt * Nx*Ny*Nz;
				SchindlerListLength++;
			}
			if (Iy > 0) {
				ShiftArray[SchindlerListLength] = Iz * (Nx*Ny) + (Iy - 1) * Nx + Ix + IsFracInt * Nx*Ny*Nz;
				SchindlerListLength++;
			}
			if (Iy < (Ny - 1)) {
				ShiftArray[SchindlerListLength] = Iz * (Nx*Ny) + (Iy + 1) * Nx + Ix + IsFracInt * Nx*Ny*Nz;
				SchindlerListLength++;
			}
			if (Iz > 0) {
				ShiftArray[SchindlerListLength] = (Iz - 1)* (Nx*Ny) + Iy * Nx + Ix + IsFracInt * Nx*Ny*Nz;
				SchindlerListLength++;
			}
			if (Iz < (Nz - 1)) {
				ShiftArray[SchindlerListLength] = (Iz + 1)* (Nx*Ny) + Iy * Nx + Ix + IsFracInt * Nx*Ny*Nz;
				SchindlerListLength++;
			}
		//}
	}

	////////////////////////batcher sort////////////////////////////
	j = SchindlerListLength;
	k = 0;
	while (j) {
		j >>= 1;
		k++;
	}
	j = 1;
	j <<= (k - 1);

	if (j == SchindlerListLength) C2UMSize = SchindlerListLength;
	else C2UMSize = j << 1;

	i = C2UMSize >> 1;
	for (p = i; p >= 1; p >>= 1) {
		r = 0;
		d = p;
		for (q = i; q >= p; q >>= 1) {
			for (k = 0; k < (C2UMSize - d); k++) {
				if ((k & p) != r) continue;
				if ((k < SchindlerListLength) && ((k + d) < SchindlerListLength) && (ShiftArray[k] > ShiftArray[k + d])) {
					itemp = ShiftArray[k];
					ShiftArray[k] = ShiftArray[k + d];
					ShiftArray[k + d] = itemp;
				}
			}
			d = q - p;
			r = p;
		}
	}
	/////////////////////////////////////////////////////////////


	/////////////////////////////remove redundant////////////////////////////////////////////
	j = 0;
	for (i = 0; i < SchindlerListLength; i++) {		
		while (((i + 1) < SchindlerListLength) && (ShiftArray[i] == ShiftArray[i + 1])) i++;
		SchindlerList[j] = ShiftArray[i];
		j++;
	}
	SchindlerListLength = j;
	///////////////////////////////////////////////////////////////////////////////////////////
}

int CreateAllocationTable(void) {
	register int j, i;

	PJac();
	i = (2 * Nx*Ny*Nz);
	j = i / MPIsize;

	if (i % MPIsize) {
		SadeqIndex = (MPIRank)*(j + 1);
		if (MPIRank == (MPIsize - 1)) NewBlockSizes =i - SadeqIndex;
		else NewBlockSizes = (j + 1);
	}
	else {
		SadeqIndex = (MPIRank)*j;
		NewBlockSizes = j;
	}

	return (preConRow[SadeqIndex + NewBlockSizes] - preConRow[SadeqIndex]);
}
