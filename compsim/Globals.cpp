#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdlib.h>
#include "Functions.h"
#include <time.h>
//#include <stdio.h>

#ifdef USE_MKL_SOLVER
FType *MKLCoefMatrix, *MKLRHSVector;
#endif

#ifdef USE_PETSC
#include <petscsnes.h>
char SolverType[256];
FType KSPRelTol, KSPAbsTol, KSPDivTol;
PetscInt MaxIters;
int NOOptionsDB;
char **OptionsDBFirst, **OptionsDBSecond;
PetscErrorCode ierr;
int MPIsize, MPIRank;
FType SSimTime;
SNESConvergedReason snesreason;
AppCtx userData;
int MatrixSize;
int *SchindlerList;
int SchindlerListLength;
FType *SolvedXs;
int *ShiftArray;
//int *OffDiagonalElements, *OnDiagonalElements;
//int LocalJacobianBlockSize;
int NewBlockSizes, SadeqIndex;
//int *LocalRowsAIJ;
//char *AllStates, *StatesCopy;
FType PCoeff;
int *PressureRow;
bool IsConstDif;
FType CDiffCoef;
FType SigmaMF;
FType MethQM, SMethQ;
char BMethQF, BMethQM;

//FType *ScalingVector, *ScalingVectorBackUp;
//int FraqColBase;
#endif

#include <iostream>
#include <fstream>
#include "MIfstream.h"
#include <math.h>





///////////////////GLOBALS//////////////////////////////////
int Nx, Ny, Nz;		//Reservoir Dimension
FType *gridDim;	//Grid size
FType ***porosity;		//Porosity Distribution
FType ****perm;		//Permeability tensor distribution
FType refP, cpor, dcpor;		//Compressibility
int PNc, UNc, Nc;		//Number of Components
unsigned char PR, SRK;		//EOS Type, boolean
FType **fluidProp;			//Fluids components properties
int Nswt, Nsgt;		//Number of saturation tables inputs
FType **swt, **sgt;			//Saturation tables
int initCond;		//Initial Condition State
FType ****sat, ****bsat;		//Saturation and its backup
FType ****fsat, ****fbsat;		//Fractured Saturation and its backup
FType ****P;		//Pressure
FType *****comp;		//Phases Composition
int refL;		//Reference layer in initial condition
int wellNO;		//Number of wells
int **welli;	//integer well data: X-Well, Y_Well, Well_Type
FType **wellf;	//Double well data: Pwf or Q
FType ****IFTran;		//Permeability average on blocks boundary, kA/h
FType ****fIFTran;		//Fracture Permeability average on blocks boundary, kA/h
FType resTemp;		//Reservoir Temperature
FType Lf;//Half of Fracture opening
FType Re_f;//Half Fractured length
//FType Re;
FType *****blockFProps;		//Calculated fluids properties in each gridblock
FType *****fblockFProps;     //Calculated fluids properties in each Fracture gridblock 
FType *****trans;		//Transmissibility
FType *****ftrans;		//Fracture Transmissibility
FType ****Wtran;		//Water Transmissibility
FType ****fWtran;		//FRacture Water Transmissibility
FType ****relPerm;		//Relative Permeability
FType ****frelPerm;		//Fracture Relative Permeability
FType watRo, watMu;	//Water properties.
					//FType **jac;		//Jacobian matrix
FType *ans;		//Answer matrix
FType ****preProp;		//Multiplication of properties in previous time step
FType ****fpreProp;		//Fracture Multiplication of properties in previous time step
FType Dt;		//Timestep size, seconds
FType *****dE;		//Phase density derivative
FType *****fdE;		//fracture Phase density derivative
char *****transS;		//upstream weighting condition
char *****ftransS;		//Fracture upstream weighting condition
FType **satJac, *satAns, *Xm, *Xms;	//Minor NR Matrixes
signed char ***phaseStat;		//Number and type of the phases in GBs, 1->L, -1->G, 0->2Ph
signed char ***fphaseStat;		//Fracture Number and type of the phases in GBs, 1->L, -1->G, 0->2Ph

						//int curSize;		//Jacobian Matrix current size
						//FType *Unk;		//Unknowns matrix
FType WOCHeight;	//WOC
FType totalTime;	//total simulation time, days input, seconds change
FType *Ki;
FType QoT, QgT, QwT;		//Well flow rates
FType cumQo, cumQg, cumQw;		//cummulative rates
FType ***blockH;		//Block Heights
int incCount = 0;		//Time Step Control Counter
FType *****bcomp;		//Composition Backup
FType *****fbcomp;		//Fracture Composition Backup
FType ****bP;		//Pressure Backup
FType ****fbP;		//Fracture Pressure Backup
FType Eps;		//Machine Epsilon
FType *CSRjac;	//Compressed Sparse Row Jacobian
int *CSRrow, *CSRcol;		//Compressed sparse row location holders
int bCSRSize;		//Compressed Sparse Row Size for a block
					//int jacIndex;		//Compressed Sparse Row Jacobian Index Holder
					//int CSRrowIndex;		//Compressed Sparse Row Index Holder
int ***pJHolder;		//Jacobian Position Holder: Parallel design
int ***fpJHolder;		//Jacobian Position Holder in fractured: Parallel design
//FType PCoeff;		//Mean P for Preconditioning
FType *preCon;		//Preconditioner
int *preConRow;		//Preconditioner row holder
int *preConIndex;		//Preconditioner index holder
						// cl_platform_id *platformIds;
						// cl_device_id deviceId;
						// cl_context context;
						// cl_program program;
						// cl_kernel xZero, xZeroF, init, f1, MatMul, f2, f3, clBuildPrecon, dot_persist_kernel, MatMulT, Pf1, Pf2, Pf3, Pf4, dot_persist_kernelF, MatMulTF, MatMulF, ClScaling, clBinPartialSort;
						// cl_command_queue commandQueue;
						//char chWellStat=0;
int totalJac, totalRow;
FType **EdgeComp; // compesion in interface between matrix & Fracture
FType *****eq_comp;//compesion in interface between matrix & Fracture
FType *f_m_relperm; //relative permeability between Fracture And Matrix
FType **f_m_drelperm; //moshtagh relative permeability between Fracture And Matrix  0 = dkrwdsw; 1 = dkrodsw; 2 = dkrodsg; 3=dkrgdsg; 4=dkrodso 
FType **bic;
FType *TStepMarker;
int TSMSize;
char repStat = 0;
FType **STcomp;
FType wellrate[50];
signed char ***bphaseStat;
signed char ***fbphaseStat;
int TSMCtrl = 0;
char maxNR = 0;
FType ****dRelPerm;
FType ****fdRelPerm;
FType ****dWtran;
FType ****fdWtran;
char gasInjStat;
//FType **fullPre;
int *preConCSRrow;
int *preConCSRcol;
double solverTimeAcc = 0;
double solverTime = 0;
clock_t BICGStart, BICGEnd;
FType ***tor, *****diffusion;
FType  *****fdiffusion;
FType *f_m_diffusion;//Diffusion between Farcture and Matrix
FType Qw, Qo, Qg;
FType wGLR;
char buildPreconFlag = -1;
//cl_mem clM, clMCSRrow, clMCSRcol;
FType ***Bift;
FType ***fBift;
FType QgCumInj = 0, QgInj;
FType QoCumProd = 0, QoProd;
FType QoProduced, QgInjected, SumQoProduced = 0, SumQgInjected = 0;
FType nesbat = 0;

long int nprocs;
int *CPUTempArray;
int *SadeqSize;
bool *WellCondition;
bool *TempWellCondition;
bool *bWellCondition;

FType Qginj, Qoprod;
FType Qginj_cum = 0;
FType Qoprod_cum = 0;
FType WellCompSadeq[3];
FType WellCompSadeq_cum[3];
FType WellCompSadeq2[3];
FType *Jac1Row;
int *Jac1Col;
FType Max_Dt;





///////////////////////////new for Fracture///////////////

FType ***fporosity;		//Fracture Porosity Distribution
FType ****fperm;		//Fracture Permeability tensor distribution
FType ****fP;		//Pressure of fracture
FType *****fcomp;		//Phases Composition of fracture
						//int FractureBase;

						//void TerM(char *str);
						//FType dSgn(FType);
						//void CLFinish(void);*
						//void ParaSolver(FType *, int *, int *, FType *, FType *);

FType *IComp, *IOComp, *IGComp;

FType Solve_Z(FType a1, FType a2, FType a3, char state) {
	FType Q, J, D, x1, x2, x3, t, t1, Z;
	FType a13, qsqrtemp;

	Q = (3 * a2 - a1*a1) / 9;
	J = (9 * a1*a2 - 27 * a3 - 2 * a1*a1*a1) / 54;
	D = Q*Q*Q + J*J;

	if (D<0) {
		qsqrtemp = 2 * sqrt(-Q);
		a13 = a1 / 3;
		t = acos(J / sqrt(-Q*Q*Q)) / 3;
		x1 = qsqrtemp*cos(t) - a13;
		x2 = qsqrtemp*cos(t + 2 * PI / 3) - a13;
		x3 = qsqrtemp*cos(t + 4 * PI / 3) - a13;
		if ((state == 'g') || (state == 'G')) {
			if ((x1 >= x2) && (x1 >= x3)) Z = x1;
			else if ((x2 >= x1) && (x2 >= x3)) Z = x2;
			else Z = x3;
		}
		else {		//Liquid state
			if ((x1 <= x2) && (x1 <= x3) && (x1>0)) Z = x1;
			else if ((x2 <= x1) && (x2 <= x3) && (x2>0)) Z = x2;
			else Z = x3;
		}
	}
	else if (D == 0) {
		if (J>0) t = pow(J, 1.0 / 3);
		else t = -pow(-J, 1.0 / 3);
		x1 = 2 * t - a1 / 3;
		x2 = -t - a1 / 3;
		if ((state == 'g') || (state == 'G')) {
			if (x2 >= x1) Z = x2;
			else Z = x1;
		}
		else {
			if ((x2 <= x1) && (x2>0)) Z = x2;
			else Z = x1;
		}
	}
	else {
		t1 = sqrt(D);
		t = J + t1;
		if (t>0) Z = pow(t, 1.0 / 3);
		else Z = -pow(-t, 1.0 / 3);
		t = J - t1;
		if (t>0) Z += pow(t, 1.0 / 3);
		else Z -= pow(-t, 1.0 / 3);
		Z -= a1 / 3;
	}

	return Z;
}

void TerM(char *str) {
	register int i, j, k, n;

#ifdef USE_MKL_SOLVER
	free(MKLCoefMatrix);
	free(MKLRHSVector);
#endif
	std::cout << str;

	for (i = 0; i < NOOptionsDB; i++) {
		free(OptionsDBFirst[i]);
	}
	free(OptionsDBFirst);
	for (i = 0; i < NOOptionsDB; i++) {
		free(OptionsDBSecond[i]);
	}
	free(OptionsDBSecond);
	//free(ScalingVector);
	//free(ScalingVectorBackUp);
	free(SolvedXs);
	//free(OffDiagonalElements);
	//free(OnDiagonalElements);
	//free(AllStates);
	//free(StatesCopy);
	free(ShiftArray);
	free(SchindlerList);
	free(gridDim);
	free(WellCondition);
	free(bWellCondition);
	free(TempWellCondition);

	free(SadeqSize);

	free(IComp);
	free(IOComp);
	free(IGComp);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) free(porosity[i][j]);
		free(porosity[i]);
	}
	free(porosity);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) free(fporosity[i][j]);
		free(fporosity[i]);
	}
	free(fporosity);


	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) free(Bift[i][j]);
		free(Bift[i]);
	}
	free(Bift);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) free(fBift[i][j]);
		free(fBift[i]);
	}
	free(fBift);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) free(tor[i][j]);
		free(tor[i]);
	}
	free(tor);

	for (i = 0; i<(Nx + 2); i++) {
		for (j = 0; j<(Ny + 2); j++) {
			for (k = 0; k<(Nz + 2); k++) free(P[i][j][k]);
			free(P[i][j]);
		}
		free(P[i]);
	}
	free(P);

	for (i = 0; i < (Nx + 2); i++) {
		for (j = 0; j < (Ny + 2); j++) {
			for (k = 0; k < (Nz + 2); k++) free(fP[i][j][k]);
			free(fP[i][j]);
		}
		free(fP[i]);
	}
	free(fP);


	for (i = 0; i<(Nx + 2); i++) {
		for (j = 0; j<(Ny + 2); j++) {
			for (k = 0; k<(Nz + 2); k++) free(bP[i][j][k]);
			free(bP[i][j]);
		}
		free(bP[i]);
	}
	free(bP);

	for (i = 0; i < (Nx + 2); i++) {
		for (j = 0; j < (Ny + 2); j++) {
			for (k = 0; k < (Nz + 2); k++) free(fbP[i][j][k]);
			free(fbP[i][j]);
		}
		free(fbP[i]);
	}
	free(fbP);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) free(perm[i][j][k]);
			free(perm[i][j]);
		}
		free(perm[i]);
	}
	free(perm);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) free(fperm[i][j][k]);
			free(fperm[i][j]);
		}
		free(fperm[i]);
	}
	free(fperm);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) free(sat[i][j][k]);
			free(sat[i][j]);
		}
		free(sat[i]);
	}
	free(sat);


	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) free(fsat[i][j][k]);
			free(fsat[i][j]);
		}
		free(fsat[i]);
	}
	free(fsat);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) free(bsat[i][j][k]);
			free(bsat[i][j]);
		}
		free(bsat[i]);
	}
	free(bsat);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) free(fbsat[i][j][k]);
			free(fbsat[i][j]);
		}
		free(fbsat[i]);
	}
	free(fbsat);

	for (i = 0; i<Nc; i++) {
		free(fluidProp[i]);
	}
	free(fluidProp);

	for (i = 0; i<Nswt; i++) {
		free(swt[i]);
	}
	free(swt);

	for (i = 0; i<Nsgt; i++) {
		free(sgt[i]);
	}
	free(sgt);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) {
				for (n = 0; n<Nc; n++) free(comp[i][j][k][n]);
				free(comp[i][j][k]);
			}
			free(comp[i][j]);
		}
		free(comp[i]);
	}
	free(comp);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) {
				for (n = 0; n < Nc; n++) free(fcomp[i][j][k][n]);
				free(fcomp[i][j][k]);
			}
			free(fcomp[i][j]);
		}
		free(fcomp[i]);
	}
	free(fcomp);

	for (i = 0; i < Nc; i++) {
		free(EdgeComp[i]);
	}
	free(EdgeComp);

	free(f_m_relperm);

	free(f_m_diffusion);


	for (i = 0; i < 3; i++) {
		free(f_m_drelperm[i]);
	}
	free(f_m_drelperm);


	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) {
				for (n = 0; n<Nc; n++) free(diffusion[i][j][k][n]);
				free(diffusion[i][j][k]);
			}
			free(diffusion[i][j]);
		}
		free(diffusion[i]);
	}
	free(diffusion);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) {
				for (n = 0; n < Nc; n++) free(fdiffusion[i][j][k][n]);
				free(fdiffusion[i][j][k]);
			}
			free(fdiffusion[i][j]);
		}
		free(fdiffusion[i]);
	}
	free(fdiffusion);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) {
				for (n = 0; n<Nc; n++) free(bcomp[i][j][k][n]);
				free(bcomp[i][j][k]);
			}
			free(bcomp[i][j]);
		}
		free(bcomp[i]);
	}
	free(bcomp);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) {
				for (n = 0; n < Nc; n++) free(fbcomp[i][j][k][n]);
				free(fbcomp[i][j][k]);
			}
			free(fbcomp[i][j]);
		}
		free(fbcomp[i]);
	}
	free(fbcomp);

	for (i = 0; i<wellNO; i++) {
		free(welli[i]);
		free(wellf[i]);
	}
	free(welli);
	free(wellf);


	//////////////////////////////////////////////////
	//////////////NON_INPUT///////////////////////////
	/////////////////////////////////////////////////
	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) free(IFTran[i][j][k]);
			free(IFTran[i][j]);
		}
		free(IFTran[i]);
	}
	free(IFTran);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) free(fIFTran[i][j][k]);
			free(fIFTran[i][j]);
		}
		free(fIFTran[i]);
	}
	free(fIFTran);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) {
				for (n = 0; n<BLOCK_F_PROPS; n++) free(blockFProps[i][j][k][n]);
				free(blockFProps[i][j][k]);
			}
			free(blockFProps[i][j]);
		}
		free(blockFProps[i]);
	}
	free(blockFProps);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) {
				for (n = 0; n < BLOCK_F_PROPS; n++) free(fblockFProps[i][j][k][n]);
				free(fblockFProps[i][j][k]);
			}
			free(fblockFProps[i][j]);
		}
		free(fblockFProps[i]);
	}
	free(fblockFProps);

	for (i = 0; i<(Nx + 1); i++) {
		for (j = 0; j<(Ny + 1); j++) {
			for (k = 0; k<(Nz + 1); k++) {
				for (n = 0; n<12; n++) {
					free(trans[i][j][k][n]);
				}
				free(trans[i][j][k]);
			}
			free(trans[i][j]);
		}
		free(trans[i]);
	}
	free(trans);

	for (i = 0; i < (Nx + 1); i++) {
		for (j = 0; j < (Ny + 1); j++) {
			for (k = 0; k < (Nz + 1); k++) {
				for (n = 0; n < 12; n++) {
					free(ftrans[i][j][k][n]);
				}
				free(ftrans[i][j][k]);
			}
			free(ftrans[i][j]);
		}
		free(ftrans[i]);
	}
	free(ftrans);

	for (i = 0; i<(Nx + 1); i++) {
		for (j = 0; j<(Ny + 1); j++) {
			for (k = 0; k<(Nz + 1); k++) {
				for (n = 0; n<3; n++) {
					free(transS[i][j][k][n]);
				}
				free(transS[i][j][k]);
			}
			free(transS[i][j]);
		}
		free(transS[i]);
	}
	free(transS);

	for (i = 0; i < (Nx + 1); i++) {
		for (j = 0; j < (Ny + 1); j++) {
			for (k = 0; k < (Nz + 1); k++) {
				for (n = 0; n < 3; n++) {
					free(ftransS[i][j][k][n]);
				}
				free(ftransS[i][j][k]);
			}
			free(ftransS[i][j]);
		}
		free(ftransS[i]);
	}
	free(ftransS);

	for (i = 0; i<(Nx + 1); i++) {
		for (j = 0; j<(Ny + 1); j++) {
			for (k = 0; k<(Nz + 1); k++) free(Wtran[i][j][k]);
			free(Wtran[i][j]);
		}
		free(Wtran[i]);
	}
	free(Wtran);

	for (i = 0; i < (Nx + 1); i++) {
		for (j = 0; j < (Ny + 1); j++) {
			for (k = 0; k < (Nz + 1); k++) free(fWtran[i][j][k]);
			free(fWtran[i][j]);
		}
		free(fWtran[i]);
	}
	free(fWtran);

	for (i = 0; i<(Nx + 1); i++) {
		for (j = 0; j<(Ny + 1); j++) {
			for (k = 0; k<(Nz + 1); k++) free(dWtran[i][j][k]);
			free(dWtran[i][j]);
		}
		free(dWtran[i]);
	}
	free(dWtran);

	for (i = 0; i < (Nx + 1); i++) {
		for (j = 0; j < (Ny + 1); j++) {
			for (k = 0; k < (Nz + 1); k++) free(fdWtran[i][j][k]);
			free(fdWtran[i][j]);
		}
		free(fdWtran[i]);
	}
	free(fdWtran);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) free(relPerm[i][j][k]);
			free(relPerm[i][j]);
		}
		free(relPerm[i]);
	}
	free(relPerm);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) free(frelPerm[i][j][k]);
			free(frelPerm[i][j]);
		}
		free(frelPerm[i]);
	}
	free(frelPerm);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) free(dRelPerm[i][j][k]);
			free(dRelPerm[i][j]);
		}
		free(dRelPerm[i]);
	}
	free(dRelPerm);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) free(fdRelPerm[i][j][k]);
			free(fdRelPerm[i][j]);
		}
		free(fdRelPerm[i]);
	}
	free(fdRelPerm);


	/*for (i=0; i<(Nx*Ny*Nz*(2*Nc+4)); i++) {
	free(preCon[i]);
	}
	free(preCon);*/

	free(CSRjac);
	free(CSRrow);
	free(CSRcol);

	free(preCon);

	free(ans);
	free(preConRow);
	free(preConIndex);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) free(pJHolder[i][j]);
		free(pJHolder[i]);
	}
	free(pJHolder);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) free(fpJHolder[i][j]);
		free(fpJHolder[i]);
	}
	free(fpJHolder);

	//free(Unk);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) free(preProp[i][j][k]);
			free(preProp[i][j]);
		}
		free(preProp[i]);
	}
	free(preProp);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) free(fpreProp[i][j][k]);
			free(fpreProp[i][j]);
		}
		free(fpreProp[i]);
	}
	free(fpreProp);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) {
			for (k = 0; k<Nz; k++) {
				for (n = 0; n<(Nc + 1); n++) free(dE[i][j][k][n]);
				free(dE[i][j][k]);
			}
			free(dE[i][j]);
		}
		free(dE[i]);
	}
	free(dE);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) {
				for (n = 0; n < (Nc + 1); n++) free(fdE[i][j][k][n]);
				free(fdE[i][j][k]);
			}
			free(fdE[i][j]);
		}
		free(fdE[i]);
	}
	free(fdE);

	////////////////////////////////////////////////////

	for (i = 0; i<(Nc + 1); i++) {
		free(satJac[i]);
	}
	free(satJac);

	free(satAns);
	free(Xm);
	free(Xms);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) free(phaseStat[i][j]);
		free(phaseStat[i]);
	}
	free(phaseStat);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) free(fphaseStat[i][j]);
		free(fphaseStat[i]);
	}
	free(fphaseStat);

	for (i = 0; i<Nx; i++) {
		for (j = 0; j<Ny; j++) free(bphaseStat[i][j]);
		free(bphaseStat[i]);
	}
	free(bphaseStat);

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) free(fbphaseStat[i][j]);
		free(fbphaseStat[i]);
	}
	free(fbphaseStat);

	free(Ki);

	for (i = 0; i<(Nx + 2); i++) {
		for (j = 0; j<(Ny + 2); j++) free(blockH[i][j]);
		free(blockH[i]);
	}
	free(blockH);

	for (i = 0; i<Nc; i++) {
		free(bic[i]);
	}
	free(bic);

	free(TStepMarker);

	for (i = 0; i<Nc; i++) {
		free(STcomp[i]);
	}
	free(STcomp);

	free(preConCSRrow);
	free(preConCSRcol);

	free(CPUTempArray);

	free(Jac1Row);
	free(Jac1Col);




	//CLFinish();

	//system("shutdown -s -f");

	exit(0);
}

FType dSgn(FType x) {
	if (x < 0) return -1.0;
	else if (x > 0) return 1.0;
	else return 0;
}


/*void CLFinish(void) {
free(platformIds);
if (xZero) clReleaseKernel(xZero);
if (xZeroF) clReleaseKernel(xZeroF);
if (init) clReleaseKernel(init);
if (f1) clReleaseKernel(f1);
if (MatMul) clReleaseKernel(MatMul);
if (MatMulT) clReleaseKernel(MatMulT);
if (MatMulTF) clReleaseKernel(MatMulTF);
if (MatMulF) clReleaseKernel(MatMulF);
if (f2) clReleaseKernel(f2);
if (f3) clReleaseKernel(f3);
if (clBuildPrecon) clReleaseKernel(clBuildPrecon);
if (dot_persist_kernel) clReleaseKernel(dot_persist_kernel);
if (dot_persist_kernelF) clReleaseKernel(dot_persist_kernelF);
if (Pf1) clReleaseKernel(Pf1);
if (Pf2) clReleaseKernel(Pf2);
if (Pf3) clReleaseKernel(Pf3);
if (Pf4) clReleaseKernel(Pf4);
if (ClScaling) clReleaseKernel(ClScaling);
if (clBinPartialSort) clReleaseKernel(clBinPartialSort);


if (program) clReleaseProgram(program);

if (context) clReleaseContext(context);
if (commandQueue) clReleaseCommandQueue(commandQueue);

if (clMCSRrow) clReleaseMemObject(clMCSRrow);
if (clMCSRcol) clReleaseMemObject(clMCSRcol);
if (clM) clReleaseMemObject(clM);
}*/

#endif