#ifndef DO_TIME_STEP_H
#define DO_TIME_STEP_H

#include <math.h>
//#include <stdio.h>
#include<iostream>
#include<string>
#include<fstream>
#include"MIfstream.h"
//#include "Minor_NR.h"
#include "Functions.h"


extern FType ****preProp;
extern FType **fluidProp;
extern FType *****comp;
extern int Nx, Ny, Nz;
extern FType *****blockFProps;
extern FType ****sat;
extern FType watRo;
extern FType refP, cpor, dcpor;
extern FType ****P;
extern int Nc;
extern FType Dt;
extern FType *gridDim;
extern int wellNO;
extern FType QoT, QgT, QwT;
extern FType cumQo, cumQg, cumQw;
extern int incCount;
extern FType *****bcomp;
extern FType ****bP;
extern int *preConRow;
//extern int *preConIndex;
extern int TSMSize;
extern FType *TStepMarker;
extern char repStat;
extern signed char ***bphaseStat;
extern signed char ***phaseStat;
extern int TSMCtrl;
extern char maxNR;
extern char buildPreconFlag;



//void CreateBackups(FType);
//void ZeroWellQ(void);
//FType MatBal(FType);
//void RestoreBackups(void);
//char TimeStepCtrl(double *);
//FType AvgPReport(void);
//FType bAvgPReport(void);
//void CreateBackups15(FType); 

void Pre_Mul(void) {
	register int i, j, k, n;
	FType porpor;

	for (i = 0; i<Nx; i++)
		for (j = 0; j<Ny; j++)
			for (k = 0; k<Nz; k++) {
				for (n = 0; n<(Nc - 1); n++) {
					porpor = porosity[i][j][k] * (1 + cpor*(P[i + 1][j + 1][k + 1][1] - refP) + dcpor*(P[i + 1][j + 1][k + 1][1] - refP)*(P[i + 1][j + 1][k + 1][1] - refP));
					preProp[i][j][k][n] = porpor*(comp[i][j][k][n][0] * blockFProps[i][j][k][RO][0] * sat[i][j][k][1] + comp[i][j][k][n][1] * blockFProps[i][j][k][RO][1] * sat[i][j][k][2]);
				}
				preProp[i][j][k][Nc - 1] = porpor*(blockFProps[i][j][k][RO][0] * sat[i][j][k][1] + blockFProps[i][j][k][RO][1] * sat[i][j][k][2]);
				preProp[i][j][k][Nc] = porpor*watRo*WAT_M_RO*sat[i][j][k][0];
			}
}

void fPre_Mul(void) {
	register int i, j, k, n;
	FType porpor;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++) {
				for (n = 0; n < (Nc - 1); n++) {
					porpor = fporosity[i][j][k] * (1 + cpor*(fP[i + 1][j + 1][k + 1][1] - refP) + dcpor*(fP[i + 1][j + 1][k + 1][1] - refP)*(fP[i + 1][j + 1][k + 1][1] - refP));
					fpreProp[i][j][k][n] = porpor*(fcomp[i][j][k][n][0] * fblockFProps[i][j][k][RO][0] * fsat[i][j][k][1] + fcomp[i][j][k][n][1] * fblockFProps[i][j][k][RO][1] * fsat[i][j][k][2]);
				}
				fpreProp[i][j][k][Nc - 1] = porpor*(fblockFProps[i][j][k][RO][0] * fsat[i][j][k][1] + fblockFProps[i][j][k][RO][1] * fsat[i][j][k][2]);
				fpreProp[i][j][k][Nc] = porpor*watRo*WAT_M_RO*fsat[i][j][k][0];
			}
}

void CreateBackups(FType simTime) {
	register int Ix, Iy, Iz, i;
	//std::ofstream fp;

	//if ((fp=fopen("test.rst","wb"))==NULL) puts("Warning: can not write in restart file!\n");
	//if (!MPIRank) fp.open("test.rst", std::ios::binary);
	//if (!MPIRank) if (!fp.is_open()) std::cout << "Warning: can not write in restart file!\n";

	for (Ix = 0; Ix<(Nx + 2); Ix++)
		for (Iy = 0; Iy<(Ny + 2); Iy++)
			for (Iz = 0; Iz<(Nz + 2); Iz++)
				for (i = 0; i<3; i++) {
					bP[Ix][Iy][Iz][i] = P[Ix][Iy][Iz][i];
					//fwrite(&P[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
					//if (!MPIRank) fp.write((char *)&P[Ix][Iy][Iz][i], sizeof(FType));
				}
	for (Ix = 0; Ix < (Nx + 2); Ix++)
		for (Iy = 0; Iy < (Ny + 2); Iy++)
			for (Iz = 0; Iz < (Nz + 2); Iz++)
				for (i = 0; i < 3; i++) {
					fbP[Ix][Iy][Iz][i] = fP[Ix][Iy][Iz][i];
					//fwrite(&P[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
					//if (!MPIRank) fp.write((char *)&fP[Ix][Iy][Iz][i], sizeof(FType));
				}
	for (Ix = 0; Ix<Nx; Ix++)
		for (Iy = 0; Iy<Ny; Iy++)
			for (Iz = 0; Iz<Nz; Iz++) {
				//bP[Ix][Iy][Iz]=P[Ix][Iy][Iz][1];
				for (i = 0; i<3; i++) {
					bsat[Ix][Iy][Iz][i] = sat[Ix][Iy][Iz][i];
					//fwrite(&sat[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
					//if (!MPIRank) fp.write((char *)&sat[Ix][Iy][Iz][i], sizeof(FType));
				}

				for (i = 0; i<Nc; i++) {
					bcomp[Ix][Iy][Iz][i][0] = comp[Ix][Iy][Iz][i][0];
					//fwrite(&comp[Ix][Iy][Iz][i][0], sizeof(FType), 1, fp);
					//if (!MPIRank) fp.write((char *)&comp[Ix][Iy][Iz][i][0], sizeof(FType));
					bcomp[Ix][Iy][Iz][i][1] = comp[Ix][Iy][Iz][i][1];
					//fwrite(&comp[Ix][Iy][Iz][i][1], sizeof(FType), 1, fp);
					//if (!MPIRank) fp.write((char *)&comp[Ix][Iy][Iz][i][1], sizeof(FType));
				}

				bphaseStat[Ix][Iy][Iz] = phaseStat[Ix][Iy][Iz];
				//fwrite(&phaseStat[Ix][Iy][Iz], sizeof(char), 1, fp);
				//if (!MPIRank) fp.write((char *)&phaseStat[Ix][Iy][Iz], sizeof(char));

			}
	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {
				//bP[Ix][Iy][Iz]=P[Ix][Iy][Iz][1];
				for (i = 0; i < 3; i++) {
					fbsat[Ix][Iy][Iz][i] = fsat[Ix][Iy][Iz][i];
					//fwrite(&sat[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
					//if (!MPIRank) fp.write((char *)&fsat[Ix][Iy][Iz][i], sizeof(FType));
				}

				for (i = 0; i < Nc; i++) {
					fbcomp[Ix][Iy][Iz][i][0] = fcomp[Ix][Iy][Iz][i][0];
					//fwrite(&comp[Ix][Iy][Iz][i][0], sizeof(FType), 1, fp);
					//if (!MPIRank) fp.write((char *)&fcomp[Ix][Iy][Iz][i][0], sizeof(FType));
					fbcomp[Ix][Iy][Iz][i][1] = fcomp[Ix][Iy][Iz][i][1];
					//fwrite(&comp[Ix][Iy][Iz][i][1], sizeof(FType), 1, fp);
					//if (!MPIRank) fp.write((char *)&fcomp[Ix][Iy][Iz][i][1], sizeof(FType));
				}

				fbphaseStat[Ix][Iy][Iz] = fphaseStat[Ix][Iy][Iz];
				//fwrite(&phaseStat[Ix][Iy][Iz], sizeof(char), 1, fp);
				//if (!MPIRank) fp.write((char *)&fphaseStat[Ix][Iy][Iz], sizeof(char));

			}
	for (i = 0; i < wellNO; i++) bWellCondition[i] = WellCondition[i];


	/*fwrite(&simTime, sizeof(FType), 1, fp);
	fwrite(&Dt, sizeof(FType), 1, fp);
	fwrite(&TSMCtrl, sizeof(int), 1, fp);

	fwrite(&SumQoProduced, sizeof(FType), 1, fp);
	fwrite(&SumQgInjected, sizeof(FType), 1, fp);

	fclose(fp);*/

	/*if (!MPIRank) fp.write((char *)&simTime, sizeof(FType));
	if (!MPIRank) fp.write((char *)&Dt, sizeof(FType));
	if (!MPIRank) fp.write((char *)&TSMCtrl, sizeof(int));

	if (!MPIRank) fp.write((char *)&SumQoProduced, sizeof(FType));
	if (!MPIRank) fp.write((char *)&SumQgInjected, sizeof(FType));
	if (!MPIRank) fp.close();*/
}

void WriteRST(FType simTime) {
	register int Ix, Iy, Iz, i;
	std::ofstream fp;

	fp.open("test.rst", std::ios::binary | std::ios::out);
	if (!fp.is_open()) std::cout << "Warning: can not write in restart file!\n";

	for (Ix = 0; Ix < (Nx + 2); Ix++)
		for (Iy = 0; Iy < (Ny + 2); Iy++)
			for (Iz = 0; Iz < (Nz + 2); Iz++)
				for (i = 0; i < 3; i++) {					
					fp.write((char *)&P[Ix][Iy][Iz][i], sizeof(FType));
				}
	for (Ix = 0; Ix < (Nx + 2); Ix++)
		for (Iy = 0; Iy < (Ny + 2); Iy++)
			for (Iz = 0; Iz < (Nz + 2); Iz++)
				for (i = 0; i < 3; i++) {					
					fp.write((char *)&fP[Ix][Iy][Iz][i], sizeof(FType));
				}
	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {				
				for (i = 0; i < 3; i++) {					
					fp.write((char *)&sat[Ix][Iy][Iz][i], sizeof(FType));
				}

				for (i = 0; i < Nc; i++) {					
					fp.write((char *)&comp[Ix][Iy][Iz][i][0], sizeof(FType));					
					fp.write((char *)&comp[Ix][Iy][Iz][i][1], sizeof(FType));
				}				
				fp.write((char *)&phaseStat[Ix][Iy][Iz], sizeof(char));

			}

	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {
				for (i = 0; i < 3; i++) {
					fp.write((char *)&fsat[Ix][Iy][Iz][i], sizeof(FType));
				}

				for (i = 0; i < Nc; i++) {
					fp.write((char *)&fcomp[Ix][Iy][Iz][i][0], sizeof(FType));
					fp.write((char *)&fcomp[Ix][Iy][Iz][i][1], sizeof(FType));
				}
				fp.write((char *)&fphaseStat[Ix][Iy][Iz], sizeof(char));

			}	

	fp.write((char *)&simTime, sizeof(FType));
	fp.write((char *)&Dt, sizeof(FType));
	fp.write((char *)&TSMCtrl, sizeof(int));

	fp.write((char *)&SMethQ, sizeof(FType));
	fp.close();
}

FType ReadRST(void) {
	register int Ix, Iy, Iz, i;
	std::ifstream fp;
	FType simTime;

	fp.open("test.rst", std::ios::binary | std::ios::in);
	if (!fp.is_open()) TerM("can not read from restart file!\n");

	for (Ix = 0; Ix < (Nx + 2); Ix++)
		for (Iy = 0; Iy < (Ny + 2); Iy++)
			for (Iz = 0; Iz < (Nz + 2); Iz++)
				for (i = 0; i < 3; i++) {
					fp.read((char *)&P[Ix][Iy][Iz][i], sizeof(FType));
				}
	for (Ix = 0; Ix < (Nx + 2); Ix++)
		for (Iy = 0; Iy < (Ny + 2); Iy++)
			for (Iz = 0; Iz < (Nz + 2); Iz++)
				for (i = 0; i < 3; i++) {
					fp.read((char *)&fP[Ix][Iy][Iz][i], sizeof(FType));
				}
	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {
				for (i = 0; i < 3; i++) {
					fp.read((char *)&sat[Ix][Iy][Iz][i], sizeof(FType));
				}

				for (i = 0; i < Nc; i++) {
					fp.read((char *)&comp[Ix][Iy][Iz][i][0], sizeof(FType));
					fp.read((char *)&comp[Ix][Iy][Iz][i][1], sizeof(FType));
				}
				fp.read((char *)&phaseStat[Ix][Iy][Iz], sizeof(char));

			}

	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {
				for (i = 0; i < 3; i++) {
					fp.read((char *)&fsat[Ix][Iy][Iz][i], sizeof(FType));
				}

				for (i = 0; i < Nc; i++) {
					fp.read((char *)&fcomp[Ix][Iy][Iz][i][0], sizeof(FType));
					fp.read((char *)&fcomp[Ix][Iy][Iz][i][1], sizeof(FType));
				}
				fp.read((char *)&fphaseStat[Ix][Iy][Iz], sizeof(char));

			}

	fp.read((char *)&simTime, sizeof(FType));
	fp.read((char *)&Dt, sizeof(FType));
	fp.read((char *)&TSMCtrl, sizeof(int));

	fp.read((char *)&SMethQ, sizeof(FType));
	fp.close();

	return simTime;
}

FType MatBal(FType simTime) {
	register int Ix, Iy, Iz;
	FType Vb, porpor, Fsh, Fsw, rV;
	static FType ihv = 0;
	static FType iwv = 0;

	cumQo += QoT*Dt;
	cumQg += QgT*Dt;
	cumQw += QwT*Dt;

	SumQgInjected += QgInjected*Dt;
	SumQoProduced += QoProduced*Dt;


	Fsh = -cumQo - cumQg - ihv;
	Fsw = -cumQw - iwv;
	for (Ix = 0; Ix<Nx; Ix++)
		for (Iy = 0; Iy<Ny; Iy++)
			for (Iz = 0; Iz<Nz; Iz++) {
				porpor = porosity[Ix][Iy][Iz] * (1 + cpor*(P[Ix + 1][Iy + 1][Iz + 1][1] - refP) + dcpor*(P[Ix + 1][Iy + 1][Iz + 1][1] - refP)*(P[Ix + 1][Iy + 1][Iz + 1][1] - refP));
				Vb = gridDim[Ix] * gridDim[Nx + Iy] * gridDim[Nx + Ny + Iz];

				Fsh += blockFProps[Ix][Iy][Iz][RO][0] * Vb*porpor*sat[Ix][Iy][Iz][1];
				Fsh += blockFProps[Ix][Iy][Iz][RO][1] * Vb*porpor*sat[Ix][Iy][Iz][2];

				Fsw += WAT_M_RO*watRo*Vb*porpor*sat[Ix][Iy][Iz][0];
			}

	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {
				porpor = fporosity[Ix][Iy][Iz] * (1 + cpor*(fP[Ix + 1][Iy + 1][Iz + 1][1] - refP) + dcpor*(fP[Ix + 1][Iy + 1][Iz + 1][1] - refP)*(fP[Ix + 1][Iy + 1][Iz + 1][1] - refP));
				Vb = gridDim[Ix] * gridDim[Nx + Iy] * gridDim[Nx + Ny + Iz];

				Fsh += fblockFProps[Ix][Iy][Iz][RO][0] * Vb*porpor*fsat[Ix][Iy][Iz][1];
				Fsh += fblockFProps[Ix][Iy][Iz][RO][1] * Vb*porpor*fsat[Ix][Iy][Iz][2];

				Fsw += WAT_M_RO*watRo*Vb*porpor*fsat[Ix][Iy][Iz][0];
			}

	if (simTime < 0) {
		ihv = Fsh;
		iwv = Fsw;
	}

	rV = 0;
	if (iwv) rV += fabs(Fsw) / iwv;
	if (ihv) rV += fabs(Fsh) / ihv;

	return rV / 2;
}

void ZeroWellQ(void) {
	QoT = 0;
	QgT = 0;
	QwT = 0;
}

void RestoreBackups(void) {
	register int Ix, Iy, Iz, i;

	for (Ix = 0; Ix<(Nx + 2); Ix++)
		for (Iy = 0; Iy<(Ny + 2); Iy++)
			for (Iz = 0; Iz<(Nz + 2); Iz++)
				for (i = 0; i<3; i++) {
					P[Ix][Iy][Iz][i] = bP[Ix][Iy][Iz][i];
					fP[Ix][Iy][Iz][i] = fbP[Ix][Iy][Iz][i];

				}

	for (Ix = 0; Ix<Nx; Ix++)
		for (Iy = 0; Iy<Ny; Iy++)
			for (Iz = 0; Iz<Nz; Iz++) {
				//P[Ix][Iy][Iz][1]=bP[Ix][Iy][Iz];
				for (i = 0; i<3; i++) {
					sat[Ix][Iy][Iz][i] = bsat[Ix][Iy][Iz][i];
					fsat[Ix][Iy][Iz][i] = fbsat[Ix][Iy][Iz][i];

				}

				for (i = 0; i<Nc; i++) {
					comp[Ix][Iy][Iz][i][0] = bcomp[Ix][Iy][Iz][i][0];
					comp[Ix][Iy][Iz][i][1] = bcomp[Ix][Iy][Iz][i][1];
					fcomp[Ix][Iy][Iz][i][0] = fbcomp[Ix][Iy][Iz][i][0];
					fcomp[Ix][Iy][Iz][i][1] = fbcomp[Ix][Iy][Iz][i][1];
				}
				phaseStat[Ix][Iy][Iz] = bphaseStat[Ix][Iy][Iz];
				fphaseStat[Ix][Iy][Iz] = fbphaseStat[Ix][Iy][Iz];

			}
	for (i = 0; i < wellNO; i++) WellCondition[i] = bWellCondition[i];


}

char TimeStepCtrl(double *SimTime) {
	register int Ix, Iy, Iz;
	char TFlag = -1;
	//static int TSMCtrl=0;
	//static char MsgNecessary=-1;

	/*for (Ix = 0; Ix<Nx; Ix++)
		for (Iy = 0; Iy<Ny; Iy++)
			for (Iz = 0; Iz<Nz; Iz++) {
				if ((fabs(P[Ix + 1][Iy + 1][Iz + 1][1] - bP[Ix + 1][Iy + 1][Iz + 1][1]))>MAXPP) {
					TFlag = 0;
					break;
				}
				//if (phaseStat[Ix][Iy][Iz]) critTest(Ix, Iy, Iz);
			}*/

	if (repStat) {
		TFlag = 0;
		repStat = 0;
	}

	/*if (Dt<MINDT) {
	//TFlag=-1;
	//if (MsgNecessary) puts("\nMINDT was hit!\n");
	//else {
	//	MsgNecessary=-1;
	//	Dt=MINDT;
	//}

	puts("\nMINDT was hit!\n");
	Dt=MINDT;
	}*/

	if (TFlag) {
		QgCumInj += QgInj*Dt;
		QoCumProd += QoProd*Dt;
		(*SimTime) += Dt;
		incCount++;
		if (incCount>NTIMESTERINC) {
			buildPreconFlag = -1;
			if (!maxNR) Dt *= TIME_STEP_INC;
			incCount = 0;
		}
		if (Dt<MINDT) Dt = MINDT;
		if (Dt>Max_Dt) Dt = Max_Dt;

		if ((((*SimTime) + Dt)>TStepMarker[TSMCtrl]) || ((fabs((*SimTime) + Dt - TStepMarker[TSMCtrl]))<TSTEPIF)) {
			Dt = TStepMarker[TSMCtrl] - (*SimTime);
			TSMCtrl++;
			//MsgNecessary=0;		
		}
		if (TSMCtrl == TSMSize) TFlag = 1;
	}
	else {
		incCount = 0;
		Dt /= TIME_STEP_DEC;
		buildPreconFlag = -1;
		if (Dt<MINDT) {
#ifndef USE_MKL_SOLVER
			PetscFinishAll();
#endif
			TerM("\nMINDT was hit!\n");
			Dt = MINDT;
			(*SimTime) += Dt;
		}
		else {
			RestoreBackups();
			PJac();
		}
		//PhaseCtrl();
	}

	//if (SimTime>500) Dt=100;
	maxNR = 0;
	repStat = 0;



	return TFlag;
}

FType AvgPReport(void) {
	register int i, j, k;
	FType porpor, Vp, sVp, meanP;

	meanP = 0;
	sVp = 0;
	for (i = 0; i<Nx; i++)
		for (j = 0; j<Ny; j++)
			for (k = 0; k<Nz; k++) {
				porpor = porosity[i][j][k] * (1 + cpor*(P[i + 1][j + 1][k + 1][1] - refP) + dcpor*(P[i + 1][j + 1][k + 1][1] - refP)*(P[i + 1][j + 1][k + 1][1] - refP));
				Vp = porpor*gridDim[i] * gridDim[Nx + j] * gridDim[Nx + Ny + k];
				sVp += Vp;
				meanP += P[i + 1][j + 1][k + 1][1] * Vp;
			}
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++) {
				porpor = fporosity[i][j][k] * (1 + cpor*(fP[i + 1][j + 1][k + 1][1] - refP) + dcpor*(fP[i + 1][j + 1][k + 1][1] - refP)*(fP[i + 1][j + 1][k + 1][1] - refP));
				Vp = porpor*gridDim[i] * gridDim[Nx + j] * gridDim[Nx + Ny + k];
				sVp += Vp;
				meanP += fP[i + 1][j + 1][k + 1][1] * Vp;
			}

	//CalcPCoeff();

	return meanP / sVp;
}

void CreateBackups15(FType simTime) {
	register int Ix, Iy, Iz, i;
	std::ofstream fp;
	static char chstat = -1;

	if (chstat) {
		//if ((fp=fopen("test15.rst","wb"))==NULL) puts("Warning: can not write in restart file!\n");
		fp.open("test15.rst", std::ios::binary);
		if (!fp.is_open()) std::cout << "Warning: can not write in restart file!\n";


		for (Ix = 0; Ix<(Nx + 2); Ix++)
			for (Iy = 0; Iy<(Ny + 2); Iy++)
				for (Iz = 0; Iz<(Nz + 2); Iz++)
					for (i = 0; i<3; i++) {
						//bP[Ix][Iy][Iz][i]=P[Ix][Iy][Iz][i];

						//fwrite(&P[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
						fp.write((char *)&P[Ix][Iy][Iz][i], sizeof(FType));
					}
		for (Ix = 0; Ix<Nx; Ix++)
			for (Iy = 0; Iy<Ny; Iy++)
				for (Iz = 0; Iz<Nz; Iz++) {
					//bP[Ix][Iy][Iz]=P[Ix][Iy][Iz][1];
					for (i = 0; i<3; i++) {
						//bsat[Ix][Iy][Iz][i]=sat[Ix][Iy][Iz][i];	
						//fwrite(&sat[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
						fp.write((char *)&sat[Ix][Iy][Iz][i], sizeof(FType));
					}

					for (i = 0; i<Nc; i++) {
						//bcomp[Ix][Iy][Iz][i][0]=comp[Ix][Iy][Iz][i][0];
						//fwrite(&comp[Ix][Iy][Iz][i][0], sizeof(FType), 1, fp);
						fp.write((char *)&comp[Ix][Iy][Iz][i][0], sizeof(FType));
						//bcomp[Ix][Iy][Iz][i][1]=comp[Ix][Iy][Iz][i][1];
						//fwrite(&comp[Ix][Iy][Iz][i][1], sizeof(FType), 1, fp);
						fp.write((char *)&comp[Ix][Iy][Iz][i][1], sizeof(FType));
					}
					//bphaseStat[Ix][Iy][Iz]=phaseStat[Ix][Iy][Iz];
					//fwrite(&phaseStat[Ix][Iy][Iz], sizeof(char), 1, fp);
					fp.write((char *)&phaseStat[Ix][Iy][Iz], sizeof(char));
				}

		//fwrite(&simTime, sizeof(FType), 1, fp);
		//fwrite(&Dt, sizeof(FType), 1, fp);
		//fwrite(&TSMCtrl, sizeof(int), 1, fp);


		//fclose(fp);
		fp.write((char *)&simTime, sizeof(FType));
		fp.write((char *)&Dt, sizeof(FType));
		fp.write((char *)&TSMCtrl, sizeof(int));

		fp.close();
		chstat = 0;
	}

}

FType bAvgPReport(void) {
	register int i, j, k;
	FType porpor, Vp, sVp, meanP;

	meanP = 0;
	sVp = 0;
	for (i = 0; i<Nx; i++)
		for (j = 0; j<Ny; j++)
			for (k = 0; k<Nz; k++) {
				porpor = porosity[i][j][k] * (1 + cpor*(bP[i + 1][j + 1][k + 1][1] - refP) + dcpor*(bP[i + 1][j + 1][k + 1][1] - refP)*(bP[i + 1][j + 1][k + 1][1] - refP));
				Vp = porpor*gridDim[i] * gridDim[Nx + j] * gridDim[Nx + Ny + k];
				sVp += Vp;
				meanP += bP[i + 1][j + 1][k + 1][1] * Vp;
			}
	return meanP / sVp;
}

#endif