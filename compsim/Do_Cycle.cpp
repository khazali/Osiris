#ifndef DO_CYCLE_H
#define DO_CYCLE_H

#include <math.h>
#include <time.h>
#include <float.h>
//#include "Globals.h"
//#include "Minor_NR.h"
//#include "Gauss_Jordan.h"
//#include "Serial_BiCGStab.h"
//#include "MKL_Solver.h"
//#include "Do_Time_Step.h"
//#include "Orthomin.h"
//#include "GMRES.h"
//#include "CL_BiCGS.h"
//#include "Numeric.h"
//#include "Pre.h"
//#include "BiCGStab_l.h"
//#include "Mass_Tran.h"
//#include "ParSol.hpp"
#include "Functions.h"

extern int Nc;
extern FType **fluidProp;
extern FType *****comp;
extern int Nx, Ny, Nz;
extern FType *****blockFProps;
extern FType *****trans;
extern FType ****relPerm;
extern int Nswt;
extern int Nsgt;
extern FType *gridDim;
extern FType watRo, watMu;
extern unsigned char PR, SRK;
//extern FType **jac;
extern FType *ans;
extern FType ****sat, ****bsat;
extern FType Dt;
extern FType *****dE;
extern FType refP, cpor, dcpor;
extern char *****transS;
extern FType ****P;
extern FType resTemp;
extern signed char ***phaseStat;
extern FType *Xm;
//extern int curSize;
//extern FType *Unk;
extern FType ****Wtran;
extern FType QoT, QgT, QwT;
extern FType cumQo, cumQg, cumQw;
extern FType ***blockH;
//extern int jacIndex;
//extern int CSRrowIndex;
extern int ***pJHolder;
//extern FType PCoeff;
extern FType *preCon;
extern int *preConRow;
extern int *preConIndex;
//extern char chWellStat;
extern int totalJac, totalRow;
extern FType **bic;
extern FType **STcomp;
extern FType wellrate[50];
extern FType ****dRelPerm;
extern FType ****dWtran;
extern char gasInjStat;
extern FType *CSRjac;
extern int *CSRrow, *CSRcol;
extern int *preConCSRrow;
extern int *preConCSRcol;
extern double solverTime, solverTimeAcc;
extern clock_t BICGStart, BICGEnd;
extern FType ***tor, *****diffusion;
extern FType Qw, Qo, Qg;
extern FType wGLR;
extern char buildPreconFlag;
extern FType ***Bift;




//void CalcSatFs(void);
//void Calc_Visco(void);
//void Calc_Transes(void);
//void Calc_Phase_P(void);
//void System_Reduce(void);
//void PhaseCtrl(void);
//FType System_Update(void);
//void CalcPhaseGrav(void);
//void DelJacRow(int, int);*
//void DelJacCol(int, int);*
//void BuildPrecon(void);
//double oilRateCh(int, int, int, double);
//double oilRateCh2(int, int, int, double, double);
//void JacGuass(void);
//void CalcIFT(void);


void CalcSatFs(void) {
	register int i, j, k, n;
	FType Krog, Krow, sdKrog, sdKrow;
	FType AA, BB, dAA, dBB;
	//FType satL, satG;
	FType Fift;

	//FType CoreySor=0;
	//FILE *fp;
	//fp=fopen("IND.txt", "wt");

	//swt[0][SAT]=0;


	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++) {
				/*relPerm[i][j][k][0]=RELPERM0;
				dRelPerm[i][j][k][0]=0;

				if (sat[i][j][k][1]<CoreySor) {
				relPerm[i][j][k][1]=RELPERM0;
				dRelPerm[i][j][k][1]=0;
				dRelPerm[i][j][k][2]=0;
				dRelPerm[i][j][k][4]=0;
				}
				else if (sat[i][j][k][1]>=1) {
				relPerm[i][j][k][1]=1;
				dRelPerm[i][j][k][1]=0;
				dRelPerm[i][j][k][2]=0;
				dRelPerm[i][j][k][4]=0;
				}
				else {
				relPerm[i][j][k][1]=pow((sat[i][j][k][1]-CoreySor)/(1-CoreySor), 3.0);
				dRelPerm[i][j][k][1]=0;
				dRelPerm[i][j][k][2]=0;
				dRelPerm[i][j][k][4]=(3.0/(1-CoreySor))*pow((sat[i][j][k][1]-CoreySor)/(1-CoreySor), 2.0);
				}
				//fprintf(fp, "%d\t%d\t%d\to=%f\t%f\n", i, j, k, relPerm[i][j][k][1], dRelPerm[i][j][k][4]);

				if (sat[i][j][k][2]<=0) {
				relPerm[i][j][k][2]=RELPERM0;
				dRelPerm[i][j][k][3]=0;
				blockFProps[i][j][k][BLOCK_PC][1]=0;
				}
				else if (sat[i][j][k][2]>=(1-CoreySor)) {
				relPerm[i][j][k][2]=1;
				dRelPerm[i][j][k][3]=0;
				blockFProps[i][j][k][BLOCK_PC][1]=3000941.83744005;
				}
				else {
				relPerm[i][j][k][2]=pow(sat[i][j][k][2]/(1-CoreySor), 1.5);

				dRelPerm[i][j][k][3]=(1.5/(1-CoreySor))*pow(sat[i][j][k][2]/(1-CoreySor), 0.5);
				blockFProps[i][j][k][BLOCK_PC][1]=103400*pow((sat[i][j][k][1]-CoreySor)/(1-CoreySor), -1.0/2.7);
				}
				blockFProps[i][j][k][BLOCK_PC][0]=0;
				//fprintf(fp, "%d\t%d\t%d\tg=%f\t%f\t%f\n", i, j, k, relPerm[i][j][k][2], dRelPerm[i][j][k][3], blockFProps[i][j][k][BLOCK_PC][1]);


				*/

				dRelPerm[i][j][k][4] = 0;			//Oil to Oil
													//Relative permeability, Baker's model & Capillary pressures
				for (n = 0; n < Nswt; n++) {
					if (sat[i][j][k][0] < swt[n][SAT]) break;
				}
				if (n == 0) {
					relPerm[i][j][k][0] = RELPERM0;
					blockFProps[i][j][k][BLOCK_PC][0] = swt[0][CAPILLARYPRESSURE];
					Krow = swt[0][KRO];

					dRelPerm[i][j][k][0] = 0;
					sdKrow = 0;
				}
				else if (n == Nswt) {
					relPerm[i][j][k][0] = swt[Nswt - 1][KR];
					blockFProps[i][j][k][BLOCK_PC][0] = swt[Nswt - 1][CAPILLARYPRESSURE];
					Krow = swt[Nswt - 1][KRO];

					dRelPerm[i][j][k][0] = 0;
					sdKrow = 0;
				}
				else {
					relPerm[i][j][k][0] = (swt[n][KR] - swt[n - 1][KR]) / (swt[n][SAT] - swt[n - 1][SAT])*(sat[i][j][k][0] - swt[n][SAT]) + swt[n][KR];
					Krow = (swt[n][KRO] - swt[n - 1][KRO]) / (swt[n][SAT] - swt[n - 1][SAT])*(sat[i][j][k][0] - swt[n][SAT]) + swt[n][KRO];
					blockFProps[i][j][k][BLOCK_PC][0] = (swt[n][CAPILLARYPRESSURE] - swt[n - 1][CAPILLARYPRESSURE]) / (swt[n][SAT] - swt[n - 1][SAT])*(sat[i][j][k][0] - swt[n][SAT]) + swt[n][CAPILLARYPRESSURE];

					dRelPerm[i][j][k][0] = (swt[n][KR] - swt[n - 1][KR]) / (swt[n][SAT] - swt[n - 1][SAT]);
					sdKrow = (swt[n][KRO] - swt[n - 1][KRO]) / (swt[n][SAT] - swt[n - 1][SAT]);
				}

				for (n = 0; n < Nsgt; n++) {
					if (sat[i][j][k][2] < sgt[n][SAT]) break;
				}
				if (n == 0) {
					relPerm[i][j][k][2] = RELPERM0;
					blockFProps[i][j][k][BLOCK_PC][1] = sgt[0][CAPILLARYPRESSURE];
					Krog = sgt[0][KRO];

					dRelPerm[i][j][k][3] = 0;
					sdKrog = 0;
				}
				else if (n == Nsgt) {
					relPerm[i][j][k][2] = sgt[Nsgt - 1][KR];
					blockFProps[i][j][k][BLOCK_PC][1] = sgt[Nsgt - 1][CAPILLARYPRESSURE];
					Krog = sgt[Nsgt - 1][KRO];

					dRelPerm[i][j][k][3] = 0;
					sdKrog = 0;
				}
				else {
					relPerm[i][j][k][2] = (sgt[n][KR] - sgt[n - 1][KR]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(sat[i][j][k][2] - sgt[n][SAT]) + sgt[n][KR];
					Krog = (sgt[n][KRO] - sgt[n - 1][KRO]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(sat[i][j][k][2] - sgt[n][SAT]) + sgt[n][KRO];
					blockFProps[i][j][k][BLOCK_PC][1] = (sgt[n][CAPILLARYPRESSURE] - sgt[n - 1][CAPILLARYPRESSURE]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(sat[i][j][k][2] - sgt[n][SAT]) + sgt[n][CAPILLARYPRESSURE];

					dRelPerm[i][j][k][3] = (sgt[n][KR] - sgt[n - 1][KR]) / (sgt[n][SAT] - sgt[n - 1][SAT]);
					sdKrog = (sgt[n][KRO] - sgt[n - 1][KRO]) / (sgt[n][SAT] - sgt[n - 1][SAT]);
				}

				//Baker Model
				//satL=bsat[i][j][k][0]-swt[0][SAT];
				//satG=bsat[i][j][k][2]-sgt[0][SAT];

				//if (bsat[i][j][k][1]) relPerm[i][j][k][1]=(satL*Krow+satG*Krog)/(satL+satG);
				//else relPerm[i][j][k][1]=0;

				//Stone II model
				//relPerm[i][j][k][1]=swt[0][KRO]*((relPerm[i][j][k][0]+Krow/swt[0][KRO])*(relPerm[i][j][k][2]+Krog/swt[0][KRO])-relPerm[i][j][k][0]-relPerm[i][j][k][2]);

				AA = relPerm[i][j][k][0] + Krow / swt[0][KRO];
				BB = relPerm[i][j][k][2] + Krog / swt[0][KRO];
				dAA = (dRelPerm[i][j][k][0] + sdKrow / swt[0][KRO]);
				dBB = (dRelPerm[i][j][k][3] + sdKrog / swt[0][KRO]);

				relPerm[i][j][k][1] = swt[0][KRO] * (AA*BB - relPerm[i][j][k][0] - relPerm[i][j][k][2]);
				dRelPerm[i][j][k][1] = swt[0][KRO] * (dAA*BB - dRelPerm[i][j][k][0]);
				dRelPerm[i][j][k][2] = swt[0][KRO] * (AA*dBB - dRelPerm[i][j][k][3]);
				////////////////////////////////////////////////
				///////////////////////////////////////////////////


				if (Bift[i][j][k] < MISCIBLEIFT) {
					Fift = pow(Bift[i][j][k] / IFTBASECASE, 1.0 / IFTPOWER);
					relPerm[i][j][k][1] = Fift * relPerm[i][j][k][1] + (1 - Fift)*sat[i][j][k][1] * (1 - swt[0][SAT]);
					relPerm[i][j][k][2] = Fift * relPerm[i][j][k][2] + (1 - Fift)*sat[i][j][k][2] * (1 - swt[0][SAT]);;

					dRelPerm[i][j][k][3] = Fift * dRelPerm[i][j][k][3] + (1 - Fift)*(1 - swt[0][SAT]);
					dRelPerm[i][j][k][1] *= Fift;
					dRelPerm[i][j][k][2] *= Fift;
					dRelPerm[i][j][k][4] = (1 - Fift)*(1 - swt[0][SAT]);
				}

				if (relPerm[i][j][k][1] < 0) relPerm[i][j][k][1] = RELPERM0;
				if (relPerm[i][j][k][1] > 1) relPerm[i][j][k][1] = 1;
				if (dRelPerm[i][j][k][1] > 0) dRelPerm[i][j][k][1] = 0;
				if (dRelPerm[i][j][k][2] > 0) dRelPerm[i][j][k][2] = 0;
			}

	//fclose(fp);
}

void fCalcSatFs(void) {
	register int i, j, k;
	//register int n,
//	FType Krog, Krow, sdKrog, sdKrow;
//	FType AA, BB, dAA, dBB;
	//FType satL, satG;
//	FType Fift;

	//FType CoreySor=0;
	//FILE *fp;
	//fp=fopen("IND.txt", "wt");

	//swt[0][SAT]=0;


	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++) {
				/*relPerm[i][j][k][0]=RELPERM0;
				dRelPerm[i][j][k][0]=0;

				if (sat[i][j][k][1]<CoreySor) {
				relPerm[i][j][k][1]=RELPERM0;
				dRelPerm[i][j][k][1]=0;
				dRelPerm[i][j][k][2]=0;
				dRelPerm[i][j][k][4]=0;
				}
				else if (sat[i][j][k][1]>=1) {
				relPerm[i][j][k][1]=1;
				dRelPerm[i][j][k][1]=0;
				dRelPerm[i][j][k][2]=0;
				dRelPerm[i][j][k][4]=0;
				}
				else {
				relPerm[i][j][k][1]=pow((sat[i][j][k][1]-CoreySor)/(1-CoreySor), 3.0);
				dRelPerm[i][j][k][1]=0;
				dRelPerm[i][j][k][2]=0;
				dRelPerm[i][j][k][4]=(3.0/(1-CoreySor))*pow((sat[i][j][k][1]-CoreySor)/(1-CoreySor), 2.0);
				}
				//fprintf(fp, "%d\t%d\t%d\to=%f\t%f\n", i, j, k, relPerm[i][j][k][1], dRelPerm[i][j][k][4]);

				if (sat[i][j][k][2]<=0) {
				relPerm[i][j][k][2]=RELPERM0;
				dRelPerm[i][j][k][3]=0;
				blockFProps[i][j][k][BLOCK_PC][1]=0;
				}
				else if (sat[i][j][k][2]>=(1-CoreySor)) {
				relPerm[i][j][k][2]=1;
				dRelPerm[i][j][k][3]=0;
				blockFProps[i][j][k][BLOCK_PC][1]=3000941.83744005;
				}
				else {
				relPerm[i][j][k][2]=pow(sat[i][j][k][2]/(1-CoreySor), 1.5);

				dRelPerm[i][j][k][3]=(1.5/(1-CoreySor))*pow(sat[i][j][k][2]/(1-CoreySor), 0.5);
				blockFProps[i][j][k][BLOCK_PC][1]=103400*pow((sat[i][j][k][1]-CoreySor)/(1-CoreySor), -1.0/2.7);
				}
				blockFProps[i][j][k][BLOCK_PC][0]=0;
				//fprintf(fp, "%d\t%d\t%d\tg=%f\t%f\t%f\n", i, j, k, relPerm[i][j][k][2], dRelPerm[i][j][k][3], blockFProps[i][j][k][BLOCK_PC][1]);




				fdRelPerm[i][j][k][4] = 1;			//Oil to Oil
				//Relative permeability, Baker's model & Capillary pressures
				for (n = 0; n<Nswt; n++) {
				if (fsat[i][j][k][0]<swt[n][SAT]) break;
				}
				if (n == 0) {
				frelPerm[i][j][k][0] = RELPERM0;
				fblockFProps[i][j][k][BLOCK_PC][0] = swt[0][CAPILLARYPRESSURE];
				Krow = swt[0][KRO];

				fdRelPerm[i][j][k][0] = 1;
				sdKrow = 0;
				}
				else if (n == Nswt) {
				frelPerm[i][j][k][0] = swt[Nswt - 1][KR];
				fblockFProps[i][j][k][BLOCK_PC][0] = swt[Nswt - 1][CAPILLARYPRESSURE];
				Krow = swt[Nswt - 1][KRO];

				fdRelPerm[i][j][k][0] = 1;
				sdKrow = 0;
				}
				else {
				frelPerm[i][j][k][0] = (swt[n][KR] - swt[n - 1][KR]) / (swt[n][SAT] - swt[n - 1][SAT])*(fsat[i][j][k][0] - swt[n][SAT]) + swt[n][KR];
				Krow = (swt[n][KRO] - swt[n - 1][KRO]) / (swt[n][SAT] - swt[n - 1][SAT])*(fsat[i][j][k][0] - swt[n][SAT]) + swt[n][KRO];
				fblockFProps[i][j][k][BLOCK_PC][0] = (swt[n][CAPILLARYPRESSURE] - swt[n - 1][CAPILLARYPRESSURE]) / (swt[n][SAT] - swt[n - 1][SAT])*(fsat[i][j][k][0] - swt[n][SAT]) + swt[n][CAPILLARYPRESSURE];

				fdRelPerm[i][j][k][0] = (swt[n][KR] - swt[n - 1][KR]) / (swt[n][SAT] - swt[n - 1][SAT]);
				sdKrow = (swt[n][KRO] - swt[n - 1][KRO]) / (swt[n][SAT] - swt[n - 1][SAT]);
				}

				for (n = 0; n<Nsgt; n++) {
				if (fsat[i][j][k][2]<sgt[n][SAT]) break;
				}
				if (n == 0) {
				frelPerm[i][j][k][2] = RELPERM0;
				fblockFProps[i][j][k][BLOCK_PC][1] = sgt[0][CAPILLARYPRESSURE];
				Krog = sgt[0][KRO];

				fdRelPerm[i][j][k][3] = 1;
				sdKrog = 0;
				}
				else if (n == Nsgt) {
				frelPerm[i][j][k][2] = sgt[Nsgt - 1][KR];
				fblockFProps[i][j][k][BLOCK_PC][1] = sgt[Nsgt - 1][CAPILLARYPRESSURE];
				Krog = sgt[Nsgt - 1][KRO];

				fdRelPerm[i][j][k][3] = 1;
				sdKrog = 0;
				}
				else {
				frelPerm[i][j][k][2] = (sgt[n][KR] - sgt[n - 1][KR]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(fsat[i][j][k][2] - sgt[n][SAT]) + sgt[n][KR];
				Krog = (sgt[n][KRO] - sgt[n - 1][KRO]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(fsat[i][j][k][2] - sgt[n][SAT]) + sgt[n][KRO];
				fblockFProps[i][j][k][BLOCK_PC][1] = (sgt[n][CAPILLARYPRESSURE] - sgt[n - 1][CAPILLARYPRESSURE]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(fsat[i][j][k][2] - sgt[n][SAT]) + sgt[n][CAPILLARYPRESSURE];

				fdRelPerm[i][j][k][3] = (sgt[n][KR] - sgt[n - 1][KR]) / (sgt[n][SAT] - sgt[n - 1][SAT]);
				sdKrog = (sgt[n][KRO] - sgt[n - 1][KRO]) / (sgt[n][SAT] - sgt[n - 1][SAT]);
				}

				//Baker Model
				//satL=bsat[i][j][k][0]-swt[0][SAT];
				//satG=bsat[i][j][k][2]-sgt[0][SAT];

				//if (bsat[i][j][k][1]) relPerm[i][j][k][1]=(satL*Krow+satG*Krog)/(satL+satG);
				//else relPerm[i][j][k][1]=0;

				//Stone II model
				//relPerm[i][j][k][1]=swt[0][KRO]*((relPerm[i][j][k][0]+Krow/swt[0][KRO])*(relPerm[i][j][k][2]+Krog/swt[0][KRO])-relPerm[i][j][k][0]-relPerm[i][j][k][2]);

				AA = frelPerm[i][j][k][0] + Krow / swt[0][KRO];
				BB = frelPerm[i][j][k][2] + Krog / swt[0][KRO];
				dAA = (fdRelPerm[i][j][k][0] + sdKrow / swt[0][KRO]);
				dBB = (fdRelPerm[i][j][k][3] + sdKrog / swt[0][KRO]);

				frelPerm[i][j][k][1] = swt[0][KRO] * (AA*BB - frelPerm[i][j][k][0] - frelPerm[i][j][k][2]);
				fdRelPerm[i][j][k][1] = swt[0][KRO] * (dAA*BB - fdRelPerm[i][j][k][0]);
				fdRelPerm[i][j][k][2] = swt[0][KRO] * (AA*dBB - fdRelPerm[i][j][k][3]);
				////////////////////////////////////////////////
				///////////////////////////////////////////////////


				if (fBift[i][j][k]<MISCIBLEIFT) {
				Fift = pow(Bift[i][j][k] / IFTBASECASE, 1.0 / IFTPOWER);
				frelPerm[i][j][k][1] = Fift*frelPerm[i][j][k][1] + (1 - Fift)*fsat[i][j][k][1] * (1 - swt[0][SAT]);
				frelPerm[i][j][k][2] = Fift*frelPerm[i][j][k][2] + (1 - Fift)*fsat[i][j][k][2] * (1 - swt[0][SAT]);;

				fdRelPerm[i][j][k][3] = Fift*fdRelPerm[i][j][k][3] + (1 - Fift)*(1 - swt[0][SAT]);
				fdRelPerm[i][j][k][1] *= Fift;
				fdRelPerm[i][j][k][2] *= Fift;
				fdRelPerm[i][j][k][4] = (1 - Fift)*(1 - swt[0][SAT]);
				}

				if (frelPerm[i][j][k][1]<0) frelPerm[i][j][k][1] = RELPERM0;
				if (frelPerm[i][j][k][1]>1) frelPerm[i][j][k][1] = 1;
				if (fdRelPerm[i][j][k][1]>0) fdRelPerm[i][j][k][1] = 1;
				if (fdRelPerm[i][j][k][2]>0) fdRelPerm[i][j][k][2] = -1;*/

				frelPerm[i][j][k][0] = fsat[i][j][k][0];//-swt[0][0];
				frelPerm[i][j][k][1] = fsat[i][j][k][1];
				frelPerm[i][j][k][2] = fsat[i][j][k][2];

				fdRelPerm[i][j][k][0] = 1;
				fdRelPerm[i][j][k][1] = 0;
				fdRelPerm[i][j][k][2] = 0;
				fdRelPerm[i][j][k][3] = 1;
				fdRelPerm[i][j][k][4] = 1;


				fblockFProps[i][j][k][BLOCK_PC][0] = 0;
				fblockFProps[i][j][k][BLOCK_PC][1] = 0;

			}

	//fclose(fp);
}
void Calc_Visco(void) {
	FType s1Ml, s2Ml, s1Mg, s2Mg, tMl, tMg, MuL, MuG, roL, roG, compTCL, compMWL, compPCL, compTCG, compMWG, compPCG;
	FType ethaL, ethaG, eAAl, eAAg, eBBl, eBBg, eCCl, eCCg;
	FType ErL, ErG, ro1L, ro1G;

	FType Ul, Wl, Al, Bl, Zl, Ug, Wg, Ag, Bg, Zg;
	FType Cal, Cbl, Cag, Cbg;
	FType dAdxl, dBdxl, dZdxl, dUdxl, dWdxl, dAdxg, dBdxg, dZdxg, dUdxg, dWdxg;
	FType dAdpl, dBdpl, dUdpl, dWdpl, dAdpg, dBdpg, dUdpg, dWdpg, dZdpl, dZdpg;

	register int i, j, k, n, Ix, Iy, Iz;
	FType tempSQR;



	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {
				////////////////////////////////////////////////////////////////////////
				//Fanar derivative


				Al = 0;
				Bl = 0;
				Ag = 0;
				Bg = 0;
				for (i = 0; i < Nc; i++) {
					Bl += comp[Ix][Iy][Iz][i][0] * fluidProp[i][EOS_B];
					Bg += comp[Ix][Iy][Iz][i][1] * fluidProp[i][EOS_B];
					for (j = 0; j < Nc; j++) {
						tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
						Al += comp[Ix][Iy][Iz][i][0] * comp[Ix][Iy][Iz][j][0] * tempSQR;
						Ag += comp[Ix][Iy][Iz][i][1] * comp[Ix][Iy][Iz][j][1] * tempSQR;
					}
				}

				Cal = Al;
				Cbl = Bl;
				Cag = Ag;
				Cbg = Bg;
				Al *= P[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*RGAS*resTemp*resTemp);
				Bl *= P[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*resTemp);
				Ag *= P[Ix + 1][Iy + 1][Iz + 1][2] / (RGAS*RGAS*resTemp*resTemp);
				Bg *= P[Ix + 1][Iy + 1][Iz + 1][2] / (RGAS*resTemp);

				dAdpl = Cal / (RGAS*RGAS*resTemp*resTemp);
				dAdpg = Cag / (RGAS*RGAS*resTemp*resTemp);
				dBdpl = Cbl / (RGAS*resTemp);
				dBdpg = Cbg / (RGAS*resTemp);


				if (SRK) {
					Ul = Bl;
					Wl = 0;
					Ug = Bg;
					Wg = 0;

					dUdpl = dBdpl;
					dWdpl = 0;
					dUdpg = dBdpg;
					dWdpg = 0;
				}
				else if (PR) {
					Ul = 2 * Bl;
					Wl = Bl;
					Ug = 2 * Bg;
					Wg = Bg;

					dUdpl = 2 * dBdpl;
					dWdpl = dBdpl;
					dUdpg = 2 * dBdpg;
					dWdpg = dBdpg;
				}

				Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl * Ul - Ul - Wl * Wl, -(Al*Bl - Bl * Wl*Wl - Wl * Wl), 'l');
				Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg * Ug - Ug - Wg * Wg, -(Ag*Bg - Bg * Wg*Wg - Wg * Wg), 'g');



				dZdpg = -(dAdpg*(Zg - Bg) + dBdpg * (-Zg * Zg - Ug * Zg - Ag + Wg * Wg) + dUdpg * (Zg*Zg - Bg * Zg - Zg) + dWdpg * (-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg * Ug - Ug - Wg * Wg);		//Modified
				dZdpl = -(dAdpl*(Zl - Bl) + dBdpl * (-Zl * Zl - Ul * Zl - Al + Wl * Wl) + dUdpl * (Zl*Zl - Bl * Zl - Zl) + dWdpl * (-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl * Ul - Ul - Wl * Wl);		//Modified


				dE[Ix][Iy][Iz][Nc][0] = (1 - P[Ix + 1][Iy + 1][Iz + 1][1] * dZdpl / Zl) / (RGAS*resTemp*Zl);
				dE[Ix][Iy][Iz][Nc][1] = (1 - P[Ix + 1][Iy + 1][Iz + 1][2] * dZdpg / Zg) / (RGAS*resTemp*Zg);


				for (k = 0; k < Nc; k++) {
					dAdxl = 0;
					dAdxg = 0;
					for (n = 0; n < Nc; n++) {
						tempSQR = sqrt(fluidProp[n][EOS_A] * fluidProp[k][EOS_A])*bic[n][k];
						dAdxl += comp[Ix][Iy][Iz][n][0] * tempSQR;
						dAdxg += comp[Ix][Iy][Iz][n][1] * tempSQR;
					}
					dAdxl *= 2 * P[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*RGAS*resTemp*resTemp);
					dAdxg *= 2 * P[Ix + 1][Iy + 1][Iz + 1][2] / (RGAS*RGAS*resTemp*resTemp);
					dBdxl = fluidProp[k][EOS_B] * P[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*resTemp);
					dBdxg = fluidProp[k][EOS_B] * P[Ix + 1][Iy + 1][Iz + 1][2] / (RGAS*resTemp);

					if (SRK) {
						dUdxl = dBdxl;
						dWdxl = 0;
						dUdxg = dBdxg;
						dWdxg = 0;
					}
					else if (PR) {
						dUdxl = 2 * dBdxl;
						dWdxl = dBdxl;
						dUdxg = 2 * dBdxg;
						dWdxg = dBdxg;
					}

					dZdxg = -(dAdxg*(Zg - Bg) + dBdxg * (-Zg * Zg - Ug * Zg - Ag + Wg * Wg) + dUdxg * (Zg*Zg - Bg * Zg - Zg) + dWdxg * (-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg * Ug - Ug - Wg * Wg);		//Modified
					dZdxl = -(dAdxl*(Zl - Bl) + dBdxl * (-Zl * Zl - Ul * Zl - Al + Wl * Wl) + dUdxl * (Zl*Zl - Bl * Zl - Zl) + dWdxl * (-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl * Ul - Ul - Wl * Wl);		//Modified


					dE[Ix][Iy][Iz][k][0] = -(P[Ix + 1][Iy + 1][Iz + 1][1] * dZdxl / (RGAS*resTemp*Zl*Zl));
					dE[Ix][Iy][Iz][k][1] = -(P[Ix + 1][Iy + 1][Iz + 1][2] * dZdxg / (RGAS*resTemp*Zg*Zg));

				}
				////////////////////////////////////////////////////////////////////////////////////

				blockFProps[Ix][Iy][Iz][RO][0] = P[Ix + 1][Iy + 1][Iz + 1][1] / (Zl*RGAS*resTemp);
				blockFProps[Ix][Iy][Iz][RO][1] = P[Ix + 1][Iy + 1][Iz + 1][2] / (Zg*RGAS*resTemp);

				s1Ml = 0;
				s2Ml = 0;
				s1Mg = 0;
				s2Mg = 0;
				roL = 0;
				roG = 0;
				compMWL = 0;
				compPCL = 0;
				compTCL = 0;
				compMWG = 0;
				compPCG = 0;
				compTCG = 0;
				for (n = 0; n < Nc; n++) {
					tempSQR = sqrt(fluidProp[n][MW]);
					tMl = comp[Ix][Iy][Iz][n][0] * tempSQR;
					tMg = comp[Ix][Iy][Iz][n][1] * tempSQR;
					s1Ml += tMl * fluidProp[n][MU];
					s2Ml += tMl;
					s1Mg += tMg * fluidProp[n][MU];
					s2Mg += tMg;
					roL += comp[Ix][Iy][Iz][n][0] * fluidProp[n][VCRIT];
					roG += comp[Ix][Iy][Iz][n][1] * fluidProp[n][VCRIT];
					compTCL += comp[Ix][Iy][Iz][n][0] * fluidProp[n][TCRIT];
					compPCL += comp[Ix][Iy][Iz][n][0] * fluidProp[n][PCRIT] / 101325;
					compMWL += comp[Ix][Iy][Iz][n][0] * fluidProp[n][MW];
					compTCG += comp[Ix][Iy][Iz][n][1] * fluidProp[n][TCRIT];
					compPCG += comp[Ix][Iy][Iz][n][1] * fluidProp[n][PCRIT] / 101325;
					compMWG += comp[Ix][Iy][Iz][n][1] * fluidProp[n][MW];
				}
				MuL = s1Ml / s2Ml;
				MuG = s1Mg / s2Mg;
				ErL = blockFProps[Ix][Iy][Iz][RO][0] * roL;
				ErG = blockFProps[Ix][Iy][Iz][RO][1] * roG;

				eAAl = pow(compTCL, 1.0 / 6);
				eAAg = pow(compTCG, 1.0 / 6);

				eBBl = sqrt(compMWL);
				eBBg = sqrt(compMWG);

				eCCl = pow(compPCL, 2.0 / 3);
				eCCg = pow(compPCG, 2.0 / 3);

				ethaL = eAAl / (eBBl*eCCl);
				ethaG = eAAg / (eBBg*eCCg);

				ro1L = 0.10230 + 0.023364*ErL + 0.058533*ErL*ErL - 0.040758*ErL*ErL*ErL + 0.0093324*ErL*ErL*ErL*ErL;
				ro1G = 0.10230 + 0.023364*ErG + 0.058533*ErG*ErG - 0.040758*ErG*ErG*ErG + 0.0093324*ErG*ErG*ErG*ErG;

				blockFProps[Ix][Iy][Iz][BMU][0] = MuL + (ro1L*ro1L*ro1L*ro1L - 1e-4) / ethaL;		//Viscosity, centipoise
				blockFProps[Ix][Iy][Iz][BMU][1] = MuG + (ro1G*ro1G*ro1G*ro1G - 1e-4) / ethaG;
			}
}


void fCalc_Visco(void) {
	FType s1Ml, s2Ml, s1Mg, s2Mg, tMl, tMg, MuL, MuG, roL, roG, compTCL, compMWL, compPCL, compTCG, compMWG, compPCG;
	FType ethaL, ethaG, eAAl, eAAg, eBBl, eBBg, eCCl, eCCg;
	FType ErL, ErG, ro1L, ro1G;

	FType Ul, Wl, Al, Bl, Zl, Ug, Wg, Ag, Bg, Zg;
	FType Cal, Cbl, Cag, Cbg;
	FType dAdxl, dBdxl, dZdxl, dUdxl, dWdxl, dAdxg, dBdxg, dZdxg, dUdxg, dWdxg;
	FType dAdpl, dBdpl, dUdpl, dWdpl, dAdpg, dBdpg, dUdpg, dWdpg, dZdpl, dZdpg;

	register int i, j, k, n, Ix, Iy, Iz;
	FType tempSQR;



	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {
				////////////////////////////////////////////////////////////////////////
				//Fanar derivative


				Al = 0;
				Bl = 0;
				Ag = 0;
				Bg = 0;
				for (i = 0; i < Nc; i++) {
					Bl += fcomp[Ix][Iy][Iz][i][0] * fluidProp[i][EOS_B];
					Bg += fcomp[Ix][Iy][Iz][i][1] * fluidProp[i][EOS_B];
					for (j = 0; j < Nc; j++) {
						tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
						Al += fcomp[Ix][Iy][Iz][i][0] * fcomp[Ix][Iy][Iz][j][0] * tempSQR;
						Ag += fcomp[Ix][Iy][Iz][i][1] * fcomp[Ix][Iy][Iz][j][1] * tempSQR;
					}
				}

				Cal = Al;
				Cbl = Bl;
				Cag = Ag;
				Cbg = Bg;
				Al *= fP[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*RGAS*resTemp*resTemp);
				Bl *= fP[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*resTemp);
				Ag *= fP[Ix + 1][Iy + 1][Iz + 1][2] / (RGAS*RGAS*resTemp*resTemp);
				Bg *= fP[Ix + 1][Iy + 1][Iz + 1][2] / (RGAS*resTemp);

				dAdpl = Cal / (RGAS*RGAS*resTemp*resTemp);
				dAdpg = Cag / (RGAS*RGAS*resTemp*resTemp);
				dBdpl = Cbl / (RGAS*resTemp);
				dBdpg = Cbg / (RGAS*resTemp);


				if (SRK) {
					Ul = Bl;
					Wl = 0;
					Ug = Bg;
					Wg = 0;

					dUdpl = dBdpl;
					dWdpl = 0;
					dUdpg = dBdpg;
					dWdpg = 0;
				}
				else if (PR) {
					Ul = 2 * Bl;
					Wl = Bl;
					Ug = 2 * Bg;
					Wg = Bg;

					dUdpl = 2 * dBdpl;
					dWdpl = dBdpl;
					dUdpg = 2 * dBdpg;
					dWdpg = dBdpg;
				}

				Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl * Ul - Ul - Wl * Wl, -(Al*Bl - Bl * Wl*Wl - Wl * Wl), 'l');
				Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg * Ug - Ug - Wg * Wg, -(Ag*Bg - Bg * Wg*Wg - Wg * Wg), 'g');



				dZdpg = -(dAdpg*(Zg - Bg) + dBdpg * (-Zg * Zg - Ug * Zg - Ag + Wg * Wg) + dUdpg * (Zg*Zg - Bg * Zg - Zg) + dWdpg * (-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg * Ug - Ug - Wg * Wg);		//Modified
				dZdpl = -(dAdpl*(Zl - Bl) + dBdpl * (-Zl * Zl - Ul * Zl - Al + Wl * Wl) + dUdpl * (Zl*Zl - Bl * Zl - Zl) + dWdpl * (-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl * Ul - Ul - Wl * Wl);		//Modified


				fdE[Ix][Iy][Iz][Nc][0] = (1 - fP[Ix + 1][Iy + 1][Iz + 1][1] * dZdpl / Zl) / (RGAS*resTemp*Zl);
				fdE[Ix][Iy][Iz][Nc][1] = (1 - fP[Ix + 1][Iy + 1][Iz + 1][2] * dZdpg / Zg) / (RGAS*resTemp*Zg);


				for (k = 0; k < Nc; k++) {
					dAdxl = 0;
					dAdxg = 0;
					for (n = 0; n < Nc; n++) {
						tempSQR = sqrt(fluidProp[n][EOS_A] * fluidProp[k][EOS_A])*bic[n][k];
						dAdxl += fcomp[Ix][Iy][Iz][n][0] * tempSQR;
						dAdxg += fcomp[Ix][Iy][Iz][n][1] * tempSQR;
					}
					dAdxl *= 2 * fP[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*RGAS*resTemp*resTemp);
					dAdxg *= 2 * fP[Ix + 1][Iy + 1][Iz + 1][2] / (RGAS*RGAS*resTemp*resTemp);
					dBdxl = fluidProp[k][EOS_B] * fP[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*resTemp);
					dBdxg = fluidProp[k][EOS_B] * fP[Ix + 1][Iy + 1][Iz + 1][2] / (RGAS*resTemp);

					if (SRK) {
						dUdxl = dBdxl;
						dWdxl = 0;
						dUdxg = dBdxg;
						dWdxg = 0;
					}
					else if (PR) {
						dUdxl = 2 * dBdxl;
						dWdxl = dBdxl;
						dUdxg = 2 * dBdxg;
						dWdxg = dBdxg;
					}

					dZdxg = -(dAdxg*(Zg - Bg) + dBdxg * (-Zg * Zg - Ug * Zg - Ag + Wg * Wg) + dUdxg * (Zg*Zg - Bg * Zg - Zg) + dWdxg * (-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg * Ug - Ug - Wg * Wg);		//Modified
					dZdxl = -(dAdxl*(Zl - Bl) + dBdxl * (-Zl * Zl - Ul * Zl - Al + Wl * Wl) + dUdxl * (Zl*Zl - Bl * Zl - Zl) + dWdxl * (-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl * Ul - Ul - Wl * Wl);		//Modified


					fdE[Ix][Iy][Iz][k][0] = -(fP[Ix + 1][Iy + 1][Iz + 1][1] * dZdxl / (RGAS*resTemp*Zl*Zl));
					fdE[Ix][Iy][Iz][k][1] = -(fP[Ix + 1][Iy + 1][Iz + 1][2] * dZdxg / (RGAS*resTemp*Zg*Zg));

				}
				////////////////////////////////////////////////////////////////////////////////////

				fblockFProps[Ix][Iy][Iz][RO][0] = fP[Ix + 1][Iy + 1][Iz + 1][1] / (Zl*RGAS*resTemp);
				fblockFProps[Ix][Iy][Iz][RO][1] = fP[Ix + 1][Iy + 1][Iz + 1][2] / (Zg*RGAS*resTemp);

				s1Ml = 0;
				s2Ml = 0;
				s1Mg = 0;
				s2Mg = 0;
				roL = 0;
				roG = 0;
				compMWL = 0;
				compPCL = 0;
				compTCL = 0;
				compMWG = 0;
				compPCG = 0;
				compTCG = 0;
				for (n = 0; n < Nc; n++) {
					tempSQR = sqrt(fluidProp[n][MW]);
					tMl = fcomp[Ix][Iy][Iz][n][0] * tempSQR;
					tMg = fcomp[Ix][Iy][Iz][n][1] * tempSQR;
					s1Ml += tMl * fluidProp[n][MU];
					s2Ml += tMl;
					s1Mg += tMg * fluidProp[n][MU];
					s2Mg += tMg;
					roL += fcomp[Ix][Iy][Iz][n][0] * fluidProp[n][VCRIT];
					roG += fcomp[Ix][Iy][Iz][n][1] * fluidProp[n][VCRIT];
					compTCL += fcomp[Ix][Iy][Iz][n][0] * fluidProp[n][TCRIT];
					compPCL += fcomp[Ix][Iy][Iz][n][0] * fluidProp[n][PCRIT] / 101325;
					compMWL += fcomp[Ix][Iy][Iz][n][0] * fluidProp[n][MW];
					compTCG += fcomp[Ix][Iy][Iz][n][1] * fluidProp[n][TCRIT];
					compPCG += fcomp[Ix][Iy][Iz][n][1] * fluidProp[n][PCRIT] / 101325;
					compMWG += fcomp[Ix][Iy][Iz][n][1] * fluidProp[n][MW];
				}
				MuL = s1Ml / s2Ml;
				MuG = s1Mg / s2Mg;
				ErL = fblockFProps[Ix][Iy][Iz][RO][0] * roL;
				ErG = fblockFProps[Ix][Iy][Iz][RO][1] * roG;

				eAAl = pow(compTCL, 1.0 / 6);
				eAAg = pow(compTCG, 1.0 / 6);

				eBBl = sqrt(compMWL);
				eBBg = sqrt(compMWG);

				eCCl = pow(compPCL, 2.0 / 3);
				eCCg = pow(compPCG, 2.0 / 3);

				ethaL = eAAl / (eBBl*eCCl);
				ethaG = eAAg / (eBBg*eCCg);

				ro1L = 0.10230 + 0.023364*ErL + 0.058533*ErL*ErL - 0.040758*ErL*ErL*ErL + 0.0093324*ErL*ErL*ErL*ErL;
				ro1G = 0.10230 + 0.023364*ErG + 0.058533*ErG*ErG - 0.040758*ErG*ErG*ErG + 0.0093324*ErG*ErG*ErG*ErG;

				fblockFProps[Ix][Iy][Iz][BMU][0] = MuL + (ro1L*ro1L*ro1L*ro1L - 1e-4) / ethaL;		//Viscosity, centipoise
				fblockFProps[Ix][Iy][Iz][BMU][1] = MuG + (ro1G*ro1G*ro1G*ro1G - 1e-4) / ethaG;
			}
}


void Calc_Transes(void) {
	register int Ix, Iy, Iz, i, j, k;
	FType mGam, pot;
	char Coef1, Coef2;
	//FType dKrw, dKrg, dKrow, dKrog;

	for (Ix = 0; Ix < (Nx); Ix++) {
		i = Ix + 1;
		for (Iy = 0; Iy < (Ny); Iy++) {
			j = Iy + 1;
			for (Iz = 0; Iz < (Nz); Iz++) {
				k = Iz + 1;

				//X-Direction
				if ((Nx > 1) && (i < Nx)) {
					pot = P[i][j][k][0] - P[i + 1][j][k][0] + G_ACC * watRo*(blockH[i][j][k] - blockH[i + 1][j][k]);
					if (pot > 0) {
						Wtran[i][j][k][0] = IFTran[Ix][Iy][Iz][0] * (watRo*WAT_M_RO / watMu)*relPerm[Ix][Iy][Iz][0];
						dWtran[i][j][k][0] = IFTran[Ix][Iy][Iz][0] * (watRo*WAT_M_RO / watMu)*dRelPerm[Ix][Iy][Iz][0];
						transS[i][j][k][2][0] = 0;
					}
					else {
						Wtran[i][j][k][0] = IFTran[Ix][Iy][Iz][0] * (watRo*WAT_M_RO / watMu)*relPerm[Ix + 1][Iy][Iz][0];
						dWtran[i][j][k][0] = IFTran[Ix][Iy][Iz][0] * (watRo*WAT_M_RO / watMu)*dRelPerm[Ix + 1][Iy][Iz][0];
						transS[i][j][k][2][0] = 1;
					}

					if (phaseStat[Ix][Iy][Iz] == (-1)) Coef1 = 0;
					else Coef1 = 1;
					if (phaseStat[Ix + 1][Iy][Iz] == (-1)) Coef2 = 0;
					else Coef2 = 1;

					mGam = (blockFProps[Ix][Iy][Iz][BGAM][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BGAM][0] * gridDim[Ix + 1] * Coef2) / (Coef1*gridDim[Ix] + gridDim[Ix + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = P[i][j][k][1] - P[i + 1][j][k][1] + mGam * (blockH[i][j][k] - blockH[i + 1][j][k]);
					if (pot > 0) {
						trans[i][j][k][0][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*relPerm[Ix][Iy][Iz][1];
						transS[i][j][k][0][0] = 0;
						trans[i][j][k][2][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*relPerm[Ix][Iy][Iz][1];

						trans[i][j][k][4][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][1];
						trans[i][j][k][6][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][1];
						trans[i][j][k][8][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][2];
						trans[i][j][k][9][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][2];

						trans[i][j][k][10][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][4];
						trans[i][j][k][11][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][4];
					}
					else {
						trans[i][j][k][0][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*relPerm[Ix + 1][Iy][Iz][1];
						transS[i][j][k][0][0] = 1;
						trans[i][j][k][2][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*relPerm[Ix + 1][Iy][Iz][1];

						trans[i][j][k][4][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix + 1][Iy][Iz][1];
						trans[i][j][k][6][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix + 1][Iy][Iz][1];
						trans[i][j][k][8][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix + 1][Iy][Iz][2];
						trans[i][j][k][9][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix + 1][Iy][Iz][2];

						trans[i][j][k][10][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix + 1][Iy][Iz][4];
						trans[i][j][k][11][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix + 1][Iy][Iz][4];
					}
					if ((!Coef1) && (!Coef2)) {
						trans[i][j][k][0][0] = 0;
						trans[i][j][k][2][0] = 0;
						trans[i][j][k][4][0] = 0;
						trans[i][j][k][6][0] = 0;
						trans[i][j][k][8][0] = 0;
						trans[i][j][k][9][0] = 0;

						//transS[i][j][k][0][0]=0;
					}

					if (phaseStat[Ix][Iy][Iz] == 1) Coef1 = 0;
					else Coef1 = 1;
					if (phaseStat[Ix + 1][Iy][Iz] == 1) Coef2 = 0;
					else Coef2 = 1;

					mGam = (blockFProps[Ix][Iy][Iz][BGAM][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BGAM][1] * gridDim[Ix + 1] * Coef2) / (gridDim[Ix] * Coef1 + gridDim[Ix + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = P[i][j][k][2] - P[i + 1][j][k][2] + mGam * (blockH[i][j][k] - blockH[i + 1][j][k]);
					if (pot > 0) {
						trans[i][j][k][1][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*relPerm[Ix][Iy][Iz][2];
						transS[i][j][k][1][0] = 0;
						trans[i][j][k][3][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*relPerm[Ix][Iy][Iz][2];

						trans[i][j][k][5][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][3];
						trans[i][j][k][7][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][3];

					}
					else {
						trans[i][j][k][1][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*relPerm[Ix + 1][Iy][Iz][2];
						transS[i][j][k][1][0] = 1;
						trans[i][j][k][3][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*relPerm[Ix + 1][Iy][Iz][2];

						trans[i][j][k][5][0] = IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix + 1][Iy][Iz][3];
						trans[i][j][k][7][0] = mGam * IFTran[Ix][Iy][Iz][0] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + blockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*dRelPerm[Ix + 1][Iy][Iz][3];
					}

					if ((!Coef1) && (!Coef2)) {
						trans[i][j][k][1][0] = 0;
						trans[i][j][k][3][0] = 0;
						trans[i][j][k][5][0] = 0;
						trans[i][j][k][7][0] = 0;

						//transS[i][j][k][1][0]=0;
					}

				}


				//Y-Direction
				if ((Ny > 1) && (j < Ny)) {
					pot = P[i][j][k][0] - P[i][j + 1][k][0] + G_ACC * watRo*(blockH[i][j][k] - blockH[i][j + 1][k]);
					if (pot > 0) {
						Wtran[i][j][k][1] = IFTran[Ix][Iy][Iz][1] * (watRo*WAT_M_RO / watMu)*relPerm[Ix][Iy][Iz][0];
						dWtran[i][j][k][1] = IFTran[Ix][Iy][Iz][1] * (watRo*WAT_M_RO / watMu)*dRelPerm[Ix][Iy][Iz][0];
						transS[i][j][k][2][1] = 0;
					}
					else {
						Wtran[i][j][k][1] = IFTran[Ix][Iy][Iz][1] * (watRo*WAT_M_RO / watMu)*relPerm[Ix][Iy + 1][Iz][0];
						dWtran[i][j][k][1] = IFTran[Ix][Iy][Iz][1] * (watRo*WAT_M_RO / watMu)*dRelPerm[Ix][Iy + 1][Iz][0];
						transS[i][j][k][2][1] = 1;
					}

					if (phaseStat[Ix][Iy][Iz] == (-1)) Coef1 = 0;
					else Coef1 = 1;
					if (phaseStat[Ix][Iy + 1][Iz] == (-1)) Coef2 = 0;
					else Coef2 = 1;

					mGam = (blockFProps[Ix][Iy][Iz][BGAM][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BGAM][0] * gridDim[Nx + Iy + 1] * Coef2) / (gridDim[Nx + Iy] * Coef1 + gridDim[Nx + Iy + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = P[i][j][k][1] - P[i][j + 1][k][1] + mGam * (blockH[i][j][k] - blockH[i][j + 1][k]);
					if (pot > 0) {
						trans[i][j][k][0][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*relPerm[Ix][Iy][Iz][1];
						transS[i][j][k][0][1] = 0;
						trans[i][j][k][2][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*relPerm[Ix][Iy][Iz][1];


						trans[i][j][k][4][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][1];
						trans[i][j][k][6][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][1];
						trans[i][j][k][8][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][2];
						trans[i][j][k][9][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][2];

						trans[i][j][k][10][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][4];
						trans[i][j][k][11][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][4];
					}
					else {
						trans[i][j][k][0][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*relPerm[Ix][Iy + 1][Iz][1];
						transS[i][j][k][0][1] = 1;
						trans[i][j][k][2][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*relPerm[Ix][Iy + 1][Iz][1];

						trans[i][j][k][4][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy + 1][Iz][1];
						trans[i][j][k][6][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy + 1][Iz][1];
						trans[i][j][k][8][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy + 1][Iz][2];
						trans[i][j][k][9][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy + 1][Iz][2];

						trans[i][j][k][10][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy + 1][Iz][4];
						trans[i][j][k][11][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy + 1][Iz][4];

					}

					if ((!Coef1) && (!Coef2)) {
						trans[i][j][k][0][1] = 0;
						trans[i][j][k][2][1] = 0;
						trans[i][j][k][4][1] = 0;
						trans[i][j][k][6][1] = 0;
						trans[i][j][k][8][1] = 0;
						trans[i][j][k][9][1] = 0;

						//transS[i][j][k][0][1]=0;
					}

					if (phaseStat[Ix][Iy][Iz] == 1) Coef1 = 0;
					else Coef1 = 1;
					if (phaseStat[Ix][Iy + 1][Iz] == 1) Coef2 = 0;
					else Coef2 = 1;

					mGam = (blockFProps[Ix][Iy][Iz][BGAM][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BGAM][1] * gridDim[Nx + Iy + 1] * Coef2) / (gridDim[Nx + Iy] * Coef1 + gridDim[Nx + Iy + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = P[i][j][k][2] - P[i][j + 1][k][2] + mGam * (blockH[i][j][k] - blockH[i][j + 1][k]);
					if (pot > 0) {
						trans[i][j][k][1][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*relPerm[Ix][Iy][Iz][2];
						transS[i][j][k][1][1] = 0;
						trans[i][j][k][3][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*relPerm[Ix][Iy][Iz][2];

						trans[i][j][k][5][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][3];
						trans[i][j][k][7][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][3];

					}
					else {
						trans[i][j][k][1][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*relPerm[Ix][Iy + 1][Iz][2];
						transS[i][j][k][1][1] = 1;
						trans[i][j][k][3][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*relPerm[Ix][Iy + 1][Iz][2];

						trans[i][j][k][5][1] = IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy + 1][Iz][3];
						trans[i][j][k][7][1] = mGam * IFTran[Ix][Iy][Iz][1] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + blockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*dRelPerm[Ix][Iy + 1][Iz][3];
					}

					if ((!Coef1) && (!Coef2)) {
						trans[i][j][k][1][1] = 0;
						trans[i][j][k][3][1] = 0;
						trans[i][j][k][5][1] = 0;
						trans[i][j][k][7][1] = 0;

						//transS[i][j][k][1][1]=0;
					}

				}


				//Z-Direction
				if ((Nz > 1) && (k < Nz)) {
					pot = P[i][j][k][0] - P[i][j][k + 1][0] + G_ACC * watRo*(blockH[i][j][k] - blockH[i][j][k + 1]);
					if (pot > 0) {
						Wtran[i][j][k][2] = IFTran[Ix][Iy][Iz][2] * (watRo*WAT_M_RO / watMu)*relPerm[Ix][Iy][Iz][0];
						dWtran[i][j][k][2] = IFTran[Ix][Iy][Iz][2] * (watRo*WAT_M_RO / watMu)*dRelPerm[Ix][Iy][Iz][0];
						transS[i][j][k][2][2] = 0;
					}
					else {
						Wtran[i][j][k][2] = IFTran[Ix][Iy][Iz][2] * (watRo*WAT_M_RO / watMu)*relPerm[Ix][Iy][Iz + 1][0];
						dWtran[i][j][k][2] = IFTran[Ix][Iy][Iz][2] * (watRo*WAT_M_RO / watMu)*dRelPerm[Ix][Iy][Iz + 1][0];
						transS[i][j][k][2][2] = 1;
					}

					if (phaseStat[Ix][Iy][Iz] == (-1)) Coef1 = 0;
					else Coef1 = 1;
					if (phaseStat[Ix][Iy][Iz + 1] == (-1)) Coef2 = 0;
					else Coef2 = 1;


					mGam = (blockFProps[Ix][Iy][Iz][BGAM][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BGAM][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (gridDim[Nx + Ny + Iz] * Coef1 + gridDim[Nx + Ny + Iz + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = P[i][j][k][1] - P[i][j][k + 1][1] + mGam * (blockH[i][j][k] - blockH[i][j][k + 1]);
					if (pot > 0) {
						trans[i][j][k][0][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*relPerm[Ix][Iy][Iz][1];
						transS[i][j][k][0][2] = 0;
						trans[i][j][k][2][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*relPerm[Ix][Iy][Iz][1];

						trans[i][j][k][4][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][1];
						trans[i][j][k][6][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][1];
						trans[i][j][k][8][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][2];
						trans[i][j][k][9][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][2];

						trans[i][j][k][10][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][4];
						trans[i][j][k][11][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][4];
					}
					else {
						trans[i][j][k][0][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*relPerm[Ix][Iy][Iz + 1][1];
						transS[i][j][k][0][2] = 1;
						trans[i][j][k][2][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*relPerm[Ix][Iy][Iz + 1][1];

						trans[i][j][k][4][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz + 1][1];
						trans[i][j][k][6][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz + 1][1];
						trans[i][j][k][8][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz + 1][2];
						trans[i][j][k][9][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz + 1][2];

						trans[i][j][k][10][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][4];
						trans[i][j][k][11][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][4];
					}
					if ((!Coef1) && (!Coef2)) {
						trans[i][j][k][0][2] = 0;
						trans[i][j][k][2][2] = 0;
						trans[i][j][k][4][2] = 0;
						trans[i][j][k][6][2] = 0;
						trans[i][j][k][8][2] = 0;
						trans[i][j][k][9][2] = 0;

						//transS[i][j][k][0][2]=0;
					}

					if (phaseStat[Ix][Iy][Iz] == 1) Coef1 = 0;
					else Coef1 = 1;
					if (phaseStat[Ix][Iy][Iz + 1] == 1) Coef2 = 0;
					else Coef2 = 1;

					mGam = (blockFProps[Ix][Iy][Iz][BGAM][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BGAM][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (gridDim[Nx + Ny + Iz] * Coef1 + gridDim[Nx + Ny + Iz + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = P[i][j][k][2] - P[i][j][k + 1][2] + mGam * (blockH[i][j][k] - blockH[i][j][k + 1]);
					if (pot > 0) {
						trans[i][j][k][1][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*relPerm[Ix][Iy][Iz][2];
						transS[i][j][k][1][2] = 0;
						trans[i][j][k][3][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*relPerm[Ix][Iy][Iz][2];

						trans[i][j][k][5][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][3];
						trans[i][j][k][7][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz][3];

					}
					else {
						trans[i][j][k][1][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*relPerm[Ix][Iy][Iz + 1][2];
						transS[i][j][k][1][2] = 1;
						trans[i][j][k][3][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*relPerm[Ix][Iy][Iz + 1][2];

						trans[i][j][k][5][2] = IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz + 1][3];
						trans[i][j][k][7][2] = mGam * IFTran[Ix][Iy][Iz][2] * (blockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (blockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + blockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*dRelPerm[Ix][Iy][Iz + 1][3];
					}
					if ((!Coef1) && (!Coef2)) {
						trans[i][j][k][1][2] = 0;
						trans[i][j][k][3][2] = 0;
						trans[i][j][k][5][2] = 0;
						trans[i][j][k][7][2] = 0;

						//transS[i][j][k][1][2]=0;
					}

				}
			}
		}
	}
}
void fCalc_Transes(void) {
	register int Ix, Iy, Iz, i, j, k;
	FType mGam, pot;
	char Coef1, Coef2;
	//FType dKrw, dKrg, dKrow, dKrog;

	for (Ix = 0; Ix < (Nx); Ix++) {
		i = Ix + 1;
		for (Iy = 0; Iy < (Ny); Iy++) {
			j = Iy + 1;
			for (Iz = 0; Iz < (Nz); Iz++) {
				k = Iz + 1;

				//X-Direction
				if ((Nx > 1) && (i < Nx)) {
					pot = fP[i][j][k][0] - fP[i + 1][j][k][0] + G_ACC * watRo*(blockH[i][j][k] - blockH[i + 1][j][k]);
					if (pot > 0) {
						fWtran[i][j][k][0] = fIFTran[Ix][Iy][Iz][0] * (watRo*WAT_M_RO / watMu)*fsat[Ix][Iy][Iz][0];
						fdWtran[i][j][k][0] = fIFTran[Ix][Iy][Iz][0] * (watRo*WAT_M_RO / watMu);
						ftransS[i][j][k][2][0] = 0;
					}
					else {
						fWtran[i][j][k][0] = fIFTran[Ix][Iy][Iz][0] * (watRo*WAT_M_RO / watMu)*fsat[Ix + 1][Iy][Iz][0];
						fdWtran[i][j][k][0] = fIFTran[Ix][Iy][Iz][0] * (watRo*WAT_M_RO / watMu);
						ftransS[i][j][k][2][0] = 1;
					}

					if (fphaseStat[Ix][Iy][Iz] == (-1)) Coef1 = 0;
					else Coef1 = 1;
					if (fphaseStat[Ix + 1][Iy][Iz] == (-1)) Coef2 = 0;
					else Coef2 = 1;

					mGam = (fblockFProps[Ix][Iy][Iz][BGAM][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BGAM][0] * gridDim[Ix + 1] * Coef2) / (Coef1*gridDim[Ix] + gridDim[Ix + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = fP[i][j][k][1] - fP[i + 1][j][k][1] + mGam * (blockH[i][j][k] - blockH[i + 1][j][k]);
					if (pot > 0) {
						ftrans[i][j][k][0][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*fsat[Ix][Iy][Iz][1];
						ftransS[i][j][k][0][0] = 0;
						ftrans[i][j][k][2][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*fsat[Ix][Iy][Iz][1];

						ftrans[i][j][k][4][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*(-1);
						ftrans[i][j][k][6][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*(-1);
						ftrans[i][j][k][8][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*(-1);
						ftrans[i][j][k][9][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*(-1);

						ftrans[i][j][k][10][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2);
						ftrans[i][j][k][11][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2);
					}
					else {
						ftrans[i][j][k][0][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*fsat[Ix + 1][Iy][Iz][1];
						ftransS[i][j][k][0][0] = 1;
						ftrans[i][j][k][2][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*fsat[Ix + 1][Iy][Iz][1];

						ftrans[i][j][k][4][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*(-1);
						ftrans[i][j][k][6][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*(-1);
						ftrans[i][j][k][8][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*(-1);
						ftrans[i][j][k][9][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2)*(-1);

						ftrans[i][j][k][10][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2);
						ftrans[i][j][k][11][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][0] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][0] * gridDim[Ix + 1] * Coef2);
					}
					if ((!Coef1) && (!Coef2)) {
						ftrans[i][j][k][0][0] = 0;
						ftrans[i][j][k][2][0] = 0;
						ftrans[i][j][k][4][0] = 0;
						ftrans[i][j][k][6][0] = 0;
						ftrans[i][j][k][8][0] = 0;
						ftrans[i][j][k][9][0] = 0;

						//ftransS[i][j][k][0][0]=0;
					}

					if (fphaseStat[Ix][Iy][Iz] == 1) Coef1 = 0;
					else Coef1 = 1;
					if (fphaseStat[Ix + 1][Iy][Iz] == 1) Coef2 = 0;
					else Coef2 = 1;

					mGam = (fblockFProps[Ix][Iy][Iz][BGAM][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BGAM][1] * gridDim[Ix + 1] * Coef2) / (gridDim[Ix] * Coef1 + gridDim[Ix + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = fP[i][j][k][2] - fP[i + 1][j][k][2] + mGam * (blockH[i][j][k] - blockH[i + 1][j][k]);
					if (pot > 0) {
						ftrans[i][j][k][1][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*fsat[Ix][Iy][Iz][2];
						ftransS[i][j][k][1][0] = 0;
						ftrans[i][j][k][3][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*fsat[Ix][Iy][Iz][2];

						ftrans[i][j][k][5][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2);
						ftrans[i][j][k][7][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2);

					}
					else {
						ftrans[i][j][k][1][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*fsat[Ix + 1][Iy][Iz][2];
						ftransS[i][j][k][1][0] = 1;
						ftrans[i][j][k][3][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2)*fsat[Ix + 1][Iy][Iz][2];

						ftrans[i][j][k][5][0] = fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2);
						ftrans[i][j][k][7][0] = mGam * fIFTran[Ix][Iy][Iz][0] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][RO][1] * gridDim[Ix + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Ix] * Coef1 + fblockFProps[Ix + 1][Iy][Iz][BMU][1] * gridDim[Ix + 1] * Coef2);
					}

					if ((!Coef1) && (!Coef2)) {
						ftrans[i][j][k][1][0] = 0;
						ftrans[i][j][k][3][0] = 0;
						ftrans[i][j][k][5][0] = 0;
						ftrans[i][j][k][7][0] = 0;

						//ftransS[i][j][k][1][0]=0;
					}

				}


				//Y-Direction
				if ((Ny > 1) && (j < Ny)) {
					pot = fP[i][j][k][0] - fP[i][j + 1][k][0] + G_ACC * watRo*(blockH[i][j][k] - blockH[i][j + 1][k]);
					if (pot > 0) {
						fWtran[i][j][k][1] = fIFTran[Ix][Iy][Iz][1] * (watRo*WAT_M_RO / watMu)*fsat[Ix][Iy][Iz][0];
						fdWtran[i][j][k][1] = fIFTran[Ix][Iy][Iz][1] * (watRo*WAT_M_RO / watMu);
						ftransS[i][j][k][2][1] = 0;
					}
					else {
						fWtran[i][j][k][1] = fIFTran[Ix][Iy][Iz][1] * (watRo*WAT_M_RO / watMu)*fsat[Ix][Iy + 1][Iz][0];
						fdWtran[i][j][k][1] = fIFTran[Ix][Iy][Iz][1] * (watRo*WAT_M_RO / watMu);
						ftransS[i][j][k][2][1] = 1;
					}

					if (fphaseStat[Ix][Iy][Iz] == (-1)) Coef1 = 0;
					else Coef1 = 1;
					if (fphaseStat[Ix][Iy + 1][Iz] == (-1)) Coef2 = 0;
					else Coef2 = 1;

					mGam = (fblockFProps[Ix][Iy][Iz][BGAM][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BGAM][0] * gridDim[Nx + Iy + 1] * Coef2) / (gridDim[Nx + Iy] * Coef1 + gridDim[Nx + Iy + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = fP[i][j][k][1] - fP[i][j + 1][k][1] + mGam * (blockH[i][j][k] - blockH[i][j + 1][k]);
					if (pot > 0) {
						ftrans[i][j][k][0][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*fsat[Ix][Iy][Iz][1];
						ftransS[i][j][k][0][1] = 0;
						ftrans[i][j][k][2][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*fsat[Ix][Iy][Iz][1];


						ftrans[i][j][k][4][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*(-1);
						ftrans[i][j][k][6][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*(-1);
						ftrans[i][j][k][8][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*(-1);
						ftrans[i][j][k][9][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*(-1);

						ftrans[i][j][k][10][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2);
						ftrans[i][j][k][11][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2);
					}
					else {
						ftrans[i][j][k][0][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*fsat[Ix][Iy + 1][Iz][1];
						ftransS[i][j][k][0][1] = 1;
						ftrans[i][j][k][2][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*fsat[Ix][Iy + 1][Iz][1];

						ftrans[i][j][k][4][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*(-1);
						ftrans[i][j][k][6][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*(-1);
						ftrans[i][j][k][8][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*(-1);
						ftrans[i][j][k][9][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2)*(-1);

						ftrans[i][j][k][10][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2);
						ftrans[i][j][k][11][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][0] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][0] * gridDim[Nx + Iy + 1] * Coef2);

					}

					if ((!Coef1) && (!Coef2)) {
						ftrans[i][j][k][0][1] = 0;
						ftrans[i][j][k][2][1] = 0;
						ftrans[i][j][k][4][1] = 0;
						ftrans[i][j][k][6][1] = 0;
						ftrans[i][j][k][8][1] = 0;
						ftrans[i][j][k][9][1] = 0;

						//ftransS[i][j][k][0][1]=0;
					}

					if (fphaseStat[Ix][Iy][Iz] == 1) Coef1 = 0;
					else Coef1 = 1;
					if (fphaseStat[Ix][Iy + 1][Iz] == 1) Coef2 = 0;
					else Coef2 = 1;

					mGam = (fblockFProps[Ix][Iy][Iz][BGAM][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BGAM][1] * gridDim[Nx + Iy + 1] * Coef2) / (gridDim[Nx + Iy] * Coef1 + gridDim[Nx + Iy + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = fP[i][j][k][2] - fP[i][j + 1][k][2] + mGam * (blockH[i][j][k] - blockH[i][j + 1][k]);
					if (pot > 0) {
						ftrans[i][j][k][1][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*fsat[Ix][Iy][Iz][2];
						ftransS[i][j][k][1][1] = 0;
						ftrans[i][j][k][3][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*fsat[Ix][Iy][Iz][2];

						ftrans[i][j][k][5][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2);
						ftrans[i][j][k][7][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2);

					}
					else {
						ftrans[i][j][k][1][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*fsat[Ix][Iy + 1][Iz][2];
						ftransS[i][j][k][1][1] = 1;
						ftrans[i][j][k][3][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2)*fsat[Ix][Iy + 1][Iz][2];

						ftrans[i][j][k][5][1] = fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2);
						ftrans[i][j][k][7][1] = mGam * fIFTran[Ix][Iy][Iz][1] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][RO][1] * gridDim[Nx + Iy + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Iy] * Coef1 + fblockFProps[Ix][Iy + 1][Iz][BMU][1] * gridDim[Nx + Iy + 1] * Coef2);
					}

					if ((!Coef1) && (!Coef2)) {
						ftrans[i][j][k][1][1] = 0;
						ftrans[i][j][k][3][1] = 0;
						ftrans[i][j][k][5][1] = 0;
						ftrans[i][j][k][7][1] = 0;

						//ftransS[i][j][k][1][1]=0;
					}

				}


				//Z-Direction
				if ((Nz > 1) && (k < Nz)) {
					pot = fP[i][j][k][0] - fP[i][j][k + 1][0] + G_ACC * watRo*(blockH[i][j][k] - blockH[i][j][k + 1]);
					if (pot > 0) {
						fWtran[i][j][k][2] = fIFTran[Ix][Iy][Iz][2] * (watRo*WAT_M_RO / watMu)*fsat[Ix][Iy][Iz][0];
						fdWtran[i][j][k][2] = fIFTran[Ix][Iy][Iz][2] * (watRo*WAT_M_RO / watMu);
						ftransS[i][j][k][2][2] = 0;
					}
					else {
						fWtran[i][j][k][2] = fIFTran[Ix][Iy][Iz][2] * (watRo*WAT_M_RO / watMu)*fsat[Ix][Iy][Iz + 1][0];
						fdWtran[i][j][k][2] = fIFTran[Ix][Iy][Iz][2] * (watRo*WAT_M_RO / watMu);
						ftransS[i][j][k][2][2] = 1;
					}

					if (fphaseStat[Ix][Iy][Iz] == (-1)) Coef1 = 0;
					else Coef1 = 1;
					if (fphaseStat[Ix][Iy][Iz + 1] == (-1)) Coef2 = 0;
					else Coef2 = 1;


					mGam = (fblockFProps[Ix][Iy][Iz][BGAM][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BGAM][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (gridDim[Nx + Ny + Iz] * Coef1 + gridDim[Nx + Ny + Iz + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = fP[i][j][k][1] - fP[i][j][k + 1][1] + mGam * (blockH[i][j][k] - blockH[i][j][k + 1]);
					if (pot > 0) {
						ftrans[i][j][k][0][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*fsat[Ix][Iy][Iz][1];
						ftransS[i][j][k][0][2] = 0;
						ftrans[i][j][k][2][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*fsat[Ix][Iy][Iz][1];

						ftrans[i][j][k][4][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*(-1);
						ftrans[i][j][k][6][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*(-1);
						ftrans[i][j][k][8][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*(-1);
						ftrans[i][j][k][9][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*(-1);

						ftrans[i][j][k][10][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2);
						ftrans[i][j][k][11][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2);
					}
					else {
						ftrans[i][j][k][0][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*fsat[Ix][Iy][Iz + 1][1];
						ftransS[i][j][k][0][2] = 1;
						ftrans[i][j][k][2][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*fsat[Ix][Iy][Iz + 1][1];

						ftrans[i][j][k][4][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*(-1);
						ftrans[i][j][k][6][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*(-1);
						ftrans[i][j][k][8][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*(-1);
						ftrans[i][j][k][9][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2)*(-1);

						ftrans[i][j][k][10][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2);
						ftrans[i][j][k][11][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][0] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][0] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][0] * gridDim[Nx + Ny + Iz + 1] * Coef2);
					}
					if ((!Coef1) && (!Coef2)) {
						ftrans[i][j][k][0][2] = 0;
						ftrans[i][j][k][2][2] = 0;
						ftrans[i][j][k][4][2] = 0;
						ftrans[i][j][k][6][2] = 0;
						ftrans[i][j][k][8][2] = 0;
						ftrans[i][j][k][9][2] = 0;

						//ftransS[i][j][k][0][2]=0;
					}

					if (fphaseStat[Ix][Iy][Iz] == 1) Coef1 = 0;
					else Coef1 = 1;
					if (fphaseStat[Ix][Iy][Iz + 1] == 1) Coef2 = 0;
					else Coef2 = 1;

					mGam = (fblockFProps[Ix][Iy][Iz][BGAM][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BGAM][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (gridDim[Nx + Ny + Iz] * Coef1 + gridDim[Nx + Ny + Iz + 1] * Coef2);
					if ((!Coef1) && (!Coef2)) mGam = 0;
					pot = fP[i][j][k][2] - fP[i][j][k + 1][2] + mGam * (blockH[i][j][k] - blockH[i][j][k + 1]);
					if (pot > 0) {
						ftrans[i][j][k][1][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*fsat[Ix][Iy][Iz][2];
						ftransS[i][j][k][1][2] = 0;
						ftrans[i][j][k][3][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*fsat[Ix][Iy][Iz][2];

						ftrans[i][j][k][5][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2);
						ftrans[i][j][k][7][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2);

					}
					else {
						ftrans[i][j][k][1][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*fsat[Ix][Iy][Iz + 1][2];
						ftransS[i][j][k][1][2] = 1;
						ftrans[i][j][k][3][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2)*fsat[Ix][Iy][Iz + 1][2];

						ftrans[i][j][k][5][2] = fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2);
						ftrans[i][j][k][7][2] = mGam * fIFTran[Ix][Iy][Iz][2] * (fblockFProps[Ix][Iy][Iz][RO][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][RO][1] * gridDim[Nx + Ny + Iz + 1] * Coef2) / (fblockFProps[Ix][Iy][Iz][BMU][1] * gridDim[Nx + Ny + Iz] * Coef1 + fblockFProps[Ix][Iy][Iz + 1][BMU][1] * gridDim[Nx + Ny + Iz + 1] * Coef2);
					}
					if ((!Coef1) && (!Coef2)) {
						ftrans[i][j][k][1][2] = 0;
						ftrans[i][j][k][3][2] = 0;
						ftrans[i][j][k][5][2] = 0;
						ftrans[i][j][k][7][2] = 0;

						//ftransS[i][j][k][1][2]=0;
					}

				}
			}
		}
	}
}
void Calc_Phase_P(void) {
	register int i, j, k;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++) {
				P[i + 1][j + 1][k + 1][0] = P[i + 1][j + 1][k + 1][1] - blockFProps[i][j][k][BLOCK_PC][0];
				if (Bift[i][j][k] < MISCIBLEIFT) P[i + 1][j + 1][k + 1][2] = P[i + 1][j + 1][k + 1][1] + blockFProps[i][j][k][BLOCK_PC][1] * Bift[i][j][k] / IFTBASECASE;
				else P[i + 1][j + 1][k + 1][2] = P[i + 1][j + 1][k + 1][1] + blockFProps[i][j][k][BLOCK_PC][1];
			}
}

void fCalc_Phase_P(void) {
	register int i, j, k;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++) {
				//fP[i + 1][j + 1][k + 1][0] = fP[i + 1][j + 1][k + 1][1] - fblockFProps[i][j][k][BLOCK_PC][0];
				fP[i + 1][j + 1][k + 1][0] = fP[i + 1][j + 1][k + 1][1];
				//if (fBift[i][j][k] < MISCIBLEIFT) fP[i + 1][j + 1][k + 1][2] = fP[i + 1][j + 1][k + 1][1] + fblockFProps[i][j][k][BLOCK_PC][1];
				//else fP[i + 1][j + 1][k + 1][2] = fP[i + 1][j + 1][k + 1][1] ;
				//fP[i + 1][j + 1][k + 1][2] = fP[i + 1][j + 1][k + 1][1] + fblockFProps[i][j][k][BLOCK_PC][1];
				fP[i + 1][j + 1][k + 1][2] = fP[i + 1][j + 1][k + 1][1];
			}
}



void System_Reduce(void) {
	/*std::ofstream testout;
	int Ix, Iy, Iz, i, nnn, j;
	testout.open("testout.txt", std::ios::out);
	i = 0;
	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				//i = Iz*(Nx*Ny) + Iy*Nx + Ix;
				testout << pJHolder[Ix][Iy][Iz] << std::endl;
				if (phaseStat[Ix][Iy][Iz]) nnn = Nc + 3;
				else nnn = 2 * Nc + 4;

				for (j = 0; j < nnn; j++) {
					testout << "           " << CSRrow[i] << std::endl;
					i++;
				}



			}

	testout << "frac\n";
	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				//i = Iz*(Nx*Ny) + Iy*Nx + Ix;
				testout << fpJHolder[Ix][Iy][Iz] << std::endl;
				if (fphaseStat[Ix][Iy][Iz]) nnn = Nc + 3;
				else nnn = 2 * Nc + 4;

				for (j = 0; j < nnn; j++) {
					testout << "           " << CSRrow[i] << std::endl;
					i++;
				}



			}

	testout.close();*/
	CSRrow[0] = totalRow;
	CSRrow[totalRow] = totalJac;


	/*preConRow[0]=0;
	preConIndex[0]=0;

	for (Iz=0; Iz<Nz; Iz++)
	for (Iy=0; Iy<Ny; Iy++)
	for (Ix=0; Ix<Nx; Ix++) {
	BLNo=(Iz*(Nx*Ny)+Iy*Nx+Ix)*(2*Nc+4);
	//printf("%d\n", BLNo);

	if (phaseStat[Ix][Iy][Iz]==1) {
	for (i=(BLNo-reRow); i<(totalSize-reRow-Nc); i++) {
	ans[i]=ans[i+Nc];
	}
	DelJacRow(BLNo-reRow, Nc);

	reRow+=Nc;



	for (i=(BLNo-reRow+2*Nc+2); i<(totalSize-reRow-1); i++) {
	ans[i]=ans[i+1];
	}
	DelJacRow(BLNo-reRow+2*Nc+2, 1);
	reRow++;

	DelJacCol(BLNo-reCol+Nc, Nc);
	reCol+=Nc;

	DelJacCol(BLNo-reCol+2*Nc+3, 1);
	reCol++;

	bSize=Nc+3;
	}

	else if (phaseStat[Ix][Iy][Iz]==-1) {
	for (i=(BLNo-reRow); i<(totalSize-reRow-Nc); i++) {
	ans[i]=ans[i+Nc];
	}
	DelJacRow(BLNo-reRow, Nc);
	reRow+=Nc;

	for (i=(BLNo-reRow+2*Nc+1); i<(totalSize-reRow-1); i++) {
	ans[i]=ans[i+1];
	}
	DelJacRow(BLNo-reRow+2*Nc+1, 1);
	reRow++;

	DelJacCol(BLNo-reCol, Nc);
	reCol+=Nc;
	DelJacCol(BLNo-reCol+2*Nc+2, 1);
	reCol++;

	bSize=Nc+3;
	}

	else bSize=2*Nc+4;

	preConRow[Iz*(Nx*Ny)+Iy*Nx+Ix+1]=preConRow[Iz*(Nx*Ny)+Iy*Nx+Ix]+bSize;
	preConIndex[Iz*(Nx*Ny)+Iy*Nx+Ix+1]=preConIndex[Iz*(Nx*Ny)+Iy*Nx+Ix]+bSize*bSize;
	}
	curSize=totalSize-reCol;	//Or totalSize-reRow

	*/
	//curSize=totalRow;
}

void PhaseCtrl(void) {
	int Ix, Iy, Iz, i;
	FType normSum;
	FType Sor;

	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {

				/*if (phaseStat[Ix][Iy][Iz] == (-1)) {
					for (i = 0; i < Nc; i++) comp[Ix][Iy][Iz][i][2] = comp[Ix][Iy][Iz][i][1];
					SuccFalsh(Ix, Iy, Iz);
				}
				else if (phaseStat[Ix][Iy][Iz] == 1) {
					for (i = 0; i < Nc; i++) comp[Ix][Iy][Iz][i][2] = comp[Ix][Iy][Iz][i][0];
					SuccFalsh(Ix, Iy, Iz);
				}*/

				//if (phaseStat[Ix][Iy][Iz]) critTest(Ix, Iy, Iz);

				//if ((phaseStat[Ix][Iy][Iz]==-1) && (P[Ix+1][Iy+1][Iz+1][2]<=Dew_Succ(Ix, Iy, Iz)) && (P[Ix+1][Iy+1][Iz+1][2]>=LowerDew_Succ(Ix, Iy, Iz))) { //2
				if ((phaseStat[Ix][Iy][Iz] == -1) && (P[Ix + 1][Iy + 1][Iz + 1][2] <= PAHSE_ENV_BOUNDARY*Dew_Succ(Ix, Iy, Iz))) { //2
					sat[Ix][Iy][Iz][1] = EPS_SAT;
					sat[Ix][Iy][Iz][2] -= EPS_SAT;
					//if (phaseStat[Ix][Iy][Iz]) buildPreconFlag = -1;
					phaseStat[Ix][Iy][Iz] = 0;
					for (i = 0; i < Nc; i++) comp[Ix][Iy][Iz][i][0] = Xm[i];

				}
				/*if ((phaseStat[Ix][Iy][Iz]==-1) && (P[Ix+1][Iy+1][Iz+1][2]>=LowerDew_Succ(Ix, Iy, Iz))) { //2
				sat[Ix][Iy][Iz][1]=EPS_SAT;
				sat[Ix][Iy][Iz][2]-=EPS_SAT;
				phaseStat[Ix][Iy][Iz]=0;
				for (i=0; i<Nc; i++) comp[Ix][Iy][Iz][i][0]=Xm[i];

				}*/
				else if ((phaseStat[Ix][Iy][Iz] == 1) && (P[Ix + 1][Iy + 1][Iz + 1][1] <= PAHSE_ENV_BOUNDARY*Bubl_Succ(Ix, Iy, Iz))) {
					sat[Ix][Iy][Iz][2] = EPS_SAT;
					sat[Ix][Iy][Iz][1] -= EPS_SAT;
					//if (phaseStat[Ix][Iy][Iz]) buildPreconFlag = -1;
					phaseStat[Ix][Iy][Iz] = 0;
					for (i = 0; i < Nc; i++) comp[Ix][Iy][Iz][i][1] = Xm[i];
				}
				else {
					if (sat[Ix][Iy][Iz][1] <= EPS_SAT*EPS_SAT_FRAC) {
						sat[Ix][Iy][Iz][1] = 0;
						sat[Ix][Iy][Iz][2] = 1 - sat[Ix][Iy][Iz][0];
						//if ((!phaseStat[Ix][Iy][Iz]) && ((P[Ix+1][Iy+1][Iz+1][2]>Dew_Succ(Ix, Iy, Iz))) || (P[Ix+1][Iy+1][Iz+1][2]<LowerDew_Succ(Ix, Iy, Iz))){		//2				
						//if ((phaseStat[Ix][Iy][Iz] == 0) || (phaseStat[Ix][Iy][Iz] == 1)) buildPreconFlag = -1;
						phaseStat[Ix][Iy][Iz] = -1;
						//}
					}
					if (sat[Ix][Iy][Iz][2] <= EPS_SAT * EPS_SAT_FRAC) {
						sat[Ix][Iy][Iz][2] = 0;
						sat[Ix][Iy][Iz][1] = 1 - sat[Ix][Iy][Iz][0];
						//if ((phaseStat[Ix][Iy][Iz] == 0) || (phaseStat[Ix][Iy][Iz] == (-1))) buildPreconFlag = -1;
						//if ((!phaseStat[Ix][Iy][Iz]) && (P[Ix+1][Iy+1][Iz+1][1]>Bubl_Succ(Ix, Iy, Iz))) {
						phaseStat[Ix][Iy][Iz] = 1;
						//}
					}
				}
				/*else if ((!phaseStat[Ix][Iy][Iz]) && (P[Ix+1][Iy+1][Iz+1][2]>Dew_Succ(Ix, Iy, Iz))) {
				phaseStat[Ix][Iy][Iz]=-1;
				}
				else if ((!phaseStat[Ix][Iy][Iz]) && (P[Ix+1][Iy+1][Iz+1][1]>Bubl_Succ(Ix, Iy, Iz))) {
				phaseStat[Ix][Iy][Iz]=1;
				}*/
				Sor = sat[Ix][Iy][Iz][1];
				//if ((gasInjStat) && (sat[Ix][Iy][Iz][1]<0.15) && (phaseStat[Ix][Iy][Iz]!=-1)) sat[Ix][Iy][Iz][1]=0.15
				normSum = 1.0 / (sat[Ix][Iy][Iz][0] + sat[Ix][Iy][Iz][1] + sat[Ix][Iy][Iz][2]);
				sat[Ix][Iy][Iz][0] *= normSum;
				sat[Ix][Iy][Iz][1] *= normSum;
				sat[Ix][Iy][Iz][2] *= normSum;
				//if ((gasInjStat) && (sat[Ix][Iy][Iz][1]<0.15) && (phaseStat[Ix][Iy][Iz]!=-1) && (Sor>=0.15)) sat[Ix][Iy][Iz][1]=0.15;
				//if ((gasInjStat) && (sat[Ix][Iy][Iz][1]<0.15) && (phaseStat[Ix][Iy][Iz]!=-1) && (Sor<0.15)) sat[Ix][Iy][Iz][1]=Sor;



				

			}



}

void fPhaseCtrl(void) {
	int Ix, Iy, Iz, i;
	FType normSum;
	FType Sor;

	for (Ix = 0; Ix < Nx; Ix++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Iz = 0; Iz < Nz; Iz++) {

				//if (fphaseStat[Ix][Iy][Iz]) critTest(Ix, Iy, Iz);

				//if ((fphaseStat[Ix][Iy][Iz]==-1) && (fP[Ix+1][Iy+1][Iz+1][2]<=Dew_Succ(Ix, Iy, Iz)) && (fP[Ix+1][Iy+1][Iz+1][2]>=LowerDew_Succ(Ix, Iy, Iz))) { //2
				if ((fphaseStat[Ix][Iy][Iz] == -1) && (fP[Ix + 1][Iy + 1][Iz + 1][2] <= PAHSE_ENV_BOUNDARY*fDew_Succ(Ix, Iy, Iz))) { //2
					fsat[Ix][Iy][Iz][1] = EPS_SAT;
					fsat[Ix][Iy][Iz][2] -= EPS_SAT;
					//if (fphaseStat[Ix][Iy][Iz]) buildPreconFlag = -1;
					fphaseStat[Ix][Iy][Iz] = 0;
					for (i = 0; i < Nc; i++) fcomp[Ix][Iy][Iz][i][0] = Xm[i];

				}
				/*if ((fphaseStat[Ix][Iy][Iz]==-1) && (fP[Ix+1][Iy+1][Iz+1][2]>=LowerDew_Succ(Ix, Iy, Iz))) { //2
				fsat[Ix][Iy][Iz][1]=EPS_SAT;
				fsat[Ix][Iy][Iz][2]-=EPS_SAT;
				fphaseStat[Ix][Iy][Iz]=0;
				for (i=0; i<Nc; i++) fcomp[Ix][Iy][Iz][i][0]=Xm[i];

				}*/
				else if ((fphaseStat[Ix][Iy][Iz] == 1) && (fP[Ix + 1][Iy + 1][Iz + 1][1] <= PAHSE_ENV_BOUNDARY*fBubl_Succ(Ix, Iy, Iz))) {
					fsat[Ix][Iy][Iz][2] = EPS_SAT;
					fsat[Ix][Iy][Iz][1] -= EPS_SAT;
					//if (fphaseStat[Ix][Iy][Iz]) buildPreconFlag = -1;
					fphaseStat[Ix][Iy][Iz] = 0;
					for (i = 0; i < Nc; i++) fcomp[Ix][Iy][Iz][i][1] = Xm[i];
				}
				else {
					if (fsat[Ix][Iy][Iz][1] <= EPS_SAT * EPS_SAT_FRAC) {
						fsat[Ix][Iy][Iz][1] = 0;
						fsat[Ix][Iy][Iz][2] = 1 - fsat[Ix][Iy][Iz][0];
						//if ((!fphaseStat[Ix][Iy][Iz]) && ((fP[Ix+1][Iy+1][Iz+1][2]>Dew_Succ(Ix, Iy, Iz))) || (fP[Ix+1][Iy+1][Iz+1][2]<LowerDew_Succ(Ix, Iy, Iz))){		//2				
						//if ((fphaseStat[Ix][Iy][Iz] == 0) || (fphaseStat[Ix][Iy][Iz] == 1)) buildPreconFlag = -1;
						fphaseStat[Ix][Iy][Iz] = -1;
						//}
					}
					if (fsat[Ix][Iy][Iz][2] <= EPS_SAT * EPS_SAT_FRAC) {
						fsat[Ix][Iy][Iz][2] = 0;
						fsat[Ix][Iy][Iz][1] = 1 - fsat[Ix][Iy][Iz][0];
						//if ((fphaseStat[Ix][Iy][Iz] == 0) || (fphaseStat[Ix][Iy][Iz] == (-1))) buildPreconFlag = -1;
						//if ((!fphaseStat[Ix][Iy][Iz]) && (fP[Ix+1][Iy+1][Iz+1][1]>Bubl_Succ(Ix, Iy, Iz))) {
						fphaseStat[Ix][Iy][Iz] = 1;
						//}
					}
				}
				/*else if ((!fphaseStat[Ix][Iy][Iz]) && (fP[Ix+1][Iy+1][Iz+1][2]>Dew_Succ(Ix, Iy, Iz))) {
				fphaseStat[Ix][Iy][Iz]=-1;
				}
				else if ((!fphaseStat[Ix][Iy][Iz]) && (fP[Ix+1][Iy+1][Iz+1][1]>Bubl_Succ(Ix, Iy, Iz))) {
				fphaseStat[Ix][Iy][Iz]=1;
				}*/
				Sor = fsat[Ix][Iy][Iz][1];
				//if ((gasInjStat) && (fsat[Ix][Iy][Iz][1]<0.15) && (fphaseStat[Ix][Iy][Iz]!=-1)) fsat[Ix][Iy][Iz][1]=0.15
				normSum = 1.0 / (fsat[Ix][Iy][Iz][0] + fsat[Ix][Iy][Iz][1] + fsat[Ix][Iy][Iz][2]);
				fsat[Ix][Iy][Iz][0] *= normSum;
				fsat[Ix][Iy][Iz][1] *= normSum;
				fsat[Ix][Iy][Iz][2] *= normSum;
				//if ((gasInjStat) && (fsat[Ix][Iy][Iz][1]<0.15) && (fphaseStat[Ix][Iy][Iz]!=-1) && (Sor>=0.15)) fsat[Ix][Iy][Iz][1]=0.15;
				//if ((gasInjStat) && (fsat[Ix][Iy][Iz][1]<0.15) && (fphaseStat[Ix][Iy][Iz]!=-1) && (Sor<0.15)) fsat[Ix][Iy][Iz][1]=Sor;

				/*if (fphaseStat[Ix][Iy][Iz] == (-1)) {
					for (i = 0; i < Nc; i++) fcomp[Ix][Iy][Iz][i][2] = fcomp[Ix][Iy][Iz][i][1];
					fSuccFalsh(Ix, Iy, Iz);
				}
				else if (fphaseStat[Ix][Iy][Iz] == 1) {
					for (i = 0; i < Nc; i++) fcomp[Ix][Iy][Iz][i][2] = fcomp[Ix][Iy][Iz][i][0];
					fSuccFalsh(Ix, Iy, Iz);
				}*/
			}

}

FType System_Update(void) {
	int Ix, Iy, Iz, i, BLNo, MSize;
	//FType NRTol, NRTol2;
	FType NRTol;
	FType normSuml, normSumg;
	FType *XNew;

	//clock_t start1;
	//FType end1;	

	//Scaling(jac, ans, curSize);
	//GaussJ(jac, ans, ans, curSize);
	//biCGStab(jac, ans, ans, curSize);
	MSize = CSRrow[0];
	MKLS1(CSRjac, CSRrow, CSRcol, ans);
	//Orthomin(jac, ans, ans, curSize);
	BICGStart = clock();
	//CLBiCGS(CSRjac, CSRrow, CSRcol, ans, preConRow, preConIndex, preConCSRrow, preConCSRcol);


	if ((XNew = (FType *)malloc(MSize * sizeof(FType))) == NULL) TerM("Can not allocate memory for XNEW");

	for (i = 0; i < wellNO; i++) WellCondition[i] = TempWellCondition[i];
	//Scaling1(CSRjac, CSRrow, CSRcol, ans);
	//BuildPrecon();
	//MRPre();
	//PrebiCGStab(CSRjac, CSRrow, CSRcol, ans);
	//Scaling(CSRjac, CSRrow, CSRcol, ans);
	//ParaSolver(CSRjac, CSRrow, CSRcol, ans, XNew);
	//MKLS1(CSRjac, CSRrow, CSRcol, ans);
	//Scaling1(CSRjac, CSRrow, CSRcol, ans);
	//if (buildPreconFlag) {
	//	BuildPrecon();
	//	MRPre();
	//	printf("\n****************With Preconditioner***********************\n");
	//}
	//PrebiCGStab(CSRjac, CSRrow, CSRcol, ans);	
	BICGEnd = clock();
	solverTime = BICGEnd - BICGStart;
	solverTimeAcc += solverTime;
	//GMRESm(CSRjac, CSRrow, CSRcol, ans, 20);
	//PreGMRESm(CSRjac, CSRrow, CSRcol, ans, 20);

	//NRTol=In_Product(ans, ans, curSize);
	//end1=(FType) (clock()-start1);
	//if (repStat) buildPreconFlag = -1;
	//else buildPreconFlag = 0;

	NRTol = 0;
	for (i = 0; i < totalRow; i++) {
		XNew[i] = ans[i];
		NRTol += XNew[i] * XNew[i];
		//if (ans[i]!=ans[i]) repStat=-1;
	}
	NRTol /= totalRow;

	//NRTol=0;
	BLNo = 0;
	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				if (phaseStat[Ix][Iy][Iz] == 1) {
					normSuml = 0;
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][0] += XNew[BLNo];
						if (comp[Ix][Iy][Iz][i][0] < EPS_COMP) comp[Ix][Iy][Iz][i][0] = 0;
						else if (comp[Ix][Iy][Iz][i][0] > 1) comp[Ix][Iy][Iz][i][0] = 1;

						normSuml += comp[Ix][Iy][Iz][i][0];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][0] /= normSuml;
					}

					if ((fabs(PCoeff*XNew[BLNo])) > (P[Ix + 1][Iy + 1][Iz + 1][1] * 0.2)) P[Ix + 1][Iy + 1][Iz + 1][1] *= 0.2*dSgn(XNew[BLNo]);
					else P[Ix + 1][Iy + 1][Iz + 1][1] += XNew[BLNo] * PCoeff;
					BLNo++;

					sat[Ix][Iy][Iz][0] += XNew[BLNo];
					BLNo++;

					sat[Ix][Iy][Iz][1] += XNew[BLNo];
					BLNo++;

					if (sat[Ix][Iy][Iz][0] < 0) sat[Ix][Iy][Iz][0] = 0;
					if (sat[Ix][Iy][Iz][0] > 1) sat[Ix][Iy][Iz][0] = 1;
					if (sat[Ix][Iy][Iz][1] > 1) sat[Ix][Iy][Iz][1] = 1;
					if (sat[Ix][Iy][Iz][1] < 0) sat[Ix][Iy][Iz][1] = 0;

					normSuml = sat[Ix][Iy][Iz][0] + sat[Ix][Iy][Iz][1];
					sat[Ix][Iy][Iz][0] /= normSuml;
					sat[Ix][Iy][Iz][1] /= normSuml;
				}

				else if (phaseStat[Ix][Iy][Iz] == -1) {
					normSumg = 0;
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][1] += XNew[BLNo];
						if (comp[Ix][Iy][Iz][i][1] < EPS_COMP) comp[Ix][Iy][Iz][i][1] = 0;
						else if (comp[Ix][Iy][Iz][i][1] > 1) comp[Ix][Iy][Iz][i][1] = 1;
						normSumg += comp[Ix][Iy][Iz][i][1];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][1] /= normSumg;
					}

					if ((fabs(PCoeff*XNew[BLNo])) > (P[Ix + 1][Iy + 1][Iz + 1][1] * 0.2)) P[Ix + 1][Iy + 1][Iz + 1][1] *= 0.2*dSgn(XNew[BLNo]);
					else P[Ix + 1][Iy + 1][Iz + 1][1] += XNew[BLNo] * PCoeff;
					BLNo++;

					sat[Ix][Iy][Iz][0] += XNew[BLNo];
					BLNo++;

					sat[Ix][Iy][Iz][2] += XNew[BLNo];
					BLNo++;

					if (sat[Ix][Iy][Iz][0] < 0) sat[Ix][Iy][Iz][0] = 0;
					if (sat[Ix][Iy][Iz][0] > 1) sat[Ix][Iy][Iz][0] = 1;
					if (sat[Ix][Iy][Iz][2] > 1) sat[Ix][Iy][Iz][2] = 1;
					if (sat[Ix][Iy][Iz][2] < 0) sat[Ix][Iy][Iz][2] = 0;


					normSuml = sat[Ix][Iy][Iz][0] + sat[Ix][Iy][Iz][2];
					sat[Ix][Iy][Iz][0] /= normSuml;
					sat[Ix][Iy][Iz][2] /= normSuml;
				}
				else if (phaseStat[Ix][Iy][Iz] == 0) {
					normSuml = 0;
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][0] += XNew[BLNo];
						if (comp[Ix][Iy][Iz][i][0] < EPS_COMP) comp[Ix][Iy][Iz][i][0] = 0;
						else if (comp[Ix][Iy][Iz][i][0] > 1) comp[Ix][Iy][Iz][i][0] = 1;
						normSuml += comp[Ix][Iy][Iz][i][0];
						BLNo++;
					}
					normSumg = 0;
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][1] += XNew[BLNo];
						if (comp[Ix][Iy][Iz][i][1] < EPS_COMP) comp[Ix][Iy][Iz][i][1] = 0;
						else if (comp[Ix][Iy][Iz][i][1] > 1) comp[Ix][Iy][Iz][i][1] = 1;
						normSumg += comp[Ix][Iy][Iz][i][1];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][0] /= normSuml;
						comp[Ix][Iy][Iz][i][1] /= normSumg;
					}

					if ((fabs(PCoeff*XNew[BLNo])) > (P[Ix + 1][Iy + 1][Iz + 1][1] * 0.2)) P[Ix + 1][Iy + 1][Iz + 1][1] *= 0.2*dSgn(XNew[BLNo]);
					else P[Ix + 1][Iy + 1][Iz + 1][1] += XNew[BLNo] * PCoeff;
					BLNo++;

					sat[Ix][Iy][Iz][0] += XNew[BLNo];
					BLNo++;

					sat[Ix][Iy][Iz][1] += XNew[BLNo];
					BLNo++;

					sat[Ix][Iy][Iz][2] += XNew[BLNo];
					BLNo++;

					if (sat[Ix][Iy][Iz][0] < 0) sat[Ix][Iy][Iz][0] = 0;
					if (sat[Ix][Iy][Iz][0] > 1) sat[Ix][Iy][Iz][0] = 1;
					if (sat[Ix][Iy][Iz][1] > 1) sat[Ix][Iy][Iz][1] = 1;
					if (sat[Ix][Iy][Iz][2] > 1) sat[Ix][Iy][Iz][2] = 1;
					if (sat[Ix][Iy][Iz][1] < 0) sat[Ix][Iy][Iz][1] = 0;
					if (sat[Ix][Iy][Iz][2] < 0) sat[Ix][Iy][Iz][2] = 0;

				}

				normSuml = sat[Ix][Iy][Iz][0] + sat[Ix][Iy][Iz][1] + sat[Ix][Iy][Iz][2];
				sat[Ix][Iy][Iz][0] /= normSuml;
				sat[Ix][Iy][Iz][1] /= normSuml;
				sat[Ix][Iy][Iz][2] /= normSuml;
			}


	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				if (fphaseStat[Ix][Iy][Iz] == 1) {
					normSuml = 0;
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][0] += XNew[BLNo];
						if (fcomp[Ix][Iy][Iz][i][0] < EPS_COMP) fcomp[Ix][Iy][Iz][i][0] = 0;
						else if (fcomp[Ix][Iy][Iz][i][0] > 1) fcomp[Ix][Iy][Iz][i][0] = 1;

						normSuml += fcomp[Ix][Iy][Iz][i][0];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][0] /= normSuml;
					}


					if ((fabs(PCoeff*XNew[BLNo])) > (fP[Ix + 1][Iy + 1][Iz + 1][1] * 0.2)) fP[Ix + 1][Iy + 1][Iz + 1][1] *= 0.2*dSgn(XNew[BLNo]);
					else fP[Ix + 1][Iy + 1][Iz + 1][1] += XNew[BLNo] * PCoeff;
					BLNo++;

					fsat[Ix][Iy][Iz][0] += XNew[BLNo];
					BLNo++;

					fsat[Ix][Iy][Iz][1] += XNew[BLNo];
					BLNo++;

					if (fsat[Ix][Iy][Iz][0] < 0) fsat[Ix][Iy][Iz][0] = 0;
					if (fsat[Ix][Iy][Iz][0] > 1) fsat[Ix][Iy][Iz][0] = 1;
					if (fsat[Ix][Iy][Iz][1] > 1) fsat[Ix][Iy][Iz][1] = 1;
					if (fsat[Ix][Iy][Iz][1] < 0) fsat[Ix][Iy][Iz][1] = 0;

					normSuml = fsat[Ix][Iy][Iz][0] + fsat[Ix][Iy][Iz][1];
					fsat[Ix][Iy][Iz][0] /= normSuml;
					fsat[Ix][Iy][Iz][1] /= normSuml;
				}

				else if (fphaseStat[Ix][Iy][Iz] == -1) {
					normSumg = 0;
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][1] += XNew[BLNo];
						if (fcomp[Ix][Iy][Iz][i][1] < EPS_COMP) fcomp[Ix][Iy][Iz][i][1] = 0;
						else if (fcomp[Ix][Iy][Iz][i][1] > 1) fcomp[Ix][Iy][Iz][i][1] = 1;
						normSumg += fcomp[Ix][Iy][Iz][i][1];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][1] /= normSumg;
					}
					if ((fabs(PCoeff*XNew[BLNo])) > (fP[Ix + 1][Iy + 1][Iz + 1][1] * 0.2)) fP[Ix + 1][Iy + 1][Iz + 1][1] *= 0.2*dSgn(XNew[BLNo]);
					else fP[Ix + 1][Iy + 1][Iz + 1][1] += XNew[BLNo] * PCoeff;
					BLNo++;

					fsat[Ix][Iy][Iz][0] += XNew[BLNo];
					BLNo++;

					fsat[Ix][Iy][Iz][2] += XNew[BLNo];
					BLNo++;

					if (fsat[Ix][Iy][Iz][0] < 0) fsat[Ix][Iy][Iz][0] = 0;
					if (fsat[Ix][Iy][Iz][0] > 1) fsat[Ix][Iy][Iz][0] = 1;
					if (fsat[Ix][Iy][Iz][2] > 1) fsat[Ix][Iy][Iz][2] = 1;
					if (fsat[Ix][Iy][Iz][2] < 0) fsat[Ix][Iy][Iz][2] = 0;

					normSuml = fsat[Ix][Iy][Iz][0] + fsat[Ix][Iy][Iz][2];
					fsat[Ix][Iy][Iz][0] /= normSuml;
					fsat[Ix][Iy][Iz][2] /= normSuml;
				}
				else if (fphaseStat[Ix][Iy][Iz] == 0) {
					normSuml = 0;
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][0] += XNew[BLNo];
						if (fcomp[Ix][Iy][Iz][i][0] < EPS_COMP) fcomp[Ix][Iy][Iz][i][0] = 0;
						else if (fcomp[Ix][Iy][Iz][i][0] > 1) fcomp[Ix][Iy][Iz][i][0] = 1;
						normSuml += fcomp[Ix][Iy][Iz][i][0];
						BLNo++;
					}
					normSumg = 0;
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][1] += XNew[BLNo];
						if (fcomp[Ix][Iy][Iz][i][1] < EPS_COMP) fcomp[Ix][Iy][Iz][i][1] = 0;
						else if (fcomp[Ix][Iy][Iz][i][1] > 1) fcomp[Ix][Iy][Iz][i][1] = 1;
						normSumg += fcomp[Ix][Iy][Iz][i][1];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][0] /= normSuml;
						fcomp[Ix][Iy][Iz][i][1] /= normSumg;
					}
					if ((fabs(PCoeff*XNew[BLNo])) > (fP[Ix + 1][Iy + 1][Iz + 1][1] * 0.2)) fP[Ix + 1][Iy + 1][Iz + 1][1] *= 0.2*dSgn(XNew[BLNo]);
					else fP[Ix + 1][Iy + 1][Iz + 1][1] += XNew[BLNo] * PCoeff;
					BLNo++;

					fsat[Ix][Iy][Iz][0] += XNew[BLNo];
					BLNo++;

					fsat[Ix][Iy][Iz][1] += XNew[BLNo];
					BLNo++;

					fsat[Ix][Iy][Iz][2] += XNew[BLNo];
					BLNo++;

					if (fsat[Ix][Iy][Iz][0] < 0) fsat[Ix][Iy][Iz][0] = 0;
					if (fsat[Ix][Iy][Iz][0] > 1) fsat[Ix][Iy][Iz][0] = 1;
					if (fsat[Ix][Iy][Iz][1] > 1) fsat[Ix][Iy][Iz][1] = 1;
					if (fsat[Ix][Iy][Iz][2] > 1) fsat[Ix][Iy][Iz][2] = 1;
					if (fsat[Ix][Iy][Iz][1] < 0) fsat[Ix][Iy][Iz][1] = 0;
					if (fsat[Ix][Iy][Iz][2] < 0) fsat[Ix][Iy][Iz][2] = 0;

					normSuml = fsat[Ix][Iy][Iz][0] + fsat[Ix][Iy][Iz][1] + fsat[Ix][Iy][Iz][2];
					fsat[Ix][Iy][Iz][0] /= normSuml;
					fsat[Ix][Iy][Iz][1] /= normSuml;
					fsat[Ix][Iy][Iz][2] /= normSuml;
				}
			}
	free(XNew);
	return NRTol;
}

void CalcPhaseGrav(void) {
	register int i, j, k, n;
	FType compMWL, compMWG;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++) {
				compMWL = 0;
				compMWG = 0;
				for (n = 0; n < Nc; n++) {
					compMWL += comp[i][j][k][n][0] * fluidProp[n][MW];
					compMWG += comp[i][j][k][n][1] * fluidProp[n][MW];
				}
				blockFProps[i][j][k][BGAM][0] = G_ACCG * blockFProps[i][j][k][RO][0] * compMWL;
				blockFProps[i][j][k][BGAM][1] = G_ACCG * blockFProps[i][j][k][RO][1] * compMWG;
			}
}

void fCalcPhaseGrav(void) {
	register int i, j, k, n;
	FType compMWL, compMWG;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++) {
				compMWL = 0;
				compMWG = 0;
				for (n = 0; n < Nc; n++) {
					compMWL += fcomp[i][j][k][n][0] * fluidProp[n][MW];
					compMWG += fcomp[i][j][k][n][1] * fluidProp[n][MW];
				}
				fblockFProps[i][j][k][BGAM][0] = G_ACCG * fblockFProps[i][j][k][RO][0] * compMWL;
				fblockFProps[i][j][k][BGAM][1] = G_ACCG * fblockFProps[i][j][k][RO][1] * compMWG;
			}
}

/*void DelJacRow(int startRow, int rowNO) {
register int i;
int shiftP, CSRrow0;

if (startRow) CSRrow0=CSRrow[startRow];
else CSRrow0=0;

shiftP=CSRrow[startRow+rowNO]-CSRrow0;

for (i=CSRrow0; i<(CSRrow[CSRrow[0]]-shiftP); i++) {
CSRjac[i]=CSRjac[i+shiftP];
CSRcol[i]=CSRcol[i+shiftP];
}

for (i=(startRow+1); i<=(CSRrow[0]-rowNO); i++) {
CSRrow[i]=CSRrow[i+rowNO]-shiftP;
}

CSRrow[0]-=rowNO;
}

void DelJacCol(int startCol, int colNO) {
register int i, j, offSet;
int colHolderSize;
int *colHolder;


if ((colHolder=(int *) malloc(colNO*(Nc+1)*sizeof(int)))==NULL) TerM("Can not allocate memory for colHolder matrix");

colHolderSize=0;
for (i=0; i<CSRrow[CSRrow[0]]; i++) {
if ((CSRcol[i]>=startCol) && (CSRcol[i]<(startCol+colNO))) {
colHolder[colHolderSize]=i;
colHolderSize++;
}
}

j=0;
for (i=1; i<=CSRrow[0]; i++) {
while (j<(colHolderSize) && (colHolder[j]<CSRrow[i])) j++;
CSRrow[i]-=j;
}

offSet=1;
for (i=colHolder[0]; i<CSRrow[CSRrow[0]]; i++) {
while ((i+offSet)==colHolder[offSet]) offSet++;

CSRjac[i]=CSRjac[i+offSet];
CSRcol[i]=CSRcol[i+offSet];
}

for (i=0; i<CSRrow[CSRrow[0]]; i++) {
if (CSRcol[i]>startCol) CSRcol[i]-=colNO;
}

free(colHolder);
}*/

/*void BuildPrecon(void) {
register int Ix, Iy, Iz, i, j, k, n;
int bLNO, temp, wn;
int bSize, preIndex, MSize;
FType **blockM, **blockP;
FType ra;


if ((blockM = (FType **)malloc((2 * Nc + 4) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for blockM matrix");
for (i = 0; i<(2 * Nc + 4); i++) {
if ((blockM[i] = (FType *)malloc((2 * Nc + 4) * sizeof(FType))) == NULL) TerM("Can not allocate memory for blockM matrix");
}
if ((blockP = (FType **)malloc((2 * Nc + 4) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for blockP matrix");
for (i = 0; i<(2 * Nc + 4); i++) {
if ((blockP[i] = (FType *)malloc((2 * Nc + 4) * sizeof(FType))) == NULL) TerM("Can not allocate memory for blockP matrix");
}


MSize = CSRrow[0];
CSRrow[0] = 0;
for (Iz = 0; Iz<Nz; Iz++)
for (Iy = 0; Iy<Ny; Iy++)
for (Ix = 0; Ix<Nx; Ix++) {
temp = Iz*(Nx*Ny) + Iy*Nx + Ix;
bLNO = preConRow[temp];
bSize = preConRow[temp + 1] - bLNO;
preIndex = preConIndex[temp];

if (phaseStat[Ix][Iy][Iz]) wn = Nc;
else wn = 2 * Nc;

for (i = 0; i<bSize; i++)
for (j = 0; j<bSize; j++) {
if (i == j) blockP[i][j] = 1;
else blockP[i][j] = 0;
blockM[i][j] = 0;
}

for (i = bLNO; i<(bLNO + bSize); i++) {
for (j = CSRrow[i]; j<CSRrow[i + 1]; j++) {
temp = CSRcol[j] - bLNO;
if ((temp<bSize) && (temp >= 0)) blockM[i - bLNO][temp] = CSRjac[j];
}
}


for (i = 0; i<bSize; i++) {
for (j = 0; j<bSize; j++) {
if (i == j) continue;

if (!blockM[i][i]) {
for (n = (i + 1); n<bSize; n++) {
if (blockM[n][i]) break;
}
if (n == bSize) continue;
for (k = 0; k<bSize; k++) {
blockM[i][k] += blockM[n][k];
blockP[i][k] += blockP[n][k];
}
}

ra = -blockM[j][i] / blockM[i][i];

for (k = 0; k<bSize; k++) {
blockM[j][k] += ra*blockM[i][k];
blockP[j][k] += ra*blockP[i][k];
}
}
}

k = 0;
for (i = 0; i<bSize; i++) {
ra = blockM[i][i];
preConCSRrow[bLNO + i] = preIndex + k;
for (j = 0; j<bSize; j++) {
preCon[preIndex + k] = blockP[i][j] / ra;
preConCSRcol[preIndex + k] = bLNO + j;
k++;
}
if (i == wn) {
if (bLNO<(MSize / 2)) {
for (n = (MSize - 1); n>(MSize - MRMAXNONZERO + bSize - 1); n--) {
preCon[preIndex + k] = 0;
preConCSRcol[preIndex + k] = n;
k++;
}
}
else {
for (n = 0; n<(MRMAXNONZERO - bSize); n++) {
preCon[preIndex + k] = 0;
preConCSRcol[preIndex + k] = n;
k++;
}
}
}
}

}

CSRrow[0] = MSize;

for (i = 0; i<(2 * Nc + 4); i++) {
free(blockM[i]);
}
free(blockM);

for (i = 0; i<(2 * Nc + 4); i++) {
free(blockP[i]);
}
free(blockP);

}*/

double oilRateCh(int Ix, int Iy, int Iz, double cGas) {
	register int i, j, k;
	FType L, F, dF, tempD;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG;
	FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType ESum;
	FType tempSQR;
	FType newD;
	FType Lb = 10;

	nesbat = cGas;

	for (i = 0; i < Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / 101325)*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / 288.705555555556));
		STcomp[i][2] = (comp[Ix][Iy][Iz][i][0] + cGas * comp[Ix][Iy][Iz][i][1]) / (1 + cGas);
		//STcomp[i][2]=comp[Ix][Iy][Iz][i][0];
	}


	k = 0;
	do {
		k++;
		L = 0.7;
		j = 0;
		do {
			F = 0;
			dF = 0;
			for (i = 0; i < Nc; i++) {
				if (Ki[i] >= 0) tempD = (1 - Ki[i]) / (L + (1 - L)*Ki[i]);
				else tempD = 1 / (L - 1);
				F += STcomp[i][2] * tempD;
				dF -= STcomp[i][2] * tempD*tempD;
			}
			tempD = F / dF;
			L -= tempD;
			j++;
		} while (((tempD*tempD) > NEWTON_TOL) && (j < MAX_FLASH_N_ITR));

		if (L > 1) {
			L = 1;
			//break;
		}
		else if (L < 0) {
			L = 0;
			//break;
		}

		ESum = fabs(Lb - L);
		Lb = L;

		for (i = 0; i < Nc; i++) {
			if (Ki[i] >= 0) {
				STcomp[i][0] = STcomp[i][2] / (L + Ki[i] * (1 - L));
				STcomp[i][1] = Ki[i] * STcomp[i][0];
			}
			else {
				STcomp[i][0] = 0;
				STcomp[i][1] = STcomp[i][2] / (1 - L);
			}
		}

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += STcomp[i][0] * fluidProp[i][EOS_B];
			Bg += STcomp[i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += STcomp[i][0] * STcomp[j][0] * tempSQR;
				Ag += STcomp[i][1] * STcomp[j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= 0.0175848235251236;
		Bl *= 42.2111625483491;
		Ag *= 0.0175848235251236;     //gas-oil
		Bg *= 42.2111625483491;


		if (SRK) {
			Ul = Bl;
			Wl = 0;
			Ug = Bg;
			Wg = 0;
			d1 = 1;
			d2 = 0;
		}
		else if (PR) {
			Ul = 2 * Bl;
			Wl = Bl;
			Ug = 2 * Bg;
			Wg = Bg;
			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl * Ul - Ul - Wl * Wl, -(Al*Bl - Bl * Wl*Wl - Wl * Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg * Ug - Ug - Wg * Wg, -(Ag*Bg - Bg * Wg*Wg - Wg * Wg), 'g');


		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += STcomp[j][0] * tempSQR;
				tg += STcomp[j][1] * tempSQR;
			}
			t1l = fluidProp[i][EOS_B] / Cbl;
			t1g = fluidProp[i][EOS_B] / Cbg;

			EEl = 2 * tl / Cal - t1l;
			EEg = 2 * tg / Cag - t1g;

			DDl = Al / (Bl*(d1 - d2));
			DDg = Ag / (Bg*(d1 - d2));

			FFl = log((Zl + d2 * Bl) / (Zl + d1 * Bl));
			FFg = log((Zg + d2 * Bg) / (Zg + d1 * Bg));

			phiL = exp(t1l*(Zl - 1) - log(Zl - Bl) + DDl * EEl*FFl);
			phiG = exp(t1g*(Zg - 1) - log(Zg - Bg) + DDg * EEg*FFg);

			fL = STcomp[i][0] * phiL;		//eliminated P
			fG = STcomp[i][1] * phiG;	//gas-oil

			if (STcomp[i][2]) {
				if (fG) {
					Ki[i] *= fL / fG;
				}
				else Ki[i] = -1;
				if (fL != fL) Ki[i] = -1;
				if (fG != fG) Ki[i] = 1;
			}
		}
	} while ((ESum > FLASH_TOL) && (k < 500));
	//////////////////////////////////////////////////////////////

	newD = 42.2111625483491 / (Zl*L);
	if (newD != newD) newD = 0;
	//if (!_finite(newD)) newD=0;


	wGLR = Zg * (1 - L) / (Zl*L);

	return newD;
}

void JacGuass(void) {
	register int Ix, Iy, Iz, i, j, k, n, xn;
	int temp, bLNO, offSet;
	int pivot, MSize;
	FType holder;

	MSize = CSRrow[0];
	CSRrow[0] = 0;
	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				temp = Iz * (Nx*Ny) + Iy * Nx + Ix;
				bLNO = preConRow[temp];

				if (phaseStat[Ix][Iy][Iz] == 1) {
					offSet = Nc - 1;
					n = 0;
					xn = 1;
				}
				else if (phaseStat[Ix][Iy][Iz] == (-1)) {
					offSet = Nc - 1;
					n = 1;
					xn = 1;
				}
				else {
					offSet = 2 * Nc - 1;
					n = 0;
					xn = 1;
				}

				pivot = CSRrow[bLNO + offSet];

				for (i = 1; i < Nc; i++) {
					k = 0;
					for (j = CSRrow[bLNO + offSet - i]; j < CSRrow[bLNO + offSet - i + 1]; j++) {
						holder = CSRjac[j];
						CSRjac[j] -= CSRjac[pivot + k] * comp[Ix][Iy][Iz][Nc - i - 1][n] * xn;
						if (machine_epsv(holder) >= fabs(CSRjac[j])) CSRjac[j] = 0;
						//if ((CSRjac[j]/holder)<1e-16) CSRjac[j]=0;
						k++;
					}
					ans[bLNO + offSet - i] -= ans[bLNO + offSet] * comp[Ix][Iy][Iz][Nc - i - 1][n] * xn;
				}
			}

	CSRrow[0] = MSize;
}

void CalcIFT(void) {
	register int Ix, Iy, Iz, n;
	FType pSum;

	for (Iz = 0; Iz < Nz; Iz++) {
		for (Iy = 0; Iy < Ny; Iy++) {
			for (Ix = 0; Ix < Nx; Ix++) {
				if (!phaseStat[Ix][Iy][Iz]) {
					pSum = 0;
					for (n = 0; n < Nc; n++) pSum += fluidProp[n][PARACHOR] * (comp[Ix][Iy][Iz][n][0] * blockFProps[Ix][Iy][Iz][RO][0] - comp[Ix][Iy][Iz][n][1] * blockFProps[Ix][Iy][Iz][RO][1]);
					Bift[Ix][Iy][Iz] = pow(pSum*1e-6, 3.88);
				}
				else Bift[Ix][Iy][Iz] = IFTBASECASE;
			}

		}
	}
}

void fCalcIFT(void) {
	register int Ix, Iy, Iz, n;
	FType pSum;

	for (Iz = 0; Iz < Nz; Iz++) {
		for (Iy = 0; Iy < Ny; Iy++) {
			for (Ix = 0; Ix < Nx; Ix++) {
				if (!fphaseStat[Ix][Iy][Iz]) {
					pSum = 0;
					for (n = 0; n < Nc; n++) pSum += fluidProp[n][PARACHOR] * (fcomp[Ix][Iy][Iz][n][0] * fblockFProps[Ix][Iy][Iz][RO][0] - fcomp[Ix][Iy][Iz][n][1] * fblockFProps[Ix][Iy][Iz][RO][1]);
					fBift[Ix][Iy][Iz] = pow(pSum*1e-6, 3.88);
				}
				else fBift[Ix][Iy][Iz] = IFTBASECASE;
			}

		}
	}
}

double oilRateCh2(int Ix, int Iy, int Iz, double flowO, double flowG) {
	register int i, j, k;
	FType L, F, dF, tempD;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG;
	FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType ESum;
	FType tempSQR;
	FType newD;
	FType Lb = 10;

	for (i = 0; i < Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / 101325)*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / 288.705555555556));
		STcomp[i][2] = (comp[Ix][Iy][Iz][i][0] * flowO + flowG * comp[Ix][Iy][Iz][i][1]) / (flowO + flowG);
		//STcomp[i][2]=comp[Ix][Iy][Iz][i][0];
	}


	k = 0;
	do {
		k++;
		L = 0.6;
		j = 0;
		do {
			F = 0;
			dF = 0;
			for (i = 0; i < Nc; i++) {
				if (Ki[i] >= 0) tempD = (1 - Ki[i]) / (L + (1 - L)*Ki[i]);
				else tempD = 1 / (L - 1);
				F += STcomp[i][2] * tempD;
				dF -= STcomp[i][2] * tempD*tempD;
			}
			tempD = F / dF;
			L -= tempD;
			j++;
		} while (((tempD*tempD) > NEWTON_TOL) && (j < MAX_FLASH_N_ITR));

		if (L > 1) {
			L = 1;
			//break;
		}
		else if (L < 0) {
			L = 0;
			//break;
		}

		ESum = fabs(Lb - L);
		Lb = L;

		for (i = 0; i < Nc; i++) {
			if (Ki[i] >= 0) {
				STcomp[i][0] = STcomp[i][2] / (L + Ki[i] * (1 - L));
				STcomp[i][1] = Ki[i] * STcomp[i][0];
			}
			else {
				STcomp[i][0] = 0;
				STcomp[i][1] = STcomp[i][2] / (1 - L);
			}
		}

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += STcomp[i][0] * fluidProp[i][EOS_B];
			Bg += STcomp[i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += STcomp[i][0] * STcomp[j][0] * tempSQR;
				Ag += STcomp[i][1] * STcomp[j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= 0.0175848235251236;
		Bl *= 42.2111625483491;
		Ag *= 0.0175848235251236;     //gas-oil
		Bg *= 42.2111625483491;


		if (SRK) {
			Ul = Bl;
			Wl = 0;
			Ug = Bg;
			Wg = 0;
			d1 = 1;
			d2 = 0;
		}
		else if (PR) {
			Ul = 2 * Bl;
			Wl = Bl;
			Ug = 2 * Bg;
			Wg = Bg;
			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl * Ul - Ul - Wl * Wl, -(Al*Bl - Bl * Wl*Wl - Wl * Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg * Ug - Ug - Wg * Wg, -(Ag*Bg - Bg * Wg*Wg - Wg * Wg), 'g');


		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += STcomp[j][0] * tempSQR;
				tg += STcomp[j][1] * tempSQR;
			}
			t1l = fluidProp[i][EOS_B] / Cbl;
			t1g = fluidProp[i][EOS_B] / Cbg;

			EEl = 2 * tl / Cal - t1l;
			EEg = 2 * tg / Cag - t1g;

			DDl = Al / (Bl*(d1 - d2));
			DDg = Ag / (Bg*(d1 - d2));

			FFl = log((Zl + d2 * Bl) / (Zl + d1 * Bl));
			FFg = log((Zg + d2 * Bg) / (Zg + d1 * Bg));

			phiL = exp(t1l*(Zl - 1) - log(Zl - Bl) + DDl * EEl*FFl);
			phiG = exp(t1g*(Zg - 1) - log(Zg - Bg) + DDg * EEg*FFg);

			fL = STcomp[i][0] * phiL;		//eliminated P
			fG = STcomp[i][1] * phiG;	//gas-oil

			if (STcomp[i][2]) {
				if (fG) {
					Ki[i] *= fL / fG;
				}
				else Ki[i] = -1;
				if (fL != fL) Ki[i] = -1;
				if (fG != fG) Ki[i] = 1;
			}
		}
	} while (ESum > FLASH_TOL);
	//////////////////////////////////////////////////////////////

	newD = 0.023690416199152*(flowO + flowG)*Zl*L;
	if (newD != newD) newD = 0;
	//if (!_finite(newD)) newD=0;


	wGLR = Zg * (1 - L) / (Zl*L);

	return newD;
}

void Scaling(FType *A, int *row, int *col, FType *b) {
	register int i, j, MSize;
	FType rowMax, a;

	MSize = row[0];
	row[0] = 0;

	for (i = 0; i < MSize; i++) {
		rowMax = 0;
		for (j = row[i]; j < row[i + 1]; j++) {
			a = fabs(A[j]);
			if (a > rowMax) rowMax = a;
		}

		if (rowMax) {
			for (j = row[i]; j < row[i + 1]; j++) {
				A[j] /= rowMax;
			}

			b[i] /= rowMax;
		}
	}

	row[0] = MSize;
}

double gasRateCh(int Ix, int Iy, int Iz, double cOil) {
	register int i, j, k;
	FType L, F, dF, tempD;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG;
	FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType ESum;
	FType tempSQR;
	FType newD;
	FType Lb = 10;

	//nesbat = cOil;

	for (i = 0; i < Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / 101325)*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / 288.705555555556));
		STcomp[i][2] = (cOil*comp[Ix][Iy][Iz][i][0] + comp[Ix][Iy][Iz][i][1]) / (1 + cOil);
		//STcomp[i][2]=comp[Ix][Iy][Iz][i][0];
	}


	k = 0;
	do {
		k++;
		L = 0.7;
		j = 0;
		do {
			F = 0;
			dF = 0;
			for (i = 0; i < Nc; i++) {
				if (Ki[i] >= 0) tempD = (1 - Ki[i]) / (L + (1 - L)*Ki[i]);
				else tempD = 1 / (L - 1);
				F += STcomp[i][2] * tempD;
				dF -= STcomp[i][2] * tempD*tempD;
			}
			tempD = F / dF;
			L -= tempD;
			j++;
		} while (((tempD*tempD) > NEWTON_TOL) && (j < MAX_FLASH_N_ITR));

		if (L > 1) {
			L = 1;
			//break;
		}
		else if (L < 0) {
			L = 0;
			//break;
		}

		ESum = fabs(Lb - L);
		Lb = L;

		for (i = 0; i < Nc; i++) {
			if (Ki[i] >= 0) {
				STcomp[i][0] = STcomp[i][2] / (L + Ki[i] * (1 - L));
				STcomp[i][1] = Ki[i] * STcomp[i][0];
			}
			else {
				STcomp[i][0] = 0;
				STcomp[i][1] = STcomp[i][2] / (1 - L);
				//STcomp[i][1] = STcomp[i][2];
			}
		}

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += STcomp[i][0] * fluidProp[i][EOS_B];
			Bg += STcomp[i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += STcomp[i][0] * STcomp[j][0] * tempSQR;
				Ag += STcomp[i][1] * STcomp[j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= 0.0175848235251236;
		Bl *= 42.2111625483491;
		Ag *= 0.0175848235251236;     //gas-oil
		Bg *= 42.2111625483491;


		if (SRK) {
			Ul = Bl;
			Wl = 0;
			Ug = Bg;
			Wg = 0;
			d1 = 1;
			d2 = 0;
		}
		else if (PR) {
			Ul = 2 * Bl;
			Wl = Bl;
			Ug = 2 * Bg;
			Wg = Bg;
			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl * Ul - Ul - Wl * Wl, -(Al*Bl - Bl * Wl*Wl - Wl * Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg * Ug - Ug - Wg * Wg, -(Ag*Bg - Bg * Wg*Wg - Wg * Wg), 'g');


		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += STcomp[j][0] * tempSQR;
				tg += STcomp[j][1] * tempSQR;
			}
			t1l = fluidProp[i][EOS_B] / Cbl;
			t1g = fluidProp[i][EOS_B] / Cbg;

			EEl = 2 * tl / Cal - t1l;
			EEg = 2 * tg / Cag - t1g;

			DDl = Al / (Bl*(d1 - d2));
			DDg = Ag / (Bg*(d1 - d2));

			FFl = log((Zl + d2 * Bl) / (Zl + d1 * Bl));
			FFg = log((Zg + d2 * Bg) / (Zg + d1 * Bg));

			phiL = exp(t1l*(Zl - 1) - log(Zl - Bl) + DDl * EEl*FFl);
			phiG = exp(t1g*(Zg - 1) - log(Zg - Bg) + DDg * EEg*FFg);

			fL = STcomp[i][0] * phiL;		//eliminated P
			fG = STcomp[i][1] * phiG;	//gas-oil

			if (STcomp[i][2]) {
				if (fG) {
					Ki[i] *= fL / fG;
				}
				else Ki[i] = -1;
				if (fL != fL) Ki[i] = -1;
				if (fG != fG) Ki[i] = 1;
			}
		}
	} while ((ESum > FLASH_TOL) && (k < 500));
	//////////////////////////////////////////////////////////////

	newD = 42.2111625483491 / (Zg*(1 - L));
	if (newD != newD) newD = 0;
	//if (!_finite(newD)) newD=0;


	//wGLR = Zg*(1 - L) / (Zl*L);

	return newD;
}

#endif