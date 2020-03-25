#ifndef INIT_CALC_H
#define INIT_CALC_H

#include <math.h>
#include <stdlib.h>
#include <malloc.h>


//#include "Globals.h"
//#include "MKL_Solver.h"
//#include "Do_Time_Step.h"
#include "Functions.h"



void InitWells(void) {
	register int i;

	for (i = 0; i < wellNO; i++) WellCondition[i] = false;
}
void IFK_Calc(void) {	//calculate kA/h
	register int i, j, k;

	for (i = 0; i<Nx; i++)
		for (j = 0; j<Ny; j++)
			for (k = 0; k<Nz; k++) {
				if ((Nx != 1) && (i != (Nx - 1)) && (porosity[i][j][k] > 1e-10) && (porosity[i+1][j][k] > 1e-10)) IFTran[i][j][k][0] = SI_C * 2 * gridDim[Nx + j] * gridDim[Nx + Ny + k] / (gridDim[i] / perm[i][j][k][0] + gridDim[i + 1] / perm[i + 1][j][k][0]);
				else IFTran[i][j][k][0] = 0;
				if ((Ny != 1) && (j != (Ny - 1)) && (porosity[i][j][k] > 1e-10) && (porosity[i][j+1][k] > 1e-10)) IFTran[i][j][k][1] = SI_C * 2 * gridDim[i] * gridDim[Nx + Ny + k] / (gridDim[Nx + j] / perm[i][j][k][1] + gridDim[Nx + j + 1] / perm[i][j + 1][k][1]);
				else IFTran[i][j][k][1] = 0;
				if ((Nz != 1) && (k != (Nz - 1)) && (porosity[i][j][k] > 1e-10) && (porosity[i][j][k+1] > 1e-10)) IFTran[i][j][k][2] = SI_C * 2 * gridDim[i] * gridDim[Nx + j] / (gridDim[Nx + Ny + k] / perm[i][j][k][2] + gridDim[Nx + Ny + k + 1] / perm[i][j][k + 1][2]);
				else IFTran[i][j][k][2] = 0;

				if ((Nx != 1) && (i != (Nx - 1)) && (fporosity[i][j][k] > 1e-10) && (fporosity[i + 1][j][k] > 1e-10)) fIFTran[i][j][k][0] = SI_C * 2 * gridDim[Nx + j] * gridDim[Nx + Ny + k] / (gridDim[i] / fperm[i][j][k][0] + gridDim[i + 1] / fperm[i + 1][j][k][0]);
				else fIFTran[i][j][k][0] = 0;
				if ((Ny != 1) && (j != (Ny - 1)) && (fporosity[i][j][k] > 1e-10) && (fporosity[i][j+1][k] > 1e-10)) fIFTran[i][j][k][1] = SI_C * 2 * gridDim[i] * gridDim[Nx + Ny + k] / (gridDim[Nx + j] / fperm[i][j][k][1] + gridDim[Nx + j + 1] / fperm[i][j + 1][k][1]);
				else fIFTran[i][j][k][1] = 0;
				if ((Nz != 1) && (k != (Nz - 1)) && (fporosity[i][j][k] > 1e-10) && (fporosity[i][j][k+1] > 1e-10)) fIFTran[i][j][k][2] = SI_C * 2 * gridDim[i] * gridDim[Nx + j] / (gridDim[Nx + Ny + k] / fperm[i][j][k][2] + gridDim[Nx + Ny + k + 1] / fperm[i][j][k + 1][2]);
				else fIFTran[i][j][k][2] = 0;

			}

	for (i = 0; i<50; i++) {
		wellrate[i] = ((FType)(i + 1)) / 50.0;
	}
}

void Tr_Mu(void) {
	register int i;
	FType gamma;
	for (i = 0; i<Nc; i++) {
		//fluidProp[i][TR]=resTemp/fluidProp[i][TCRIT];
		gamma = pow(fluidProp[i][TCRIT], (1.0 / 6)) / (sqrt(fluidProp[i][MW])*pow((fluidProp[i][PCRIT] / 101325), (2.0 / 3)));
		/*if (fluidProp[i][TR]>1.5) {
		fluidProp[i][MU]=(1.778e-4)*pow((4.58*fluidProp[i][TR]-1.67), 0.625)/gamma;
		}
		else {
		fluidProp[i][MU]=(3.4e-4)*pow(fluidProp[i][TR], 0.94)/gamma;
		}*/

		fluidProp[i][MU] = (4.61*pow(fluidProp[i][TR], 0.618) - 2.04*exp(-0.449*fluidProp[i][TR]) + 1.94*exp(-4.058*fluidProp[i][TR]) + 0.1)*1e-4 / gamma;

	}
}

void EOS_Init(void) {
	register int i;
	FType m, ac;

	for (i = 0; i<Nc; i++) {
		fluidProp[i][COLDIA] = 0.1866*pow(fluidProp[i][VCRIT], 1.0 / 3.0)*pow(fluidProp[i][ZCRIT], -1.2);
		fluidProp[i][EPDIA] = 65.3*fluidProp[i][TCRIT] * pow(fluidProp[i][ZCRIT], 3.6);

		fluidProp[i][TR] = resTemp / fluidProp[i][TCRIT];
		fluidProp[i][VCRIT] /= 1000000.0;	
		
		//fluidProp[i][COLDIA] = (2.3551 - 0.087*fluidProp[i][AC])*pow((fluidProp[i][TCRIT] * 101325 / fluidProp[i][PCRIT]), (1.0 / 3.0));
		//fluidProp[i][EPDIA] = (0.7915 + 0.1963*fluidProp[i][AC])*fluidProp[i][TCRIT];
	}

	if (SRK) {
		for (i = 0; i<Nc; i++) {
			ac = 0.42747*RGAS*RGAS*fluidProp[i][TCRIT] * fluidProp[i][TCRIT] / fluidProp[i][PCRIT];
			m = 0.48508 + 1.55171*fluidProp[i][AC] - 0.15613*fluidProp[i][AC] * fluidProp[i][AC];
			fluidProp[i][EOS_A] = ac*(1 + m*(1 - sqrt(fluidProp[i][TR])))*(1 + m*(1 - sqrt(fluidProp[i][TR])));
			fluidProp[i][EOS_B] = 0.08664*RGAS*fluidProp[i][TCRIT] / fluidProp[i][PCRIT];
		}
	}
	else if (PR) {
		for (i = 0; i<Nc; i++) {
			ac = 0.457235*RGAS*RGAS*fluidProp[i][TCRIT] * fluidProp[i][TCRIT] / fluidProp[i][PCRIT];
			m = 0.3796 + 1.485*fluidProp[i][AC] - 0.1644*fluidProp[i][AC] * fluidProp[i][AC] + 0.01667*fluidProp[i][AC] * fluidProp[i][AC] * fluidProp[i][AC];
			fluidProp[i][EOS_A] = ac*(1 + m*(1 - sqrt(fluidProp[i][TR])))*(1 + m*(1 - sqrt(fluidProp[i][TR])));
			fluidProp[i][EOS_B] = 0.077796*RGAS*fluidProp[i][TCRIT] / fluidProp[i][PCRIT];
		}
	}
}

void Edge_Trans(void) {
	register int Ix, Iy, Iz;

	for (Iy = 0; Iy<(Ny + 1); Iy++)
		for (Iz = 0; Iz<(Nz + 1); Iz++) {
			trans[0][Iy][Iz][0][0] = 0;
			trans[0][Iy][Iz][1][0] = 0;

			trans[Nx][Iy][Iz][0][0] = 0;
			trans[Nx][Iy][Iz][1][0] = 0;

			Wtran[0][Iy][Iz][0] = 0;
			Wtran[Nx][Iy][Iz][0] = 0;

			dWtran[0][Iy][Iz][0] = 0;
			dWtran[Nx][Iy][Iz][0] = 0;

			transS[0][Iy][Iz][0][0] = 1;
			transS[0][Iy][Iz][1][0] = 1;
			transS[0][Iy][Iz][2][0] = 1;


			transS[Nx][Iy][Iz][0][0] = 0;
			transS[Nx][Iy][Iz][1][0] = 0;
			transS[Nx][Iy][Iz][2][0] = 0;

			trans[0][Iy][Iz][2][0] = 0;
			trans[0][Iy][Iz][3][0] = 0;

			trans[Nx][Iy][Iz][2][0] = 0;
			trans[Nx][Iy][Iz][3][0] = 0;
			//////////////////////////
			trans[0][Iy][Iz][4][0] = 0;
			trans[0][Iy][Iz][5][0] = 0;

			trans[Nx][Iy][Iz][4][0] = 0;
			trans[Nx][Iy][Iz][5][0] = 0;

			trans[0][Iy][Iz][6][0] = 0;
			trans[0][Iy][Iz][7][0] = 0;

			trans[Nx][Iy][Iz][6][0] = 0;
			trans[Nx][Iy][Iz][7][0] = 0;

			trans[0][Iy][Iz][8][0] = 0;
			trans[0][Iy][Iz][9][0] = 0;

			trans[Nx][Iy][Iz][8][0] = 0;
			trans[Nx][Iy][Iz][9][0] = 0;

			trans[0][Iy][Iz][10][0] = 0;
			trans[0][Iy][Iz][11][0] = 0;

			trans[Nx][Iy][Iz][10][0] = 0;
			trans[Nx][Iy][Iz][11][0] = 0;
		}

	for (Ix = 0; Ix<(Nx + 1); Ix++)
		for (Iz = 0; Iz<(Nz + 1); Iz++) {
			trans[Ix][0][Iz][0][1] = 0;
			trans[Ix][0][Iz][1][1] = 0;

			trans[Ix][Ny][Iz][0][1] = 0;
			trans[Ix][Ny][Iz][1][1] = 0;

			Wtran[Ix][0][Iz][1] = 0;
			Wtran[Ix][Ny][Iz][1] = 0;

			dWtran[Ix][0][Iz][1] = 0;
			dWtran[Ix][Ny][Iz][1] = 0;


			transS[Ix][0][Iz][0][1] = 1;
			transS[Ix][0][Iz][1][1] = 1;
			transS[Ix][0][Iz][2][1] = 1;

			transS[Ix][Ny][Iz][0][1] = 0;
			transS[Ix][Ny][Iz][1][1] = 0;
			transS[Ix][Ny][Iz][2][1] = 0;

			trans[Ix][0][Iz][2][1] = 0;
			trans[Ix][0][Iz][3][1] = 0;

			trans[Ix][Ny][Iz][2][1] = 0;
			trans[Ix][Ny][Iz][3][1] = 0;
			/////////////////////////////////
			trans[Ix][0][Iz][4][1] = 0;
			trans[Ix][0][Iz][5][1] = 0;

			trans[Ix][Ny][Iz][4][1] = 0;
			trans[Ix][Ny][Iz][5][1] = 0;

			trans[Ix][0][Iz][6][1] = 0;
			trans[Ix][0][Iz][7][1] = 0;

			trans[Ix][Ny][Iz][6][1] = 0;
			trans[Ix][Ny][Iz][7][1] = 0;

			trans[Ix][0][Iz][8][1] = 0;
			trans[Ix][0][Iz][9][1] = 0;

			trans[Ix][Ny][Iz][8][1] = 0;
			trans[Ix][Ny][Iz][9][1] = 0;

			trans[Ix][0][Iz][10][1] = 0;
			trans[Ix][0][Iz][11][1] = 0;

			trans[Ix][Ny][Iz][10][1] = 0;
			trans[Ix][Ny][Iz][11][1] = 0;
		}

	for (Ix = 0; Ix<(Nx + 1); Ix++)
		for (Iy = 0; Iy<(Ny + 1); Iy++) {
			trans[Ix][Iy][0][0][2] = 0;
			trans[Ix][Iy][0][1][2] = 0;

			trans[Ix][Iy][Nz][0][2] = 0;
			trans[Ix][Iy][Nz][1][2] = 0;

			Wtran[Ix][Iy][0][2] = 0;
			Wtran[Ix][Iy][Nz][2] = 0;

			dWtran[Ix][Iy][0][2] = 0;
			dWtran[Ix][Iy][Nz][2] = 0;


			transS[Ix][Iy][0][0][2] = 1;
			transS[Ix][Iy][0][1][2] = 1;
			transS[Ix][Iy][0][2][2] = 1;

			transS[Ix][Iy][Nz][0][2] = 0;
			transS[Ix][Iy][Nz][1][2] = 0;
			transS[Ix][Iy][Nz][2][2] = 0;

			trans[Ix][Iy][0][2][2] = 0;
			trans[Ix][Iy][0][3][2] = 0;

			trans[Ix][Iy][Nz][2][2] = 0;
			trans[Ix][Iy][Nz][3][2] = 0;
			///////////////////////////////
			trans[Ix][Iy][0][4][2] = 0;
			trans[Ix][Iy][0][5][2] = 0;

			trans[Ix][Iy][Nz][4][2] = 0;
			trans[Ix][Iy][Nz][5][2] = 0;

			trans[Ix][Iy][0][6][2] = 0;
			trans[Ix][Iy][0][7][2] = 0;

			trans[Ix][Iy][Nz][6][2] = 0;
			trans[Ix][Iy][Nz][7][2] = 0;

			trans[Ix][Iy][0][8][2] = 0;
			trans[Ix][Iy][0][9][2] = 0;

			trans[Ix][Iy][Nz][8][2] = 0;
			trans[Ix][Iy][Nz][9][2] = 0;

			trans[Ix][Iy][0][10][2] = 0;
			trans[Ix][Iy][0][11][2] = 0;

			trans[Ix][Iy][Nz][10][2] = 0;
			trans[Ix][Iy][Nz][11][2] = 0;
		}
	////////////////////////////////////////////////
	for (Iy = 1; Iy<(Ny + 1); Iy++)
		for (Iz = 1; Iz<(Nz + 1); Iz++) {
			P[0][Iy][Iz][0] = 0;
			P[0][Iy][Iz][1] = 0;
			P[0][Iy][Iz][2] = 0;

			P[Nx + 1][Iy][Iz][0] = 0;
			P[Nx + 1][Iy][Iz][1] = 0;
			P[Nx + 1][Iy][Iz][2] = 0;

			blockH[0][Iy][Iz] = 0;
			blockH[Nx + 1][Iy][Iz] = 0;
		}

	for (Ix = 1; Ix<(Nx + 1); Ix++)
		for (Iz = 1; Iz<(Nz + 1); Iz++) {
			P[Ix][0][Iz][0] = 0;
			P[Ix][0][Iz][1] = 0;
			P[Ix][0][Iz][2] = 0;

			P[Ix][Ny + 1][Iz][0] = 0;
			P[Ix][Ny + 1][Iz][1] = 0;
			P[Ix][Ny + 1][Iz][2] = 0;


			blockH[Ix][0][Iz] = 0;
			blockH[Ix][Ny + 1][Iz] = 0;
		}

	for (Ix = 1; Ix<(Nx + 1); Ix++)
		for (Iy = 1; Iy<(Ny + 1); Iy++) {
			P[Ix][Iy][0][0] = 0;
			P[Ix][Iy][0][1] = 0;
			P[Ix][Iy][0][2] = 0;

			P[Ix][Iy][Nz + 1][0] = 0;
			P[Ix][Iy][Nz + 1][1] = 0;
			P[Ix][Iy][Nz + 1][2] = 0;


			blockH[Ix][Iy][0] = 0;
			blockH[Ix][Iy][Nz + 1] = 0;
		}
}


void Edge_fTrans(void) {
	register int Ix, Iy, Iz;

	for (Iy = 0; Iy < (Ny + 1); Iy++)
		for (Iz = 0; Iz < (Nz + 1); Iz++) {
			ftrans[0][Iy][Iz][0][0] = 0;
			ftrans[0][Iy][Iz][1][0] = 0;

			ftrans[Nx][Iy][Iz][0][0] = 0;
			ftrans[Nx][Iy][Iz][1][0] = 0;

			fWtran[0][Iy][Iz][0] = 0;
			fWtran[Nx][Iy][Iz][0] = 0;

			fdWtran[0][Iy][Iz][0] = 0;
			fdWtran[Nx][Iy][Iz][0] = 0;

			ftransS[0][Iy][Iz][0][0] = 1;
			ftransS[0][Iy][Iz][1][0] = 1;
			ftransS[0][Iy][Iz][2][0] = 1;


			ftransS[Nx][Iy][Iz][0][0] = 0;
			ftransS[Nx][Iy][Iz][1][0] = 0;
			ftransS[Nx][Iy][Iz][2][0] = 0;

			ftrans[0][Iy][Iz][2][0] = 0;
			ftrans[0][Iy][Iz][3][0] = 0;

			ftrans[Nx][Iy][Iz][2][0] = 0;
			ftrans[Nx][Iy][Iz][3][0] = 0;
			//////////////////////////
			ftrans[0][Iy][Iz][4][0] = 0;
			ftrans[0][Iy][Iz][5][0] = 0;

			ftrans[Nx][Iy][Iz][4][0] = 0;
			ftrans[Nx][Iy][Iz][5][0] = 0;

			ftrans[0][Iy][Iz][6][0] = 0;
			ftrans[0][Iy][Iz][7][0] = 0;

			ftrans[Nx][Iy][Iz][6][0] = 0;
			ftrans[Nx][Iy][Iz][7][0] = 0;

			ftrans[0][Iy][Iz][8][0] = 0;
			ftrans[0][Iy][Iz][9][0] = 0;

			ftrans[Nx][Iy][Iz][8][0] = 0;
			ftrans[Nx][Iy][Iz][9][0] = 0;

			ftrans[0][Iy][Iz][10][0] = 0;
			ftrans[0][Iy][Iz][11][0] = 0;

			ftrans[Nx][Iy][Iz][10][0] = 0;
			ftrans[Nx][Iy][Iz][11][0] = 0;
		}

	for (Ix = 0; Ix < (Nx + 1); Ix++)
		for (Iz = 0; Iz < (Nz + 1); Iz++) {
			ftrans[Ix][0][Iz][0][1] = 0;
			ftrans[Ix][0][Iz][1][1] = 0;

			ftrans[Ix][Ny][Iz][0][1] = 0;
			ftrans[Ix][Ny][Iz][1][1] = 0;

			fWtran[Ix][0][Iz][1] = 0;
			fWtran[Ix][Ny][Iz][1] = 0;

			fdWtran[Ix][0][Iz][1] = 0;
			fdWtran[Ix][Ny][Iz][1] = 0;


			ftransS[Ix][0][Iz][0][1] = 1;
			ftransS[Ix][0][Iz][1][1] = 1;
			ftransS[Ix][0][Iz][2][1] = 1;

			ftransS[Ix][Ny][Iz][0][1] = 0;
			ftransS[Ix][Ny][Iz][1][1] = 0;
			ftransS[Ix][Ny][Iz][2][1] = 0;

			ftrans[Ix][0][Iz][2][1] = 0;
			ftrans[Ix][0][Iz][3][1] = 0;

			ftrans[Ix][Ny][Iz][2][1] = 0;
			ftrans[Ix][Ny][Iz][3][1] = 0;
			/////////////////////////////////
			ftrans[Ix][0][Iz][4][1] = 0;
			ftrans[Ix][0][Iz][5][1] = 0;

			ftrans[Ix][Ny][Iz][4][1] = 0;
			ftrans[Ix][Ny][Iz][5][1] = 0;

			ftrans[Ix][0][Iz][6][1] = 0;
			ftrans[Ix][0][Iz][7][1] = 0;

			ftrans[Ix][Ny][Iz][6][1] = 0;
			ftrans[Ix][Ny][Iz][7][1] = 0;

			ftrans[Ix][0][Iz][8][1] = 0;
			ftrans[Ix][0][Iz][9][1] = 0;

			ftrans[Ix][Ny][Iz][8][1] = 0;
			ftrans[Ix][Ny][Iz][9][1] = 0;

			ftrans[Ix][0][Iz][10][1] = 0;
			ftrans[Ix][0][Iz][11][1] = 0;

			ftrans[Ix][Ny][Iz][10][1] = 0;
			ftrans[Ix][Ny][Iz][11][1] = 0;
		}

	for (Ix = 0; Ix < (Nx + 1); Ix++)
		for (Iy = 0; Iy < (Ny + 1); Iy++) {
			ftrans[Ix][Iy][0][0][2] = 0;
			ftrans[Ix][Iy][0][1][2] = 0;

			ftrans[Ix][Iy][Nz][0][2] = 0;
			ftrans[Ix][Iy][Nz][1][2] = 0;

			fWtran[Ix][Iy][0][2] = 0;
			fWtran[Ix][Iy][Nz][2] = 0;

			fdWtran[Ix][Iy][0][2] = 0;
			fdWtran[Ix][Iy][Nz][2] = 0;


			ftransS[Ix][Iy][0][0][2] = 1;
			ftransS[Ix][Iy][0][1][2] = 1;
			ftransS[Ix][Iy][0][2][2] = 1;

			ftransS[Ix][Iy][Nz][0][2] = 0;
			ftransS[Ix][Iy][Nz][1][2] = 0;
			ftransS[Ix][Iy][Nz][2][2] = 0;

			ftrans[Ix][Iy][0][2][2] = 0;
			ftrans[Ix][Iy][0][3][2] = 0;

			ftrans[Ix][Iy][Nz][2][2] = 0;
			ftrans[Ix][Iy][Nz][3][2] = 0;
			///////////////////////////////
			ftrans[Ix][Iy][0][4][2] = 0;
			ftrans[Ix][Iy][0][5][2] = 0;

			ftrans[Ix][Iy][Nz][4][2] = 0;
			ftrans[Ix][Iy][Nz][5][2] = 0;

			ftrans[Ix][Iy][0][6][2] = 0;
			ftrans[Ix][Iy][0][7][2] = 0;

			ftrans[Ix][Iy][Nz][6][2] = 0;
			ftrans[Ix][Iy][Nz][7][2] = 0;

			ftrans[Ix][Iy][0][8][2] = 0;
			ftrans[Ix][Iy][0][9][2] = 0;

			ftrans[Ix][Iy][Nz][8][2] = 0;
			ftrans[Ix][Iy][Nz][9][2] = 0;

			ftrans[Ix][Iy][0][10][2] = 0;
			ftrans[Ix][Iy][0][11][2] = 0;

			ftrans[Ix][Iy][Nz][10][2] = 0;
			ftrans[Ix][Iy][Nz][11][2] = 0;
		}
	////////////////////////////////////////////////
	for (Iy = 1; Iy < (Ny + 1); Iy++)
		for (Iz = 1; Iz < (Nz + 1); Iz++) {
			fP[0][Iy][Iz][0] = 0;
			fP[0][Iy][Iz][1] = 0;
			fP[0][Iy][Iz][2] = 0;

			fP[Nx + 1][Iy][Iz][0] = 0;
			fP[Nx + 1][Iy][Iz][1] = 0;
			fP[Nx + 1][Iy][Iz][2] = 0;

			blockH[0][Iy][Iz] = 0;
			blockH[Nx + 1][Iy][Iz] = 0;
		}

	for (Ix = 1; Ix < (Nx + 1); Ix++)
		for (Iz = 1; Iz < (Nz + 1); Iz++) {
			fP[Ix][0][Iz][0] = 0;
			fP[Ix][0][Iz][1] = 0;
			fP[Ix][0][Iz][2] = 0;

			fP[Ix][Ny + 1][Iz][0] = 0;
			fP[Ix][Ny + 1][Iz][1] = 0;
			fP[Ix][Ny + 1][Iz][2] = 0;


			blockH[Ix][0][Iz] = 0;
			blockH[Ix][Ny + 1][Iz] = 0;
		}

	for (Ix = 1; Ix < (Nx + 1); Ix++)
		for (Iy = 1; Iy < (Ny + 1); Iy++) {
			fP[Ix][Iy][0][0] = 0;
			fP[Ix][Iy][0][1] = 0;
			fP[Ix][Iy][0][2] = 0;

			fP[Ix][Iy][Nz + 1][0] = 0;
			fP[Ix][Iy][Nz + 1][1] = 0;
			fP[Ix][Iy][Nz + 1][2] = 0;


			blockH[Ix][Iy][0] = 0;
			blockH[Ix][Iy][Nz + 1] = 0;
		}
}


void CPlus_Props(void) {
	int i;
	FType Tb, API, Tc, Pc, Zc;

	for (i = PNc; i<Nc; i++) {
		Tb = fluidProp[i][TB] * 5 / 9 - 459.67;
		API = 141.5 / fluidProp[i][SG] - 131.5;

		//Cavett Correlation
		Tc = 768.071 + 1.7134*Tb - 1.0834e-3*Tb*Tb + 3.889e-7*Tb*Tb*Tb - 8.9213e-3*Tb*API + 5.3095e-6*Tb*Tb*API + 3.2712e-8*Tb*Tb*API*API;
		fluidProp[i][TCRIT] = (9 / 5)*Tc;

		Pc = 2.829 + 9.412e-4*Tb - 3.0475e-6*Tb*Tb + 1.5184e-9*Tb*Tb*Tb - 2.0876e-5*Tb*API + 1.1048e-8*Tb*Tb*API - 4.827e-8*Tb*API*API + 1.395e-10*Tb*Tb*API*API;
		fluidProp[i][PCRIT] = 6897.549353301566*exp(Pc*LN10);

		//Edmister Correlation
		fluidProp[i][AC] = (3.0 / 7.0)*log(fluidProp[i][PCRIT] / 101325) / (LN10*(fluidProp[i][TCRIT] / fluidProp[i][TB] - 1)) - 1;

		//Pitzer correlation
		Zc = 0.2901 - 0.0879*fluidProp[i][AC];


		fluidProp[i][VCRIT] = Zc*RGAS*fluidProp[i][TCRIT] / fluidProp[i][PCRIT];
	}

}
void WaterSat(void) {
	register int i, j, k, n;
	FType watP;

	for (i = 0; i<Nx; i++)
		for (j = 0; j<Ny; j++)
			for (k = 0; k<Nz; k++) {
				watP = G_ACC*(blockH[i + 1][j + 1][k + 1] - WOCHeight)*watRo;


				for (n = 0; n<Nswt; n++) {
					if (watP>swt[n][CAPILLARYPRESSURE]) break;
				}
				if (n == 0) {
					sat[i][j][k][0] = swt[0][SAT];
				}
				else if (n == Nswt) {
					sat[i][j][k][0] = 1;
				}
				else {
					sat[i][j][k][0] = (swt[n][SAT] - swt[n - 1][SAT]) / (swt[n][CAPILLARYPRESSURE] - swt[n - 1][CAPILLARYPRESSURE])*(watP - swt[n][CAPILLARYPRESSURE]) + swt[n][SAT];
				}
			}
}
void SuccFalsh(int Ix, int Iy, int Iz) {
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

	for (i = 0; i<Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / P[Ix + 1][Iy + 1][Iz + 1][1])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));
	}

	k = 0;
	do {
		L = 0.6;
		j = 0;
		do {
			F = 0;
			dF = 0;
			for (i = 0; i<Nc; i++) {
				if (Ki[i] >= 0) tempD = (1 - Ki[i]) / (L + (1 - L)*Ki[i]);
				else tempD = 1 / (L - 1);
				F += comp[Ix][Iy][Iz][i][2] * tempD;
				dF -= comp[Ix][Iy][Iz][i][2] * tempD*tempD;
			}
			tempD = F / dF;
			L -= tempD;
			j++;
		} while (((tempD*tempD)>NEWTON_TOL) && (j<MAX_FLASH_N_ITR));

		if (L>1) {
			L = 1;
			//break;
		}
		else if (L<0) {
			L = 0;
			//break;
		}

		ESum = fabs(Lb - L);
		Lb = L;

		for (i = 0; i<Nc; i++) {
			if (Ki[i] >= 0) {
				comp[Ix][Iy][Iz][i][0] = comp[Ix][Iy][Iz][i][2] / (L + Ki[i] * (1 - L));
				comp[Ix][Iy][Iz][i][1] = Ki[i] * comp[Ix][Iy][Iz][i][0];
			}
			else {
				comp[Ix][Iy][Iz][i][0] = 0;
				comp[Ix][Iy][Iz][i][1] = comp[Ix][Iy][Iz][i][2] / (1 - L);
			}
		}

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i<Nc; i++) {
			Bl += comp[Ix][Iy][Iz][i][0] * fluidProp[i][EOS_B];
			Bg += comp[Ix][Iy][Iz][i][1] * fluidProp[i][EOS_B];
			for (j = 0; j<Nc; j++) {
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
		Ag *= P[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*RGAS*resTemp*resTemp);		//gas-oil
		Bg *= P[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*resTemp);


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

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl*Ul - Ul - Wl*Wl, -(Al*Bl - Bl*Wl*Wl - Wl*Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg*Ug - Ug - Wg*Wg, -(Ag*Bg - Bg*Wg*Wg - Wg*Wg), 'g');


		for (i = 0; i<Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j<Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += comp[Ix][Iy][Iz][j][0] * tempSQR;
				tg += comp[Ix][Iy][Iz][j][1] * tempSQR;
			}
			t1l = fluidProp[i][EOS_B] / Cbl;
			t1g = fluidProp[i][EOS_B] / Cbg;

			EEl = 2 * tl / Cal - t1l;
			EEg = 2 * tg / Cag - t1g;

			DDl = Al / (Bl*(d1 - d2));
			DDg = Ag / (Bg*(d1 - d2));

			FFl = log((Zl + d2*Bl) / (Zl + d1*Bl));
			FFg = log((Zg + d2*Bg) / (Zg + d1*Bg));

			phiL = exp(t1l*(Zl - 1) - log(Zl - Bl) + DDl*EEl*FFl);
			phiG = exp(t1g*(Zg - 1) - log(Zg - Bg) + DDg*EEg*FFg);

			fL = comp[Ix][Iy][Iz][i][0] * phiL;		//eliminated P
			fG = comp[Ix][Iy][Iz][i][1] * phiG;	//gas-oil

			if (comp[Ix][Iy][Iz][i][2]) {
				if (fG) {
					Ki[i] *= fL / fG;
				}
				else Ki[i] = -1;
				if (fL != fL) Ki[i] = -1;
				if (fG != fG) Ki[i] = 1;
			}
		}
		k++;
	} while ((ESum>FLASH_TOL) && (k<500));
	//////////////////////////////////////////////////////////////

	if (fabs(L)<DZERO) {
		phaseStat[Ix][Iy][Iz] = -1;
		//L=0;
	}
	else if (fabs(1 - L)<DZERO) {
		phaseStat[Ix][Iy][Iz] = 1;
		//L=1;
	}
	else {
		phaseStat[Ix][Iy][Iz] = 0;
	}


	newD = Zl*L / (Zl*L + Zg*(1 - L));

	sat[Ix][Iy][Iz][1] = newD*(1 - sat[Ix][Iy][Iz][0]);
	sat[Ix][Iy][Iz][2] = (1 - newD)*(1 - sat[Ix][Iy][Iz][0]);

	blockFProps[Ix][Iy][Iz][RO][0] = P[Ix + 1][Iy + 1][Iz + 1][1] / (Zl*RGAS*resTemp);
	blockFProps[Ix][Iy][Iz][RO][1] = P[Ix + 1][Iy + 1][Iz + 1][1] / (Zg*RGAS*resTemp);	//gas-oil

}

void fSuccFalsh(int Ix, int Iy, int Iz) {
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
		Ki[i] = (fluidProp[i][PCRIT] / fP[Ix + 1][Iy + 1][Iz + 1][1])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));
	}

	k = 0;
	do {
		L = 0.6;
		j = 0;
		do {
			F = 0;
			dF = 0;
			for (i = 0; i < Nc; i++) {
				if (Ki[i] >= 0) tempD = (1 - Ki[i]) / (L + (1 - L)*Ki[i]);
				else tempD = 1 / (L - 1);
				F += fcomp[Ix][Iy][Iz][i][2] * tempD;
				dF -= fcomp[Ix][Iy][Iz][i][2] * tempD*tempD;
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
				fcomp[Ix][Iy][Iz][i][0] = fcomp[Ix][Iy][Iz][i][2] / (L + Ki[i] * (1 - L));
				fcomp[Ix][Iy][Iz][i][1] = Ki[i] * fcomp[Ix][Iy][Iz][i][0];
			}
			else {
				fcomp[Ix][Iy][Iz][i][0] = 0;
				fcomp[Ix][Iy][Iz][i][1] = fcomp[Ix][Iy][Iz][i][2] / (1 - L);
			}
		}

		/////////////////////////////////////////////////////////////

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
		Ag *= fP[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*RGAS*resTemp*resTemp);		//gas-oil
		Bg *= fP[Ix + 1][Iy + 1][Iz + 1][1] / (RGAS*resTemp);


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

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl*Ul - Ul - Wl*Wl, -(Al*Bl - Bl*Wl*Wl - Wl*Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg*Ug - Ug - Wg*Wg, -(Ag*Bg - Bg*Wg*Wg - Wg*Wg), 'g');


		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += fcomp[Ix][Iy][Iz][j][0] * tempSQR;
				tg += fcomp[Ix][Iy][Iz][j][1] * tempSQR;
			}
			t1l = fluidProp[i][EOS_B] / Cbl;
			t1g = fluidProp[i][EOS_B] / Cbg;

			EEl = 2 * tl / Cal - t1l;
			EEg = 2 * tg / Cag - t1g;

			DDl = Al / (Bl*(d1 - d2));
			DDg = Ag / (Bg*(d1 - d2));

			FFl = log((Zl + d2*Bl) / (Zl + d1*Bl));
			FFg = log((Zg + d2*Bg) / (Zg + d1*Bg));

			phiL = exp(t1l*(Zl - 1) - log(Zl - Bl) + DDl*EEl*FFl);
			phiG = exp(t1g*(Zg - 1) - log(Zg - Bg) + DDg*EEg*FFg);

			fL = fcomp[Ix][Iy][Iz][i][0] * phiL;		//eliminated P
			fG = fcomp[Ix][Iy][Iz][i][1] * phiG;	//gas-oil

			if (fcomp[Ix][Iy][Iz][i][2]) {
				if (fG) {
					Ki[i] *= fL / fG;
				}
				else Ki[i] = -1;
				if (fL != fL) Ki[i] = -1;
				if (fG != fG) Ki[i] = 1;
			}
		}
		k++;
	} while ((ESum > FLASH_TOL) && (k < 500));
	//////////////////////////////////////////////////////////////

	if (fabs(L) < DZERO) {
		fphaseStat[Ix][Iy][Iz] = -1;
		//L=0;
	}
	else if (fabs(1 - L) < DZERO) {
		fphaseStat[Ix][Iy][Iz] = 1;
		//L=1;
	}
	else {
		fphaseStat[Ix][Iy][Iz] = 0;
	}


	newD = Zl*L / (Zl*L + Zg*(1 - L));

	fsat[Ix][Iy][Iz][1] = newD*(1 - fsat[Ix][Iy][Iz][0]);
	fsat[Ix][Iy][Iz][2] = (1 - newD)*(1 - fsat[Ix][Iy][Iz][0]);

	fblockFProps[Ix][Iy][Iz][RO][0] = fP[Ix + 1][Iy + 1][Iz + 1][1] / (Zl*RGAS*resTemp);
	fblockFProps[Ix][Iy][Iz][RO][1] = fP[Ix + 1][Iy + 1][Iz + 1][1] / (Zg*RGAS*resTemp);	//gas-oil

}


FType Pcgo(int Ix, int Iy, int Iz) {
	register int n;
	FType rd;

	for (n = 0; n<Nsgt; n++) {
		if (sat[Ix][Iy][Iz][2]<sgt[n][SAT]) break;
	}
	if (n == 0) {
		rd = 0;
	}
	else if (n == Nsgt) {
		rd = sgt[Nsgt - 1][CAPILLARYPRESSURE];
	}
	else {
		rd = (sgt[n][CAPILLARYPRESSURE] - sgt[n - 1][CAPILLARYPRESSURE]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(sat[Ix][Iy][Iz][2] - sgt[n][SAT]) + sgt[n][CAPILLARYPRESSURE];
	}

	return rd;
}

FType Pcow(int Ix, int Iy, int Iz) {
	register int n;
	FType rd;

	for (n = 0; n<Nswt; n++) {
		if (sat[Ix][Iy][Iz][0]<swt[n][SAT]) break;
	}
	if (n == 0) {
		rd = swt[0][CAPILLARYPRESSURE];
	}
	else if (n == Nswt) {
		rd = 0;
	}
	else {
		rd = (swt[n][CAPILLARYPRESSURE] - swt[n - 1][CAPILLARYPRESSURE]) / (swt[n][SAT] - swt[n - 1][SAT])*(sat[Ix][Iy][Iz][0] - swt[n][SAT]) + swt[n][CAPILLARYPRESSURE];
	}

	return rd;
}
void AllFlash(void) {
	register int Ix, Iy, Iz;
	for (Ix = 0; Ix<Nx; Ix++)
		for (Iy = 0; Iy<Ny; Iy++)
			for (Iz = 0; Iz<Nz; Iz++) {
				SuccFalsh(Ix, Iy, Iz);
				fSuccFalsh(Ix, Iy, Iz);
				blockFProps[Ix][Iy][Iz][BLOCK_PC][1] = Pcgo(Ix, Iy, Iz);
				blockFProps[Ix][Iy][Iz][BLOCK_PC][0] = Pcow(Ix, Iy, Iz);
				//fblockFProps[Ix][Iy][Iz][BLOCK_PC][1] = Pcgo(Ix, Iy, Iz);
				//fblockFProps[Ix][Iy][Iz][BLOCK_PC][0] = Pcow(Ix, Iy, Iz);
			}

}

/*FType ShapeFactor(int Ix, int Iy, int Iz) {    //Do_cycle 
	return (4 * ((1 / (gridDim[Ix] * gridDim[Ix])) + (1 / (gridDim[Nx + Iy] * gridDim[Nx + Iy])) + (1 / (gridDim[Nx + Ny + Iz] * gridDim[Nx + Ny + Iz]))));
	//return 0;
}*/


FType ShapeFactor_1(int Ix, int Iy, int Iz) {    //Waren & Rooth 
	if (IsConstDif && (!CDiffCoef)) return 0;
	else return (4 * SigmaMF*((1 / (gridDim[Ix] * gridDim[Ix])) + (1 / (gridDim[Nx + Iy] * gridDim[Nx + Iy])) + (1 / (gridDim[Nx + Ny + Iz] * gridDim[Nx + Ny + Iz]))));
}

FType ShapeFactor_2(int Ix, int Iy, int Iz) {    // permeability With Shapefactor 
	return (4 * SigmaMF*((perm[Ix][Iy][Iz][0] / (gridDim[Ix] * gridDim[Ix])) + (perm[Ix][Iy][Iz][1] / (gridDim[Nx + Iy] * gridDim[Nx + Iy])) + (perm[Ix][Iy][Iz][2] / (gridDim[Nx + Ny + Iz] * gridDim[Nx + Ny + Iz]))));
	//return 0;
}

FType ShapeFactor_3(int Ix, int Iy, int Iz) {    // Shapefactor for Equilibrium mole Fraction
	if (IsConstDif && (!CDiffCoef)) return 0;
	else return (4* SigmaMF * ((1 / gridDim[Ix]) + (1 / gridDim[Nx + Iy] ) + (1 / gridDim[Nx + Ny + Iz] )) / Lf); 
}

void QInitValue(void) {
	cumQo = 0;
	cumQg = 0;
	cumQw = 0;

	incCount = 0;
}

void CalcBlockHeight(void) {
	register int i, j, k;
	FType height;

	for (i = 0; i<Nx; i++)
		for (j = 0; j<Ny; j++)
			for (k = (Nz - 1); k >= 0; k--) {
				if (k == (Nz - 1)) {
					height = 0.5*gridDim[Nx + Ny + Nz - 1];
				}
				else {
					height += 0.5*(gridDim[Nx + Ny + k] + gridDim[Nx + Ny + k + 1]);
				}

				blockH[i + 1][j + 1][k + 1] = height;
			}
}



void PJac(void) {
	register int Ix, Iy, Iz;
	int blockN, bSize;
	preConRow[0] = 0;
	preConIndex[0] = 0;
	PressureRow[0] = 0;

	blockN = 18 * Nc*Nc + 35 * Nc + 19;				//n*(2*n+1)+8*n*(2*n+4)+8*2+2*n+3
	totalRow = 2 * Nx*Ny*Nz*(2 * Nc + 4);
	

	totalJac = 0;
	for (Iz = 0; Iz<Nz; Iz++)
		for (Iy = 0; Iy<Ny; Iy++)
			for (Ix = 0; Ix<Nx; Ix++) {
				pJHolder[Ix][Iy][Iz] = totalJac;
				totalJac += blockN;

				if (Ix == 0) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;		//n*(2*n+4)+2
				if (Ix == (Nx - 1)) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (Iy == 0) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (Iy == (Ny - 1)) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (Iz == 0) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (Iz == (Nz - 1)) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (phaseStat[Ix][Iy][Iz] == 1) {
					//ierr = PetscBarrier(NULL);
					totalJac -= 3 * Nc*Nc + 3 * Nc + 1;		//n*(2*n+1)+n*(n+1)+n+1
					totalRow -= Nc + 1;
					bSize = Nc + 3;
					PressureRow[Iz*(Nx*Ny) + Iy * Nx + Ix] += Nc;					
				}
				else if (phaseStat[Ix][Iy][Iz] == -1) {
					//ierr = PetscBarrier(NULL);
					totalJac -= 3 * Nc*Nc + 4 * Nc + 1;		////n*(2*n+1)+n*(n+2)+n+1
					totalRow -= Nc + 1;
					bSize = Nc + 3;
					PressureRow[Iz*(Nx*Ny) + Iy * Nx + Ix] += Nc;				
				}
				else {
					//ierr = PetscBarrier(NULL);
					bSize = 2 * Nc + 4;
					PressureRow[Iz*(Nx*Ny) + Iy * Nx + Ix] += 2 * Nc;					
				}

				if (fphaseStat[Ix][Iy][Iz]) {
					totalJac -= Nc*Nc+Nc;
				}

				SadeqSize[Iz*(Nx*Ny) + Iy*Nx + Ix] = bSize;
				preConRow[Iz*(Nx*Ny) + Iy*Nx + Ix + 1] = preConRow[Iz*(Nx*Ny) + Iy*Nx + Ix] + bSize;
				PressureRow[Iz*(Nx*Ny) + Iy * Nx + Ix + 1] = preConRow[Iz*(Nx*Ny) + Iy * Nx + Ix] + bSize;
				preConIndex[Iz*(Nx*Ny) + Iy*Nx + Ix + 1] = preConIndex[Iz*(Nx*Ny) + Iy*Nx + Ix] + bSize*(bSize - 1) + MRMAXNONZERO;

			}
	preConCSRrow[preConRow[Nx*Ny*Nz]] = preConIndex[Nx*Ny*Nz];
	//FractureBase = preConRow[Nx*Ny*Nz];
	//FraqColBase = totalJac;

	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				fpJHolder[Ix][Iy][Iz] = totalJac;
				totalJac += blockN;

				if (Ix == 0) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (Ix == (Nx - 1)) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (Iy == 0) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (Iy == (Ny - 1)) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (Iz == 0) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (Iz == (Nz - 1)) totalJac -= 2 * Nc*Nc + 4 * Nc + 2;
				if (fphaseStat[Ix][Iy][Iz] == 1) {
					//ierr = PetscBarrier(NULL);
					totalJac -= 3 * Nc*Nc + 3 * Nc + 1;
					totalRow -= Nc + 1;
					bSize = Nc + 3;
					PressureRow[Nx*Ny*Nz + Iz * (Nx*Ny) + Iy * Nx + Ix] += Nc;					
				}
				else if (fphaseStat[Ix][Iy][Iz] == -1) {
					//ierr = PetscBarrier(NULL);
					totalJac -= 3 * Nc*Nc + 4 * Nc + 1;
					totalRow -= Nc + 1;
					bSize = Nc + 3;
					PressureRow[Nx*Ny*Nz + Iz * (Nx*Ny) + Iy * Nx + Ix] += Nc;					
				}
				else {
					//ierr = PetscBarrier(NULL);
					bSize = 2 * Nc + 4;
					PressureRow[Nx*Ny*Nz + Iz * (Nx*Ny) + Iy * Nx + Ix] += 2 * Nc;					
				}

				if (phaseStat[Ix][Iy][Iz]) {
					totalJac -= Nc*Nc + Nc;
				}

				SadeqSize[Nx*Ny*Nz + Iz*(Nx*Ny) + Iy*Nx + Ix] = bSize;
				preConRow[Nx*Ny*Nz + Iz*(Nx*Ny) + Iy*Nx + Ix + 1] = preConRow[Nx*Ny*Nz + Iz*(Nx*Ny) + Iy*Nx + Ix] + bSize;
				PressureRow[Nx*Ny*Nz + Iz * (Nx*Ny) + Iy * Nx + Ix + 1] = preConRow[Nx*Ny*Nz + Iz * (Nx*Ny) + Iy * Nx + Ix] + bSize;
				preConIndex[Nx*Ny*Nz + Iz*(Nx*Ny) + Iy*Nx + Ix + 1] = preConIndex[Nx*Ny*Nz + Iz*(Nx*Ny) + Iy*Nx + Ix] + bSize*(bSize - 1) + MRMAXNONZERO;

			}
	preConCSRrow[preConRow[2 * Nx*Ny*Nz]] = preConIndex[2 * Nx*Ny*Nz];
	//totalJac=jacIndex;
	MatrixSize = preConRow[2 * Nx*Ny*Nz];
}

void CalcPCoeff(void) {
	register int i, j, k;
	PCoeff = 0;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++) {				
				if (P[i + 1][j + 1][k + 1][1] > PCoeff) PCoeff = P[i + 1][j + 1][k + 1][1];
				if (fP[i + 1][j + 1][k + 1][1] > PCoeff) PCoeff = fP[i + 1][j + 1][k + 1][1];
			}
	
	//PCoeff = 5e6;
}


void ManageTSMarker(void) {
	register int i, j, k, n;
	FType temp;

	n = 2 * wellNO;
	//Sort
	for (i = 0; i<n; i++)
		for (j = (i + 1); j <= n; j++) if (TStepMarker[i]>TStepMarker[j]) {
			temp = TStepMarker[i];
			TStepMarker[i] = TStepMarker[j];
			TStepMarker[j] = temp;
		}

	j = 0;
	for (i = 0; i<n; i++) {
		if ((fabs(TStepMarker[i - j] - TStepMarker[i - j + 1])<TSTEPTOL) || (!TStepMarker[i - j])) {
			for (k = (i - j); k<(n - j); k++) TStepMarker[k] = TStepMarker[k + 1];
			j++;
		}
	}

	TSMSize = n - j + 1;

	if (Dt >= TStepMarker[0]) Dt = TStepMarker[0] / 2.1;
}



void FlashInjection(FType *IComposition, FType IPresure, FType *IOComposition, FType *IGComposition, FType &OilFraction, FType &GasFraction, FType &OilDensity, FType &GasDensity, FType &OilViscosity, FType &GasViscosity) {
	register int i, j, k, n;
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
	FType Lb = 10;
	FType s1Ml, s2Ml, s1Mg, s2Mg, tMl, tMg, MuL, MuG, roL, roG, compTCL, compMWL, compPCL, compTCG, compMWG, compPCG;
	FType ethaL, ethaG, eAAl, eAAg, eBBl, eBBg, eCCl, eCCg;
	FType ErL, ErG, ro1L, ro1G;

	for (i = 0; i < Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / IPresure)*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));
	}

	k = 0;
	do {
		L = 0.6;
		j = 0;
		do {
			F = 0;
			dF = 0;
			for (i = 0; i < Nc; i++) {
				if (Ki[i] >= 0) tempD = (1 - Ki[i]) / (L + (1 - L)*Ki[i]);
				else tempD = 1 / (L - 1);
				F += IComposition[i] * tempD;
				dF -= IComposition[i] * tempD*tempD;
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
				IOComposition[i] = IComposition[i] / (L + Ki[i] * (1 - L));
				IGComposition[i] = Ki[i] * IOComposition[i];
			}
			else {
				IOComposition[i] = 0;
				IGComposition[i] = IGComposition[i] / (1 - L);
			}
		}

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += IOComposition[i] * fluidProp[i][EOS_B];
			Bg += IGComposition[i] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += IOComposition[i] * IOComposition[j] * tempSQR;
				Ag += IGComposition[i] * IGComposition[j] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= IPresure / (RGAS*RGAS*resTemp*resTemp);
		Bl *= IPresure / (RGAS*resTemp);
		Ag *= IPresure / (RGAS*RGAS*resTemp*resTemp);		//gas-oil
		Bg *= IPresure / (RGAS*resTemp);


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

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl*Ul - Ul - Wl*Wl, -(Al*Bl - Bl*Wl*Wl - Wl*Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg*Ug - Ug - Wg*Wg, -(Ag*Bg - Bg*Wg*Wg - Wg*Wg), 'g');


		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += IOComposition[j] * tempSQR;
				tg += IGComposition[j] * tempSQR;
			}
			t1l = fluidProp[i][EOS_B] / Cbl;
			t1g = fluidProp[i][EOS_B] / Cbg;

			EEl = 2 * tl / Cal - t1l;
			EEg = 2 * tg / Cag - t1g;

			DDl = Al / (Bl*(d1 - d2));
			DDg = Ag / (Bg*(d1 - d2));

			FFl = log((Zl + d2*Bl) / (Zl + d1*Bl));
			FFg = log((Zg + d2*Bg) / (Zg + d1*Bg));

			phiL = exp(t1l*(Zl - 1) - log(Zl - Bl) + DDl*EEl*FFl);
			phiG = exp(t1g*(Zg - 1) - log(Zg - Bg) + DDg*EEg*FFg);

			fL = IOComposition[i] * phiL;		//eliminated P
			fG = IGComposition[i] * phiG;	//gas-oil

			if (IComposition[i]) {
				if (fG) {
					Ki[i] *= fL / fG;
				}
				else Ki[i] = -1;
				if (fL != fL) Ki[i] = -1;
				if (fG != fG) Ki[i] = 1;
			}
		}
		k++;
	} while ((ESum > FLASH_TOL) && (k < 500));
	//////////////////////////////////////////////////////////////

	if (fabs(L) < DZERO) {
		OilFraction = 0;
		GasFraction = 1;
	}
	else if (fabs(1 - L) < DZERO) {
		OilFraction = 1;
		GasFraction = 0;
	}
	else {
		OilFraction = L;
		GasFraction = 1-L;
	}

	OilDensity = IPresure / (Zl*RGAS*resTemp);
	GasDensity = IPresure / (Zg*RGAS*resTemp);	//gas-oil
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
		tMl = IOComposition[n] * tempSQR;
		tMg = IGComposition[n] * tempSQR;
		s1Ml += tMl*fluidProp[n][MU];
		s2Ml += tMl;
		s1Mg += tMg*fluidProp[n][MU];
		s2Mg += tMg;
		roL += IOComposition[n] * fluidProp[n][VCRIT];
		roG += IGComposition[n] * fluidProp[n][VCRIT];
		compTCL += IOComposition[n] * fluidProp[n][TCRIT];
		compPCL += IOComposition[n] * fluidProp[n][PCRIT] / 101325;
		compMWL += IOComposition[n] * fluidProp[n][MW];
		compTCG += IGComposition[n] * fluidProp[n][TCRIT];
		compPCG += IGComposition[n] * fluidProp[n][PCRIT] / 101325;
		compMWG += IGComposition[n] * fluidProp[n][MW];
	}
	MuL = s1Ml / s2Ml;
	MuG = s1Mg / s2Mg;
	ErL = OilDensity * roL;
	ErG = GasDensity * roG;

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

	OilViscosity = MuL + (ro1L*ro1L*ro1L*ro1L - 1e-4) / ethaL;		//Viscosity, centipoise
	GasViscosity = MuG + (ro1G*ro1G*ro1G*ro1G - 1e-4) / ethaG;
}

#endif