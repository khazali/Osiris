#include "Functions.h"

void MixedCalcIFT(void) {
	register int iix, ppy, n;
	int Ix, Iy, Iz, qxtemp;
	FType pSum;

	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
			if (!phaseStat[Ix][Iy][Iz]) {
				pSum = 0;
				for (n = 0; n < Nc; n++) pSum += fluidProp[n][PARACHOR] * (comp[Ix][Iy][Iz][n][0] * blockFProps[Ix][Iy][Iz][RO][0] - comp[Ix][Iy][Iz][n][1] * blockFProps[Ix][Iy][Iz][RO][1]);
				Bift[Ix][Iy][Iz] = pow(pSum*1e-6, 3.88);
			}
			else Bift[Ix][Iy][Iz] = IFTBASECASE;



			if (!fphaseStat[Ix][Iy][Iz]) {
				pSum = 0;
				for (n = 0; n < Nc; n++) pSum += fluidProp[n][PARACHOR] * (fcomp[Ix][Iy][Iz][n][0] * fblockFProps[Ix][Iy][Iz][RO][0] - fcomp[Ix][Iy][Iz][n][1] * fblockFProps[Ix][Iy][Iz][RO][1]);
				fBift[Ix][Iy][Iz] = pow(pSum*1e-6, 3.88);
			}
			else fBift[Ix][Iy][Iz] = IFTBASECASE;
		}
}

/*void MixedCalcEquilibrium(void) {
	register int iix, ppy;
	int Ix, Iy, Iz, qxtemp;
	
	for (iix = 0; iix < SchindlerListLength; iix++) {
		ppy = SchindlerList[iix];
		if (ppy < Nx*Ny*Nz) {
			//Matrix System
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;

			Equilibrium_Concentrations(Ix, Iy, Iz);
		}
		else {
			//Fracture System
			ppy -= Nx * Ny*Nz;
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;


			Equilibrium_Concentrations(Ix, Iy, Iz);
		}
	}
}*/

void MixedCalcGrav(void) {
	register int iix, ppy, n;
	int Ix, Iy, Iz, qxtemp;
	FType compMWL, compMWG;

	for (iix = 0; iix < SchindlerListLength; iix++) {
		ppy = SchindlerList[iix];
		if (ppy < Nx*Ny*Nz) {
			//Matrix System
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;


			compMWL = 0;
			compMWG = 0;
			for (n = 0; n < Nc; n++) {
				compMWL += comp[Ix][Iy][Iz][n][0] * fluidProp[n][MW];
				compMWG += comp[Ix][Iy][Iz][n][1] * fluidProp[n][MW];
			}
			blockFProps[Ix][Iy][Iz][BGAM][0] = G_ACCG * blockFProps[Ix][Iy][Iz][RO][0] * compMWL;
			blockFProps[Ix][Iy][Iz][BGAM][1] = G_ACCG * blockFProps[Ix][Iy][Iz][RO][1] * compMWG;
		}
		else {
			//Fracture System
			ppy -= Nx * Ny*Nz;
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;



			compMWL = 0;
			compMWG = 0;
			for (n = 0; n < Nc; n++) {
				compMWL += fcomp[Ix][Iy][Iz][n][0] * fluidProp[n][MW];
				compMWG += fcomp[Ix][Iy][Iz][n][1] * fluidProp[n][MW];
			}
			fblockFProps[Ix][Iy][Iz][BGAM][0] = G_ACCG * fblockFProps[Ix][Iy][Iz][RO][0] * compMWL;
			fblockFProps[Ix][Iy][Iz][BGAM][1] = G_ACCG * fblockFProps[Ix][Iy][Iz][RO][1] * compMWG;
		}
	}
}

void MixedCalcDiffusion(void) {
	register int iix, ppy;
	int Ix, Iy, Iz, qxtemp;
	register int i, j, k;
	FType roRD, pDij;
	FType sSum1, sSum2, Tij, Sij, Oij, difTemp;

	for (iix = 0; iix < SchindlerListLength; iix++) {
		ppy = SchindlerList[iix];
		if (ppy < Nx*Ny*Nz) {
			//Matrix System
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;


			for (i = 0; i < 2; i++) {				
				sSum1 = 0;
				sSum2 = 0;
				for (j = 0; j < Nc; j++) {
					sSum1 += comp[Ix][Iy][Iz][j][i] * pow(fluidProp[j][VCRIT], (5.0 / 3.0));
					sSum2 += comp[Ix][Iy][Iz][j][i] * pow(fluidProp[j][VCRIT], (2.0 / 3.0));
				}
				sSum1 = blockFProps[Ix][Iy][Iz][RO][i] * sSum1 / sSum2;
				roRD = (0.99589 + 0.096016*sSum1 - 0.22035*sSum1*sSum1 + 0.032874*sSum1*sSum1*sSum1) / blockFProps[Ix][Iy][Iz][RO][i];

				for (j = 0; j < Nc; j++) {
					sSum1 = 0;
					for (k = 0; k < Nc; k++) {
						if (j == k) continue;
						Sij = (fluidProp[j][COLDIA] + fluidProp[k][COLDIA]) / 2;
						Tij = resTemp / sqrt(fluidProp[j][EPDIA] * fluidProp[k][EPDIA]);
						Oij = 1.06306*pow(Tij, -0.1561) + 0.193*exp(-0.47635*Tij) + 1.03587*exp(-1.52996*Tij) + 1.76474*exp(-3.89411*Tij);
						pDij = 188.2922475*sqrt(resTemp*(1.0 / fluidProp[j][MW] + 1.0 / fluidProp[k][MW])) / (RGAS*Oij*Sij*Sij);
						sSum1 += comp[Ix][Iy][Iz][k][i] / (pDij*roRD);
					}
					difTemp = (1 - comp[Ix][Iy][Iz][j][i]) / (sSum1 * 10000);
					if (comp[Ix][Iy][Iz][j][i] == 1) {
						diffusion[Ix][Iy][Iz][j][i] = roRD * pDij;
					}
					else if ((difTemp < 0) || (!sSum1)) diffusion[Ix][Iy][Iz][j][i] = 0;
					else diffusion[Ix][Iy][Iz][j][i] = DIFF_INTERNAL_COEFF*difTemp;
					if (IsConstDif) diffusion[Ix][Iy][Iz][j][i] = CDiffCoef;					
				}
			}
		}
		else {
			//Fracture System
			ppy -= Nx * Ny*Nz;
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;


			for (i = 0; i < 2; i++) {
				sSum1 = 0;
				sSum2 = 0;
				for (j = 0; j < Nc; j++) {
					sSum1 += fcomp[Ix][Iy][Iz][j][i] * pow(fluidProp[j][VCRIT], (5.0 / 3.0));
					sSum2 += fcomp[Ix][Iy][Iz][j][i] * pow(fluidProp[j][VCRIT], (2.0 / 3.0));
				}
				sSum1 = fblockFProps[Ix][Iy][Iz][RO][i] * sSum1 / sSum2;
				roRD = (0.99589 + 0.096016*sSum1 - 0.22035*sSum1*sSum1 + 0.032874*sSum1*sSum1*sSum1) / fblockFProps[Ix][Iy][Iz][RO][i];

				for (j = 0; j < Nc; j++) {
					sSum1 = 0;
					for (k = 0; k < Nc; k++) {
						if (j == k) continue;
						Sij = (fluidProp[j][COLDIA] + fluidProp[k][COLDIA]) / 2;
						Tij = resTemp / sqrt(fluidProp[j][EPDIA] * fluidProp[k][EPDIA]);
						Oij = 1.06306*pow(Tij, -0.1561) + 0.193*exp(-0.47635*Tij) + 1.03587*exp(-1.52996*Tij) + 1.76474*exp(-3.89411*Tij);
						pDij = 188.2922475*sqrt(resTemp*(1.0 / fluidProp[j][MW] + 1.0 / fluidProp[k][MW])) / (RGAS*Oij*Sij*Sij);
						sSum1 += fcomp[Ix][Iy][Iz][k][i] / (pDij*roRD);
					}
					difTemp = (1 - fcomp[Ix][Iy][Iz][j][i]) / (sSum1 * 10000);
					if (fcomp[Ix][Iy][Iz][j][i] == 1) {
						fdiffusion[Ix][Iy][Iz][j][i] = roRD * pDij;
					}
					else if ((difTemp < 0) || (!sSum1)) fdiffusion[Ix][Iy][Iz][j][i] = 0;
					else fdiffusion[Ix][Iy][Iz][j][i] = DIFF_INTERNAL_COEFF*difTemp;

					if (IsConstDif) fdiffusion[Ix][Iy][Iz][j][i] = CDiffCoef;						
				}
			}
		}
	}
}

void MixedCalcSatFs(void) {
	register int iix, ppy, n;
	int Ix, Iy, Iz, qxtemp;
	
	FType Krog, Krow, sdKrog, sdKrow;
	FType AA, BB, dAA, dBB;
	//FType satL, satG;
	FType Fift;


	for (iix = 0; iix < SchindlerListLength; iix++) {
		ppy = SchindlerList[iix];
		if (ppy < Nx*Ny*Nz) {
			//Matrix System
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;


			/*relPerm[Ix][Iy][Iz][0]=RELPERM0;
				dRelPerm[Ix][Iy][Iz][0]=0;

				if (sat[Ix][Iy][Iz][1]<CoreySor) {
				relPerm[Ix][Iy][Iz][1]=RELPERM0;
				dRelPerm[Ix][Iy][Iz][1]=0;
				dRelPerm[Ix][Iy][Iz][2]=0;
				dRelPerm[Ix][Iy][Iz][4]=0;
				}
				else if (sat[Ix][Iy][Iz][1]>=1) {
				relPerm[Ix][Iy][Iz][1]=1;
				dRelPerm[Ix][Iy][Iz][1]=0;
				dRelPerm[Ix][Iy][Iz][2]=0;
				dRelPerm[Ix][Iy][Iz][4]=0;
				}
				else {
				relPerm[Ix][Iy][Iz][1]=pow((sat[Ix][Iy][Iz][1]-CoreySor)/(1-CoreySor), 3.0);
				dRelPerm[Ix][Iy][Iz][1]=0;
				dRelPerm[Ix][Iy][Iz][2]=0;
				dRelPerm[Ix][Iy][Iz][4]=(3.0/(1-CoreySor))*pow((sat[Ix][Iy][Iz][1]-CoreySor)/(1-CoreySor), 2.0);
				}
				//fprintf(fp, "%d\t%d\t%d\to=%f\t%f\n", i, j, k, relPerm[Ix][Iy][Iz][1], dRelPerm[Ix][Iy][Iz][4]);

				if (sat[Ix][Iy][Iz][2]<=0) {
				relPerm[Ix][Iy][Iz][2]=RELPERM0;
				dRelPerm[Ix][Iy][Iz][3]=0;
				blockFProps[Ix][Iy][Iz][BLOCK_PC][1]=0;
				}
				else if (sat[Ix][Iy][Iz][2]>=(1-CoreySor)) {
				relPerm[Ix][Iy][Iz][2]=1;
				dRelPerm[Ix][Iy][Iz][3]=0;
				blockFProps[Ix][Iy][Iz][BLOCK_PC][1]=3000941.83744005;
				}
				else {
				relPerm[Ix][Iy][Iz][2]=pow(sat[Ix][Iy][Iz][2]/(1-CoreySor), 1.5);

				dRelPerm[Ix][Iy][Iz][3]=(1.5/(1-CoreySor))*pow(sat[Ix][Iy][Iz][2]/(1-CoreySor), 0.5);
				blockFProps[Ix][Iy][Iz][BLOCK_PC][1]=103400*pow((sat[Ix][Iy][Iz][1]-CoreySor)/(1-CoreySor), -1.0/2.7);
				}
				blockFProps[Ix][Iy][Iz][BLOCK_PC][0]=0;
				//fprintf(fp, "%d\t%d\t%d\tg=%f\t%f\t%f\n", i, j, k, relPerm[Ix][Iy][Iz][2], dRelPerm[Ix][Iy][Iz][3], blockFProps[Ix][Iy][Iz][BLOCK_PC][1]);


				*/

			dRelPerm[Ix][Iy][Iz][4] = 0;			//Oil to Oil
												//Relative permeability, Baker's model & Capillary pressures
			for (n = 0; n < Nswt; n++) {
				if (sat[Ix][Iy][Iz][0] < swt[n][SAT]) break;
			}
			if (n == 0) {
				relPerm[Ix][Iy][Iz][0] = RELPERM0;
				//blockFProps[Ix][Iy][Iz][BLOCK_PC][0] = swt[0][CAPILLARYPRESSURE];
				Krow = swt[0][KRO];

				dRelPerm[Ix][Iy][Iz][0] = 0;
				sdKrow = 0;
			}
			else if (n == Nswt) {
				relPerm[Ix][Iy][Iz][0] = swt[Nswt - 1][KR];
				//blockFProps[Ix][Iy][Iz][BLOCK_PC][0] = swt[Nswt - 1][CAPILLARYPRESSURE];
				Krow = swt[Nswt - 1][KRO];

				dRelPerm[Ix][Iy][Iz][0] = 0;
				sdKrow = 0;
			}
			else {
				relPerm[Ix][Iy][Iz][0] = (swt[n][KR] - swt[n - 1][KR]) / (swt[n][SAT] - swt[n - 1][SAT])*(sat[Ix][Iy][Iz][0] - swt[n][SAT]) + swt[n][KR];
				Krow = (swt[n][KRO] - swt[n - 1][KRO]) / (swt[n][SAT] - swt[n - 1][SAT])*(sat[Ix][Iy][Iz][0] - swt[n][SAT]) + swt[n][KRO];
				//blockFProps[Ix][Iy][Iz][BLOCK_PC][0] = (swt[n][CAPILLARYPRESSURE] - swt[n - 1][CAPILLARYPRESSURE]) / (swt[n][SAT] - swt[n - 1][SAT])*(sat[Ix][Iy][Iz][0] - swt[n][SAT]) + swt[n][CAPILLARYPRESSURE];

				dRelPerm[Ix][Iy][Iz][0] = (swt[n][KR] - swt[n - 1][KR]) / (swt[n][SAT] - swt[n - 1][SAT]);
				sdKrow = (swt[n][KRO] - swt[n - 1][KRO]) / (swt[n][SAT] - swt[n - 1][SAT]);
			}

			for (n = 0; n < Nsgt; n++) {
				if (sat[Ix][Iy][Iz][2] < sgt[n][SAT]) break;
			}
			if (n == 0) {
				relPerm[Ix][Iy][Iz][2] = RELPERM0;
				//blockFProps[Ix][Iy][Iz][BLOCK_PC][1] = sgt[0][CAPILLARYPRESSURE];
				Krog = sgt[0][KRO];

				dRelPerm[Ix][Iy][Iz][3] = 0;
				sdKrog = 0;
			}
			else if (n == Nsgt) {
				relPerm[Ix][Iy][Iz][2] = sgt[Nsgt - 1][KR];
				//blockFProps[Ix][Iy][Iz][BLOCK_PC][1] = sgt[Nsgt - 1][CAPILLARYPRESSURE];
				Krog = sgt[Nsgt - 1][KRO];

				dRelPerm[Ix][Iy][Iz][3] = 0;
				sdKrog = 0;
			}
			else {
				relPerm[Ix][Iy][Iz][2] = (sgt[n][KR] - sgt[n - 1][KR]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(sat[Ix][Iy][Iz][2] - sgt[n][SAT]) + sgt[n][KR];
				Krog = (sgt[n][KRO] - sgt[n - 1][KRO]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(sat[Ix][Iy][Iz][2] - sgt[n][SAT]) + sgt[n][KRO];
				//blockFProps[Ix][Iy][Iz][BLOCK_PC][1] = (sgt[n][CAPILLARYPRESSURE] - sgt[n - 1][CAPILLARYPRESSURE]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(sat[Ix][Iy][Iz][2] - sgt[n][SAT]) + sgt[n][CAPILLARYPRESSURE];

				dRelPerm[Ix][Iy][Iz][3] = (sgt[n][KR] - sgt[n - 1][KR]) / (sgt[n][SAT] - sgt[n - 1][SAT]);
				sdKrog = (sgt[n][KRO] - sgt[n - 1][KRO]) / (sgt[n][SAT] - sgt[n - 1][SAT]);
			}

			//Baker Model
			//satL=bsat[Ix][Iy][Iz][0]-swt[0][SAT];
			//satG=bsat[Ix][Iy][Iz][2]-sgt[0][SAT];

			//if (bsat[Ix][Iy][Iz][1]) relPerm[Ix][Iy][Iz][1]=(satL*Krow+satG*Krog)/(satL+satG);
			//else relPerm[Ix][Iy][Iz][1]=0;

			//Stone II model
			//relPerm[Ix][Iy][Iz][1]=swt[0][KRO]*((relPerm[Ix][Iy][Iz][0]+Krow/swt[0][KRO])*(relPerm[Ix][Iy][Iz][2]+Krog/swt[0][KRO])-relPerm[Ix][Iy][Iz][0]-relPerm[Ix][Iy][Iz][2]);

			AA = relPerm[Ix][Iy][Iz][0] + Krow / swt[0][KRO];
			BB = relPerm[Ix][Iy][Iz][2] + Krog / swt[0][KRO];
			dAA = (dRelPerm[Ix][Iy][Iz][0] + sdKrow / swt[0][KRO]);
			dBB = (dRelPerm[Ix][Iy][Iz][3] + sdKrog / swt[0][KRO]);

			relPerm[Ix][Iy][Iz][1] = swt[0][KRO] * (AA*BB - relPerm[Ix][Iy][Iz][0] - relPerm[Ix][Iy][Iz][2]);
			dRelPerm[Ix][Iy][Iz][1] = swt[0][KRO] * (dAA*BB - dRelPerm[Ix][Iy][Iz][0]);
			dRelPerm[Ix][Iy][Iz][2] = swt[0][KRO] * (AA*dBB - dRelPerm[Ix][Iy][Iz][3]);
			////////////////////////////////////////////////
			///////////////////////////////////////////////////


			if (Bift[Ix][Iy][Iz] < MISCIBLEIFT) {
				Fift = pow(Bift[Ix][Iy][Iz] / IFTBASECASE, 1.0 / IFTPOWER);
				relPerm[Ix][Iy][Iz][1] = Fift * relPerm[Ix][Iy][Iz][1] + (1 - Fift)*sat[Ix][Iy][Iz][1] * (1 - swt[0][SAT]);
				relPerm[Ix][Iy][Iz][2] = Fift * relPerm[Ix][Iy][Iz][2] + (1 - Fift)*sat[Ix][Iy][Iz][2] * (1 - swt[0][SAT]);;

				dRelPerm[Ix][Iy][Iz][3] = Fift * dRelPerm[Ix][Iy][Iz][3] + (1 - Fift)*(1 - swt[0][SAT]);
				dRelPerm[Ix][Iy][Iz][1] *= Fift;
				dRelPerm[Ix][Iy][Iz][2] *= Fift;
				dRelPerm[Ix][Iy][Iz][4] = (1 - Fift)*(1 - swt[0][SAT]);
			}

			if (relPerm[Ix][Iy][Iz][1] < 0) relPerm[Ix][Iy][Iz][1] = RELPERM0;
			if (relPerm[Ix][Iy][Iz][1] > 1) relPerm[Ix][Iy][Iz][1] = 1;
			if (dRelPerm[Ix][Iy][Iz][1] > 0) dRelPerm[Ix][Iy][Iz][1] = 0;
			if (dRelPerm[Ix][Iy][Iz][2] > 0) dRelPerm[Ix][Iy][Iz][2] = 0;
		}
		else {
			//Fracture System
			ppy -= Nx * Ny*Nz;
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;


			/*relPerm[Ix][Iy][Iz][0]=RELPERM0;
				dRelPerm[Ix][Iy][Iz][0]=0;

				if (sat[Ix][Iy][Iz][1]<CoreySor) {
				relPerm[Ix][Iy][Iz][1]=RELPERM0;
				dRelPerm[Ix][Iy][Iz][1]=0;
				dRelPerm[Ix][Iy][Iz][2]=0;
				dRelPerm[Ix][Iy][Iz][4]=0;
				}
				else if (sat[Ix][Iy][Iz][1]>=1) {
				relPerm[Ix][Iy][Iz][1]=1;
				dRelPerm[Ix][Iy][Iz][1]=0;
				dRelPerm[Ix][Iy][Iz][2]=0;
				dRelPerm[Ix][Iy][Iz][4]=0;
				}
				else {
				relPerm[Ix][Iy][Iz][1]=pow((sat[Ix][Iy][Iz][1]-CoreySor)/(1-CoreySor), 3.0);
				dRelPerm[Ix][Iy][Iz][1]=0;
				dRelPerm[Ix][Iy][Iz][2]=0;
				dRelPerm[Ix][Iy][Iz][4]=(3.0/(1-CoreySor))*pow((sat[Ix][Iy][Iz][1]-CoreySor)/(1-CoreySor), 2.0);
				}
				//fprintf(fp, "%d\t%d\t%d\to=%f\t%f\n", i, j, k, relPerm[Ix][Iy][Iz][1], dRelPerm[Ix][Iy][Iz][4]);

				if (sat[Ix][Iy][Iz][2]<=0) {
				relPerm[Ix][Iy][Iz][2]=RELPERM0;
				dRelPerm[Ix][Iy][Iz][3]=0;
				blockFProps[Ix][Iy][Iz][BLOCK_PC][1]=0;
				}
				else if (sat[Ix][Iy][Iz][2]>=(1-CoreySor)) {
				relPerm[Ix][Iy][Iz][2]=1;
				dRelPerm[Ix][Iy][Iz][3]=0;
				blockFProps[Ix][Iy][Iz][BLOCK_PC][1]=3000941.83744005;
				}
				else {
				relPerm[Ix][Iy][Iz][2]=pow(sat[Ix][Iy][Iz][2]/(1-CoreySor), 1.5);

				dRelPerm[Ix][Iy][Iz][3]=(1.5/(1-CoreySor))*pow(sat[Ix][Iy][Iz][2]/(1-CoreySor), 0.5);
				blockFProps[Ix][Iy][Iz][BLOCK_PC][1]=103400*pow((sat[Ix][Iy][Iz][1]-CoreySor)/(1-CoreySor), -1.0/2.7);
				}
				blockFProps[Ix][Iy][Iz][BLOCK_PC][0]=0;
				//fprintf(fp, "%d\t%d\t%d\tg=%f\t%f\t%f\n", i, j, k, relPerm[Ix][Iy][Iz][2], dRelPerm[Ix][Iy][Iz][3], blockFProps[Ix][Iy][Iz][BLOCK_PC][1]);




				fdRelPerm[Ix][Iy][Iz][4] = 1;			//Oil to Oil
				//Relative permeability, Baker's model & Capillary pressures
				for (n = 0; n<Nswt; n++) {
				if (fsat[Ix][Iy][Iz][0]<swt[n][SAT]) break;
				}
				if (n == 0) {
				frelPerm[Ix][Iy][Iz][0] = RELPERM0;
				fblockFProps[Ix][Iy][Iz][BLOCK_PC][0] = swt[0][CAPILLARYPRESSURE];
				Krow = swt[0][KRO];

				fdRelPerm[Ix][Iy][Iz][0] = 1;
				sdKrow = 0;
				}
				else if (n == Nswt) {
				frelPerm[Ix][Iy][Iz][0] = swt[Nswt - 1][KR];
				fblockFProps[Ix][Iy][Iz][BLOCK_PC][0] = swt[Nswt - 1][CAPILLARYPRESSURE];
				Krow = swt[Nswt - 1][KRO];

				fdRelPerm[Ix][Iy][Iz][0] = 1;
				sdKrow = 0;
				}
				else {
				frelPerm[Ix][Iy][Iz][0] = (swt[n][KR] - swt[n - 1][KR]) / (swt[n][SAT] - swt[n - 1][SAT])*(fsat[Ix][Iy][Iz][0] - swt[n][SAT]) + swt[n][KR];
				Krow = (swt[n][KRO] - swt[n - 1][KRO]) / (swt[n][SAT] - swt[n - 1][SAT])*(fsat[Ix][Iy][Iz][0] - swt[n][SAT]) + swt[n][KRO];
				fblockFProps[Ix][Iy][Iz][BLOCK_PC][0] = (swt[n][CAPILLARYPRESSURE] - swt[n - 1][CAPILLARYPRESSURE]) / (swt[n][SAT] - swt[n - 1][SAT])*(fsat[Ix][Iy][Iz][0] - swt[n][SAT]) + swt[n][CAPILLARYPRESSURE];

				fdRelPerm[Ix][Iy][Iz][0] = (swt[n][KR] - swt[n - 1][KR]) / (swt[n][SAT] - swt[n - 1][SAT]);
				sdKrow = (swt[n][KRO] - swt[n - 1][KRO]) / (swt[n][SAT] - swt[n - 1][SAT]);
				}

				for (n = 0; n<Nsgt; n++) {
				if (fsat[Ix][Iy][Iz][2]<sgt[n][SAT]) break;
				}
				if (n == 0) {
				frelPerm[Ix][Iy][Iz][2] = RELPERM0;
				fblockFProps[Ix][Iy][Iz][BLOCK_PC][1] = sgt[0][CAPILLARYPRESSURE];
				Krog = sgt[0][KRO];

				fdRelPerm[Ix][Iy][Iz][3] = 1;
				sdKrog = 0;
				}
				else if (n == Nsgt) {
				frelPerm[Ix][Iy][Iz][2] = sgt[Nsgt - 1][KR];
				fblockFProps[Ix][Iy][Iz][BLOCK_PC][1] = sgt[Nsgt - 1][CAPILLARYPRESSURE];
				Krog = sgt[Nsgt - 1][KRO];

				fdRelPerm[Ix][Iy][Iz][3] = 1;
				sdKrog = 0;
				}
				else {
				frelPerm[Ix][Iy][Iz][2] = (sgt[n][KR] - sgt[n - 1][KR]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(fsat[Ix][Iy][Iz][2] - sgt[n][SAT]) + sgt[n][KR];
				Krog = (sgt[n][KRO] - sgt[n - 1][KRO]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(fsat[Ix][Iy][Iz][2] - sgt[n][SAT]) + sgt[n][KRO];
				fblockFProps[Ix][Iy][Iz][BLOCK_PC][1] = (sgt[n][CAPILLARYPRESSURE] - sgt[n - 1][CAPILLARYPRESSURE]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(fsat[Ix][Iy][Iz][2] - sgt[n][SAT]) + sgt[n][CAPILLARYPRESSURE];

				fdRelPerm[Ix][Iy][Iz][3] = (sgt[n][KR] - sgt[n - 1][KR]) / (sgt[n][SAT] - sgt[n - 1][SAT]);
				sdKrog = (sgt[n][KRO] - sgt[n - 1][KRO]) / (sgt[n][SAT] - sgt[n - 1][SAT]);
				}

				//Baker Model
				//satL=bsat[Ix][Iy][Iz][0]-swt[0][SAT];
				//satG=bsat[Ix][Iy][Iz][2]-sgt[0][SAT];

				//if (bsat[Ix][Iy][Iz][1]) relPerm[Ix][Iy][Iz][1]=(satL*Krow+satG*Krog)/(satL+satG);
				//else relPerm[Ix][Iy][Iz][1]=0;

				//Stone II model
				//relPerm[Ix][Iy][Iz][1]=swt[0][KRO]*((relPerm[Ix][Iy][Iz][0]+Krow/swt[0][KRO])*(relPerm[Ix][Iy][Iz][2]+Krog/swt[0][KRO])-relPerm[Ix][Iy][Iz][0]-relPerm[Ix][Iy][Iz][2]);

				AA = frelPerm[Ix][Iy][Iz][0] + Krow / swt[0][KRO];
				BB = frelPerm[Ix][Iy][Iz][2] + Krog / swt[0][KRO];
				dAA = (fdRelPerm[Ix][Iy][Iz][0] + sdKrow / swt[0][KRO]);
				dBB = (fdRelPerm[Ix][Iy][Iz][3] + sdKrog / swt[0][KRO]);

				frelPerm[Ix][Iy][Iz][1] = swt[0][KRO] * (AA*BB - frelPerm[Ix][Iy][Iz][0] - frelPerm[Ix][Iy][Iz][2]);
				fdRelPerm[Ix][Iy][Iz][1] = swt[0][KRO] * (dAA*BB - fdRelPerm[Ix][Iy][Iz][0]);
				fdRelPerm[Ix][Iy][Iz][2] = swt[0][KRO] * (AA*dBB - fdRelPerm[Ix][Iy][Iz][3]);
				////////////////////////////////////////////////
				///////////////////////////////////////////////////


				if (fBift[Ix][Iy][Iz]<MISCIBLEIFT) {
				Fift = pow(Bift[Ix][Iy][Iz] / IFTBASECASE, 1.0 / IFTPOWER);
				frelPerm[Ix][Iy][Iz][1] = Fift*frelPerm[Ix][Iy][Iz][1] + (1 - Fift)*fsat[Ix][Iy][Iz][1] * (1 - swt[0][SAT]);
				frelPerm[Ix][Iy][Iz][2] = Fift*frelPerm[Ix][Iy][Iz][2] + (1 - Fift)*fsat[Ix][Iy][Iz][2] * (1 - swt[0][SAT]);;

				fdRelPerm[Ix][Iy][Iz][3] = Fift*fdRelPerm[Ix][Iy][Iz][3] + (1 - Fift)*(1 - swt[0][SAT]);
				fdRelPerm[Ix][Iy][Iz][1] *= Fift;
				fdRelPerm[Ix][Iy][Iz][2] *= Fift;
				fdRelPerm[Ix][Iy][Iz][4] = (1 - Fift)*(1 - swt[0][SAT]);
				}

				if (frelPerm[Ix][Iy][Iz][1]<0) frelPerm[Ix][Iy][Iz][1] = RELPERM0;
				if (frelPerm[Ix][Iy][Iz][1]>1) frelPerm[Ix][Iy][Iz][1] = 1;
				if (fdRelPerm[Ix][Iy][Iz][1]>0) fdRelPerm[Ix][Iy][Iz][1] = 1;
				if (fdRelPerm[Ix][Iy][Iz][2]>0) fdRelPerm[Ix][Iy][Iz][2] = -1;*/

			frelPerm[Ix][Iy][Iz][0] = fsat[Ix][Iy][Iz][0];//-swt[0][0];
			frelPerm[Ix][Iy][Iz][1] = fsat[Ix][Iy][Iz][1];
			frelPerm[Ix][Iy][Iz][2] = fsat[Ix][Iy][Iz][2];

			fdRelPerm[Ix][Iy][Iz][0] = 1;
			fdRelPerm[Ix][Iy][Iz][1] = 0;
			fdRelPerm[Ix][Iy][Iz][2] = 0;
			fdRelPerm[Ix][Iy][Iz][3] = 1;
			fdRelPerm[Ix][Iy][Iz][4] = 1;


			//fblockFProps[Ix][Iy][Iz][BLOCK_PC][0] = 0;
			//fblockFProps[Ix][Iy][Iz][BLOCK_PC][1] = 0;

		}

		for (Iz = 0; Iz < Nz; Iz++)
			for (Iy = 0; Iy < Ny; Iy++)
				for (Ix = 0; Ix < Nx; Ix++) {
					fblockFProps[Ix][Iy][Iz][BLOCK_PC][0] = 0;
					fblockFProps[Ix][Iy][Iz][BLOCK_PC][1] = 0;




					
					for (n = 0; n < Nswt; n++) {
						if (sat[Ix][Iy][Iz][0] < swt[n][SAT]) break;
					}
					if (n == 0) {						
						blockFProps[Ix][Iy][Iz][BLOCK_PC][0] = swt[0][CAPILLARYPRESSURE];						
					}
					else if (n == Nswt) {						
						blockFProps[Ix][Iy][Iz][BLOCK_PC][0] = swt[Nswt - 1][CAPILLARYPRESSURE];						
					}
					else {						
						blockFProps[Ix][Iy][Iz][BLOCK_PC][0] = (swt[n][CAPILLARYPRESSURE] - swt[n - 1][CAPILLARYPRESSURE]) / (swt[n][SAT] - swt[n - 1][SAT])*(sat[Ix][Iy][Iz][0] - swt[n][SAT]) + swt[n][CAPILLARYPRESSURE];
					}

					for (n = 0; n < Nsgt; n++) {
						if (sat[Ix][Iy][Iz][2] < sgt[n][SAT]) break;
					}
					if (n == 0) {						
						blockFProps[Ix][Iy][Iz][BLOCK_PC][1] = sgt[0][CAPILLARYPRESSURE];						
					}
					else if (n == Nsgt) {						
						blockFProps[Ix][Iy][Iz][BLOCK_PC][1] = sgt[Nsgt - 1][CAPILLARYPRESSURE];						
					}
					else {						
						blockFProps[Ix][Iy][Iz][BLOCK_PC][1] = (sgt[n][CAPILLARYPRESSURE] - sgt[n - 1][CAPILLARYPRESSURE]) / (sgt[n][SAT] - sgt[n - 1][SAT])*(sat[Ix][Iy][Iz][2] - sgt[n][SAT]) + sgt[n][CAPILLARYPRESSURE];
					}
				}

	}
}

void MixedCalcP(void) {
	register int iix, ppy;
	int Ix, Iy, Iz, qxtemp;
	

	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {


			P[Ix + 1][Iy + 1][Iz + 1][0] = P[Ix + 1][Iy + 1][Iz + 1][1] - blockFProps[Ix][Iy][Iz][BLOCK_PC][0];
			if (Bift[Ix][Iy][Iz] < MISCIBLEIFT) P[Ix + 1][Iy + 1][Iz + 1][2] = P[Ix + 1][Iy + 1][Iz + 1][1] + blockFProps[Ix][Iy][Iz][BLOCK_PC][1] * Bift[Ix][Iy][Iz] / IFTBASECASE;
			else P[Ix + 1][Iy + 1][Iz + 1][2] = P[Ix + 1][Iy + 1][Iz + 1][1] + blockFProps[Ix][Iy][Iz][BLOCK_PC][1];
		


			//fP[Ix + 1][Iy + 1][Iz + 1][0] = fP[Ix + 1][Iy + 1][Iz + 1][1] - fblockFProps[Ix][Iy][Iz][BLOCK_PC][0];
			fP[Ix + 1][Iy + 1][Iz + 1][0] = fP[Ix + 1][Iy + 1][Iz + 1][1];
			//if (fBift[Ix][Iy][Iz] < MISCIBLEIFT) fP[Ix + 1][Iy + 1][Iz + 1][2] = fP[Ix + 1][Iy + 1][Iz + 1][1] + fblockFProps[Ix][Iy][Iz][BLOCK_PC][1];
			//else fP[Ix + 1][Iy + 1][Iz + 1][2] = fP[Ix + 1][Iy + 1][Iz + 1][1] ;
			//fP[Ix + 1][Iy + 1][Iz + 1][2] = fP[Ix + 1][Iy + 1][Iz + 1][1] + fblockFProps[Ix][Iy][Iz][BLOCK_PC][1];
			fP[Ix + 1][Iy + 1][Iz + 1][2] = fP[Ix + 1][Iy + 1][Iz + 1][1];
		}



	
}

void MixedCalcVisco(void) {
	register int iix, ppy, n;
	int Ix, Iy, Iz, qxtemp;
	FType s1Ml, s2Ml, s1Mg, s2Mg, tMl, tMg, MuL, MuG, roL, roG, compTCL, compMWL, compPCL, compTCG, compMWG, compPCG;
	FType ethaL, ethaG, eAAl, eAAg, eBBl, eBBg, eCCl, eCCg;
	FType ErL, ErG, ro1L, ro1G;

	FType Ul, Wl, Al, Bl, Zl, Ug, Wg, Ag, Bg, Zg;
	FType Cal, Cbl, Cag, Cbg;
	FType dAdxl, dBdxl, dZdxl, dUdxl, dWdxl, dAdxg, dBdxg, dZdxg, dUdxg, dWdxg;
	FType dAdpl, dBdpl, dUdpl, dWdpl, dAdpg, dBdpg, dUdpg, dWdpg, dZdpl, dZdpg;

	register int i, j, k;
	FType tempSQR;

	for (iix = 0; iix < SchindlerListLength; iix++) {
		ppy = SchindlerList[iix];
		if (ppy < Nx*Ny*Nz) {
			//Matrix System
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;


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
		else {
			//Fracture System
			ppy -= Nx * Ny*Nz;
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;


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
}

void MixedCalcTranses(void) {
	register int iix, ppy;
	int Ix, Iy, Iz, qxtemp;
	register int i, j, k;
	FType mGam, pot;
	char Coef1, Coef2;

	for (iix = 0; iix < SchindlerListLength; iix++) {
		ppy = SchindlerList[iix];
		if (ppy < Nx*Ny*Nz) {
			//Matrix System
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;

			i = Ix + 1;
			j = Iy + 1;
			k = Iz + 1;


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
		else {
			//Fracture System
			ppy -= Nx * Ny*Nz;
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;

			i = Ix + 1;
			j = Iy + 1;
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

void MixedPreMul(void) {
	register int iix, ppy, n;
	int Ix, Iy, Iz, qxtemp;
	FType porpor;

	for (iix = 0; iix < SchindlerListLength; iix++) {
		ppy = SchindlerList[iix];
		if (ppy < Nx*Ny*Nz) {
			//Matrix System
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;

			for (n = 0; n < (Nc - 1); n++) {
				porpor = porosity[Ix][Iy][Iz] * (1 + cpor * (P[Ix + 1][Iy + 1][Iz + 1][1] - refP) + dcpor * (P[Ix + 1][Iy + 1][Iz + 1][1] - refP)*(P[Ix + 1][Iy + 1][Iz + 1][1] - refP));
				preProp[Ix][Iy][Iz][n] = porpor * (comp[Ix][Iy][Iz][n][0] * blockFProps[Ix][Iy][Iz][RO][0] * sat[Ix][Iy][Iz][1] + comp[Ix][Iy][Iz][n][1] * blockFProps[Ix][Iy][Iz][RO][1] * sat[Ix][Iy][Iz][2]);
			}
			preProp[Ix][Iy][Iz][Nc - 1] = porpor * (blockFProps[Ix][Iy][Iz][RO][0] * sat[Ix][Iy][Iz][1] + blockFProps[Ix][Iy][Iz][RO][1] * sat[Ix][Iy][Iz][2]);
			preProp[Ix][Iy][Iz][Nc] = porpor * watRo*WAT_M_RO*sat[Ix][Iy][Iz][0];
		}
		else {
			//Fracture System
			ppy -= Nx * Ny*Nz;
			Iz = ppy / (Nx*Ny);
			qxtemp = ppy % (Nx*Ny);
			Iy = qxtemp / Nx;
			Ix = qxtemp % Nx;


			for (n = 0; n < (Nc - 1); n++) {
				porpor = fporosity[Ix][Iy][Iz] * (1 + cpor * (fP[Ix + 1][Iy + 1][Iz + 1][1] - refP) + dcpor * (fP[Ix + 1][Iy + 1][Iz + 1][1] - refP)*(fP[Ix + 1][Iy + 1][Iz + 1][1] - refP));
				fpreProp[Ix][Iy][Iz][n] = porpor * (fcomp[Ix][Iy][Iz][n][0] * fblockFProps[Ix][Iy][Iz][RO][0] * fsat[Ix][Iy][Iz][1] + fcomp[Ix][Iy][Iz][n][1] * fblockFProps[Ix][Iy][Iz][RO][1] * fsat[Ix][Iy][Iz][2]);
			}
			fpreProp[Ix][Iy][Iz][Nc - 1] = porpor * (fblockFProps[Ix][Iy][Iz][RO][0] * fsat[Ix][Iy][Iz][1] + fblockFProps[Ix][Iy][Iz][RO][1] * fsat[Ix][Iy][Iz][2]);
			fpreProp[Ix][Iy][Iz][Nc] = porpor * watRo*WAT_M_RO*fsat[Ix][Iy][Iz][0];
		}
	}
}