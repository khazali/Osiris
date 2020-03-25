#ifndef MASS_TRAN_H
#define MASS_TRAN_H

//#include "Globals.h"
#include <math.h>
#include "Functions.h"

extern FType *****diffusion;
extern FType *****comp;
extern FType **fluidProp;
extern int Nx, Ny, Nz, Nc;
extern FType *****blockFProps;
extern FType resTemp;
extern FType *****dE;
extern FType ****P;

//FType DifDiffusion(int, int, int, int, int, int);
//FType DiffVisco(int, int, int, int, int);


void CalcDiffusion(void) {
	register int Ix, Iy, Iz, i, j, k;
	FType roRD, pDij;
	FType sSum1, sSum2, Tij, Sij, Oij, difTemp;

	for (Iz = 0; Iz<Nz; Iz++) {
		for (Iy = 0; Iy<Ny; Iy++) {
			for (Ix = 0; Ix<Nx; Ix++) {
				for (i = 0; i<2; i++) {
					sSum1 = 0;
					sSum2 = 0;
					for (j = 0; j<Nc; j++) {
						sSum1 += comp[Ix][Iy][Iz][j][i] * pow(fluidProp[j][VCRIT], (5.0 / 3.0));
						sSum2 += comp[Ix][Iy][Iz][j][i] * pow(fluidProp[j][VCRIT], (2.0 / 3.0));
					}
					sSum1 = blockFProps[Ix][Iy][Iz][RO][i] * sSum1 / sSum2;
					roRD = (0.99589 + 0.096016*sSum1 - 0.22035*sSum1*sSum1 + 0.032874*sSum1*sSum1*sSum1) / blockFProps[Ix][Iy][Iz][RO][i];

					for (j = 0; j<Nc; j++) {
						sSum1 = 0;
						for (k = 0; k<Nc; k++) {
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
						else if ((difTemp<0) || (!sSum1)) diffusion[Ix][Iy][Iz][j][i] = 0;
						else diffusion[Ix][Iy][Iz][j][i]= DIFF_INTERNAL_COEFF*difTemp;
						if (IsConstDif) diffusion[Ix][Iy][Iz][j][i] = CDiffCoef;
					}
				}

			}
		}
	}
}

void fCalcDiffusion(void) {
	register int Ix, Iy, Iz, i, j, k;
	FType roRD, pDij;
	FType sSum1, sSum2, Tij, Sij, Oij, difTemp;

	for (Iz = 0; Iz < Nz; Iz++) {
		for (Iy = 0; Iy < Ny; Iy++) {
			for (Ix = 0; Ix < Nx; Ix++) {
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
						else fdiffusion[Ix][Iy][Iz][j][i]= DIFF_INTERNAL_COEFF*difTemp;
						if (IsConstDif) fdiffusion[Ix][Iy][Iz][j][i] = CDiffCoef;
					}
				}

			}
		}
	}
}

FType DifDiffusion(int Ix, int Iy, int Iz, int Component, int dComponent, int Phase) {
	register int j, k;
	FType sSum1, sSum2, Erl;
	FType dErdx, AA, dAAdx;
	FType dErdp, dAAdp;
	FType roRD, pDij;
	FType Tij, Sij, Oij;
	FType dDiffdx, dDiff_1dx, finDif;
	FType dDiffdp, dDiff_1dp, Dij;

	if (dComponent<Nc) {
		sSum1 = 0;
		sSum2 = 0;
		for (k = 0; k<Nc; k++) {
			sSum1 += comp[Ix][Iy][Iz][k][Phase] * pow(fluidProp[k][VCRIT], (5.0 / 3.0));
			sSum2 += comp[Ix][Iy][Iz][k][Phase] * pow(fluidProp[k][VCRIT], (2.0 / 3.0));
		}
		Erl = blockFProps[Ix][Iy][Iz][RO][Phase] * sSum1 / sSum2;
		dErdx = dE[Ix][Iy][Iz][dComponent][Phase] * sSum1 / sSum2 + blockFProps[Ix][Iy][Iz][RO][Phase] * (pow(fluidProp[dComponent][VCRIT], (5.0 / 3.0))*sSum2 - pow(fluidProp[dComponent][VCRIT], (2.0 / 3.0))*sSum1) / (sSum2*sSum2);
		AA = 0.99589 + 0.096016*Erl - 0.22035*Erl*Erl + 0.032874*Erl*Erl*Erl;
		dAAdx = 0.096016*dErdx - 0.4407*Erl*dErdx + 0.098622*Erl*Erl*dErdx;
		roRD = AA / blockFProps[Ix][Iy][Iz][RO][Phase];

		if (Component == dComponent) {
			sSum1 = 0;
			sSum2 = 0;
			for (j = 0; j<Nc; j++) {
				if (j == Component) continue;
				Sij = (fluidProp[j][COLDIA] + fluidProp[Component][COLDIA]) / 2;
				Tij = resTemp / sqrt(fluidProp[j][EPDIA] * fluidProp[Component][EPDIA]);
				Oij = 1.06306*pow(Tij, -0.1561) + 0.193*exp(-0.47635*Tij) + 1.03587*exp(-1.52996*Tij) + 1.76474*exp(-3.89411*Tij);
				pDij = 188.2922475*sqrt(resTemp*(1.0 / fluidProp[j][MW] + 1.0 / fluidProp[Component][MW])) / (RGAS*Oij*Sij*Sij);
				Dij = pDij*roRD;

				dDiffdx = pDij*(dAAdx*blockFProps[Ix][Iy][Iz][RO][Phase] - dE[Ix][Iy][Iz][dComponent][Phase] * AA) / (blockFProps[Ix][Iy][Iz][RO][Phase] * blockFProps[Ix][Iy][Iz][RO][Phase]);
				dDiff_1dx = -dDiffdx / (Dij*Dij);
				sSum1 += comp[Ix][Iy][Iz][j][Phase] * dDiff_1dx;
				sSum2 += comp[Ix][Iy][Iz][j][Phase] / Dij;
			}
			if ((AA<0) || (!sSum2)) finDif = 0;
			else finDif = -(sSum2 + sSum1*(1 - comp[Ix][Iy][Iz][Component][Phase])) / (sSum2*sSum2);
		}
		else {
			sSum1 = 0;
			sSum2 = 0;
			for (j = 0; j<Nc; j++) {
				if (j == Component) continue;
				Sij = (fluidProp[j][COLDIA] + fluidProp[Component][COLDIA]) / 2;
				Tij = resTemp / sqrt(fluidProp[j][EPDIA] * fluidProp[Component][EPDIA]);
				Oij = 1.06306*pow(Tij, -0.1561) + 0.193*exp(-0.47635*Tij) + 1.03587*exp(-1.52996*Tij) + 1.76474*exp(-3.89411*Tij);
				pDij = 188.2922475*sqrt(resTemp*(1.0 / fluidProp[j][MW] + 1.0 / fluidProp[Component][MW])) / (RGAS*Oij*Sij*Sij);
				Dij = pDij*roRD;

				if (j == dComponent) sSum1 += 1.0 / (pDij*roRD);

				dDiffdx = pDij*(dAAdx*blockFProps[Ix][Iy][Iz][RO][Phase] - dE[Ix][Iy][Iz][dComponent][Phase] * AA) / (blockFProps[Ix][Iy][Iz][RO][Phase] * blockFProps[Ix][Iy][Iz][RO][Phase]);
				dDiff_1dx = -dDiffdx / (Dij*Dij);
				sSum1 += comp[Ix][Iy][Iz][j][Phase] * dDiff_1dx;
				sSum2 += comp[Ix][Iy][Iz][j][Phase] / Dij;
			}
			if ((AA<0) || (!sSum2)) finDif = 0;
			else finDif = -sSum1*(1 - comp[Ix][Iy][Iz][Component][Phase]) / (sSum2*sSum2);
		}
	}
	else {
		sSum1 = 0;
		sSum2 = 0;
		for (k = 0; k<Nc; k++) {
			sSum1 += comp[Ix][Iy][Iz][k][Phase] * pow(fluidProp[k][VCRIT], (5.0 / 3.0));
			sSum2 += comp[Ix][Iy][Iz][k][Phase] * pow(fluidProp[k][VCRIT], (2.0 / 3.0));
		}
		Erl = blockFProps[Ix][Iy][Iz][RO][Phase] * sSum1 / sSum2;
		dErdp = dE[Ix][Iy][Iz][dComponent][Phase] * sSum1 / sSum2;
		AA = 0.99589 + 0.096016*Erl - 0.22035*Erl*Erl + 0.032874*Erl*Erl*Erl;
		dAAdp = 0.096016*dErdp - 0.4407*Erl*dErdp + 0.098622*Erl*Erl*dErdp;
		roRD = AA / blockFProps[Ix][Iy][Iz][RO][Phase];

		sSum1 = 0;
		sSum2 = 0;
		for (j = 0; j<Nc; j++) {
			if (j == Component) continue;
			Sij = (fluidProp[j][COLDIA] + fluidProp[Component][COLDIA]) / 2;
			Tij = resTemp / sqrt(fluidProp[j][EPDIA] * fluidProp[Component][EPDIA]);
			Oij = 1.06306*pow(Tij, -0.1561) + 0.193*exp(-0.47635*Tij) + 1.03587*exp(-1.52996*Tij) + 1.76474*exp(-3.89411*Tij);
			pDij = 188.2922475*sqrt(resTemp*(1.0 / fluidProp[j][MW] + 1.0 / fluidProp[Component][MW])) / (RGAS*Oij*Sij*Sij);
			Dij = pDij*roRD;

			dDiffdp = pDij*(dAAdp*blockFProps[Ix][Iy][Iz][RO][Phase] - dE[Ix][Iy][Iz][dComponent][Phase] * AA) / (blockFProps[Ix][Iy][Iz][RO][Phase] * blockFProps[Ix][Iy][Iz][RO][Phase]);
			dDiff_1dp = -dDiffdp / (Dij*Dij);
			sSum1 += comp[Ix][Iy][Iz][j][Phase] * dDiff_1dp;
			sSum2 += comp[Ix][Iy][Iz][j][Phase] / Dij;
		}
		if ((AA<0) || (!sSum2)) finDif = 0;
		else finDif = -sSum1*(1 - comp[Ix][Iy][Iz][Component][Phase]) / (sSum2*sSum2);
	}

	if (IsConstDif) return 0;
	else return DIFF_INTERNAL_COEFF*finDif;
}

FType fDifDiffusion(int Ix, int Iy, int Iz, int Component, int dComponent, int Phase) {
	register int j, k;
	FType sSum1, sSum2, Erl;
	FType dErdx, AA, dAAdx;
	FType dErdp, dAAdp;
	FType roRD, pDij;
	FType Tij, Sij, Oij;
	FType dDiffdx, dDiff_1dx, finDif;
	FType dDiffdp, dDiff_1dp, Dij;

	if (dComponent < Nc) {
		sSum1 = 0;
		sSum2 = 0;
		for (k = 0; k < Nc; k++) {
			sSum1 += fcomp[Ix][Iy][Iz][k][Phase] * pow(fluidProp[k][VCRIT], (5.0 / 3.0));
			sSum2 += fcomp[Ix][Iy][Iz][k][Phase] * pow(fluidProp[k][VCRIT], (2.0 / 3.0));
		}
		Erl = fblockFProps[Ix][Iy][Iz][RO][Phase] * sSum1 / sSum2;
		dErdx = fdE[Ix][Iy][Iz][dComponent][Phase] * sSum1 / sSum2 + fblockFProps[Ix][Iy][Iz][RO][Phase] * (pow(fluidProp[dComponent][VCRIT], (5.0 / 3.0))*sSum2 - pow(fluidProp[dComponent][VCRIT], (2.0 / 3.0))*sSum1) / (sSum2*sSum2);
		AA = 0.99589 + 0.096016*Erl - 0.22035*Erl*Erl + 0.032874*Erl*Erl*Erl;
		dAAdx = 0.096016*dErdx - 0.4407*Erl*dErdx + 0.098622*Erl*Erl*dErdx;
		roRD = AA / fblockFProps[Ix][Iy][Iz][RO][Phase];

		if (Component == dComponent) {
			sSum1 = 0;
			sSum2 = 0;
			for (j = 0; j < Nc; j++) {
				if (j == Component) continue;
				Sij = (fluidProp[j][COLDIA] + fluidProp[Component][COLDIA]) / 2;
				Tij = resTemp / sqrt(fluidProp[j][EPDIA] * fluidProp[Component][EPDIA]);
				Oij = 1.06306*pow(Tij, -0.1561) + 0.193*exp(-0.47635*Tij) + 1.03587*exp(-1.52996*Tij) + 1.76474*exp(-3.89411*Tij);
				pDij = 188.2922475*sqrt(resTemp*(1.0 / fluidProp[j][MW] + 1.0 / fluidProp[Component][MW])) / (RGAS*Oij*Sij*Sij);
				Dij = pDij*roRD;

				dDiffdx = pDij*(dAAdx*fblockFProps[Ix][Iy][Iz][RO][Phase] - fdE[Ix][Iy][Iz][dComponent][Phase] * AA) / (fblockFProps[Ix][Iy][Iz][RO][Phase] * fblockFProps[Ix][Iy][Iz][RO][Phase]);
				dDiff_1dx = -dDiffdx / (Dij*Dij);
				sSum1 += fcomp[Ix][Iy][Iz][j][Phase] * dDiff_1dx;
				sSum2 += fcomp[Ix][Iy][Iz][j][Phase] / Dij;
			}
			if ((AA < 0) || (!sSum2)) finDif = 0;
			else finDif = -(sSum2 + sSum1*(1 - fcomp[Ix][Iy][Iz][Component][Phase])) / (sSum2*sSum2);
		}
		else {
			sSum1 = 0;
			sSum2 = 0;
			for (j = 0; j < Nc; j++) {
				if (j == Component) continue;
				Sij = (fluidProp[j][COLDIA] + fluidProp[Component][COLDIA]) / 2;
				Tij = resTemp / sqrt(fluidProp[j][EPDIA] * fluidProp[Component][EPDIA]);
				Oij = 1.06306*pow(Tij, -0.1561) + 0.193*exp(-0.47635*Tij) + 1.03587*exp(-1.52996*Tij) + 1.76474*exp(-3.89411*Tij);
				pDij = 188.2922475*sqrt(resTemp*(1.0 / fluidProp[j][MW] + 1.0 / fluidProp[Component][MW])) / (RGAS*Oij*Sij*Sij);
				Dij = pDij*roRD;

				if (j == dComponent) sSum1 += 1.0 / (pDij*roRD);

				dDiffdx = pDij*(dAAdx*fblockFProps[Ix][Iy][Iz][RO][Phase] - fdE[Ix][Iy][Iz][dComponent][Phase] * AA) / (fblockFProps[Ix][Iy][Iz][RO][Phase] * fblockFProps[Ix][Iy][Iz][RO][Phase]);
				dDiff_1dx = -dDiffdx / (Dij*Dij);
				sSum1 += fcomp[Ix][Iy][Iz][j][Phase] * dDiff_1dx;
				sSum2 += fcomp[Ix][Iy][Iz][j][Phase] / Dij;
			}
			if ((AA < 0) || (!sSum2)) finDif = 0;
			else finDif = -sSum1*(1 - fcomp[Ix][Iy][Iz][Component][Phase]) / (sSum2*sSum2);
		}
	}
	else {
		sSum1 = 0;
		sSum2 = 0;
		for (k = 0; k < Nc; k++) {
			sSum1 += fcomp[Ix][Iy][Iz][k][Phase] * pow(fluidProp[k][VCRIT], (5.0 / 3.0));
			sSum2 += fcomp[Ix][Iy][Iz][k][Phase] * pow(fluidProp[k][VCRIT], (2.0 / 3.0));
		}
		Erl = fblockFProps[Ix][Iy][Iz][RO][Phase] * sSum1 / sSum2;
		dErdp = fdE[Ix][Iy][Iz][dComponent][Phase] * sSum1 / sSum2;
		AA = 0.99589 + 0.096016*Erl - 0.22035*Erl*Erl + 0.032874*Erl*Erl*Erl;
		dAAdp = 0.096016*dErdp - 0.4407*Erl*dErdp + 0.098622*Erl*Erl*dErdp;
		roRD = AA / fblockFProps[Ix][Iy][Iz][RO][Phase];

		sSum1 = 0;
		sSum2 = 0;
		for (j = 0; j < Nc; j++) {
			if (j == Component) continue;
			Sij = (fluidProp[j][COLDIA] + fluidProp[Component][COLDIA]) / 2;
			Tij = resTemp / sqrt(fluidProp[j][EPDIA] * fluidProp[Component][EPDIA]);
			Oij = 1.06306*pow(Tij, -0.1561) + 0.193*exp(-0.47635*Tij) + 1.03587*exp(-1.52996*Tij) + 1.76474*exp(-3.89411*Tij);
			pDij = 188.2922475*sqrt(resTemp*(1.0 / fluidProp[j][MW] + 1.0 / fluidProp[Component][MW])) / (RGAS*Oij*Sij*Sij);
			Dij = pDij*roRD;

			dDiffdp = pDij*(dAAdp*fblockFProps[Ix][Iy][Iz][RO][Phase] - fdE[Ix][Iy][Iz][dComponent][Phase] * AA) / (fblockFProps[Ix][Iy][Iz][RO][Phase] * fblockFProps[Ix][Iy][Iz][RO][Phase]);
			dDiff_1dp = -dDiffdp / (Dij*Dij);
			sSum1 += fcomp[Ix][Iy][Iz][j][Phase] * dDiff_1dp;
			sSum2 += fcomp[Ix][Iy][Iz][j][Phase] / Dij;
		}
		if ((AA < 0) || (!sSum2)) finDif = 0;
		else finDif = -sSum1*(1 - fcomp[Ix][Iy][Iz][Component][Phase]) / (sSum2*sSum2);
	}

	if (IsConstDif) return 0;
	else return DIFF_INTERNAL_COEFF*finDif;
}

FType DiffVisco(int Ix, int Iy, int Iz, int Component, int Phase) {
	FType s1Ml, s2Ml, tMl, MuL, roL, compTCL, compMWL, compPCL;
	FType ethaL, eAAl, eBBl, eCCl;
	FType ErL, ro1L;
	register int n;
	FType dMuSdx, dEthadx, dErdx, dMMdx, dMuFdx;

	s1Ml = 0;
	s2Ml = 0;
	roL = 0;
	compMWL = 0;
	compPCL = 0;
	compTCL = 0;
	for (n = 0; n<Nc; n++) {
		tMl = comp[Ix][Iy][Iz][n][Phase] * sqrt(fluidProp[n][MW]);
		s1Ml += tMl*fluidProp[n][MU];
		s2Ml += tMl;
		roL += comp[Ix][Iy][Iz][n][Phase] * fluidProp[n][VCRIT];
		compTCL += comp[Ix][Iy][Iz][n][Phase] * fluidProp[n][TCRIT];
		compPCL += comp[Ix][Iy][Iz][n][Phase] * fluidProp[n][PCRIT] / 101325;
		compMWL += comp[Ix][Iy][Iz][n][Phase] * fluidProp[n][MW];
	}
	MuL = s1Ml / s2Ml;
	ErL = blockFProps[Ix][Iy][Iz][RO][Phase] * roL;
	eAAl = pow(compTCL, 1.0 / 6);
	eBBl = sqrt(compMWL);
	eCCl = pow(compPCL, 2.0 / 3);
	ethaL = eAAl / (eBBl*eCCl);
	ro1L = 0.10230 + 0.023364*ErL + 0.058533*ErL*ErL - 0.040758*ErL*ErL*ErL + 0.0093324*ErL*ErL*ErL*ErL;

	if (Component<Nc) {
		dMuSdx = (fluidProp[Component][MU] * sqrt(fluidProp[Component][MW])*s2Ml - sqrt(fluidProp[Component][MW])*s1Ml) / (s2Ml*s2Ml);
		dEthadx = ((1.0 / 6.0)*pow(compTCL, -5.0 / 6.0)*fluidProp[Component][TCRIT] * eBBl*eCCl - eAAl*(0.5*fluidProp[Component][MW] * eCCl / sqrt(compMWL) + (2.0 / 3.0)*pow(compPCL, -1.0 / 3.0)*fluidProp[Component][PCRIT] * eBBl / 101325)) / (eBBl*eBBl*eCCl*eCCl);
		dErdx = dE[Ix][Iy][Iz][Component][Phase] * roL + blockFProps[Ix][Iy][Iz][RO][Phase] * fluidProp[Component][VCRIT];
	}
	else {
		dMuSdx = 0;
		dEthadx = 0;
		dErdx = dE[Ix][Iy][Iz][Component][Phase] * roL;
	}

	dMMdx = 4 * ro1L*ro1L*ro1L*(0.023364*dErdx + 0.117066*dErdx*ErL - 0.122274*dErdx*ErL*ErL + 0.0373296*dErdx*ErL*ErL*ErL);
	dMuFdx = (dMMdx*ethaL - (ro1L*ro1L*ro1L*ro1L - 1e-4)*dEthadx) / (ethaL*ethaL) + dMuSdx;

	return dMuFdx;
}

FType fDiffVisco(int Ix, int Iy, int Iz, int Component, int Phase) {
	FType s1Ml, s2Ml, tMl, MuL, roL, compTCL, compMWL, compPCL;
	FType ethaL, eAAl, eBBl, eCCl;
	FType ErL, ro1L;
	register int n;
	FType dMuSdx, dEthadx, dErdx, dMMdx, dMuFdx;

	s1Ml = 0;
	s2Ml = 0;
	roL = 0;
	compMWL = 0;
	compPCL = 0;
	compTCL = 0;
	for (n = 0; n < Nc; n++) {
		tMl = fcomp[Ix][Iy][Iz][n][Phase] * sqrt(fluidProp[n][MW]);
		s1Ml += tMl*fluidProp[n][MU];
		s2Ml += tMl;
		roL += fcomp[Ix][Iy][Iz][n][Phase] * fluidProp[n][VCRIT];
		compTCL += fcomp[Ix][Iy][Iz][n][Phase] * fluidProp[n][TCRIT];
		compPCL += fcomp[Ix][Iy][Iz][n][Phase] * fluidProp[n][PCRIT] / 101325;
		compMWL += fcomp[Ix][Iy][Iz][n][Phase] * fluidProp[n][MW];
	}
	MuL = s1Ml / s2Ml;
	ErL = fblockFProps[Ix][Iy][Iz][RO][Phase] * roL;
	eAAl = pow(compTCL, 1.0 / 6);
	eBBl = sqrt(compMWL);
	eCCl = pow(compPCL, 2.0 / 3);
	ethaL = eAAl / (eBBl*eCCl);
	ro1L = 0.10230 + 0.023364*ErL + 0.058533*ErL*ErL - 0.040758*ErL*ErL*ErL + 0.0093324*ErL*ErL*ErL*ErL;

	if (Component < Nc) {
		dMuSdx = (fluidProp[Component][MU] * sqrt(fluidProp[Component][MW])*s2Ml - sqrt(fluidProp[Component][MW])*s1Ml) / (s2Ml*s2Ml);
		dEthadx = ((1.0 / 6.0)*pow(compTCL, -5.0 / 6.0)*fluidProp[Component][TCRIT] * eBBl*eCCl - eAAl*(0.5*fluidProp[Component][MW] * eCCl / sqrt(compMWL) + (2.0 / 3.0)*pow(compPCL, -1.0 / 3.0)*fluidProp[Component][PCRIT] * eBBl / 101325)) / (eBBl*eBBl*eCCl*eCCl);
		dErdx = fdE[Ix][Iy][Iz][Component][Phase] * roL + fblockFProps[Ix][Iy][Iz][RO][Phase] * fluidProp[Component][VCRIT];
	}
	else {
		dMuSdx = 0;
		dEthadx = 0;
		dErdx = fdE[Ix][Iy][Iz][Component][Phase] * roL;
	}

	dMMdx = 4 * ro1L*ro1L*ro1L*(0.023364*dErdx + 0.117066*dErdx*ErL - 0.122274*dErdx*ErL*ErL + 0.0373296*dErdx*ErL*ErL*ErL);
	dMuFdx = (dMMdx*ethaL - (ro1L*ro1L*ro1L*ro1L - 1e-4)*dEthadx) / (ethaL*ethaL) + dMuSdx;

	return dMuFdx;
}

int RotatingOne(int n) {
	switch (n % 4) {
	case 0:
		return 1;
	case 1:
		return 0;
	case 2:
		return -1;
	case 3:
		return 0;
	}

	return 0;
}

FType CalcDiffIntegral(int Ix, int Iy, int Iz, int jComponent, int DirectionInt, char wPhase) {
	int Ixc, Iyc, Izc, sgn;
	FType sum, tol, a, v, me, ex, L, x1, x2, D, w1, w2, DeltaC, sWs, T, DeltaX;
	register int n;
	
	Ixc = Ix;
	Iyc = Iy;
	Izc = Iz;
	switch (DirectionInt) {
	case 1:
		Ixc++;
		w1 = gridDim[Ix];
		w2 = gridDim[Ixc];		
		break;
	case -1:
		Ixc--;
		w1 = gridDim[Ix];
		w2 = gridDim[Ixc];
		break;
	case 2:
		Iyc++;
		w1 = gridDim[Nx + Iy];
		w2 = gridDim[Nx + Iyc];
		break;
	case -2:
		Iyc--;
		w1 = gridDim[Nx + Iy];
		w2 = gridDim[Nx + Iyc];
		break;
	case 3:
		Izc++;
		w1 = gridDim[Nx + Ny + Iz];
		w2 = gridDim[Nx + Ny + Izc];
		break;
	case -3:
		Izc--;
		w1 = gridDim[Nx + Ny + Iz];
		w2 = gridDim[Nx + Ny + Izc];
		break;
	default:
		TerM("Error in new diffusion function!");
	}

	x1 = comp[Ix][Iy][Iz][jComponent][wPhase];
	x2 = comp[Ixc][Iyc][Izc][jComponent][wPhase];

	if (DirectionInt > 0) DeltaX = x2 - x1;
	else DeltaX = x1 - x2;

	if (DeltaX > 0) sgn = 1;
	else if (DeltaX < 0) sgn = -1;
	else sgn = 0;

	sWs = w1 + w2;

	T = (tor[Ix][Iy][Iz] * w1 + tor[Ixc][Iyc][Izc] * w2) / sWs;
	D = (diffusion[Ix][Iy][Iz][jComponent][wPhase] * w1 + diffusion[Ixc][Iyc][Izc][jComponent][wPhase] * w2) / (sWs*T);
	sum = 0;	
	n = 1;
	do {
		v = 2 * (x2 - x1)*(PI*n*RotatingOne(n)) / ((PI*n)*(PI*n));
		ex = exp(-D * Dt*(n*PI / L)*(n*PI / L));
		me = 2 * (1 - RotatingOne(n)) / (n*PI);
		a = v * ex*me;
		sum = sum + a;

		if (a) tol = fabs(a);
		else tol = 1;
		
		n = n + 1;
	} while (tol > 1e-60);
	sum += x2 / 4 - x1 / 4;

	DeltaC = fabs((0.5*sum*w1*sWs) / (D*Dt));

	//if (DeltaC < fabs(DeltaX)) DeltaC = fabs(DeltaX);
	DeltaC *= sgn;

	return DeltaC;
}

FType fCalcDiffIntegral(int Ix, int Iy, int Iz, int jComponent, int DirectionInt, char wPhase) {
	int Ixc, Iyc, Izc, sgn;
	FType sum, tol, a, v, me, ex, L, x1, x2, D, w1, w2, DeltaC, sWs, DeltaX;
	register int n;

	Ixc = Ix;
	Iyc = Iy;
	Izc = Iz;
	switch (DirectionInt) {
	case 1:
		Ixc++;
		w1 = gridDim[Ix];
		w2 = gridDim[Ixc];
		break;
	case -1:
		Ixc--;
		w1 = gridDim[Ix];
		w2 = gridDim[Ixc];
		break;
	case 2:
		Iyc++;
		w1 = gridDim[Nx + Iy];
		w2 = gridDim[Nx + Iyc];
		break;
	case -2:
		Iyc--;
		w1 = gridDim[Nx + Iy];
		w2 = gridDim[Nx + Iyc];
		break;
	case 3:
		Izc++;
		w1 = gridDim[Nx + Ny + Iz];
		w2 = gridDim[Nx + Ny + Izc];
		break;
	case -3:
		Izc--;
		w1 = gridDim[Nx + Ny + Iz];
		w2 = gridDim[Nx + Ny + Izc];
		break;
	default:
		TerM("Error in new diffusion function!");
	}

	x1 = fcomp[Ix][Iy][Iz][jComponent][wPhase];
	x2 = fcomp[Ixc][Iyc][Izc][jComponent][wPhase];

	if (DirectionInt > 0) DeltaX = x2 - x1;
	else DeltaX = x1 - x2;

	if (DeltaX > 0) sgn = 1;
	else if (DeltaX < 0) sgn = -1;
	else sgn = 0;

	sWs = w1 + w2;

	D = (fdiffusion[Ix][Iy][Iz][jComponent][wPhase] * w1 + fdiffusion[Ixc][Iyc][Izc][jComponent][wPhase] * w2) / sWs;
	sum = 0;
	n = 1;
	do {
		v = 2 * (x2 - x1)*(PI*n*RotatingOne(n)) / ((PI*n)*(PI*n));
		ex = exp(-D * Dt*(n*PI / L)*(n*PI / L));
		me = 2 * (1 - RotatingOne(n)) / (n*PI);
		a = v * ex*me;
		sum = sum + a;

		if (a) tol = fabs(a);
		else tol = 1;

		n = n + 1;
	} while (tol > 1e-60);
	sum += x2 / 4 - x1 / 4;

	DeltaC = fabs((0.5*sum*w1*sWs) / (D*Dt));

	//if (DeltaC < fabs(DeltaX)) DeltaC = fabs(DeltaX);
	DeltaC *= sgn;

	return DeltaC;
}


#endif