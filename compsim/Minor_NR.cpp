#ifndef MINOR_NR_H
#define MINOR_NR_H

//#include "Globals.h"
//#include "Do_Cycle.h" 
#include <math.h>
//#include <stdio.h>
//#include "Serial_BiCGStab.h"
//#include "Gauss_Jordan.h"
//#include <MKL_Solver.h>
#include "Functions.h"

extern int Nc;
extern FType **fluidProp;
extern FType *****comp;
extern FType resTemp;
extern FType LF;
extern FType ****P;
extern unsigned char PR, SRK;
extern FType **satJac, *satAns, *Xm, *Xms;
extern FType *Ki;
extern FType **bic;
extern char gasInjStat;
extern char buildPreconFlag;


//FType Dew_NR(int, int, int);
//FType Bubl_Succ(int , int, int);
//FType Dew_Succ(int, int, int);
//FType LowerDew_Succ(int, int, int);
//void critTest(int, int, int);

FType Bubl_NR(int Ix, int Iy, int Iz) {
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG, dAdxg, dBdxg, dZdxg;
	FType dAAdxg, dBBdxg, DDl, DDg, EEl, EEg, FFl, FFg;
	FType dFFdxg, dEEdxg, dDDdxg, dPhidxg;
	FType dAdpg, dBdpg, dZdpg;
	FType dAAdpg, dBBdpg, dDDdpg, dFFdpg, dPhidpg;
	FType dUdpg, dWdpg;
	FType dUdxg, dWdxg;
	register int i, j, k, n;

	FType NRKi;
	FType SXm;
	FType XmTol;

	FType dAdpl, dBdpl, dZdpl;
	FType dAAdpl, dBBdpl, dDDdpl, dFFdpl, dPhidpl;
	FType dUdpl, dWdpl;
	FType tempSQR;
	FType mSum;


	Xm[Nc] = P[Ix + 1][Iy + 1][Iz + 1][1];
	for (i = 0; i < Nc; i++) {
		NRKi = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));		//Wilson Equation
		Xm[i] = NRKi*comp[Ix][Iy][Iz][i][0];
	}

	do {
		SXm = 0;
		for (i = 0; i < Nc; i++) {
			if (Xm[i] < 0) Xm[i] = 0;
			else if (Xm[i] > 1) Xm[i] = 1;
			SXm += Xm[i];
			satJac[Nc][i] = 1;
		}

		if (Xm[Nc] < 0) Xm[Nc] = PZERO;

		satAns[Nc] = -SXm + 1;
		satJac[Nc][Nc] = 0;



		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += comp[Ix][Iy][Iz][i][0] * fluidProp[i][EOS_B];
			Bg += Xm[i] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += comp[Ix][Iy][Iz][i][0] * comp[Ix][Iy][Iz][j][0] * tempSQR;
				Ag += Xm[i] * Xm[j] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bg *= Xm[Nc] / (RGAS*resTemp);


		dAdpg = Cag / (RGAS*RGAS*resTemp*resTemp);
		dBdpg = Cbg / (RGAS*resTemp);
		dAdpl = Cal / (RGAS*RGAS*resTemp*resTemp);
		dBdpl = Cbl / (RGAS*resTemp);



		if (SRK) {
			Ul = Bl;
			Wl = 0;
			Ug = Bg;
			Wg = 0;

			dUdpg = dBdpg;
			dWdpg = 0;
			dUdpl = dBdpl;
			dWdpl = 0;

			d1 = 1;
			d2 = 0;
		}
		else if (PR) {
			Ul = 2 * Bl;
			Wl = Bl;
			Ug = 2 * Bg;
			Wg = Bg;

			dUdpg = 2 * dBdpg;
			dWdpg = dBdpg;
			dUdpl = 2 * dBdpl;
			dWdpl = dBdpl;

			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl*Ul - Ul - Wl*Wl, -(Al*Bl - Bl*Wl*Wl - Wl*Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg*Ug - Ug - Wg*Wg, -(Ag*Bg - Bg*Wg*Wg - Wg*Wg), 'g');


		dZdpg = -(dAdpg*(Zg - Bg) + dBdpg*(-Zg*Zg - Ug*Zg - Ag + Wg*Wg) + dUdpg*(Zg*Zg - Bg*Zg - Zg) + dWdpg*(-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg*Ug - Ug - Wg*Wg);		//Modified
		dZdpl = -(dAdpl*(Zl - Bl) + dBdpl*(-Zl*Zl - Ul*Zl - Al + Wl*Wl) + dUdpl*(Zl*Zl - Bl*Zl - Zl) + dWdpl*(-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl*Ul - Ul - Wl*Wl);		//Modified


		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += comp[Ix][Iy][Iz][j][0] * tempSQR;
				tg += Xm[j] * tempSQR;
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

			satAns[i] = -comp[Ix][Iy][Iz][i][0] * phiL + Xm[i] * phiG;

			dAAdpg = fluidProp[i][EOS_B] * dZdpg / Cbg;
			dAAdpl = fluidProp[i][EOS_B] * dZdpl / Cbl;

			dBBdpg = (dBdpg - dZdpg) / (Zg - Bg);
			dBBdpl = (dBdpl - dZdpl) / (Zl - Bl);

			dDDdpg = (Bg*dAdpg - Ag*dBdpg) / (Bg*Bg*(d1 - d2));
			dDDdpl = (Bl*dAdpl - Al*dBdpl) / (Bl*Bl*(d1 - d2));

			dFFdpg = (dZdpg + d2*dBdpg) / (Zg + d2*Bg) - (dZdpg + d1*dBdpg) / (Zg + d1*Bg);
			dFFdpl = (dZdpl + d2*dBdpl) / (Zl + d2*Bl) - (dZdpl + d1*dBdpl) / (Zl + d1*Bl);

			dPhidpg = phiG*(dAAdpg + dBBdpg + EEg*(dDDdpg*FFg + DDg*dFFdpg));
			dPhidpl = phiL*(dAAdpl + dBBdpl + EEl*(dDDdpl*FFl + DDl*dFFdpl));

			satJac[i][Nc] = comp[Ix][Iy][Iz][i][0] * dPhidpl - Xm[i] * dPhidpg;		//Modified

																					//i: Row of Jac
																					//k: Col of Jac
			for (k = 0; k < Nc; k++) {
				dAdxg = 0;
				for (n = 0; n < Nc; n++) {		//SSFF
					dAdxg += Xm[n] * sqrt(fluidProp[n][EOS_A] * fluidProp[k][EOS_A])*bic[n][k];
					//This loop can be reduced
				}
				dAdxg *= 2 * Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
				dBdxg = fluidProp[k][EOS_B] * Xm[Nc] / (RGAS*resTemp);

				if (SRK) {
					dUdxg = dBdxg;
					dWdxg = 0;
				}
				else if (PR) {
					dUdxg = 2 * dBdxg;
					dWdxg = dBdxg;
				}


				dZdxg = -(dAdxg*(Zg - Bg) + dBdxg*(-Zg*Zg - Ug*Zg - Ag + Wg*Wg) + dUdxg*(Zg*Zg - Bg*Zg - Zg) + dWdxg*(-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg*Ug - Ug - Wg*Wg);		//Modified

				dAAdxg = -fluidProp[k][EOS_B] * fluidProp[i][EOS_B] * (Zg - 1) / (Cbg*Cbg) + fluidProp[i][EOS_B] * dZdxg / Cbg;

				dBBdxg = (dBdxg - dZdxg) / (Zg - Bg);

				dDDdxg = (dAdxg*Bg - dBdxg*Ag) / (Bg*Bg*(d1 - d2));

				dFFdxg = (dZdxg + d2*dBdxg) / (Zg + d2*Bg) - (dZdxg + d1*dBdxg) / (Zg + d1*Bg);

				dEEdxg = 2 * (bic[i][k] * sqrt(fluidProp[i][EOS_A] * fluidProp[k][EOS_A]) / Cag - dAdxg*RGAS*RGAS*resTemp*resTemp*tg / (Xm[Nc] * Cag*Cag)) + fluidProp[i][EOS_B] * fluidProp[k][EOS_B] / (Cbg*Cbg);		//modified

				dPhidxg = phiG*(dAAdxg + dBBdxg + dDDdxg*EEg*FFg + DDg*dEEdxg*FFg + DDg*EEg*dFFdxg);

				if (i == k) {
					satJac[i][k] = -phiG - Xm[i] * dPhidxg;
				}
				else {
					satJac[i][k] = -Xm[i] * dPhidxg;
				}


			}


		}



		//Scaling(satJac, satAns, Nc+1);
		//biCGStab(satJac, satAns, Xms, Nc+1);
		//GaussJ(satJac, satAns, Xms, Nc+1);

		MKLS(satJac, satAns, Xms, Nc + 1);             //Return later Sadeq

		
		XmTol = 0;
		for (i = 0; i < (Nc + 1); i++) {
			Xm[i] += Xms[i];
			XmTol += fabs(Xms[i] / Xm[i]);			
		}
		XmTol /= (Nc + 1);

		

		//XmTol=In_Product(Xms, Xms, Nc+1);

	} while ((XmTol > XMTOL) && (fabs(Xms[Nc]) > PRESSTOL));

	mSum = 0;
	for (i = 0; i < Nc; i++) {
		mSum += fabs(1 - Xm[i] / comp[Ix][Iy][Iz][i][0]);
	}

	return Xm[Nc];
}

FType Dew_NR(int Ix, int Iy, int Iz) {
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType dAdxl, dBdxl, dZdxl, phiL, phiG;
	FType dAAdxl, dBBdxl, DDl, DDg, EEl, EEg, FFl, FFg;
	FType dFFdxl, dEEdxl, dDDdxl, dPhidxl;
	FType dAdpl, dBdpl, dZdpl;
	FType dAAdpl, dBBdpl, dDDdpl, dFFdpl, dPhidpl;
	FType dUdpl, dWdpl;
	FType dUdxl, dWdxl;
	register int i, j, k, n;

	FType NRKi;
	FType SXm;
	FType XmTol;

	FType dAdpg, dBdpg, dZdpg;
	FType dAAdpg, dBBdpg, dDDdpg, dFFdpg, dPhidpg;
	FType dUdpg, dWdpg;
	FType tempSQR;



	Xm[Nc] = P[Ix + 1][Iy + 1][Iz + 1][2];
	for (i = 0; i < Nc; i++) {
		NRKi = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));		//Wilson Equation
		Xm[i] = comp[Ix][Iy][Iz][i][1] / NRKi;
	}


	do {
		SXm = 0;
		for (i = 0; i < Nc; i++) {
			if (Xm[i] < 0) Xm[i] = 0;
			else if (Xm[i] > 1) Xm[i] = 1;
			SXm += Xm[i];
			satJac[Nc][i] = 1;
		}

		if (Xm[Nc] < 0) Xm[Nc] = PZERO;

		satAns[Nc] = -SXm + 1;
		satJac[Nc][Nc] = 0;


		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += Xm[i] * fluidProp[i][EOS_B];
			Bg += comp[Ix][Iy][Iz][i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += Xm[i] * Xm[j] * tempSQR;
				Ag += comp[Ix][Iy][Iz][i][1] * comp[Ix][Iy][Iz][j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bg *= Xm[Nc] / (RGAS*resTemp);

		dAdpl = Cal / (RGAS*RGAS*resTemp*resTemp);
		dBdpl = Cbl / (RGAS*resTemp);
		dAdpg = Cag / (RGAS*RGAS*resTemp*resTemp);
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

			d1 = 1;
			d2 = 0;
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

			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl*Ul - Ul - Wl*Wl, -(Al*Bl - Bl*Wl*Wl - Wl*Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg*Ug - Ug - Wg*Wg, -(Ag*Bg - Bg*Wg*Wg - Wg*Wg), 'g');



		dZdpg = -(dAdpg*(Zg - Bg) + dBdpg*(-Zg*Zg - Ug*Zg - Ag + Wg*Wg) + dUdpg*(Zg*Zg - Bg*Zg - Zg) + dWdpg*(-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg*Ug - Ug - Wg*Wg);		//Modified
		dZdpl = -(dAdpl*(Zl - Bl) + dBdpl*(-Zl*Zl - Ul*Zl - Al + Wl*Wl) + dUdpl*(Zl*Zl - Bl*Zl - Zl) + dWdpl*(-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl*Ul - Ul - Wl*Wl);		//Modified

		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += Xm[j] * tempSQR;
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

			satAns[i] = -Xm[i] * phiL + comp[Ix][Iy][Iz][i][1] * phiG;

			dAAdpl = fluidProp[i][EOS_B] * dZdpl / Cbl;
			dAAdpg = fluidProp[i][EOS_B] * dZdpg / Cbg;

			dBBdpl = (dBdpl - dZdpl) / (Zl - Bl);
			dBBdpg = (dBdpg - dZdpg) / (Zg - Bg);

			dDDdpl = (Bl*dAdpl - Al*dBdpl) / (Bl*Bl*(d1 - d2));
			dDDdpg = (Bg*dAdpg - Ag*dBdpg) / (Bg*Bg*(d1 - d2));

			dFFdpl = (dZdpl + d2*dBdpl) / (Zl + d2*Bl) - (dZdpl + d1*dBdpl) / (Zl + d1*Bl);
			dFFdpg = (dZdpg + d2*dBdpg) / (Zg + d2*Bg) - (dZdpg + d1*dBdpg) / (Zg + d1*Bg);

			dPhidpl = phiL*(dAAdpl + dBBdpl + EEl*(dDDdpl*FFl + DDl*dFFdpl));
			dPhidpg = phiG*(dAAdpg + dBBdpg + EEg*(dDDdpg*FFg + DDg*dFFdpg));

			satJac[i][Nc] = Xm[i] * dPhidpl - comp[Ix][Iy][Iz][i][1] * dPhidpg;


			//i: Row of Jac
			//k: Col of Jac
			for (k = 0; k < Nc; k++) {
				dAdxl = 0;
				for (n = 0; n < Nc; n++) {		//SSFF
					dAdxl += Xm[n] * sqrt(fluidProp[n][EOS_A] * fluidProp[k][EOS_A])*bic[n][k];
					//This loop can be reduced
				}
				dAdxl *= 2 * Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
				dBdxl = fluidProp[k][EOS_B] * Xm[Nc] / (RGAS*resTemp);

				if (SRK) {
					dUdxl = dBdxl;
					dWdxl = 0;
				}
				else if (PR) {
					dUdxl = 2 * dBdxl;
					dWdxl = dBdxl;
				}



				dZdxl = -(dAdxl*(Zl - Bl) + dBdxl*(-Zl*Zl - Ul*Zl - Al + Wl*Wl) + dUdxl*(Zl*Zl - Bl*Zl - Zl) + dWdxl*(-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl*Ul - Ul - Wl*Wl);		//Modified

				dAAdxl = -fluidProp[k][EOS_B] * fluidProp[i][EOS_B] * (Zl - 1) / (Cbl*Cbl) + fluidProp[i][EOS_B] * dZdxl / Cbl;

				dBBdxl = (dBdxl - dZdxl) / (Zl - Bl);

				dDDdxl = (dAdxl*Bl - dBdxl*Al) / (Bl*Bl*(d1 - d2));

				dFFdxl = (dZdxl + d2*dBdxl) / (Zl + d2*Bl) - (dZdxl + d1*dBdxl) / (Zl + d1*Bl);

				dEEdxl = 2 * (bic[i][k] * sqrt(fluidProp[i][EOS_A] * fluidProp[k][EOS_A]) / Cal - dAdxl*RGAS*RGAS*resTemp*resTemp*tl / (Xm[Nc] * Cal*Cal)) + fluidProp[i][EOS_B] * fluidProp[k][EOS_B] / (Cbl*Cbl);

				dPhidxl = phiL*(dAAdxl + dBBdxl + dDDdxl*EEl*FFl + DDl*dEEdxl*FFl + DDl*EEl*dFFdxl);

				if (i == k) {
					satJac[i][k] = phiL + Xm[i] * dPhidxl;
				}
				else {
					satJac[i][k] = Xm[i] * dPhidxl;
				}


			}
		}

		//Scaling(satJac, satAns, Nc+1);
		//biCGStab(satJac, satAns, Xms, Nc+1);
		//GaussJ(satJac, satAns, Xms, Nc+1);
		MKLS(satJac, satAns, Xms, Nc + 1);         //Sadeq


		XmTol = 0;
		for (i = 0; i < (Nc + 1); i++) {
			Xm[i] += Xms[i];
			XmTol += fabs(Xms[i] / Xm[i]);
		}
		XmTol /= (Nc + 1);

		//XmTol=In_Product(Xms, Xms, Nc+1);

	} while ((XmTol > XMTOL) && (fabs(Xms[Nc]) > PRESSTOL));

	return Xm[Nc];
}

FType Bubl_Succ(int Ix, int Iy, int Iz) {
	register int nIter;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	//FType phiL, phiG;
	//FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType SXm;
	FType tempSQR;
	FType sumK;

	//FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	//FType Ug, Wg;
	//FType tl, tg, d1, d2, t1l, t1g;
	//FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG, dAdxg, dBdxg, dZdxg;
	FType dAAdxg, dBBdxg, DDl, DDg, EEl, EEg, FFl, FFg;
	FType dFFdxg, dEEdxg, dDDdxg, dPhidxg;
	FType dAdpg, dBdpg, dZdpg;
	FType dAAdpg, dBBdpg, dDDdpg, dFFdpg, dPhidpg;
	FType dUdpg, dWdpg;
	FType dUdxg, dWdxg;
	register int i, j, k, n;

	FType NRKi;
	//FType SXm;
	FType XmTol;

	FType dAdpl, dBdpl, dZdpl;
	FType dAAdpl, dBBdpl, dDDdpl, dFFdpl, dPhidpl;
	FType dUdpl, dWdpl;
	//FType tempSQR;
	//FType mSum;
	char returnVal = 0;
	FType xSum;


	Xm[Nc] = P[Ix + 1][Iy + 1][Iz + 1][1];
	for (i = 0; i < Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));		//Wilson Equation
	}

	nIter = 0;
	do {
		xSum = 0;
		for (i = 0; i < Nc; i++) {
			Xm[i] = comp[Ix][Iy][Iz][i][1] / Ki[i];
			if (Xm[i] < 0) Xm[i] = 0;
			else if (Xm[i] > 1) Xm[i] = 1;
			xSum += Xm[i];
		}
		if (Xm[Nc] < 0) Xm[Nc] = PZERO;
		for (i = 0; i < Nc; i++) Xm[i] /= xSum;

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += comp[Ix][Iy][Iz][i][0] * fluidProp[i][EOS_B];
			Bg += Xm[i] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += comp[Ix][Iy][Iz][i][0] * comp[Ix][Iy][Iz][j][0] * tempSQR;
				Ag += Xm[i] * Xm[j] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);		//gas-oil
		Bg *= Xm[Nc] / (RGAS*resTemp);


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


		SXm = 0;
		k = 0;
		sumK = 0;
		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += comp[Ix][Iy][Iz][j][0] * tempSQR;
				tg += Xm[j] * tempSQR;
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

			fL = comp[Ix][Iy][Iz][i][0] * phiL;
			fG = Xm[i] * phiG;	//gas-oil

			if (comp[Ix][Iy][Iz][i][0]) {
				Ki[i] *= fL / fG;
				SXm += Ki[i] * comp[Ix][Iy][Iz][i][0];
				k++;
				sumK += fabs(1 - Ki[i]);
			}
		}

		Xm[Nc] *= SXm;
		nIter++;
		if (nIter > MAX_MNR_ITER) {
			returnVal = -1;
			break;
		}
	} while ((fabs(SXm - 1) / k) > SXMTOL);
	//////////////////////////////////////////////////////////////
	
	if ((sumK / k) < TRIV_K) returnVal = -1;
	if ((returnVal != (-1)) && ((Xm[Nc] - Xm[Nc]) == 0)) return Xm[Nc];	
	
	

	////////////////////////////////////////////////////////////////////////////////
	


	Xm[Nc] = P[Ix + 1][Iy + 1][Iz + 1][1];
	for (i = 0; i < Nc; i++) {
		NRKi = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));		//Wilson Equation
		Xm[i] = NRKi*comp[Ix][Iy][Iz][i][0];
	}

	nIter = 0;
	do {
		SXm = 0;
		for (i = 0; i < Nc; i++) {			
			SXm += Xm[i];
			satJac[Nc][i] = 1;
		}
		satAns[Nc] = -SXm + 1;
		satJac[Nc][Nc] = 0;



		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += comp[Ix][Iy][Iz][i][0] * fluidProp[i][EOS_B];
			Bg += Xm[i] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += comp[Ix][Iy][Iz][i][0] * comp[Ix][Iy][Iz][j][0] * tempSQR;
				Ag += Xm[i] * Xm[j] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bg *= Xm[Nc] / (RGAS*resTemp);


		dAdpg = Cag / (RGAS*RGAS*resTemp*resTemp);
		dBdpg = Cbg / (RGAS*resTemp);
		dAdpl = Cal / (RGAS*RGAS*resTemp*resTemp);
		dBdpl = Cbl / (RGAS*resTemp);



		if (SRK) {
			Ul = Bl;
			Wl = 0;
			Ug = Bg;
			Wg = 0;

			dUdpg = dBdpg;
			dWdpg = 0;
			dUdpl = dBdpl;
			dWdpl = 0;

			d1 = 1;
			d2 = 0;
		}
		else if (PR) {
			Ul = 2 * Bl;
			Wl = Bl;
			Ug = 2 * Bg;
			Wg = Bg;

			dUdpg = 2 * dBdpg;
			dWdpg = dBdpg;
			dUdpl = 2 * dBdpl;
			dWdpl = dBdpl;

			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl*Ul - Ul - Wl*Wl, -(Al*Bl - Bl*Wl*Wl - Wl*Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg*Ug - Ug - Wg*Wg, -(Ag*Bg - Bg*Wg*Wg - Wg*Wg), 'g');


		dZdpg = -(dAdpg*(Zg - Bg) + dBdpg*(-Zg*Zg - Ug*Zg - Ag + Wg*Wg) + dUdpg*(Zg*Zg - Bg*Zg - Zg) + dWdpg*(-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg*Ug - Ug - Wg*Wg);		//Modified
		dZdpl = -(dAdpl*(Zl - Bl) + dBdpl*(-Zl*Zl - Ul*Zl - Al + Wl*Wl) + dUdpl*(Zl*Zl - Bl*Zl - Zl) + dWdpl*(-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl*Ul - Ul - Wl*Wl);		//Modified


		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += comp[Ix][Iy][Iz][j][0] * tempSQR;
				tg += Xm[j] * tempSQR;
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

			satAns[i] = -comp[Ix][Iy][Iz][i][0] * phiL + Xm[i] * phiG;

			dAAdpg = fluidProp[i][EOS_B] * dZdpg / Cbg;
			dAAdpl = fluidProp[i][EOS_B] * dZdpl / Cbl;

			dBBdpg = (dBdpg - dZdpg) / (Zg - Bg);
			dBBdpl = (dBdpl - dZdpl) / (Zl - Bl);

			dDDdpg = (Bg*dAdpg - Ag*dBdpg) / (Bg*Bg*(d1 - d2));
			dDDdpl = (Bl*dAdpl - Al*dBdpl) / (Bl*Bl*(d1 - d2));

			dFFdpg = (dZdpg + d2*dBdpg) / (Zg + d2*Bg) - (dZdpg + d1*dBdpg) / (Zg + d1*Bg);
			dFFdpl = (dZdpl + d2*dBdpl) / (Zl + d2*Bl) - (dZdpl + d1*dBdpl) / (Zl + d1*Bl);

			dPhidpg = phiG*(dAAdpg + dBBdpg + EEg*(dDDdpg*FFg + DDg*dFFdpg));
			dPhidpl = phiL*(dAAdpl + dBBdpl + EEl*(dDDdpl*FFl + DDl*dFFdpl));

			satJac[i][Nc] = comp[Ix][Iy][Iz][i][0] * dPhidpl - Xm[i] * dPhidpg;		//Modified

																					//i: Row of Jac
																					//k: Col of Jac
			for (k = 0; k < Nc; k++) {
				dAdxg = 0;
				for (n = 0; n < Nc; n++) {		//SSFF
					dAdxg += Xm[n] * sqrt(fluidProp[n][EOS_A] * fluidProp[k][EOS_A])*bic[n][k];
					//This loop can be reduced
				}
				dAdxg *= 2 * Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
				dBdxg = fluidProp[k][EOS_B] * Xm[Nc] / (RGAS*resTemp);

				if (SRK) {
					dUdxg = dBdxg;
					dWdxg = 0;
				}
				else if (PR) {
					dUdxg = 2 * dBdxg;
					dWdxg = dBdxg;
				}


				dZdxg = -(dAdxg*(Zg - Bg) + dBdxg*(-Zg*Zg - Ug*Zg - Ag + Wg*Wg) + dUdxg*(Zg*Zg - Bg*Zg - Zg) + dWdxg*(-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg*Ug - Ug - Wg*Wg);		//Modified

				dAAdxg = -fluidProp[k][EOS_B] * fluidProp[i][EOS_B] * (Zg - 1) / (Cbg*Cbg) + fluidProp[i][EOS_B] * dZdxg / Cbg;

				dBBdxg = (dBdxg - dZdxg) / (Zg - Bg);

				dDDdxg = (dAdxg*Bg - dBdxg*Ag) / (Bg*Bg*(d1 - d2));

				dFFdxg = (dZdxg + d2*dBdxg) / (Zg + d2*Bg) - (dZdxg + d1*dBdxg) / (Zg + d1*Bg);

				dEEdxg = 2 * (bic[i][k] * sqrt(fluidProp[i][EOS_A] * fluidProp[k][EOS_A]) / Cag - dAdxg*RGAS*RGAS*resTemp*resTemp*tg / (Xm[Nc] * Cag*Cag)) + fluidProp[i][EOS_B] * fluidProp[k][EOS_B] / (Cbg*Cbg);		//modified

				dPhidxg = phiG*(dAAdxg + dBBdxg + dDDdxg*EEg*FFg + DDg*dEEdxg*FFg + DDg*EEg*dFFdxg);

				if (i == k) {
					satJac[i][k] = -phiG - Xm[i] * dPhidxg;
				}
				else {
					satJac[i][k] = -Xm[i] * dPhidxg;
				}


			}


		}



		//Scaling(satJac, satAns, Nc+1);
		//biCGStab(satJac, satAns, Xms, Nc+1);
		//GaussJ(satJac, satAns, Xms, Nc+1);

		MKLS(satJac, satAns, Xms, Nc + 1);             //Return later Sadeq


		XmTol = 0;
		xSum = 0;
		for (i = 0; i < (Nc + 1); i++) {
			Xm[i] += Xms[i];
			if (i < Nc) {
				if (Xm[i] < 0) Xm[i] = 0;
				else if (Xm[i] > 1) Xm[i] = 1;
				xSum += Xm[i];
			}
			else {
				if (Xm[Nc] < 0) Xm[Nc] = PZERO;
			}
			if (Xm[i]) XmTol += fabs(Xms[i] / Xm[i]);
			else XmTol += fabs(Xms[i]);
		}
		XmTol /= (Nc + 1);		
		for (i = 0; i < Nc; i++) Xm[i] /= xSum;

		nIter++;
		if (nIter > MAX_MNR_ITER) {
			return -1;			
		}

		//XmTol=In_Product(Xms, Xms, Nc+1);

	} while (XmTol > XMTOL);// && (fabs(Xms[Nc]) > PRESSTOL));

	sumK = 0;
	k = 0;
	for (i = 0; i < Nc; i++) {
		if (comp[Ix][Iy][Iz][i][0]) {
			sumK += fabs(1 - Xm[i] / comp[Ix][Iy][Iz][i][0]);
			k++;
		}		
	}
	if ((sumK / k) < TRIV_K) return -1;
	else return Xm[Nc];
}

FType Dew_Succ(int Ix, int Iy, int Iz) {
	register int i, j, k, nIter;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	//FType phiL, phiG;
	//FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType SXm;
	//FType tempSQR;
	FType sumK;

	//FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	//FType Ug, Wg;
	//FType tl, tg, d1, d2, t1l, t1g;
	//FType Cal, Cbl, Cag, Cbg;
	FType dAdxl, dBdxl, dZdxl, phiL, phiG;
	FType dAAdxl, dBBdxl, DDl, DDg, EEl, EEg, FFl, FFg;
	FType dFFdxl, dEEdxl, dDDdxl, dPhidxl;
	FType dAdpl, dBdpl, dZdpl;
	FType dAAdpl, dBBdpl, dDDdpl, dFFdpl, dPhidpl;
	FType dUdpl, dWdpl;
	FType dUdxl, dWdxl;
	register int n;

	FType NRKi;
	//FType SXm;
	FType XmTol;

	FType dAdpg, dBdpg, dZdpg;
	FType dAAdpg, dBBdpg, dDDdpg, dFFdpg, dPhidpg;
	FType dUdpg, dWdpg;
	FType tempSQR;
	char returnVal = 0;
	FType xSum;


	Xm[Nc] = P[Ix + 1][Iy + 1][Iz + 1][2];
	for (i = 0; i < Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));
	}

	nIter = 0;
	do {
		xSum = 0;
		for (i = 0; i < Nc; i++) {
			Xm[i] = comp[Ix][Iy][Iz][i][1] / Ki[i];
			if (Xm[i] < 0) Xm[i] = 0;
			else if (Xm[i] > 1) Xm[i] = 1;
			xSum += Xm[i];
		}
		if (Xm[Nc] < 0) Xm[Nc] = PZERO;
		for (i = 0; i < Nc; i++) Xm[i] /= xSum;

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += Xm[i] * fluidProp[i][EOS_B];
			Bg += comp[Ix][Iy][Iz][i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];;
				Al += Xm[i] * Xm[j] * tempSQR;
				Ag += comp[Ix][Iy][Iz][i][1] * comp[Ix][Iy][Iz][j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);		//gas-oil
		Bg *= Xm[Nc] / (RGAS*resTemp);


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


		SXm = 0;
		sumK = 0;
		k = 0;
		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += Xm[j] * tempSQR;
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

			fL = Xm[i] * phiL;
			fG = comp[Ix][Iy][Iz][i][1] * phiG;	//gas-oil

			if (comp[Ix][Iy][Iz][i][1]) {
				Ki[i] *= fL / fG;
				SXm += comp[Ix][Iy][Iz][i][1] / Ki[i];
				sumK += fabs(1 - Ki[i]);
				k++;
			}
		}

		Xm[Nc] *= SXm;
		nIter++;
		if (nIter > MAX_MNR_ITER) {
			returnVal = -1;
			break;
		}
	} while ((fabs(SXm - 1) / k) > SXMTOL);

	if ((sumK / k) < TRIV_K) returnVal = -1;
	if ((returnVal!=(-1)) && ((Xm[Nc] - Xm[Nc]) == 0)) return Xm[Nc];

	//////////////////////////////////////////////////////////////////////////////////////////////
	



	Xm[Nc] = P[Ix + 1][Iy + 1][Iz + 1][2];
	for (i = 0; i < Nc; i++) {
		NRKi = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));		//Wilson Equation
		Xm[i] = comp[Ix][Iy][Iz][i][1] / NRKi;
	}

	nIter = 0;
	do {
		SXm = 0;
		for (i = 0; i < Nc; i++) {
			//if (Xm[i] < 0) Xm[i] = 0;
			//else if (Xm[i] > 1) Xm[i] = 1;
			SXm += Xm[i];
			satJac[Nc][i] = 1;
		}

		//if (Xm[Nc] < 0) Xm[Nc] = PZERO;

		satAns[Nc] = -SXm + 1;
		satJac[Nc][Nc] = 0;


		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += Xm[i] * fluidProp[i][EOS_B];
			Bg += comp[Ix][Iy][Iz][i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += Xm[i] * Xm[j] * tempSQR;
				Ag += comp[Ix][Iy][Iz][i][1] * comp[Ix][Iy][Iz][j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bg *= Xm[Nc] / (RGAS*resTemp);

		dAdpl = Cal / (RGAS*RGAS*resTemp*resTemp);
		dBdpl = Cbl / (RGAS*resTemp);
		dAdpg = Cag / (RGAS*RGAS*resTemp*resTemp);
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

			d1 = 1;
			d2 = 0;
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

			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl*Ul - Ul - Wl*Wl, -(Al*Bl - Bl*Wl*Wl - Wl*Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg*Ug - Ug - Wg*Wg, -(Ag*Bg - Bg*Wg*Wg - Wg*Wg), 'g');



		dZdpg = -(dAdpg*(Zg - Bg) + dBdpg*(-Zg*Zg - Ug*Zg - Ag + Wg*Wg) + dUdpg*(Zg*Zg - Bg*Zg - Zg) + dWdpg*(-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg*Ug - Ug - Wg*Wg);		//Modified
		dZdpl = -(dAdpl*(Zl - Bl) + dBdpl*(-Zl*Zl - Ul*Zl - Al + Wl*Wl) + dUdpl*(Zl*Zl - Bl*Zl - Zl) + dWdpl*(-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl*Ul - Ul - Wl*Wl);		//Modified

		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += Xm[j] * tempSQR;
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

			satAns[i] = -Xm[i] * phiL + comp[Ix][Iy][Iz][i][1] * phiG;

			dAAdpl = fluidProp[i][EOS_B] * dZdpl / Cbl;
			dAAdpg = fluidProp[i][EOS_B] * dZdpg / Cbg;

			dBBdpl = (dBdpl - dZdpl) / (Zl - Bl);
			dBBdpg = (dBdpg - dZdpg) / (Zg - Bg);

			dDDdpl = (Bl*dAdpl - Al*dBdpl) / (Bl*Bl*(d1 - d2));
			dDDdpg = (Bg*dAdpg - Ag*dBdpg) / (Bg*Bg*(d1 - d2));

			dFFdpl = (dZdpl + d2*dBdpl) / (Zl + d2*Bl) - (dZdpl + d1*dBdpl) / (Zl + d1*Bl);
			dFFdpg = (dZdpg + d2*dBdpg) / (Zg + d2*Bg) - (dZdpg + d1*dBdpg) / (Zg + d1*Bg);

			dPhidpl = phiL*(dAAdpl + dBBdpl + EEl*(dDDdpl*FFl + DDl*dFFdpl));
			dPhidpg = phiG*(dAAdpg + dBBdpg + EEg*(dDDdpg*FFg + DDg*dFFdpg));

			satJac[i][Nc] = Xm[i] * dPhidpl - comp[Ix][Iy][Iz][i][1] * dPhidpg;


			//i: Row of Jac
			//k: Col of Jac
			for (k = 0; k < Nc; k++) {
				dAdxl = 0;
				for (n = 0; n < Nc; n++) {		//SSFF
					dAdxl += Xm[n] * sqrt(fluidProp[n][EOS_A] * fluidProp[k][EOS_A])*bic[n][k];
					//This loop can be reduced
				}
				dAdxl *= 2 * Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
				dBdxl = fluidProp[k][EOS_B] * Xm[Nc] / (RGAS*resTemp);

				if (SRK) {
					dUdxl = dBdxl;
					dWdxl = 0;
				}
				else if (PR) {
					dUdxl = 2 * dBdxl;
					dWdxl = dBdxl;
				}



				dZdxl = -(dAdxl*(Zl - Bl) + dBdxl*(-Zl*Zl - Ul*Zl - Al + Wl*Wl) + dUdxl*(Zl*Zl - Bl*Zl - Zl) + dWdxl*(-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl*Ul - Ul - Wl*Wl);		//Modified

				dAAdxl = -fluidProp[k][EOS_B] * fluidProp[i][EOS_B] * (Zl - 1) / (Cbl*Cbl) + fluidProp[i][EOS_B] * dZdxl / Cbl;

				dBBdxl = (dBdxl - dZdxl) / (Zl - Bl);

				dDDdxl = (dAdxl*Bl - dBdxl*Al) / (Bl*Bl*(d1 - d2));

				dFFdxl = (dZdxl + d2*dBdxl) / (Zl + d2*Bl) - (dZdxl + d1*dBdxl) / (Zl + d1*Bl);

				dEEdxl = 2 * (bic[i][k] * sqrt(fluidProp[i][EOS_A] * fluidProp[k][EOS_A]) / Cal - dAdxl*RGAS*RGAS*resTemp*resTemp*tl / (Xm[Nc] * Cal*Cal)) + fluidProp[i][EOS_B] * fluidProp[k][EOS_B] / (Cbl*Cbl);

				dPhidxl = phiL*(dAAdxl + dBBdxl + dDDdxl*EEl*FFl + DDl*dEEdxl*FFl + DDl*EEl*dFFdxl);

				if (i == k) {
					satJac[i][k] = phiL + Xm[i] * dPhidxl;
				}
				else {
					satJac[i][k] = Xm[i] * dPhidxl;
				}


			}
		}

		//Scaling(satJac, satAns, Nc+1);
		//biCGStab(satJac, satAns, Xms, Nc+1);
		//GaussJ(satJac, satAns, Xms, Nc+1);
		MKLS(satJac, satAns, Xms, Nc + 1);         //Sadeq


		
		XmTol = 0;
		xSum = 0;
		for (i = 0; i < (Nc + 1); i++) {
			Xm[i] += Xms[i];
			if (i < Nc) {
				if (Xm[i] < 0) Xm[i] = 0;
				else if (Xm[i] > 1) Xm[i] = 1;
				xSum += Xm[i];
			}
			else {
				if (Xm[Nc] < 0) Xm[Nc] = PZERO;
			}
			if (Xm[i]) XmTol += fabs(Xms[i] / Xm[i]);
			else XmTol += fabs(Xms[i]);
		}
		XmTol /= (Nc + 1);
		for (i = 0; i < Nc; i++) Xm[i] /= xSum;

		nIter++;
		if (nIter > MAX_MNR_ITER) {
			return -1;			
		}

	} while (XmTol > XMTOL);// && (fabs(Xms[Nc]) > PRESSTOL));

	sumK = 0;
	k = 0;
	for (i = 0; i < Nc; i++) {
		if (comp[Ix][Iy][Iz][i][1]) {
			sumK += fabs(1 - comp[Ix][Iy][Iz][i][1] / Xm[i]);
			k++;
		}
	}
	if ((sumK / k) < TRIV_K) return -1;
	else return Xm[Nc];
}

FType fBubl_Succ(int Ix, int Iy, int Iz) {
	register int nIter;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	//FType phiL, phiG;
	//FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType SXm;
	FType tempSQR;
	FType sumK;

	//FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	//FType Ug, Wg;
	//FType tl, tg, d1, d2, t1l, t1g;
	//FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG, dAdxg, dBdxg, dZdxg;
	FType dAAdxg, dBBdxg, DDl, DDg, EEl, EEg, FFl, FFg;
	FType dFFdxg, dEEdxg, dDDdxg, dPhidxg;
	FType dAdpg, dBdpg, dZdpg;
	FType dAAdpg, dBBdpg, dDDdpg, dFFdpg, dPhidpg;
	FType dUdpg, dWdpg;
	FType dUdxg, dWdxg;
	register int i, j, k, n;

	FType NRKi;
	//FType SXm;
	FType XmTol;

	FType dAdpl, dBdpl, dZdpl;
	FType dAAdpl, dBBdpl, dDDdpl, dFFdpl, dPhidpl;
	FType dUdpl, dWdpl;
	//FType tempSQR;
	//FType mSum;
	char returnVal = 0;
	FType xSum;


	Xm[Nc] = fP[Ix + 1][Iy + 1][Iz + 1][1];
	for (i = 0; i < Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));		//Wilson Equation
	}

	nIter = 0;
	do {
		xSum = 0;
		for (i = 0; i < Nc; i++) {
			Xm[i] = fcomp[Ix][Iy][Iz][i][1] / Ki[i];
			if (Xm[i] < 0) Xm[i] = 0;
			else if (Xm[i] > 1) Xm[i] = 1;
			xSum += Xm[i];
		}
		if (Xm[Nc] < 0) Xm[Nc] = PZERO;
		for (i = 0; i < Nc; i++) Xm[i] /= xSum;

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += fcomp[Ix][Iy][Iz][i][0] * fluidProp[i][EOS_B];
			Bg += Xm[i] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += fcomp[Ix][Iy][Iz][i][0] * fcomp[Ix][Iy][Iz][j][0] * tempSQR;
				Ag += Xm[i] * Xm[j] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);		//gas-oil
		Bg *= Xm[Nc] / (RGAS*resTemp);


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


		SXm = 0;
		k = 0;
		sumK = 0;
		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += fcomp[Ix][Iy][Iz][j][0] * tempSQR;
				tg += Xm[j] * tempSQR;
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

			fL = fcomp[Ix][Iy][Iz][i][0] * phiL;
			fG = Xm[i] * phiG;	//gas-oil

			if (fcomp[Ix][Iy][Iz][i][0]) {
				Ki[i] *= fL / fG;
				SXm += Ki[i] * fcomp[Ix][Iy][Iz][i][0];
				k++;
				sumK += fabs(1 - Ki[i]);
			}
		}

		Xm[Nc] *= SXm;
		nIter++;
		if (nIter > MAX_MNR_ITER) {
			returnVal = -1;
			break;
		}

	} while ((fabs(SXm - 1) / k) > SXMTOL);
	//////////////////////////////////////////////////////////////
	if ((sumK / k) < TRIV_K) returnVal = -1;
	if ((returnVal != (-1)) && ((Xm[Nc] - Xm[Nc]) == 0)) return Xm[Nc];



	////////////////////////////////////////////////////////////////////////////////



	Xm[Nc] = fP[Ix + 1][Iy + 1][Iz + 1][1];
	for (i = 0; i < Nc; i++) {
		NRKi = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));		//Wilson Equation
		Xm[i] = NRKi*fcomp[Ix][Iy][Iz][i][0];
	}

	nIter = 0;
	do {
		SXm = 0;
		for (i = 0; i < Nc; i++) {
			//if (Xm[i] < 0) Xm[i] = 0;
			//else if (Xm[i] > 1) Xm[i] = 1;
			SXm += Xm[i];
			satJac[Nc][i] = 1;
		}

		//if (Xm[Nc] < 0) Xm[Nc] = PZERO;

		satAns[Nc] = -SXm + 1;
		satJac[Nc][Nc] = 0;



		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += fcomp[Ix][Iy][Iz][i][0] * fluidProp[i][EOS_B];
			Bg += Xm[i] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += fcomp[Ix][Iy][Iz][i][0] * fcomp[Ix][Iy][Iz][j][0] * tempSQR;
				Ag += Xm[i] * Xm[j] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bg *= Xm[Nc] / (RGAS*resTemp);


		dAdpg = Cag / (RGAS*RGAS*resTemp*resTemp);
		dBdpg = Cbg / (RGAS*resTemp);
		dAdpl = Cal / (RGAS*RGAS*resTemp*resTemp);
		dBdpl = Cbl / (RGAS*resTemp);



		if (SRK) {
			Ul = Bl;
			Wl = 0;
			Ug = Bg;
			Wg = 0;

			dUdpg = dBdpg;
			dWdpg = 0;
			dUdpl = dBdpl;
			dWdpl = 0;

			d1 = 1;
			d2 = 0;
		}
		else if (PR) {
			Ul = 2 * Bl;
			Wl = Bl;
			Ug = 2 * Bg;
			Wg = Bg;

			dUdpg = 2 * dBdpg;
			dWdpg = dBdpg;
			dUdpl = 2 * dBdpl;
			dWdpl = dBdpl;

			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl*Ul - Ul - Wl*Wl, -(Al*Bl - Bl*Wl*Wl - Wl*Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg*Ug - Ug - Wg*Wg, -(Ag*Bg - Bg*Wg*Wg - Wg*Wg), 'g');


		dZdpg = -(dAdpg*(Zg - Bg) + dBdpg*(-Zg*Zg - Ug*Zg - Ag + Wg*Wg) + dUdpg*(Zg*Zg - Bg*Zg - Zg) + dWdpg*(-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg*Ug - Ug - Wg*Wg);		//Modified
		dZdpl = -(dAdpl*(Zl - Bl) + dBdpl*(-Zl*Zl - Ul*Zl - Al + Wl*Wl) + dUdpl*(Zl*Zl - Bl*Zl - Zl) + dWdpl*(-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl*Ul - Ul - Wl*Wl);		//Modified


		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += fcomp[Ix][Iy][Iz][j][0] * tempSQR;
				tg += Xm[j] * tempSQR;
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

			satAns[i] = -fcomp[Ix][Iy][Iz][i][0] * phiL + Xm[i] * phiG;

			dAAdpg = fluidProp[i][EOS_B] * dZdpg / Cbg;
			dAAdpl = fluidProp[i][EOS_B] * dZdpl / Cbl;

			dBBdpg = (dBdpg - dZdpg) / (Zg - Bg);
			dBBdpl = (dBdpl - dZdpl) / (Zl - Bl);

			dDDdpg = (Bg*dAdpg - Ag*dBdpg) / (Bg*Bg*(d1 - d2));
			dDDdpl = (Bl*dAdpl - Al*dBdpl) / (Bl*Bl*(d1 - d2));

			dFFdpg = (dZdpg + d2*dBdpg) / (Zg + d2*Bg) - (dZdpg + d1*dBdpg) / (Zg + d1*Bg);
			dFFdpl = (dZdpl + d2*dBdpl) / (Zl + d2*Bl) - (dZdpl + d1*dBdpl) / (Zl + d1*Bl);

			dPhidpg = phiG*(dAAdpg + dBBdpg + EEg*(dDDdpg*FFg + DDg*dFFdpg));
			dPhidpl = phiL*(dAAdpl + dBBdpl + EEl*(dDDdpl*FFl + DDl*dFFdpl));

			satJac[i][Nc] = fcomp[Ix][Iy][Iz][i][0] * dPhidpl - Xm[i] * dPhidpg;		//Modified

																					//i: Row of Jac
																					//k: Col of Jac
			for (k = 0; k < Nc; k++) {
				dAdxg = 0;
				for (n = 0; n < Nc; n++) {		//SSFF
					dAdxg += Xm[n] * sqrt(fluidProp[n][EOS_A] * fluidProp[k][EOS_A])*bic[n][k];
					//This loop can be reduced
				}
				dAdxg *= 2 * Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
				dBdxg = fluidProp[k][EOS_B] * Xm[Nc] / (RGAS*resTemp);

				if (SRK) {
					dUdxg = dBdxg;
					dWdxg = 0;
				}
				else if (PR) {
					dUdxg = 2 * dBdxg;
					dWdxg = dBdxg;
				}


				dZdxg = -(dAdxg*(Zg - Bg) + dBdxg*(-Zg*Zg - Ug*Zg - Ag + Wg*Wg) + dUdxg*(Zg*Zg - Bg*Zg - Zg) + dWdxg*(-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg*Ug - Ug - Wg*Wg);		//Modified

				dAAdxg = -fluidProp[k][EOS_B] * fluidProp[i][EOS_B] * (Zg - 1) / (Cbg*Cbg) + fluidProp[i][EOS_B] * dZdxg / Cbg;

				dBBdxg = (dBdxg - dZdxg) / (Zg - Bg);

				dDDdxg = (dAdxg*Bg - dBdxg*Ag) / (Bg*Bg*(d1 - d2));

				dFFdxg = (dZdxg + d2*dBdxg) / (Zg + d2*Bg) - (dZdxg + d1*dBdxg) / (Zg + d1*Bg);

				dEEdxg = 2 * (bic[i][k] * sqrt(fluidProp[i][EOS_A] * fluidProp[k][EOS_A]) / Cag - dAdxg*RGAS*RGAS*resTemp*resTemp*tg / (Xm[Nc] * Cag*Cag)) + fluidProp[i][EOS_B] * fluidProp[k][EOS_B] / (Cbg*Cbg);		//modified

				dPhidxg = phiG*(dAAdxg + dBBdxg + dDDdxg*EEg*FFg + DDg*dEEdxg*FFg + DDg*EEg*dFFdxg);

				if (i == k) {
					satJac[i][k] = -phiG - Xm[i] * dPhidxg;
				}
				else {
					satJac[i][k] = -Xm[i] * dPhidxg;
				}


			}


		}



		//Scaling(satJac, satAns, Nc+1);
		//biCGStab(satJac, satAns, Xms, Nc+1);
		//GaussJ(satJac, satAns, Xms, Nc+1);

		MKLS(satJac, satAns, Xms, Nc + 1);             //Return later Sadeq


		XmTol = 0;
		xSum = 0;
		for (i = 0; i < (Nc + 1); i++) {
			Xm[i] += Xms[i];
			if (i < Nc) {
				if (Xm[i] < 0) Xm[i] = 0;
				else if (Xm[i] > 1) Xm[i] = 1;
				xSum += Xm[i];
			}
			else {
				if (Xm[Nc] < 0) Xm[Nc] = PZERO;
			}
			if (Xm[i]) XmTol += fabs(Xms[i] / Xm[i]);
			else XmTol += fabs(Xms[i]);
		}
		XmTol /= (Nc + 1);
		for (i = 0; i < Nc; i++) Xm[i] /= xSum;

		nIter++;
		if (nIter > MAX_MNR_ITER) {
			return -1;
		}

		//XmTol=In_Product(Xms, Xms, Nc+1);

	} while (XmTol > XMTOL);// && (fabs(Xms[Nc]) > PRESSTOL));

	sumK = 0;
	k = 0;
	for (i = 0; i < Nc; i++) {
		if (fcomp[Ix][Iy][Iz][i][0]) {
			sumK += fabs(1 - Xm[i] / fcomp[Ix][Iy][Iz][i][0]);
			k++;
		}
	}
	if ((sumK / k) < TRIV_K) return -1;
	else return Xm[Nc];
}

FType fDew_Succ(int Ix, int Iy, int Iz) {
	register int i, j, k, nIter;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	//FType phiL, phiG;
	//FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType SXm;
	//FType tempSQR;
	FType sumK;

	//FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	//FType Ug, Wg;
	//FType tl, tg, d1, d2, t1l, t1g;
	//FType Cal, Cbl, Cag, Cbg;
	FType dAdxl, dBdxl, dZdxl, phiL, phiG;
	FType dAAdxl, dBBdxl, DDl, DDg, EEl, EEg, FFl, FFg;
	FType dFFdxl, dEEdxl, dDDdxl, dPhidxl;
	FType dAdpl, dBdpl, dZdpl;
	FType dAAdpl, dBBdpl, dDDdpl, dFFdpl, dPhidpl;
	FType dUdpl, dWdpl;
	FType dUdxl, dWdxl;
	register int n;

	FType NRKi;
	//FType SXm;
	FType XmTol;

	FType dAdpg, dBdpg, dZdpg;
	FType dAAdpg, dBBdpg, dDDdpg, dFFdpg, dPhidpg;
	FType dUdpg, dWdpg;
	FType tempSQR;
	char returnVal = 0;
	FType xSum;


	Xm[Nc] = fP[Ix + 1][Iy + 1][Iz + 1][2];
	for (i = 0; i < Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));
	}

	nIter = 0;
	do {
		xSum = 0;
		for (i = 0; i < Nc; i++) {
			Xm[i] = fcomp[Ix][Iy][Iz][i][1] / Ki[i];
			if (Xm[i] < 0) Xm[i] = 0;
			else if (Xm[i] > 1) Xm[i] = 1;
			xSum += Xm[i];
		}
		if (Xm[Nc] < 0) Xm[Nc] = PZERO;
		for (i = 0; i < Nc; i++) Xm[i] /= xSum;

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += Xm[i] * fluidProp[i][EOS_B];
			Bg += fcomp[Ix][Iy][Iz][i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];;
				Al += Xm[i] * Xm[j] * tempSQR;
				Ag += fcomp[Ix][Iy][Iz][i][1] * fcomp[Ix][Iy][Iz][j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);		//gas-oil
		Bg *= Xm[Nc] / (RGAS*resTemp);


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


		SXm = 0;
		sumK = 0;
		k = 0;
		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += Xm[j] * tempSQR;
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

			fL = Xm[i] * phiL;
			fG = fcomp[Ix][Iy][Iz][i][1] * phiG;	//gas-oil

			if (fcomp[Ix][Iy][Iz][i][1]) {
				Ki[i] *= fL / fG;
				SXm += fcomp[Ix][Iy][Iz][i][1] / Ki[i];
				sumK += fabs(1 - Ki[i]);
				k++;
			}
		}

		Xm[Nc] *= SXm;
		nIter++;
		if (nIter > MAX_MNR_ITER) {
			returnVal = -1;
			break;
		}

	} while ((fabs(SXm - 1) / k) > SXMTOL);

	if ((sumK / k) < TRIV_K) returnVal = -1;
	if ((returnVal != (-1)) && ((Xm[Nc] - Xm[Nc]) == 0)) return Xm[Nc];

	//////////////////////////////////////////////////////////////////////////////////////////////




	Xm[Nc] = fP[Ix + 1][Iy + 1][Iz + 1][2];
	for (i = 0; i < Nc; i++) {
		NRKi = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));		//Wilson Equation
		Xm[i] = fcomp[Ix][Iy][Iz][i][1] / NRKi;
	}

	nIter = 0;
	do {
		SXm = 0;
		for (i = 0; i < Nc; i++) {
			//if (Xm[i] < 0) Xm[i] = 0;
			//else if (Xm[i] > 1) Xm[i] = 1;
			SXm += Xm[i];
			satJac[Nc][i] = 1;
		}

		//if (Xm[Nc] < 0) Xm[Nc] = PZERO;

		satAns[Nc] = -SXm + 1;
		satJac[Nc][Nc] = 0;


		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += Xm[i] * fluidProp[i][EOS_B];
			Bg += fcomp[Ix][Iy][Iz][i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += Xm[i] * Xm[j] * tempSQR;
				Ag += fcomp[Ix][Iy][Iz][i][1] * fcomp[Ix][Iy][Iz][j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bg *= Xm[Nc] / (RGAS*resTemp);

		dAdpl = Cal / (RGAS*RGAS*resTemp*resTemp);
		dBdpl = Cbl / (RGAS*resTemp);
		dAdpg = Cag / (RGAS*RGAS*resTemp*resTemp);
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

			d1 = 1;
			d2 = 0;
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

			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}

		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl*Ul - Ul - Wl*Wl, -(Al*Bl - Bl*Wl*Wl - Wl*Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg*Ug - Ug - Wg*Wg, -(Ag*Bg - Bg*Wg*Wg - Wg*Wg), 'g');



		dZdpg = -(dAdpg*(Zg - Bg) + dBdpg*(-Zg*Zg - Ug*Zg - Ag + Wg*Wg) + dUdpg*(Zg*Zg - Bg*Zg - Zg) + dWdpg*(-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg*Ug - Ug - Wg*Wg);		//Modified
		dZdpl = -(dAdpl*(Zl - Bl) + dBdpl*(-Zl*Zl - Ul*Zl - Al + Wl*Wl) + dUdpl*(Zl*Zl - Bl*Zl - Zl) + dWdpl*(-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl*Ul - Ul - Wl*Wl);		//Modified

		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += Xm[j] * tempSQR;
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

			satAns[i] = -Xm[i] * phiL + fcomp[Ix][Iy][Iz][i][1] * phiG;

			dAAdpl = fluidProp[i][EOS_B] * dZdpl / Cbl;
			dAAdpg = fluidProp[i][EOS_B] * dZdpg / Cbg;

			dBBdpl = (dBdpl - dZdpl) / (Zl - Bl);
			dBBdpg = (dBdpg - dZdpg) / (Zg - Bg);

			dDDdpl = (Bl*dAdpl - Al*dBdpl) / (Bl*Bl*(d1 - d2));
			dDDdpg = (Bg*dAdpg - Ag*dBdpg) / (Bg*Bg*(d1 - d2));

			dFFdpl = (dZdpl + d2*dBdpl) / (Zl + d2*Bl) - (dZdpl + d1*dBdpl) / (Zl + d1*Bl);
			dFFdpg = (dZdpg + d2*dBdpg) / (Zg + d2*Bg) - (dZdpg + d1*dBdpg) / (Zg + d1*Bg);

			dPhidpl = phiL*(dAAdpl + dBBdpl + EEl*(dDDdpl*FFl + DDl*dFFdpl));
			dPhidpg = phiG*(dAAdpg + dBBdpg + EEg*(dDDdpg*FFg + DDg*dFFdpg));

			satJac[i][Nc] = Xm[i] * dPhidpl - fcomp[Ix][Iy][Iz][i][1] * dPhidpg;


			//i: Row of Jac
			//k: Col of Jac
			for (k = 0; k < Nc; k++) {
				dAdxl = 0;
				for (n = 0; n < Nc; n++) {		//SSFF
					dAdxl += Xm[n] * sqrt(fluidProp[n][EOS_A] * fluidProp[k][EOS_A])*bic[n][k];
					//This loop can be reduced
				}
				dAdxl *= 2 * Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
				dBdxl = fluidProp[k][EOS_B] * Xm[Nc] / (RGAS*resTemp);

				if (SRK) {
					dUdxl = dBdxl;
					dWdxl = 0;
				}
				else if (PR) {
					dUdxl = 2 * dBdxl;
					dWdxl = dBdxl;
				}



				dZdxl = -(dAdxl*(Zl - Bl) + dBdxl*(-Zl*Zl - Ul*Zl - Al + Wl*Wl) + dUdxl*(Zl*Zl - Bl*Zl - Zl) + dWdxl*(-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl*Ul - Ul - Wl*Wl);		//Modified

				dAAdxl = -fluidProp[k][EOS_B] * fluidProp[i][EOS_B] * (Zl - 1) / (Cbl*Cbl) + fluidProp[i][EOS_B] * dZdxl / Cbl;

				dBBdxl = (dBdxl - dZdxl) / (Zl - Bl);

				dDDdxl = (dAdxl*Bl - dBdxl*Al) / (Bl*Bl*(d1 - d2));

				dFFdxl = (dZdxl + d2*dBdxl) / (Zl + d2*Bl) - (dZdxl + d1*dBdxl) / (Zl + d1*Bl);

				dEEdxl = 2 * (bic[i][k] * sqrt(fluidProp[i][EOS_A] * fluidProp[k][EOS_A]) / Cal - dAdxl*RGAS*RGAS*resTemp*resTemp*tl / (Xm[Nc] * Cal*Cal)) + fluidProp[i][EOS_B] * fluidProp[k][EOS_B] / (Cbl*Cbl);

				dPhidxl = phiL*(dAAdxl + dBBdxl + dDDdxl*EEl*FFl + DDl*dEEdxl*FFl + DDl*EEl*dFFdxl);

				if (i == k) {
					satJac[i][k] = phiL + Xm[i] * dPhidxl;
				}
				else {
					satJac[i][k] = Xm[i] * dPhidxl;
				}


			}
		}

		//Scaling(satJac, satAns, Nc+1);
		//biCGStab(satJac, satAns, Xms, Nc+1);
		//GaussJ(satJac, satAns, Xms, Nc+1);
		MKLS(satJac, satAns, Xms, Nc + 1);         //Sadeq


		XmTol = 0;
		xSum = 0;
		for (i = 0; i < (Nc + 1); i++) {
			Xm[i] += Xms[i];
			if (i < Nc) {
				if (Xm[i] < 0) Xm[i] = 0;
				else if (Xm[i] > 1) Xm[i] = 1;
				xSum += Xm[i];
			}
			else {
				if (Xm[Nc] < 0) Xm[Nc] = PZERO;
			}
			if (Xm[i]) XmTol += fabs(Xms[i] / Xm[i]);
			else XmTol += fabs(Xms[i]);
		}
		XmTol /= (Nc + 1);
		for (i = 0; i < Nc; i++) Xm[i] /= xSum;

		//XmTol=In_Product(Xms, Xms, Nc+1);

		nIter++;
		if (nIter > MAX_MNR_ITER) {
			return -1;
		}

	} while (XmTol > XMTOL);// && (fabs(Xms[Nc]) > PRESSTOL));

	sumK = 0;
	k = 0;
	for (i = 0; i < Nc; i++) {
		if (fcomp[Ix][Iy][Iz][i][1]) {
			sumK += fabs(1 - fcomp[Ix][Iy][Iz][i][1] / Xm[i]);
			k++;
		}
	}
	if ((sumK / k) < TRIV_K) return -1;
	else return Xm[Nc];

}

FType LowerDew_Succ(int Ix, int Iy, int Iz) {
	register int i, j, k, nIter;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG;
	FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType SXm;
	FType tempSQR;
	FType sumK;

	Xm[Nc] = P[Ix + 1][Iy + 1][Iz + 1][2];

	for (i = 0; i < Nc; i++) {
		Ki[i] = (fluidProp[i][PCRIT] / Xm[Nc])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));
	}

	nIter = 0;
	do {
		for (i = 0; i < Nc; i++) {
			Xm[i] = comp[Ix][Iy][Iz][i][1] / Ki[i];
			if (Xm[i] < 0) Xm[i] = 0;
			else if (Xm[i] > 1) Xm[i] = 1;
		}
		if (Xm[Nc] < 0) Xm[Nc] = 1;

		/////////////////////////////////////////////////////////////

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += Xm[i] * fluidProp[i][EOS_B];
			Bg += comp[Ix][Iy][Iz][i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += Xm[i] * Xm[j] * tempSQR;
				Ag += comp[Ix][Iy][Iz][i][1] * comp[Ix][Iy][Iz][j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);
		Bl *= Xm[Nc] / (RGAS*resTemp);
		Ag *= Xm[Nc] / (RGAS*RGAS*resTemp*resTemp);		//gas-oil
		Bg *= Xm[Nc] / (RGAS*resTemp);


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


		SXm = 0;
		k = 0;
		sumK = 0;
		for (i = 0; i < Nc; i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += Xm[j] * tempSQR;
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

			fL = Xm[i] * phiL;
			fG = comp[Ix][Iy][Iz][i][1] * phiG;	//gas-oil

			if (comp[Ix][Iy][Iz][i][1]) {
				Ki[i] *= fL / fG;
				SXm += comp[Ix][Iy][Iz][i][1] / Ki[i];
				sumK += fabs(1 - Ki[i]);
				k++;
			}
		}

		Xm[Nc] /= SXm;
		nIter++;
		if (nIter > MAX_MNR_ITER) break;
	} while ((fabs(SXm - 1) / k) > SXMTOL);

	if ((sumK / k) < TRIV_K) Xm[Nc] = 1e100;
	return Xm[Nc];
}

void critTest(int Ix, int Iy, int Iz) {
	register int i, j, nIter;
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
	char chStat;

	for (i = 0; i < Nc; i++) {
		chStat = phaseStat[Ix][Iy][Iz];
		Ki[i] = (fluidProp[i][PCRIT] / P[Ix + 1][Iy + 1][Iz + 1][1])*exp(5.37*(1 + fluidProp[i][AC])*(1 - fluidProp[i][TCRIT] / resTemp));
		if (chStat == 1) STcomp[i][2] = comp[Ix][Iy][Iz][i][0];
		if (chStat == (-1)) STcomp[i][2] = comp[Ix][Iy][Iz][i][1];
	}

	nIter = 0;
	do {
		L = 0.5;
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

			FFl = log((Zl + d2*Bl) / (Zl + d1*Bl));
			FFg = log((Zg + d2*Bg) / (Zg + d1*Bg));

			phiL = exp(t1l*(Zl - 1) - log(Zl - Bl) + DDl*EEl*FFl);
			phiG = exp(t1g*(Zg - 1) - log(Zg - Bg) + DDg*EEg*FFg);

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
		nIter++;
		if (nIter > MAX_MNR_ITER) break;
	} while (ESum > FLASH_TOL);
	//////////////////////////////////////////////////////////////

	if ((fabs(L) < DZERO) && (chStat == 1)) {
		if ((phaseStat[Ix][Iy][Iz] == 0) || (phaseStat[Ix][Iy][Iz] == 1)) buildPreconFlag = -1;
		phaseStat[Ix][Iy][Iz] = -1;
		for (i = 0; i < Nc; i++) {
			//comp[Ix][Iy][Iz][i][0]=comp[Ix][Iy][Iz][i][1];
			comp[Ix][Iy][Iz][i][1] = STcomp[i][2];
		}
		sat[Ix][Iy][Iz][2] = sat[Ix][Iy][Iz][1];
		sat[Ix][Iy][Iz][1] = 0;
	}
	else if ((fabs(1 - L) < DZERO) && (chStat == (-1))) {
		if ((phaseStat[Ix][Iy][Iz] == 0) || (phaseStat[Ix][Iy][Iz] == (-1))) buildPreconFlag = -1;
		phaseStat[Ix][Iy][Iz] = 1;
		for (i = 0; i < Nc; i++) {
			comp[Ix][Iy][Iz][i][0] = STcomp[i][2];

		}
		sat[Ix][Iy][Iz][1] = sat[Ix][Iy][Iz][2];
		sat[Ix][Iy][Iz][2] = 0;
	}
}

void Equilibrium_Concentrations(int Ix, int Iy, int Iz) {
	bool fflag;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG, dAdxg, dBdxg, dZdxg;
	FType dAAdxg, dBBdxg, DDl, DDg, EEl, EEg, FFl, FFg;
	FType dFFdxg, dEEdxg, dDDdxg, dPhidxg;
	FType dAdpg, dBdpg, dZdpg;
	//FType dAAdpg, dBBdpg, dDDdpg, dFFdpg, dPhidpg;
	FType dUdpg, dWdpg;
	FType dUdxg, dWdxg;
	register int i, j, k, n;

	FType mLf;
	//FType SXm;
	//FType XmTol;

	FType dAdpl, dBdpl, dZdpl;
	//FType dAAdpl, dBBdpl, dDDdpl, dFFdpl, dPhidpl;
	FType dUdpl, dWdpl;
	FType tempSQR;
	FType dAAdxl, dBBdxl, dFFdxl, dEEdxl, dDDdxl, dPhidxl;
	FType dAdxl, dBdxl, dZdxl;
	FType dUdxl, dWdxl;

	FType LCoefM, LCoefF, GCoefM, GCoefF, GasP, OilP, MeanP;
	FType Lm, Norm2, porpor;
	FType **MinorJac, *MinorAns, *XXs;
	FType sumx, sumy;
	register unsigned int nIter = 0;
	//FType MinorAnsC;


	//FType fL, fG, LV, F, dF, tempD, ESum, LL;	
	//FType Lb = 10;
	//FType **CEdgeComp;

	/*if ((CEdgeComp = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for EdgeComp");
	for (i = 0; i < Nc; i++) {
		if ((CEdgeComp[i] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for EdgeComp");
	}*/
	if ((MinorJac = (FType **)malloc(2*(Nc) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for MinorJac");
	for (i = 0; i < (2*Nc); i++) {
		if ((MinorJac[i] = (FType *)malloc(2*(Nc) * sizeof(FType))) == NULL) TerM("Can not allocate memory for MinorJac");
		for (j = 0; j < (2*Nc); j++) MinorJac[i][j] = 0;
	}
	if ((MinorAns = (FType *)malloc((2*Nc) * sizeof(FType))) == NULL) TerM("Can not allocate memory for MinorAns");
	if ((XXs = (FType *)malloc((2*Nc) * sizeof(FType))) == NULL) TerM("Can not allocate memory for XXs");

	if (phaseStat[Ix][Iy][Iz] != -1) LCoefM = 0.5;
	else LCoefM = 0;
	if (fphaseStat[Ix][Iy][Iz] != -1) LCoefF = 0.5;
	else LCoefF = 0;
	if (phaseStat[Ix][Iy][Iz] != 1) GCoefM = 0.5;
	else GCoefM = 0;
	if (fphaseStat[Ix][Iy][Iz] != 1) GCoefF = 0.5;
	else GCoefF = 0;

	for (i = 0; i < Nc; i++) {
		if (phaseStat[Ix][Iy][Iz] == -1 && fphaseStat[Ix][Iy][Iz] == -1) EdgeComp[i][0] = 0;
		else EdgeComp[i][0] = (LCoefM*comp[Ix][Iy][Iz][i][0] + LCoefF*fcomp[Ix][Iy][Iz][i][0]) / (LCoefM + LCoefF);
		if (phaseStat[Ix][Iy][Iz] == 1 && fphaseStat[Ix][Iy][Iz] == 1) EdgeComp[i][1] = 0;
		else EdgeComp[i][1] = (GCoefM*comp[Ix][Iy][Iz][i][1] + GCoefF*fcomp[Ix][Iy][Iz][i][1]) / (GCoefM + GCoefF);
		//EdgeComp[i][0] = 0.5;
		//EdgeComp[i][1] = 0.5;
	}
	sumx = 0;
	sumy = 0;
	for (i = 0; i < (Nc); i++) {		
		if (EdgeComp[i][0] < 0) EdgeComp[i][0] = 0;
		if (EdgeComp[i][0] > 1) EdgeComp[i][0] = 1;		
		sumx += EdgeComp[i][0];		
		if (EdgeComp[i][1] > 1) EdgeComp[i][1] = 1;
		else if (EdgeComp[i][1] < 0) EdgeComp[i][1] = 0;		
		sumy += EdgeComp[i][1];
	}
	/*
	if (phaseStat[Ix][Iy][Iz]==(-1)) {
		GasP = P[Ix + 1][Iy + 1][Iz + 1][2];
		OilP = fP[Ix + 1][Iy + 1][Iz + 1][1];
	}
	else if (phaseStat[Ix][Iy][Iz] == 0) {
		if (fphaseStat[Ix][Iy][Iz] == (-1)) {
			GasP = fP[Ix + 1][Iy + 1][Iz + 1][2];
			OilP = P[Ix + 1][Iy + 1][Iz + 1][1];
		}
		else if (fphaseStat[Ix][Iy][Iz] == (1)) {
			GasP = P[Ix + 1][Iy + 1][Iz + 1][2];
			OilP = fP[Ix + 1][Iy + 1][Iz + 1][1];
		}
	}
	else {
		GasP = fP[Ix + 1][Iy + 1][Iz + 1][2];
		OilP = P[Ix + 1][Iy + 1][Iz + 1][1];
	}*/
	
	MeanP = (P[Ix + 1][Iy + 1][Iz + 1][1] + fP[Ix + 1][Iy + 1][Iz + 1][1]) / 2;
	porpor = porosity[Ix][Iy][Iz] * (1 + cpor * (P[Ix + 1][Iy + 1][Iz + 1][1] - refP) + dcpor * (P[Ix + 1][Iy + 1][Iz + 1][1] - refP)*(P[Ix + 1][Iy + 1][Iz + 1][1] - refP));

	OilP = MeanP;
	GasP = MeanP;


	Lm = ShapeFactor_1(Ix, Iy, Iz) *porpor / tor[Ix][Iy][Iz];
	mLf =ShapeFactor_3(Ix, Iy, Iz);

	fflag = false;
	do {

		for (i = 0; i < (Nc); i++) {
			EdgeComp[i][0] /= sumx;
			EdgeComp[i][1] /= sumy;
		}

		Al = 0;
		Bl = 0;
		Ag = 0;
		Bg = 0;
		for (i = 0; i < Nc; i++) {
			Bl += EdgeComp[i][0] * fluidProp[i][EOS_B];
			Bg += EdgeComp[i][1] * fluidProp[i][EOS_B];
			for (j = 0; j < Nc; j++) {
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				Al += EdgeComp[i][0] * EdgeComp[j][0] * tempSQR;
				Ag += EdgeComp[i][1] * EdgeComp[j][1] * tempSQR;
			}
		}

		Cal = Al;
		Cbl = Bl;
		Cag = Ag;
		Cbg = Bg;
		Al *= OilP / (RGAS*RGAS*resTemp*resTemp);
		Bl *= OilP / (RGAS*resTemp);
		Ag *= GasP / (RGAS*RGAS*resTemp*resTemp);
		Bg *= GasP / (RGAS*resTemp);

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

			d1 = 1;
			d2 = 0;
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

			d1 = 1 + SQRT2;
			d2 = 1 - SQRT2;
		}


		Zl = Solve_Z(-(1 + Bl - Ul), Al - Bl * Ul - Ul - Wl * Wl, -(Al*Bl - Bl * Wl*Wl - Wl * Wl), 'l');
		Zg = Solve_Z(-(1 + Bg - Ug), Ag - Bg * Ug - Ug - Wg * Wg, -(Ag*Bg - Bg * Wg*Wg - Wg * Wg), 'g');

		dZdpg = -(dAdpg*(Zg - Bg) + dBdpg * (-Zg * Zg - Ug * Zg - Ag + Wg * Wg) + dUdpg * (Zg*Zg - Bg * Zg - Zg) + dWdpg * (-2 * Wg*Zg + 2 * Bg*Wg + 2 * Wg)) / (3 * Zg*Zg - 2 * Zg*(1 + Bg - Ug) + Ag - Bg * Ug - Ug - Wg * Wg);		//Modified
		dZdpl = -(dAdpl*(Zl - Bl) + dBdpl * (-Zl * Zl - Ul * Zl - Al + Wl * Wl) + dUdpl * (Zl*Zl - Bl * Zl - Zl) + dWdpl * (-2 * Wl*Zl + 2 * Bl*Wl + 2 * Wl)) / (3 * Zl*Zl - 2 * Zl*(1 + Bl - Ul) + Al - Bl * Ul - Ul - Wl * Wl);		//Modified


		for (i = 0; i < (Nc - 1); i++) {
			tl = 0;
			tg = 0;
			for (j = 0; j < Nc; j++) {	//FFSS
				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[j][EOS_A])*bic[i][j];
				tl += EdgeComp[j][0] * tempSQR;
				tg += EdgeComp[j][1] * tempSQR;
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

			MinorAns[i] = -EdgeComp[i][0] * OilP * phiL + EdgeComp[i][1] * GasP * phiG;

			//i: Row of Jac
			//k: Col of Jac
			for (k = 0; k < Nc; k++) {		//for all x
				dAdxl = 0;
				dAdxg = 0;
				for (n = 0; n < Nc; n++) {		//SSFF
					tempSQR = sqrt(fluidProp[k][EOS_A] * fluidProp[n][EOS_A])*bic[k][n];
					dAdxl += EdgeComp[n][0] * tempSQR;
					dAdxg += EdgeComp[n][1] * tempSQR;
					//This loop can be reduced
				}
				dAdxl *= 2 * OilP / (RGAS*RGAS*resTemp*resTemp);
				dAdxg *= 2 * GasP / (RGAS*RGAS*resTemp*resTemp);
				dBdxl = fluidProp[k][EOS_B] * OilP / (RGAS*resTemp);
				dBdxg = fluidProp[k][EOS_B] * GasP / (RGAS*resTemp);

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

				dAAdxl = -fluidProp[k][EOS_B] * fluidProp[i][EOS_B] * (Zl - 1) / (Cbl*Cbl) + fluidProp[i][EOS_B] * dZdxl / Cbl;
				dAAdxg = -fluidProp[k][EOS_B] * fluidProp[i][EOS_B] * (Zg - 1) / (Cbg*Cbg) + fluidProp[i][EOS_B] * dZdxg / Cbg;

				dBBdxl = (dBdxl - dZdxl) / (Zl - Bl);
				dBBdxg = (dBdxg - dZdxg) / (Zg - Bg);

				dDDdxl = (dAdxl*Bl - dBdxl * Al) / (Bl*Bl*(d1 - d2));
				dDDdxg = (dAdxg*Bg - dBdxg * Ag) / (Bg*Bg*(d1 - d2));

				dFFdxl = (dZdxl + d2 * dBdxl) / (Zl + d2 * Bl) - (dZdxl + d1 * dBdxl) / (Zl + d1 * Bl);
				dFFdxg = (dZdxg + d2 * dBdxg) / (Zg + d2 * Bg) - (dZdxg + d1 * dBdxg) / (Zg + d1 * Bg);

				tempSQR = sqrt(fluidProp[i][EOS_A] * fluidProp[k][EOS_A])*bic[i][k];
				dEEdxl = 2 * (tempSQR / Cal - dAdxl * RGAS*RGAS*resTemp*resTemp*tl / (OilP * Cal*Cal)) + fluidProp[i][EOS_B] * fluidProp[k][EOS_B] / (Cbl*Cbl);
				dEEdxg = 2 * (tempSQR / Cag - dAdxg * RGAS*RGAS*resTemp*resTemp*tg / (GasP * Cag*Cag)) + fluidProp[i][EOS_B] * fluidProp[k][EOS_B] / (Cbg*Cbg);

				dPhidxl = phiL * (dAAdxl + dBBdxl + dDDdxl * EEl*FFl + DDl * dEEdxl*FFl + DDl * EEl*dFFdxl);
				dPhidxg = phiG * (dAAdxg + dBBdxg + dDDdxg * EEg*FFg + DDg * dEEdxg*FFg + DDg * EEg*dFFdxg);

				if (i == k) {
					MinorJac[i][k] = OilP * phiL + OilP * EdgeComp[i][0] * dPhidxl;

					MinorJac[i][k + Nc] = -(GasP * phiG + GasP * EdgeComp[i][1] * dPhidxg);
				}
				else {
					MinorJac[i][k] = OilP * EdgeComp[i][0] * dPhidxl;

					MinorJac[i][k + Nc] = -GasP * EdgeComp[i][1] * dPhidxg;
				}
			}
		}

		if (((phaseStat[Ix][Iy][Iz] == 1) && (fphaseStat[Ix][Iy][Iz] != 1)) || ((phaseStat[Ix][Iy][Iz] == 0) && (fphaseStat[Ix][Iy][Iz] == (-1)))) {
			for (i = 0; i < (Nc - 1); i++) {
				MinorAns[Nc + i - 1] = sat[Ix][Iy][Iz][1] * diffusion[Ix][Iy][Iz][i][0] * blockFProps[Ix][Iy][Iz][RO][0] * (comp[Ix][Iy][Iz][i][0] - EdgeComp[i][0])*Lm + fsat[Ix][Iy][Iz][2] * fdiffusion[Ix][Iy][Iz][i][1] * fblockFProps[Ix][Iy][Iz][RO][1] * (fcomp[Ix][Iy][Iz][i][1] - EdgeComp[i][1])*mLf;
				MinorJac[Nc + i - 1][i] = sat[Ix][Iy][Iz][1] * diffusion[Ix][Iy][Iz][i][0] * blockFProps[Ix][Iy][Iz][RO][0] * Lm;
				MinorJac[Nc + i - 1][i + Nc] = fsat[Ix][Iy][Iz][2] * fdiffusion[Ix][Iy][Iz][i][1] * fblockFProps[Ix][Iy][Iz][RO][1] * mLf;

			}
		}
		else if (((phaseStat[Ix][Iy][Iz] == (-1)) && (fphaseStat[Ix][Iy][Iz] != (-1))) || ((phaseStat[Ix][Iy][Iz] == 0) && (fphaseStat[Ix][Iy][Iz] == 1))) {
			for (i = 0; i < (Nc - 1); i++) {
				MinorAns[Nc + i - 1] = sat[Ix][Iy][Iz][2] * diffusion[Ix][Iy][Iz][i][1] * blockFProps[Ix][Iy][Iz][RO][1] * (comp[Ix][Iy][Iz][i][1] - EdgeComp[i][1])*Lm + fsat[Ix][Iy][Iz][1] * fdiffusion[Ix][Iy][Iz][i][0] * fblockFProps[Ix][Iy][Iz][RO][0] * (fcomp[Ix][Iy][Iz][i][0] - EdgeComp[i][0])*mLf;
				MinorJac[Nc + i - 1][i] = fsat[Ix][Iy][Iz][1] * fdiffusion[Ix][Iy][Iz][i][0] * fblockFProps[Ix][Iy][Iz][RO][0] * mLf;
				MinorJac[Nc + i - 1][i + Nc] = sat[Ix][Iy][Iz][2] * diffusion[Ix][Iy][Iz][i][1] * blockFProps[Ix][Iy][Iz][RO][1] * Lm;
			}
		}
		/*if (((phaseStat[Ix][Iy][Iz] == 1) && (fphaseStat[Ix][Iy][Iz] != 1)) || ((phaseStat[Ix][Iy][Iz] == 0) && (fphaseStat[Ix][Iy][Iz] == (-1)))) {
			i = Nc - 1;
				MinorAnsC = sat[Ix][Iy][Iz][1] * diffusion[Ix][Iy][Iz][i][0] * (comp[Ix][Iy][Iz][i][0] - EdgeComp[i][0]) / Lm + fsat[Ix][Iy][Iz][2] * fdiffusion[Ix][Iy][Iz][i][1] * (fcomp[Ix][Iy][Iz][i][1] - EdgeComp[i][1]) / mLf;
				//MinorJac[Nc + i][i] = sat[Ix][Iy][Iz][1] * diffusion[Ix][Iy][Iz][i][0] / Lm;
				//MinorJac[Nc + i][i + Nc] = fsat[Ix][Iy][Iz][2] * fdiffusion[Ix][Iy][Iz][i][1] / mLf;


		}
		else if (((phaseStat[Ix][Iy][Iz] == (-1)) && (fphaseStat[Ix][Iy][Iz] != (-1))) || ((phaseStat[Ix][Iy][Iz] == 0) && (fphaseStat[Ix][Iy][Iz] == 1))) {
			i = Nc - 1;
				MinorAnsC = sat[Ix][Iy][Iz][2] * diffusion[Ix][Iy][Iz][i][1] * (comp[Ix][Iy][Iz][i][1] - EdgeComp[i][1]) / Lm + fsat[Ix][Iy][Iz][1] * fdiffusion[Ix][Iy][Iz][i][0] * (fcomp[Ix][Iy][Iz][i][0] - EdgeComp[i][0]) / mLf;
				//MinorJac[Nc + i][i] = fsat[Ix][Iy][Iz][1] * fdiffusion[Ix][Iy][Iz][i][0] / mLf;
				//MinorJac[Nc + i][i + Nc] = sat[Ix][Iy][Iz][2] * diffusion[Ix][Iy][Iz][i][1] / Lm;

		}*/
		sumx = 0;
		sumy = 0;
		for (i = 0; i < Nc; i++) {
			MinorJac[2 * Nc - 2][i] = 1;
			MinorJac[2 * Nc - 1][i + Nc] = 1;
			sumx += EdgeComp[i][0];
			sumy += EdgeComp[i][1];
		}
		MinorAns[2 * Nc - 2] = 1 - sumx;
		MinorAns[2 * Nc - 1] = 1 - sumy;


		MKLS(MinorJac, MinorAns, XXs, 2 * Nc);


		Norm2 = 0;
		sumx = 0;
		sumy = 0;
		for (i = 0; i < (Nc); i++) {
			EdgeComp[i][0] += XXs[i];
			if (EdgeComp[i][0] < 0) EdgeComp[i][0] = 0;
			else if (EdgeComp[i][0] > 1) EdgeComp[i][0] = 1;

			if (EdgeComp[i][0]) Norm2 += fabs(XXs[i] / EdgeComp[i][0]);
			else Norm2 += fabs(XXs[i]);
			//Norm2 += XXs[i] * XXs[i];
			sumx += EdgeComp[i][0];
			EdgeComp[i][1] += XXs[i + Nc];
			if (EdgeComp[i][1] > 1) EdgeComp[i][1] = 1;
			else if (EdgeComp[i][1] < 0) EdgeComp[i][1] = 0;

			if (EdgeComp[i][1]) Norm2 += fabs(XXs[i + Nc] / EdgeComp[i][1]);
			else Norm2 += fabs(XXs[i + Nc]);
			//Norm2 += XXs[i + Nc] * XXs[i + Nc];
			sumy += EdgeComp[i][1];
		}
		//if (sumx) Norm2 /= (2 * Nc*sumx);
		Norm2 /= (2 * Nc);
		//for (i = 0; i < (Nc); i++) {
		//	EdgeComp[i][0] /= sumx;
		//	EdgeComp[i][1] /= sumy;				
		//}
		nIter++;
		if (nIter > MAX_EQL_ITER) {

			//fflag = true;
			break;
		}



	} while (((Norm2 - Norm2) != 0) || (Norm2 > XMTOLEDGE));

	for (i = 0; i < (Nc); i++) {
		EdgeComp[i][0] /= sumx;
		EdgeComp[i][1] /= sumy;
	}
		
	for (i = 0; i < Nc; i++) {
		if (EdgeComp[i][0] != EdgeComp[i][0]) {
			fflag = true;
			break;
		}
		if (EdgeComp[i][1] != EdgeComp[i][1]) {
			fflag = true;
			break;
		}
	}
	if (fflag) {

		for (i = 0; i < Nc; i++) {
			if (phaseStat[Ix][Iy][Iz] == -1 && fphaseStat[Ix][Iy][Iz] == -1) EdgeComp[i][0] = 0;
			else EdgeComp[i][0] = (LCoefM*comp[Ix][Iy][Iz][i][0] + LCoefF * fcomp[Ix][Iy][Iz][i][0]) / (LCoefM + LCoefF);
			if (phaseStat[Ix][Iy][Iz] == 1 && fphaseStat[Ix][Iy][Iz] == 1) EdgeComp[i][1] = 0;
			else EdgeComp[i][1] = (GCoefM*comp[Ix][Iy][Iz][i][1] + GCoefF * fcomp[Ix][Iy][Iz][i][1]) / (GCoefM + GCoefF);
		}

		std::cout << ("Equilibrium conds failed\n");
	}

	for (i = 0; i < 2 * Nc; i++) {
		free(MinorJac[i]);
	}
	free(MinorJac);

	free(MinorAns);
	free(XXs);
		//return 0;
}

FType Calc_Re(int Ix, int Iy, int Iz, int i) {

	FType Fraction_well;// fraction of well associated with the  well block
	FType *TransWell;//interface transmissibility factor between Block i and the well block
	FType *aj;//distance from well to its Image Well j
	FType **rij;//distance from Grid Point i to Well j or its image
	FType *ri_1;////distance from Grid Point i to Well 
	FType bSum = 0;  // sum of transmissibility
	FType X_posi, Y_posi;// position of well in block a & b
	FType iZarb = 1, jZarb = 1;
	int  j;
	FType Re;

	if ((rij = (FType **)malloc(5 * sizeof(FType *))) == NULL) TerM("Can not allocate memory for rij");
	for (int j = 0; j < 5; j++) {
		if ((rij[j] = (FType *)malloc(9 * sizeof(FType))) == NULL) TerM("Can not allocate memory for rij");
	}
	if ((TransWell = (FType *)malloc((9) * sizeof(FType))) == NULL) TerM("Can not allocate memory for TransWell");
	if ((aj = (FType *)malloc((5) * sizeof(FType))) == NULL) TerM("Can not allocate memory for aj");
	if ((ri_1 = (FType *)malloc((9) * sizeof(FType))) == NULL) TerM("Can not allocate memory for ri_1");


	X_posi = 0;// gridDim[Ix] / 2;
	Y_posi = 0;// gridDim[Nx + Iy] / 2;

	ri_1[1] = sqrt((X_posi*X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
	ri_1[2] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (Y_posi* Y_posi));
	ri_1[3] = sqrt((X_posi*X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
	ri_1[4] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (Y_posi* Y_posi));
	ri_1[5] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
	ri_1[6] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
	ri_1[7] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
	ri_1[8] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));

	TransWell[5] = 0;
	TransWell[6] = 0;
	TransWell[7] = 0;
	TransWell[8] = 0;
	TransWell[1] = (gridDim[Ix] / gridDim[Nx + Iy]) - (TransWell[5] + TransWell[6]);
	TransWell[2] = (gridDim[Nx + Iy] / gridDim[Ix]) - (TransWell[6] + TransWell[7]);
	TransWell[3] = (gridDim[Ix] / gridDim[Nx + Iy]) - (TransWell[7] + TransWell[8]);
	TransWell[4] = (gridDim[Nx + Iy] / gridDim[Ix]) - (TransWell[5] + TransWell[8]);



	//i = 0;
	if ((Ix == welli[i][1]) && (Iy == welli[i][2]) && (Iz >= welli[i][3]) && (Iz <= welli[i][4])) {

		if ((0 < welli[i][1]) && (welli[i][1] < (Nx - 1)) && (0 < welli[i][2]) && (welli[i][2] < (Ny - 1))) {//Interior Block
			Fraction_well = 1;
			for (j = 1; j < 5; j++)
				iZarb *= pow(ri_1[j], TransWell[j]);
			for (j = 1; j < 9; j++)
				bSum += TransWell[j];
			bSum = 1 / bSum;
			Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);
		}


		if ((Iy == Ny - 1) && (Ix != 0 && Ix != Nx - 1)) {
			if (Y_posi > (-gridDim[Nx + Iy] / 2) && Y_posi <= (gridDim[Nx + Iy] / 2)) {//Block on one Reservoir Boundary Iy == Ny - 1
				Fraction_well = 1;
				aj[2] = gridDim[Nx + Iy] + 2 * Y_posi;
				rij[2][2] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[2][7] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));
				rij[2][3] = sqrt((X_posi*X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));
				rij[2][8] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));
				rij[2][4] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[2][1] = 1;
				rij[2][5] = 1;
				rij[2][6] = 1;

				ri_1[5] = 1;
				ri_1[1] = 1;
				ri_1[6] = 1;

				TransWell[5] = 0;
				TransWell[1] = 0;
				TransWell[6] = 0;

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j])*pow((rij[2][j] / aj[2]), TransWell[j]);

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);
			}


			if (Y_posi == (-gridDim[Nx + Iy] / 2)) {//Well on Reservoir Boundary Iy == Ny - 1

				Fraction_well = 0.5;
				ri_1[5] = 1;
				ri_1[1] = 1;
				ri_1[6] = 1;

				TransWell[5] = 0;
				TransWell[1] = 0;
				TransWell[6] = 0;

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j]);
				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;
				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);
			}
		}


		if ((Iy == 0) && (Ix != 0 && Ix != Nx - 1)) {
			if (Y_posi < (gridDim[Nx + Iy] / 2) && Y_posi >= (-gridDim[Nx + Iy] / 2)) {//Block on one Reservoir Boundary Iy == 0
				Fraction_well = 1;
				aj[2] = gridDim[Nx + Iy] - 2 * Y_posi;
				rij[2][2] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[2][6] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));
				rij[2][1] = sqrt((X_posi*X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));
				rij[2][5] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));
				rij[2][4] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[2][3] = 1;
				rij[2][7] = 1;
				rij[2][8] = 1;

				ri_1[3] = 1;
				ri_1[7] = 1;
				ri_1[8] = 1;

				TransWell[3] = 0;
				TransWell[7] = 0;
				TransWell[8] = 0;

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j])*pow((rij[2][j] / aj[2]), TransWell[j]);



				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);

			}
			if (Y_posi == (gridDim[Nx + Iy] / 2)) {//Well on Reservoir Boundary Iy == 0

				Fraction_well = 0.5;
				ri_1[3] = 1;
				ri_1[7] = 1;
				ri_1[8] = 1;

				TransWell[3] = 0;
				TransWell[7] = 0;
				TransWell[8] = 0;

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j]);
				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;
				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);


			}
		}

		if ((Ix == Nx - 1) && (Iy != 0 && Iy != Ny - 1)) {
			if (X_posi < (gridDim[Ix] / 2) && X_posi >= (-gridDim[Ix] / 2)) {//Block on one Reservoir Boundary Ix == Nx - 1
				Fraction_well = 1;
				aj[2] = gridDim[Ix] - 2 * X_posi;
				rij[2][1] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[2][5] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[2][4] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (Y_posi*Y_posi));
				rij[2][8] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[2][3] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[2][2] = 1;
				rij[2][7] = 1;
				rij[2][6] = 1;

				ri_1[2] = 1;
				ri_1[7] = 1;
				ri_1[6] = 1;

				TransWell[2] = 0;
				TransWell[7] = 0;
				TransWell[6] = 0;

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j])*pow((rij[2][j] / aj[2]), TransWell[j]);

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);

			}
			if (X_posi == (gridDim[Ix] / 2)) {//Well on Reservoir Boundary Ix == Nx - 1
				Fraction_well = 0.5;
				ri_1[2] = 1;
				ri_1[7] = 1;
				ri_1[6] = 1;

				TransWell[2] = 0;
				TransWell[7] = 0;
				TransWell[6] = 0;

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j]);
				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;
				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);

			}
		}


		if ((Ix == 0) && (Iy != 0 && Iy != Ny - 1)) {
			if (X_posi > (-gridDim[Ix] / 2) && X_posi <= (gridDim[Ix] / 2)) {//Block on one Reservoir Boundary Ix == 0
				Fraction_well = 1;
				aj[2] = gridDim[Ix] + 2 * X_posi;
				rij[2][1] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[2][7] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[2][2] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (Y_posi*Y_posi));
				rij[2][6] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[2][3] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[2][8] = 1;
				rij[2][4] = 1;
				rij[2][5] = 1;

				ri_1[8] = 1;
				ri_1[4] = 1;
				ri_1[5] = 1;

				TransWell[8] = 0;
				TransWell[4] = 0;
				TransWell[5] = 0;

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j])*pow((rij[2][j] / aj[2]), TransWell[j]);

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);

			}
			if (X_posi == (-gridDim[Ix] / 2)) {//Well on Reservoir Boundary Ix == 0
				Fraction_well = 0.5;
				ri_1[8] = 1;
				ri_1[4] = 1;
				ri_1[5] = 1;

				TransWell[8] = 0;
				TransWell[4] = 0;
				TransWell[5] = 0;

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j]);
				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;
				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);
			}
		}

		if (Ix == 0 && Iy == Ny - 1) {//Corner Block
			if (X_posi > (-gridDim[Ix] / 2) && X_posi <= (gridDim[Ix] / 2) && Y_posi > (-gridDim[Nx + Iy] / 2) && Y_posi <= (gridDim[Nx + Iy] / 2)) {
				Fraction_well = 1;
				aj[2] = gridDim[Nx + Iy] + 2 * Y_posi;
				aj[4] = gridDim[Ix] + 2 * X_posi;
				aj[3] = sqrt((gridDim[Nx + Iy] + 2 * Y_posi)* (gridDim[Nx + Iy] + 2 * Y_posi) + (gridDim[Ix] + 2 * X_posi)*(gridDim[Ix] + 2 * X_posi));

				rij[2][2] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[2][7] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));
				rij[2][3] = sqrt((X_posi*X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));

				rij[3][2] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[3][3] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));
				rij[3][7] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));

				rij[4][2] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (Y_posi * Y_posi));
				rij[4][3] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[4][7] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));


				for (j = 1; j < 9; j++)
					if (j != 2 && j != 3 && j != 7) {
						TransWell[j] = 0;
						ri_1[j] = 1;
						rij[2][j] = 1;
						rij[3][j] = 1;
						rij[4][j] = 1;
					}

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				iZarb = 1;
				for (int t = 1; t < 5; t++) {
					jZarb = 1;
					for (j = 2; j < 5; j++)
						jZarb *= pow((rij[j][t] / aj[j]), TransWell[t]);
					iZarb *= pow(ri_1[t], TransWell[t])*jZarb;
				}

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);



			}
			if (X_posi == (-gridDim[Ix] / 2) || Y_posi == (-gridDim[Nx + Iy] / 2)) {
				if (X_posi == (-gridDim[Ix] / 2) && Y_posi == (-gridDim[Nx + Iy] / 2)) Fraction_well = 0.25;
				else Fraction_well = 0.5;

				for (j = 1; j < 9; j++)
					if (j != 2 && j != 3 && j != 7) {
						TransWell[j] = 0;
						ri_1[j] = 1;
					}

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j]);

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);


			}
		}

		if (Ix == 0 && Iy == 0) {//Corner Block
			if (X_posi > (-gridDim[Ix] / 2) && X_posi <= (gridDim[Ix] / 2) && Y_posi >= (-gridDim[Nx + Iy] / 2) && Y_posi < (gridDim[Nx + Iy] / 2)) {
				Fraction_well = 1;
				aj[2] = gridDim[Nx + Iy] - 2 * Y_posi;
				aj[4] = gridDim[Ix] + 2 * X_posi;
				aj[3] = sqrt((gridDim[Nx + Iy] - 2 * Y_posi)* (gridDim[Nx + Iy] - 2 * Y_posi) + (gridDim[Ix] + 2 * X_posi)*(gridDim[Ix] + 2 * X_posi));

				rij[2][2] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[2][6] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));
				rij[2][1] = sqrt((X_posi*X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));

				rij[3][2] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[3][1] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));
				rij[3][6] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));

				rij[4][2] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (Y_posi * Y_posi));
				rij[4][1] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[4][6] = sqrt((2 * gridDim[Ix] + X_posi)*(2 * gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));


				for (j = 1; j < 9; j++)
					if (j != 1 && j != 2 && j != 6) {
						TransWell[j] = 0;
						ri_1[j] = 1;
						rij[2][j] = 1;
						rij[3][j] = 1;
						rij[4][j] = 1;
					}

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				iZarb = 1;
				for (int t = 1; t < 5; t++) {
					jZarb = 1;
					for (j = 2; j < 5; j++)
						jZarb *= pow((rij[j][t] / aj[j]), TransWell[t]);
					iZarb *= pow(ri_1[t], TransWell[t])*jZarb;
				}

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);


			}
			if (X_posi == (-gridDim[Ix] / 2) || Y_posi == (gridDim[Nx + Iy] / 2)) {
				if (X_posi == (-gridDim[Ix] / 2) && Y_posi == (gridDim[Nx + Iy] / 2)) Fraction_well = 0.25;
				else Fraction_well = 0.5;

				for (j = 1; j < 9; j++)
					if (j != 2 && j != 1 && j != 6) {
						TransWell[j] = 0;
						ri_1[j] = 1;
					}

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j]);

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);

			}
		}

		if (Ix == Nx - 1 && Iy == Ny - 1) {//Corner Block
			if (X_posi >= (-gridDim[Ix] / 2) && X_posi < (gridDim[Ix] / 2) && Y_posi >(-gridDim[Nx + Iy] / 2) && Y_posi <= (gridDim[Nx + Iy] / 2)) {
				Fraction_well = 1;
				aj[2] = gridDim[Ix] - 2 * X_posi;
				aj[4] = gridDim[Nx + Iy] + 2 * Y_posi;
				aj[3] = sqrt((gridDim[Ix] - 2 * X_posi)* (gridDim[Ix] - 2 * X_posi) + (gridDim[Nx + Iy] + 2 * Y_posi)*(gridDim[Nx + Iy] + 2 * Y_posi));

				rij[2][4] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (Y_posi*Y_posi));
				rij[2][8] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[2][3] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));

				rij[3][4] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[3][3] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));
				rij[3][8] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));

				rij[4][3] = sqrt((2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi) + (X_posi * X_posi));
				rij[4][4] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[4][8] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (2 * gridDim[Nx + Iy] + Y_posi)*(2 * gridDim[Nx + Iy] + Y_posi));


				for (j = 1; j < 9; j++)
					if (j != 3 && j != 4 && j != 8) {
						TransWell[j] = 0;
						ri_1[j] = 1;
						rij[2][j] = 1;
						rij[3][j] = 1;
						rij[4][j] = 1;
					}

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				iZarb = 1;
				for (int t = 1; t < 5; t++) {
					jZarb = 1;
					for (j = 2; j < 5; j++)
						jZarb *= pow((rij[j][t] / aj[j]), TransWell[t]);
					iZarb *= pow(ri_1[t], TransWell[t])*jZarb;
				}

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);

			}
			if (X_posi == (gridDim[Ix] / 2) || Y_posi == (-gridDim[Nx + Iy] / 2)) {
				if (X_posi == (gridDim[Ix] / 2) && Y_posi == (-gridDim[Nx + Iy] / 2)) Fraction_well = 0.25;
				else Fraction_well = 0.5;

				for (j = 1; j < 9; j++)
					if (j != 4 && j != 8 && j != 3) {
						TransWell[j] = 0;
						ri_1[j] = 1;
					}

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j]);

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);

			}
		}

		if (Ix == Nx - 1 && Iy == 0) {//Corner Block
			if (X_posi >= (-gridDim[Ix] / 2) && X_posi < (gridDim[Ix] / 2) && Y_posi >= (-gridDim[Nx + Iy] / 2) && Y_posi < (gridDim[Nx + Iy] / 2)) {
				Fraction_well = 1;
				aj[2] = gridDim[Ix] - 2 * X_posi;
				aj[4] = gridDim[Nx + Iy] - 2 * Y_posi;
				aj[3] = sqrt((gridDim[Ix] - 2 * X_posi)* (gridDim[Ix] - 2 * X_posi) + (gridDim[Nx + Iy] - 2 * Y_posi)*(gridDim[Nx + Iy] - 2 * Y_posi));

				rij[2][1] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[2][5] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] + Y_posi)*(gridDim[Nx + Iy] + Y_posi));
				rij[2][4] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (Y_posi*Y_posi));

				rij[3][4] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[3][1] = sqrt((gridDim[Ix] - X_posi)*(gridDim[Ix] - X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));
				rij[3][5] = sqrt((2 * gridDim[Ix] - X_posi)*(2 * gridDim[Ix] - X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));

				rij[4][1] = sqrt((2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi) + (X_posi * X_posi));
				rij[4][4] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (gridDim[Nx + Iy] - Y_posi)*(gridDim[Nx + Iy] - Y_posi));
				rij[4][5] = sqrt((gridDim[Ix] + X_posi)*(gridDim[Ix] + X_posi) + (2 * gridDim[Nx + Iy] - Y_posi)*(2 * gridDim[Nx + Iy] - Y_posi));


				for (j = 1; j < 9; j++)
					if (j != 1 && j != 4 && j != 5) {
						TransWell[j] = 0;
						ri_1[j] = 1;
						rij[2][j] = 1;
						rij[3][j] = 1;
						rij[4][j] = 1;
					}

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				iZarb = 1;
				for (int t = 1; t < 5; t++) {
					jZarb = 1;
					for (j = 2; j < 5; j++)
						jZarb *= pow((rij[j][t] / aj[j]), TransWell[t]);
					iZarb *= pow(ri_1[t], TransWell[t])*jZarb;
				}

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);

			}
			if (X_posi == (gridDim[Ix] / 2) || Y_posi == (gridDim[Nx + Iy] / 2)) {
				if (X_posi == (gridDim[Ix] / 2) && Y_posi == (gridDim[Nx + Iy] / 2)) Fraction_well = 0.25;
				else Fraction_well = 0.5;

				for (j = 1; j < 9; j++)
					if (j != 4 && j != 1 && j != 5) {
						TransWell[j] = 0;
						ri_1[j] = 1;
					}

				for (j = 1; j < 5; j++)
					iZarb *= pow(ri_1[j], TransWell[j]);

				for (j = 1; j < 9; j++)
					bSum += TransWell[j];
				bSum = 1 / bSum;

				Re = pow(exp(-2 * PI*Fraction_well)*iZarb, bSum);
			}
		}
	}



	for (i = 0; i < 5; i++) {
		free(rij[i]);
	}
	free(rij);
	free(aj);
	free(TransWell);
	free(ri_1);

	return Re;
}

#endif