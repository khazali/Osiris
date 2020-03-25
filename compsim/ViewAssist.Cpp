#ifndef VIEWASSIST_H
#define VIEWASSIST_H

//#include <stdio.h>
#include <iostream>
#include <fstream>
#include "MIfstream.h"
#include "Functions.h"
#include <iomanip>


extern FType cumQo, cumQg, cumQw;


//void PrintAll(std::ofstream &fp, FType);

void Print2D(FType **mat2D, FType *mat1D, int matSize, std::ofstream &fp) {
	register int i, j;

	for (i = 0; i<matSize; i++) {
		for (j = 0; j<matSize; j++) {
			//fprintf(fp, "%e ", mat2D[i][j]);
			fp << mat2D[i][j];
		}
		//fprintf(fp, "; ");
		fp << ";";
	}

	//fprintf(fp, "\n\n\n");
	fp << "\n\n\n";

	for (i = 0; i<matSize; i++) {
		//fprintf(fp, "%e; ", mat1D[i]);
		fp << mat1D[i];
	}

}

void PrintAll(std::ofstream &fp, FType t) {
	register int Ix, Iy, Iz, i;
	static int CHours = 0;
	//FType ssum;

	//if (t>=(3600.0*CHours)) {
	for (i = 0; i<80; i++) fp << "-";
	//fprintf(fp, "%f\n", t/3600);
	//fprintf(fp, "\nX\tY\tZ\tP\t\t\tSw\t\t\tSo\t\t\tSg\t\t\t");
	fp << "\n" << t / 3600;
	fp << "\nX\tY\tZ\tP\t\t\tSw\t\t\tSo\t\t\tSg\t\t\t";
	//for (i=0; i<Nc; i++) fprintf(fp, "Xo[%d]\t\t\t", i);
	for (i = 0; i < Nc; i++) fp << "Xo[" << i << "]\t\t\t";
	//for (i=0; i<Nc; i++) fprintf(fp, "Yg[%d]\t\t\t", i);
	for (i = 0; i<Nc; i++) fp << "Yg[" << i << "]\t\t\t";

	fp << "\nM\n";
	for (Iz = 0; Iz<Nz; Iz++)
		for (Iy = 0; Iy<Ny; Iy++)
			for (Ix = 0; Ix<Nx; Ix++) {
				//fprintf(fp, "\n%d\t%d\t%d\t%f\t%f\t%f\t%f\t", Ix, Iy, Iz, P[Ix+1][Iy+1][Iz+1][1], sat[Ix][Iy][Iz][0], sat[Ix][Iy][Iz][1], sat[Ix][Iy][Iz][2]);
				fp << std::setprecision(20) << "\n" << Ix << "\t" << Iy << "\t" << Iz << "\t" << P[Ix + 1][Iy + 1][Iz + 1][1] << "\t" << sat[Ix][Iy][Iz][0] << "\t" << sat[Ix][Iy][Iz][1] << "\t" << sat[Ix][Iy][Iz][2];
				//for (i=0; i<Nc; i++) fprintf(fp, "%f\t", comp[Ix][Iy][Iz][i][0]);
				//for (i=0; i<Nc; i++) fprintf(fp, "%f\t", comp[Ix][Iy][Iz][i][1]);
				for (i = 0; i < Nc; i++) fp << std::setprecision(20) << "\t" << comp[Ix][Iy][Iz][i][0] << "\t";
				for (i = 0; i<Nc; i++) fp << std::setprecision(20) << comp[Ix][Iy][Iz][i][1] << "\t";
				//ssum = 0;
				//for (i = 0; i < Nc; i++) ssum+= (sat[Ix][Iy][Iz][1] * comp[Ix][Iy][Iz][i][0] * blockFProps[Ix][Iy][Iz][RO][0] + sat[Ix][Iy][Iz][2] * comp[Ix][Iy][Iz][i][1] * blockFProps[Ix][Iy][Iz][RO][1]) / (1 - sat[Ix][Iy][Iz][0]);
				//for (i = 0; i < Nc; i++) fp << (sat[Ix][Iy][Iz][1] * comp[Ix][Iy][Iz][i][0] * blockFProps[Ix][Iy][Iz][RO][0] + sat[Ix][Iy][Iz][2] * comp[Ix][Iy][Iz][i][1] * blockFProps[Ix][Iy][Iz][RO][1]) / ((1 - sat[Ix][Iy][Iz][0])*ssum) << "\t";
			}
	fp << "\nF\n";
	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				//fprintf(fp, "\n%d\t%d\t%d\t%f\t%f\t%f\t%f\t", Ix, Iy, Iz, P[Ix+1][Iy+1][Iz+1][1], sat[Ix][Iy][Iz][0], sat[Ix][Iy][Iz][1], sat[Ix][Iy][Iz][2]);
				fp << std::setprecision(20) << "\n" << Ix << "\t" << Iy << "\t" << Iz << "\t" << fP[Ix + 1][Iy + 1][Iz + 1][1] << "\t" << fsat[Ix][Iy][Iz][0] << "\t" << fsat[Ix][Iy][Iz][1] << "\t" << fsat[Ix][Iy][Iz][2];
				//for (i=0; i<Nc; i++) fprintf(fp, "%f\t", comp[Ix][Iy][Iz][i][0]);
				//for (i=0; i<Nc; i++) fprintf(fp, "%f\t", comp[Ix][Iy][Iz][i][1]);
				for (i = 0; i < Nc; i++) fp << std::setprecision(20) << "\t" << fcomp[Ix][Iy][Iz][i][0] << "\t";
				for (i = 0; i < Nc; i++) fp << std::setprecision(20) << fcomp[Ix][Iy][Iz][i][1] << "\t";
				//ssum = 0;
				//for (i = 0; i < Nc; i++) ssum+= (fsat[Ix][Iy][Iz][1] * fcomp[Ix][Iy][Iz][i][0] * fblockFProps[Ix][Iy][Iz][RO][0] + fsat[Ix][Iy][Iz][2] * fcomp[Ix][Iy][Iz][i][1] * fblockFProps[Ix][Iy][Iz][RO][1]) / (1 - fsat[Ix][Iy][Iz][0]);
				//for (i = 0; i < Nc; i++) fp << (fsat[Ix][Iy][Iz][1] * fcomp[Ix][Iy][Iz][i][0] * fblockFProps[Ix][Iy][Iz][RO][0] + fsat[Ix][Iy][Iz][2] * fcomp[Ix][Iy][Iz][i][1] * fblockFProps[Ix][Iy][Iz][RO][1]) / ((1 - fsat[Ix][Iy][Iz][0])*ssum) << "\t";

			}
	fp << "\n\n\n";
	fp.flush();
	CHours++;
	//}
}

void PrintAll2(int Ix, int Iy, int Iz, std::ofstream &fp, FType t) {
	register int i;
	FType ssum;
	static int CHours = 0;

	//if (t >= (3600.0*CHours)) {

		fp <<t/3600<<"\t"<< P[Ix + 1][Iy + 1][Iz + 1][1] << "\t" << sat[Ix][Iy][Iz][0] << "\t" << sat[Ix][Iy][Iz][1] << "\t" << sat[Ix][Iy][Iz][2];
		for (i = 0; i < Nc; i++) fp << "\t" << comp[Ix][Iy][Iz][i][0] << "\t";
		for (i = 0; i < Nc; i++) fp << comp[Ix][Iy][Iz][i][1] << "\t";
		ssum = 0;
		for (i = 0; i < Nc; i++) ssum += (sat[Ix][Iy][Iz][1] * comp[Ix][Iy][Iz][i][0] * blockFProps[Ix][Iy][Iz][RO][0] + sat[Ix][Iy][Iz][2] * comp[Ix][Iy][Iz][i][1] * blockFProps[Ix][Iy][Iz][RO][1]) / (1 - sat[Ix][Iy][Iz][0]);
		for (i = 0; i < Nc; i++) fp << (sat[Ix][Iy][Iz][1] * comp[Ix][Iy][Iz][i][0] * blockFProps[Ix][Iy][Iz][RO][0] + sat[Ix][Iy][Iz][2] * comp[Ix][Iy][Iz][i][1] * blockFProps[Ix][Iy][Iz][RO][1]) / ((1 - sat[Ix][Iy][Iz][0])*ssum) << "\t";




		fp << "\t" << fP[Ix + 1][Iy + 1][Iz + 1][1] << "\t" << fsat[Ix][Iy][Iz][0] << "\t" << fsat[Ix][Iy][Iz][1] << "\t" << fsat[Ix][Iy][Iz][2];
		for (i = 0; i < Nc; i++) fp << "\t" << fcomp[Ix][Iy][Iz][i][0] << "\t";
		for (i = 0; i < Nc; i++) fp << fcomp[Ix][Iy][Iz][i][1] << "\t";
		ssum = 0;
		for (i = 0; i < Nc; i++) ssum += (fsat[Ix][Iy][Iz][1] * fcomp[Ix][Iy][Iz][i][0] * fblockFProps[Ix][Iy][Iz][RO][0] + fsat[Ix][Iy][Iz][2] * fcomp[Ix][Iy][Iz][i][1] * fblockFProps[Ix][Iy][Iz][RO][1]) / (1 - fsat[Ix][Iy][Iz][0]);
		for (i = 0; i < Nc; i++) fp << (fsat[Ix][Iy][Iz][1] * fcomp[Ix][Iy][Iz][i][0] * fblockFProps[Ix][Iy][Iz][RO][0] + fsat[Ix][Iy][Iz][2] * fcomp[Ix][Iy][Iz][i][1] * fblockFProps[Ix][Iy][Iz][RO][1]) / ((1 - fsat[Ix][Iy][Iz][0])*ssum) << "\t";

		fp << std::endl;
		fp.flush();
		CHours++;
	//}
	
}


#endif