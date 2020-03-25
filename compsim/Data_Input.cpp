//#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include "MIfstream.h"
#include <fstream>
#include<string>
#include "Functions.h"



void Data_Input(MIfstream &fp) {
	char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
	register int i, j, k, n;
	FType tempL;
	int eclint;

	//Dimension
	if (!fp.FileSearch("GRID")) TerM("No GRID keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect GRID keyword format in the input file!");
	Nx = atoi(str);
	if (!fp.ReadWord(str)) TerM("Incorrect GRID keyword format in the input file!");
	Ny = atoi(str);
	if (!fp.ReadWord(str)) TerM("Incorrect GRID keyword format in the input file!");
	Nz = atoi(str);

	//Number of Components
	if (!fp.FileSearch("NC")) TerM("No NC keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect NC keyword format in the input file!");
	PNc = atoi(str);
	if (!fp.ReadWord(str)) TerM("Incorrect NC keyword format in the input file!");
	UNc = atoi(str);
	Nc = PNc + UNc;


	//Saturation Tables
	if (!fp.FileSearch("SWT")) TerM("No SWT keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect SWT keyword format in the input file!");
	Nswt = atoi(str);

	if (!fp.FileSearch("SGT")) TerM("No SGT keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect SGT keyword format in the input file!");
	Nsgt = atoi(str);

	//Wells
	if (!fp.FileSearch("WELLS")) TerM("No WELLS keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect WELLS keyword format in the input file!");
	wellNO = atoi(str);

	if (fp.FileSearch("DBS")) {
		if (!fp.ReadWord(str)) TerM("Incorrect DBS keyword format in the input file!");
		NOOptionsDB = atoi(str);
	}
	else {
		NOOptionsDB = 0;
	}


	Allocation();


	//Reservoir Temperature
	if (!fp.FileSearch("RESTEMP")) TerM("No RESTEMP keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect RESTEMP keyword format in the input file!");
	resTemp = atof(str) + 273.15;


	//Water Properties
	if (!fp.FileSearch("WATERPROPS")) TerM("No WATERPROPS keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect WATERPROPS keyword format in the input file!");
	watRo = atof(str);
	if (!fp.ReadWord(str)) TerM("Incorrect WATERPROPS keyword format in the input file!");
	watMu = atof(str);


	//////////block sizes (finite difference)//////////////
	if (!fp.FileSearch("DI")) TerM("No DI keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect DI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i = 0; i<Nx; i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect DI keyword format in the input file!");
		gridDim[i] = atof(str1);
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect DI keyword format in the input file!");
		tempL = atof(str1);
		for (i = 0; i<Nx; i++) gridDim[i] = tempL;
	}
	else {
		TerM("Incorrect DI keyword format in the input file!");
	}

	if (!fp.FileSearch("DJ")) TerM("No DJ keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect DJ keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i = Nx; i<(Nx + Ny); i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect DJ keyword format in the input file!");
		gridDim[i] = atof(str1);
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect DJ keyword format in the input file!");
		tempL = atof(str1);
		for (i = Nx; i<(Nx + Ny); i++) gridDim[i] = tempL;
	}
	else {
		TerM("Incorrect DJ keyword format in the input file!");
	}

	if (!fp.FileSearch("DK")) TerM("No DK keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect DK keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i = (Nx + Ny); i<(Nx + Ny + Nz); i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect DK keyword format in the input file!");
		gridDim[i] = atof(str1);
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect DK keyword format in the input file!");
		tempL = atof(str1);
		for (i = (Nx + Ny); i<(Nx + Ny + Nz); i++) gridDim[i] = tempL;
	}
	else {
		TerM("Incorrect DK keyword format in the input file!");
	}

	CalcBlockHeight();


	//Tortuosity
	if (!fp.FileSearch("TOR")) {
		if (!MPIRank) std::cout << "All tortuosity values reset to unity!\n";
		for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) tor[i][j][k] = 1;
	}
	else {
		if (!fp.ReadWord(str)) TerM("Incorrect TOR keyword format in the input file!");
		if (!strcmp(str, "VAR")) for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect TOR keyword format in the input file!");
			tor[i][j][k] = atof(str1);
		}
		else if (!strcmp(str, "CON")) {
			if (!fp.ReadWord(str1)) TerM("Incorrect TOR keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) tor[i][j][k] = tempL;
		}
		else if (!strcmp(str, "IVAR")) {
			for (i = 0; i<Nx; i++) {
				if (!fp.ReadWord(str1)) TerM("Incorrect TOR keyword format in the input file!");
				tempL = atof(str1);
				for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) tor[i][j][k] = tempL;
			}
		}
		else if (!strcmp(str, "JVAR")) {
			for (j = 0; i<Ny; j++) {
				if (!fp.ReadWord(str1)) TerM("Incorrect TOR keyword format in the input file!");
				tempL = atof(str1);
				for (k = 0; k<Nz; k++) for (i = 0; i<Nx; i++) tor[i][j][k] = tempL;
			}
		}
		else if (!strcmp(str, "KVAR")) {
			for (k = 0; k<Nz; k++) {
				if (!fp.ReadWord(str1)) TerM("Incorrect TOR keyword format in the input file!");
				tempL = atof(str1);
				for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) tor[i][j][k] = tempL;
			}
		}
		else {
			TerM("Incorrect TOR keyword format in the input file!");
		}
	}


	//POROSITY
	if (!fp.FileSearch("POR")) TerM("No POR keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect POR keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect POR keyword format in the input file!");
		porosity[i][j][k] = atof(str1);
		if (!porosity[i][j][k]) porosity[i][j][k] = 1e-11;
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect POR keyword format in the input file!");
		tempL = atof(str1);
		for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) porosity[i][j][k] = tempL;
	}
	else if (!strcmp(str, "IVAR")) {
		for (i = 0; i<Nx; i++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect POR keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) porosity[i][j][k] = tempL;
		}
	}
	else if (!strcmp(str, "JVAR")) {
		for (j = 0; i<Ny; j++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect POR keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k<Nz; k++) for (i = 0; i<Nx; i++) porosity[i][j][k] = tempL;
		}
	}
	else if (!strcmp(str, "KVAR")) {
		for (k = 0; k<Nz; k++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect POR keyword format in the input file!");
			tempL = atof(str1);
			for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) porosity[i][j][k] = tempL;
		}
	}
	else {
		TerM("Incorrect POR keyword format in the input file!");
	}

	if (fp.FileSearch("CDIF")) {
		IsConstDif = true;
		if (!fp.ReadWord(str)) TerM("Incorrect fracture CDIF keyword format in the input file!");
		CDiffCoef = atof(str);
	}
	else {
		IsConstDif = false;
	}

	if (fp.FileSearch("SIGMA")) {		
		if (!fp.ReadWord(str)) TerM("Incorrect fracture SIGMA keyword format in the input file!");
		SigmaMF = atof(str);
	}
	else {
		SigmaMF = 1;
	}
	
	if (fp.FileSearch("DTMAX")) {		
		if (!fp.ReadWord(str)) TerM("Incorrect fracture DTMAX keyword format in the input file!");
		Max_Dt = atof(str);
	}
	else {
		Max_Dt = 864000;
	}

	if (fp.FileSearch("ABSTOL")) {
		if (!fp.ReadWord(str)) TerM("Incorrect fracture ABSTOL keyword format in the input file!");
		KSPAbsTol = atof(str);
	}
	else {
		KSPAbsTol = 1e-50;
	}
	if (fp.FileSearch("RELTOL")) {
		if (!fp.ReadWord(str)) TerM("Incorrect fracture RelTol keyword format in the input file!");
		KSPRelTol = atof(str);
	}
	else {
		KSPRelTol = 1e-18;
	}
	if (fp.FileSearch("DIVTOL")) {
		if (!fp.ReadWord(str)) TerM("Incorrect fracture DIVTOL keyword format in the input file!");
		KSPDivTol = atof(str);
	}
	else {
		KSPDivTol = 1e3;
	}
	if (fp.FileSearch("MAXITERS")) {
		if (!fp.ReadWord(str)) TerM("Incorrect fracture MAXITERS keyword format in the input file!");
		MaxIters = atoi(str);
	}
	else {
		MaxIters = 5000;
	}
	if (fp.FileSearch("SOLVERTYPE")) {
		if (!fp.ReadWord(str)) TerM("Incorrect fracture SOLVERTYPE keyword format in the input file!");
		strcpy(SolverType, str);
	}
	else {
		strcpy(SolverType, "bcgs");
	}

	if (fp.FileSearch("DBS")) {
		fp.ReadWord(str);
		for (i = 0; i < NOOptionsDB; i++) {
			if (!fp.ReadWord(str)) TerM("Incorrect DBS options keyword format in the input file!");
			strcpy(&(OptionsDBFirst[i][0]), str);
			if (!fp.ReadWord(str)) TerM("Incorrect DBS options keyword format in the input file!");
			strcpy(&(OptionsDBSecond[i][0]), str);
		}
	}

	//if (!fp.ReadWord(str)) TerM("Incorrect fracture POR keyword format in the input file!");
	//Fracture POROSITY
	if (!fp.FileSearch("FPOR")) TerM("No fracture POR keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect fracture POR keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect fracture POR keyword format in the input file!");
		fporosity[i][j][k] = atof(str1);
		if (!fporosity[i][j][k]) fporosity[i][j][k] = 1e-11;
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect fracture POR keyword format in the input file!");
		tempL = atof(str1);
		for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fporosity[i][j][k] = tempL;
	}
	else if (!strcmp(str, "IVAR")) {
		for (i = 0; i < Nx; i++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect fracture POR keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) fporosity[i][j][k] = tempL;
		}
	}
	else if (!strcmp(str, "JVAR")) {
		for (j = 0; i < Ny; j++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect fracture POR keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k < Nz; k++) for (i = 0; i < Nx; i++) fporosity[i][j][k] = tempL;
		}
	}
	else if (!strcmp(str, "KVAR")) {
		for (k = 0; k < Nz; k++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect fracture POR keyword format in the input file!");
			tempL = atof(str1);
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fporosity[i][j][k] = tempL;
		}
	}
	else {
		TerM("Incorrect fracture POR keyword format in the input file!");
	}

	//Half of Fracture opening
	if (!fp.FileSearch("Lf")) TerM("No Lf(Half of Fracture opening) keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect Lf(Half of Fracture opening) keyword format in the input file!");
	Lf = atof(str);

	// Half of Fracture length
		if (!fp.FileSearch("Re_f")) TerM("No Re_f (Half  Fracture length) keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect Re_f (Half  Fracture length) keyword format in the input file!");
	Re_f = atof(str);

	//PERMEABILITY
	if (!fp.FileSearch("PERMI")) TerM("No PERMI keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect PERMI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect PERMI keyword format in the input file!");
		perm[i][j][k][0] = atof(str1);
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect PERMI keyword format in the input file!");
		tempL = atof(str1);
		for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) perm[i][j][k][0] = tempL;
	}
	else if (!strcmp(str, "IVAR")) {
		for (i = 0; i<Nx; i++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect PERMI keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) perm[i][j][k][0] = tempL;
		}
	}
	else if (!strcmp(str, "JVAR")) {
		for (j = 0; i<Ny; j++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect PERMI keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k<Nz; k++) for (i = 0; i<Nx; i++) perm[i][j][k][0] = tempL;
		}
	}
	else if (!strcmp(str, "KVAR")) {
		for (k = 0; k<Nz; k++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect PERMI keyword format in the input file!");
			tempL = atof(str1);
			for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) perm[i][j][k][0] = tempL;
		}
	}
	else {
		TerM("Incorrect PERMI keyword format in the input file!");
	}

	if (!fp.FileSearch("PERMJ")) TerM("No PERMJ keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect PERMJ keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect PERMI keyword format in the input file!");
		perm[i][j][k][1] = atof(str1);
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect PERJ keyword format in the input file!");
		tempL = atof(str1);
		for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) perm[i][j][k][1] = tempL;
	}
	else if (!strcmp(str, "IVAR")) {
		for (i = 0; i<Nx; i++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect PERMJ keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) perm[i][j][k][1] = tempL;
		}
	}
	else if (!strcmp(str, "JVAR")) {
		for (j = 0; i<Ny; j++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect PERMJ keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k<Nz; k++) for (i = 0; i<Nx; i++) perm[i][j][k][1] = tempL;
		}
	}
	else if (!strcmp(str, "KVAR")) {
		for (k = 0; k<Nz; k++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect PERMJ keyword format in the input file!");
			tempL = atof(str1);
			for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) perm[i][j][k][1] = tempL;
		}
	}
	else {
		TerM("Incorrect PERMJ keyword format in the input file!");
	}

	if (!fp.FileSearch("PERMK")) TerM("No PERMK keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect PERMK keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect PERMK keyword format in the input file!");
		perm[i][j][k][2] = atof(str1);
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect PERMK keyword format in the input file!");
		tempL = atof(str1);
		for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) perm[i][j][k][2] = tempL;
	}
	else if (!strcmp(str, "IVAR")) {
		for (i = 0; i<Nx; i++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect PERMK keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k<Nz; k++) for (j = 0; j<Ny; j++) perm[i][j][k][2] = tempL;
		}
	}
	else if (!strcmp(str, "JVAR")) {
		for (j = 0; i<Ny; j++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect PERMK keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k<Nz; k++) for (i = 0; i<Nx; i++) perm[i][j][k][2] = tempL;
		}
	}
	else if (!strcmp(str, "KVAR")) {
		for (k = 0; k<Nz; k++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect PERMK keyword format in the input file!");
			tempL = atof(str1);
			for (j = 0; j<Ny; j++) for (i = 0; i<Nx; i++) perm[i][j][k][2] = tempL;
		}
	}
	else {
		TerM("Incorrect PERMK keyword format in the input file!");
	}

	//Fracture PERMEABILITY
	if (!fp.FileSearch("FPERMI")) TerM("No Fracture PERMI keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect Fracture PERMI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMI keyword format in the input file!");
		fperm[i][j][k][0] = atof(str1);
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMI keyword format in the input file!");
		tempL = atof(str1);
		for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fperm[i][j][k][0] = tempL;
	}
	else if (!strcmp(str, "IVAR")) {
		for (i = 0; i < Nx; i++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMI keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) fperm[i][j][k][0] = tempL;
		}
	}
	else if (!strcmp(str, "JVAR")) {
		for (j = 0; i < Ny; j++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMI keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k < Nz; k++) for (i = 0; i < Nx; i++) fperm[i][j][k][0] = tempL;
		}
	}
	else if (!strcmp(str, "KVAR")) {
		for (k = 0; k < Nz; k++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMI keyword format in the input file!");
			tempL = atof(str1);
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fperm[i][j][k][0] = tempL;
		}
	}
	else {
		TerM("Incorrect Fracture PERMI keyword format in the input file!");
	}

	if (!fp.FileSearch("FPERMJ")) TerM("No Fracture PERMJ keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect Fracture PERMJ keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMI keyword format in the input file!");
		fperm[i][j][k][1] = atof(str1);
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERJ keyword format in the input file!");
		tempL = atof(str1);
		for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fperm[i][j][k][1] = tempL;
	}
	else if (!strcmp(str, "IVAR")) {
		for (i = 0; i < Nx; i++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMJ keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) fperm[i][j][k][1] = tempL;
		}
	}
	else if (!strcmp(str, "JVAR")) {
		for (j = 0; i < Ny; j++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMJ keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k < Nz; k++) for (i = 0; i < Nx; i++) fperm[i][j][k][1] = tempL;
		}
	}
	else if (!strcmp(str, "KVAR")) {
		for (k = 0; k < Nz; k++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMJ keyword format in the input file!");
			tempL = atof(str1);
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fperm[i][j][k][1] = tempL;
		}
	}
	else {
		TerM("Incorrect Fracture PERMJ keyword format in the input file!");
	}

	if (!fp.FileSearch("FPERMK")) TerM("No Fracture PERMK keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect Fracture PERMK keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
		if (!fp.ReadWord(str1)) TerM("Incorrect FracturePERMK keyword format in the input file!");
		fperm[i][j][k][2] = atof(str1);
	}
	else if (!strcmp(str, "CON")) {
		if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMK keyword format in the input file!");
		tempL = atof(str1);
		for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fperm[i][j][k][2] = tempL;
	}
	else if (!strcmp(str, "IVAR")) {
		for (i = 0; i < Nx; i++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMK keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) fperm[i][j][k][2] = tempL;
		}
	}
	else if (!strcmp(str, "JVAR")) {
		for (j = 0; i < Ny; j++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMK keyword format in the input file!");
			tempL = atof(str1);
			for (k = 0; k < Nz; k++) for (i = 0; i < Nx; i++) fperm[i][j][k][2] = tempL;
		}
	}
	else if (!strcmp(str, "KVAR")) {
		for (k = 0; k < Nz; k++) {
			if (!fp.ReadWord(str1)) TerM("Incorrect Fracture PERMK keyword format in the input file!");
			tempL = atof(str1);
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fperm[i][j][k][2] = tempL;
		}
	}
	else {
		TerM("Incorrect Fracture PERMK keyword format in the input file!");
	}

	//Compressibility
	if (!fp.FileSearch("REFP")) {
		if (!MPIRank) std::cout << "Warning: No REFP keyword in the input file!\n";
		refP = 0;
	}
	else {
		if (!fp.ReadWord(str)) TerM("Incorrect REFP keyword format in the input file!");
		refP = atof(str);
	}
	if (!fp.FileSearch("CPOR")) {
		if (!MPIRank) std::cout << "Warning: No CPOR keyword in the input file!\n";
		cpor = 0;
	}
	else {
		if (!fp.ReadWord(str)) TerM("Incorrect CPOR keyword format in the input file!");
		cpor = atof(str);
	}
	if (!fp.FileSearch("DCPOR")) {
		if (!MPIRank) std::cout << "Warning: No DCPOR keyword in the input file!\n";
		dcpor = 0;
	}
	else {
		if (!fp.ReadWord(str)) TerM("Incorrect DCPOR keyword format in the input file!");
		dcpor = atof(str);
	}

	//EOS Type
	if (!fp.FileSearch("MODEL")) TerM("No EOS MODEL keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect EOS MODEL keyword format in the input file!");
	if (!strcmp(str, "PR")) {
		PR = 1;
		SRK = 0;
	}
	else if (!strcmp(str, "SRK")) {
		PR = 0;
		SRK = 1;
	}
	else {
		TerM("Incorrect EOS keyword format in the input file!");
	}

	if (UNc) {
		if (!fp.FileSearch("SG")) {
			if (!MPIRank) std::cout << "Warning: No SG keyword in the input file!\n";
		}
		else for (i = PNc; i<Nc; i++) {
			if (!fp.ReadWord(str)) TerM("Incorrect SG keyword format in the input file!");
			fluidProp[i][SG] = atof(str);
		}

		if (!fp.FileSearch("TB")) {
			if (!MPIRank) std::cout << "Warning: No TB keyword in the input file!\n";
		}
		else for (i = PNc; i<Nc; i++) {
			if (!fp.ReadWord(str)) TerM("Incorrect TB keyword format in the input file!");
			fluidProp[i][TB] = atof(str);
		}

		if (!fp.FileSearch("MW")) {
			if (!MPIRank) std::cout << "Warning: No MW keyword in the input file!\n";
		}
		else for (i = PNc; i<Nc; i++) {
			if (!fp.ReadWord(str)) TerM("Incorrect MW keyword format in the input file!");
			fluidProp[i][MW] = atof(str);
		}

		if (!fp.FileSearch("AC")) {
			if (!MPIRank) std::cout << "Warning: No AC keyword in the input file!\n";
		}
		else for (i = PNc; i<Nc; i++) {
			if (!fp.ReadWord(str)) TerM("Incorrect AC keyword format in the input file!");
			fluidProp[i][AC] = atof(str);
		}

		if (!fp.FileSearch("PCRIT")) {
			if (!MPIRank) std::cout << "Warning: No PCRIT keyword in the input file!\n";
		}
		else for (i = PNc; i<Nc; i++) {
			if (!fp.ReadWord(str)) TerM("Incorrect PCRIT keyword format in the input file!");
			fluidProp[i][PCRIT] = atof(str);
		}

		if (!fp.FileSearch("TCRIT")) {
			if (!MPIRank) std::cout << "Warning: No TCRIT keyword in the input file!\n";
		}
		else for (i = PNc; i<Nc; i++) {
			if (!fp.ReadWord(str)) TerM("Incorrect TCRIT keyword format in the input file!");
			fluidProp[i][TCRIT] = atof(str);
		}

		if (!fp.FileSearch("VCRIT")) {
			if (!MPIRank) std::cout << "Warning: No VCRIT keyword in the input file!\n";
		}
		else for (i = PNc; i<Nc; i++) {
			if (!fp.ReadWord(str)) TerM("Incorrect VCRIT keyword format in the input file!");
			fluidProp[i][VCRIT] = atof(str);
		}

		if (!fp.FileSearch("ZCRIT")) {
			if (!MPIRank) std::cout << "Warning: No ZCRIT keyword in the input file!\n";
		}
		else for (i = PNc; i < Nc; i++) {
			if (!fp.ReadWord(str)) TerM("Incorrect ZCRIT keyword format in the input file!");
			fluidProp[i][ZCRIT] = atof(str);
		}

	}
	CPlus_Props();

	if (!fp.FileSearch("COMPNAME")) TerM("No COMPNAME keyword in the input file!");
	for (i = 0; i<PNc; i++) {
		if (!fp.ReadWord(str)) TerM("Incorrect COMPNAME keyword format in the input file!");
		if (!strcmp(str, "C1")) {
			fluidProp[i][MW] = 16.043;
			fluidProp[i][TCRIT] = 190.56;
			fluidProp[i][PCRIT] = 4599000;
			fluidProp[i][VCRIT] = 98.6;
			fluidProp[i][AC] = 0.008;	//edited
			fluidProp[i][PARACHOR] = 74.05;
			fluidProp[i][ZCRIT] = 0.2862;
		}
		else if (!strcmp(str, "C2")) {
			fluidProp[i][MW] = 30.070;
			fluidProp[i][TCRIT] = 305.32;
			fluidProp[i][PCRIT] = 4872000;
			fluidProp[i][VCRIT] = 145.5;
			fluidProp[i][AC] = 0.0995;
			fluidProp[i][PARACHOR] = 112.91;
			fluidProp[i][ZCRIT] = 0.2793;
		}
		else if (!strcmp(str, "C3")) {
			fluidProp[i][MW] = 44.096;
			fluidProp[i][TCRIT] = 369.83;
			fluidProp[i][PCRIT] = 4248000;
			fluidProp[i][VCRIT] = 200;
			fluidProp[i][AC] = 0.1523;
			fluidProp[i][PARACHOR] = 154.03;
			fluidProp[i][ZCRIT] = 0.2763;
		}
		else if (!strcmp(str, "iC4")) {
			fluidProp[i][MW] = 58.123;
			fluidProp[i][TCRIT] = 408.14;
			fluidProp[i][PCRIT] = 3648000;
			fluidProp[i][VCRIT] = 262.7;
			fluidProp[i][AC] = 0.1770;
			fluidProp[i][PARACHOR] = 185.32;
			fluidProp[i][ZCRIT] = 0.2824;
		}
		else if (!strcmp(str, "nC4")) {
			fluidProp[i][MW] = 58.123;
			fluidProp[i][TCRIT] = 425.12;
			fluidProp[i][PCRIT] = 3796000;
			fluidProp[i][VCRIT] = 255;
			fluidProp[i][AC] = 0.2002;
			fluidProp[i][PARACHOR] = 193.90;
			fluidProp[i][ZCRIT] = 0.2739;
		}
		else if (!strcmp(str, "iC5")) {
			fluidProp[i][MW] = 72.150;
			fluidProp[i][TCRIT] = 460.43;
			fluidProp[i][PCRIT] = 3381000;
			fluidProp[i][VCRIT] = 305.8;
			fluidProp[i][AC] = 0.2275;
			fluidProp[i][PARACHOR] = 229.37;
			fluidProp[i][ZCRIT] = 0.2701;
		}
		else if (!strcmp(str, "nC5")) {
			fluidProp[i][MW] = 72.150;
			fluidProp[i][TCRIT] = 469.7;
			fluidProp[i][PCRIT] = 3370000;
			fluidProp[i][VCRIT] = 313;
			fluidProp[i][AC] = 0.2515;
			fluidProp[i][PARACHOR] = 236.00;
			fluidProp[i][ZCRIT] = 0.2701;
		}
		else if (!strcmp(str, "nC6")) {
			fluidProp[i][MW] = 86.177;
			fluidProp[i][TCRIT] = 507.6;
			fluidProp[i][PCRIT] = 3025000;
			fluidProp[i][VCRIT] = 371;
			fluidProp[i][AC] = 0.3013;
			fluidProp[i][PARACHOR] = 276.71;
			fluidProp[i][ZCRIT] = 0.2659;
		}
		else if (!strcmp(str, "nC7")) {
			fluidProp[i][MW] = 100.204;
			fluidProp[i][TCRIT] = 540.2;
			fluidProp[i][PCRIT] = 2740000;
			fluidProp[i][VCRIT] = 428;
			fluidProp[i][AC] = 0.3495;
			fluidProp[i][PARACHOR] = 318.44;
			fluidProp[i][ZCRIT] = 0.2611;
		}
		else if (!strcmp(str, "nC8")) {
			fluidProp[i][MW] = 114.231;
			fluidProp[i][TCRIT] = 568.7;
			fluidProp[i][PCRIT] = 2490000;
			fluidProp[i][VCRIT] = 486;
			fluidProp[i][AC] = 0.3996;
			fluidProp[i][PARACHOR] = 359.33;
			fluidProp[i][ZCRIT] = 0.2559;
		}
		else if (!strcmp(str, "nC9")) {
			fluidProp[i][MW] = 128.258;
			fluidProp[i][TCRIT] = 594.6;
			fluidProp[i][PCRIT] = 2290000;
			fluidProp[i][VCRIT] = 544;
			fluidProp[i][AC] = 0.4435;
			fluidProp[i][PARACHOR] = 399.57;
			fluidProp[i][ZCRIT] = 0.2520;
		}
		else if (!strcmp(str, "nC10")) {
			fluidProp[i][MW] = 142.285;
			fluidProp[i][TCRIT] = 617.7;
			fluidProp[i][PCRIT] = 2110000;
			fluidProp[i][VCRIT] = 600;
			fluidProp[i][AC] = 0.4923;
			fluidProp[i][PARACHOR] = 440.69;
			fluidProp[i][ZCRIT] = 0.2465;
		}
		else if (!strcmp(str, "nC11")) {
			fluidProp[i][MW] = 156.312;
			fluidProp[i][TCRIT] = 639;
			fluidProp[i][PCRIT] = 1949000;
			fluidProp[i][VCRIT] = 659;
			fluidProp[i][AC] = 0.5303;
			fluidProp[i][PARACHOR] = 482.00;
			fluidProp[i][ZCRIT] = 0.2419;
		}
		else if (!strcmp(str, "nC12")) {
			fluidProp[i][MW] = 170.338;
			fluidProp[i][TCRIT] = 658;
			fluidProp[i][PCRIT] = 1820000;
			fluidProp[i][VCRIT] = 716;
			fluidProp[i][AC] = 0.5764;
			fluidProp[i][PARACHOR] = 522.26;
			fluidProp[i][ZCRIT] = 0.2382;
		}
		else if (!strcmp(str, "nC13")) {
			fluidProp[i][MW] = 184.365;
			fluidProp[i][TCRIT] = 675;
			fluidProp[i][PCRIT] = 1680000;
			fluidProp[i][VCRIT] = 775;
			fluidProp[i][AC] = 0.6174;
			fluidProp[i][PARACHOR] = 536.77;
			fluidProp[i][ZCRIT] = 0.2320;
		}
		else if (!strcmp(str, "nC14")) {
			fluidProp[i][MW] = 198.392;
			fluidProp[i][TCRIT] = 693;
			fluidProp[i][PCRIT] = 1570000;
			fluidProp[i][VCRIT] = 830;
			fluidProp[i][AC] = 0.6430;
			fluidProp[i][PARACHOR] = 606.05;
			fluidProp[i][ZCRIT] = 0.2262;
		}
		else if (!strcmp(str, "nC15")) {
			fluidProp[i][MW] = 212.419;
			fluidProp[i][TCRIT] = 708;
			fluidProp[i][PCRIT] = 1480000;
			fluidProp[i][VCRIT] = 889;
			fluidProp[i][AC] = 0.6863;
			fluidProp[i][PARACHOR] = 647.43;
			fluidProp[i][ZCRIT] = 0.2235;
		}
		else if (!strcmp(str, "nC16")) {
			fluidProp[i][MW] = 226.446;
			fluidProp[i][TCRIT] = 723;
			fluidProp[i][PCRIT] = 1400000;
			fluidProp[i][VCRIT] = 944;
			fluidProp[i][AC] = 0.7174;
			fluidProp[i][PARACHOR] = 688.50;
			fluidProp[i][ZCRIT] = 0.2199;
		}
		else if (!strcmp(str, "nC17")) {
			fluidProp[i][MW] = 240.473;
			fluidProp[i][TCRIT] = 736;
			fluidProp[i][PCRIT] = 1340000;
			fluidProp[i][VCRIT] = 1000;
			fluidProp[i][AC] = 0.7697;
			fluidProp[i][PARACHOR] = 730.05;
			fluidProp[i][ZCRIT] = 0.2190;
		}
		else if (!strcmp(str, "nC18")) {
			fluidProp[i][MW] = 254.5;
			fluidProp[i][TCRIT] = 747;
			fluidProp[i][PCRIT] = 1270000;
			fluidProp[i][VCRIT] = 1060;
			fluidProp[i][AC] = 0.8114;
			fluidProp[i][PARACHOR] = 771.95;
			fluidProp[i][ZCRIT] = 0.2168;
		}
		else if (!strcmp(str, "nC19")) {
			fluidProp[i][MW] = 268.527;
			fluidProp[i][TCRIT] = 758;
			fluidProp[i][PCRIT] = 1210000;
			fluidProp[i][VCRIT] = 1120;
			fluidProp[i][AC] = 0.8522;
			fluidProp[i][PARACHOR] = 813.85;
			fluidProp[i][ZCRIT] = 0.2150;
		}
		else if (!strcmp(str, "nC20")) {
			fluidProp[i][MW] = 282.553;
			fluidProp[i][TCRIT] = 768;
			fluidProp[i][PCRIT] = 1160000;
			fluidProp[i][VCRIT] = 1170;
			fluidProp[i][AC] = 0.9069;
			fluidProp[i][PARACHOR] = 853.67;
			fluidProp[i][ZCRIT] = 0.2126;
		}
		else if (!strcmp(str, "nC21")) {
			fluidProp[i][MW] = 296.580;
			fluidProp[i][TCRIT] = 781.7;
			fluidProp[i][PCRIT] = 1147000;
			fluidProp[i][VCRIT] = 1198;
			fluidProp[i][AC] = 0.9220;
			fluidProp[i][PARACHOR] = 897.64;
			fluidProp[i][ZCRIT] = 0.2114;
		}
		else if (!strcmp(str, "nC22")) {
			fluidProp[i][MW] = 310.610;
			fluidProp[i][TCRIT] = 791.8;
			fluidProp[i][PCRIT] = 1101000;
			fluidProp[i][VCRIT] = 1253;
			fluidProp[i][AC] = 0.9550;
			fluidProp[i][PARACHOR] = 939.55;
			fluidProp[i][ZCRIT] = 0.2095;
		}
		else if (!strcmp(str, "nC23")) {
			fluidProp[i][MW] = 324.630;
			fluidProp[i][TCRIT] = 801.3;
			fluidProp[i][PCRIT] = 1059000;
			fluidProp[i][VCRIT] = 1307;
			fluidProp[i][AC] = 0.9890;
			fluidProp[i][PARACHOR] = 981.43;
			fluidProp[i][ZCRIT] = 0.2078;
		}
		else if (!strcmp(str, "nC24")) {
			fluidProp[i][MW] = 338.680;
			fluidProp[i][TCRIT] = 810.4;
			fluidProp[i][PCRIT] = 1019000;
			fluidProp[i][VCRIT] = 1362;
			fluidProp[i][AC] = 1.0190;
			fluidProp[i][PARACHOR] = 1023.40;
			fluidProp[i][ZCRIT] = 0.2061;
		}
		else if (!strcmp(str, "CO2")) {
			fluidProp[i][MW] = 44.010;
			fluidProp[i][TCRIT] = 304.19;
			fluidProp[i][PCRIT] = 7382000;
			fluidProp[i][VCRIT] = 94;
			fluidProp[i][AC] = 0.2276;
			fluidProp[i][PARACHOR] = 82.00;
			fluidProp[i][ZCRIT] = 0.2744;
		}
		else if (!strcmp(str, "O2")) {
			fluidProp[i][MW] = 31.999;
			fluidProp[i][TCRIT] = 154.58;
			fluidProp[i][PCRIT] = 5043000;
			fluidProp[i][VCRIT] = 73.4;
			fluidProp[i][AC] = 0.0218;
			fluidProp[i][ZCRIT] = 0.2880;
		}
		else if (!strcmp(str, "N2")) {
			fluidProp[i][MW] = 28.014;
			fluidProp[i][TCRIT] = 126.1;
			fluidProp[i][PCRIT] = 3394000;
			fluidProp[i][VCRIT] = 90.1;
			fluidProp[i][AC] = 0.0403;
			fluidProp[i][PARACHOR] = 61.12;
			fluidProp[i][ZCRIT] = 0.2917;
		}
		else if (!strcmp(str, "H2S")) {
			fluidProp[i][MW] = 34.082;
			fluidProp[i][TCRIT] = 373.53;
			fluidProp[i][PCRIT] = 8963000;
			fluidProp[i][VCRIT] = 98.5;
			fluidProp[i][AC] = 0.0827;
			fluidProp[i][PARACHOR] = 85.50;
			fluidProp[i][ZCRIT] = 0.2843;
		}
		else if (!strcmp(str, "SO2")) {
			fluidProp[i][MW] = 64.065;
			fluidProp[i][TCRIT] = 430.75;
			fluidProp[i][PCRIT] = 7884000;
			fluidProp[i][VCRIT] = 122;
			fluidProp[i][AC] = 0.2451;
			fluidProp[i][ZCRIT] = 0.2686;
		}
		else if (!strcmp(str, "H2")) {
			fluidProp[i][MW] = 2.016;
			fluidProp[i][TCRIT] = 33.18;
			fluidProp[i][PCRIT] = 1313000;
			fluidProp[i][VCRIT] = 64.2;
			fluidProp[i][AC] = 0.2150;
			fluidProp[i][ZCRIT] = 0.3053;
		}
		else if (!strcmp(str, "H2O")) {
			fluidProp[i][MW] = 18.015;
			fluidProp[i][TCRIT] = 647.13;
			fluidProp[i][PCRIT] = 22055000;
			fluidProp[i][VCRIT] = 56;
			fluidProp[i][AC] = 0.3449;
			fluidProp[i][ZCRIT] = 0.2294;
		}
		else {
			TerM("Unknown Component!");
		}
	}
	EOS_Init();

	fp.FileSearch("SWT");
	fp.ReadWord(str);
	for (i = 0; i<Nswt; i++) for (j = 0; j<SAT_TABLE; j++) {
		if (!fp.ReadWord(str)) TerM("Incorrect SWT keyword format in the input file!");
		tempL = atof(str);
		if (tempL) swt[i][j] = tempL;
		else swt[i][j] = RELPERM0;
		//tempL=atof(str);
	}


	fp.FileSearch("SGT");
	fp.ReadWord(str);
	for (i = 0; i<Nsgt; i++) for (j = 0; j<SAT_TABLE; j++) {
		if (!fp.ReadWord(str)) TerM("Incorrect SGT keyword format in the input file!");
		tempL = atof(str);
		if (tempL) sgt[i][j] = tempL;
		else sgt[i][j] = RELPERM0;
		//sgt[i][j]=atof(str);
	}

	for (i = 0; i<Nc; i++)
		for (j = 0; j<Nc; j++) bic[i][j] = 1;
	if (!fp.FileSearch("BIC")) {
		if (!MPIRank) std::cout << "Warning: No BIC keyword in the input file!\n";
	}
	else {
		for (i = 1; i<Nc; i++)
			for (j = 0; j<i; j++) {
				if (!fp.ReadWord(str)) TerM("Incorrect BIC keyword format in the input file!");
				tempL = 1 - atof(str);
				bic[i][j] = tempL;
				bic[j][i] = tempL;
			}
	}

	if (!fp.FileSearch("INITCOND")) TerM("No INITCOND keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect INITCOND keyword format in the input file!");
	initCond = atoi(str);
	switch (initCond) {
	case 0:		//all
		if (!fp.FileSearch("IPRESS")) TerM("No IPRESS keyword in the input file!");		//initial pressure
		eclint = 0;
		for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
			if (!eclint) {
				if (!fp.ReadWord(str)) TerM("Incorrect IPRESS keyword format in the input file!");
				ECLStar(str, &eclint, &tempL);
			}
			P[i + 1][j + 1][k + 1][1] = tempL;
			eclint--;
		}

		if (!fp.FileSearch("FIPRESS")) TerM("No Fracture FIPRESS keyword in the input file!");		//initial Fracture pressure
		eclint = 0;
		for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
			if (!eclint) {
				if (!fp.ReadWord(str)) TerM("Incorrect Fracture IPRESS keyword format in the input file!");
				ECLStar(str, &eclint, &tempL);
			}
			fP[i + 1][j + 1][k + 1][1] = tempL;
			eclint--;
		}

		if (!fp.FileSearch("IWS")) TerM("No IWS keyword in the input file!");		//Initial water saturation
		eclint = 0;
		for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
			if (!eclint) {
				if (!fp.ReadWord(str)) TerM("Incorrect IWS keyword format in the input file!");
				ECLStar(str, &eclint, &tempL);
			}
			sat[i][j][k][0] = tempL;
			eclint--;
		}

		if (!fp.FileSearch("FIWS")) TerM("No FIWS keyword in the input file!");		//Fractured Initial water saturation
		eclint = 0;
		for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
			if (!eclint) {
				if (!fp.ReadWord(str)) TerM("Incorrect FIWS keyword format in the input file!");
				ECLStar(str, &eclint, &tempL);
			}
			fsat[i][j][k][0] = tempL;
			eclint--;
		}

		if (!fp.FileSearch("IGC")) TerM("No IGC keyword in the input file!");		//Initial global composition
		eclint = 0;
		for (n = 0; n < Nc; n++) for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
			if (!eclint) {
				if (!fp.ReadWord(str)) TerM("Incorrect IGC keyword format in the input file!");
				ECLStar(str, &eclint, &tempL);
			}
			comp[i][j][k][n][2] = tempL;
			eclint--;
		}

		if (!fp.FileSearch("FIGC")) TerM("No Fracture IGC keyword in the input file!");		//Initial Fracture global composition
		eclint = 0;
		for (n = 0; n < Nc; n++) for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) {
			if (!eclint) {
				if (!fp.ReadWord(str)) TerM("Incorrect Fracture IGC keyword format in the input file!");
				ECLStar(str, &eclint, &tempL);
			}
			fcomp[i][j][k][n][2] = tempL;
			eclint--;
		}
		AllFlash();
		break;

	case 1:		//same as 0 but depth variation only
		if (!fp.FileSearch("IPRESS")) TerM("No IPRESS keyword in the input file!");		//initial pressure
		for (k = 0; k < Nz; k++) {
			if (!fp.ReadWord(str)) TerM("Incorrect IPRESS keyword format in the input file!");
			tempL = atof(str);
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) P[i + 1][j + 1][k + 1][1] = tempL;
		}

		if (!fp.FileSearch("FIPRESS")) TerM("No FIPRESS keyword in the input file!");		//Fracture initial pressure
		for (k = 0; k < Nz; k++) {
			if (!fp.ReadWord(str)) TerM("Incorrect FIPRESS keyword format in the input file!");
			tempL = atof(str);
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fP[i + 1][j + 1][k + 1][1] = tempL;
		}

		if (!fp.FileSearch("IWS")) TerM("No IWS keyword in the input file!");		//Initial water saturation
		for (k = 0; k < Nz; k++) {
			if (!fp.ReadWord(str)) TerM("Incorrect IWS keyword format in the input file!");
			tempL = atof(str);
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) sat[i][j][k][0] = atof(str);
		}

		if (!fp.FileSearch("FIWS")) TerM("No FIWS keyword in the input file!");		//Fracture Initial water saturation
		for (k = 0; k < Nz; k++) {
			if (!fp.ReadWord(str)) TerM("Incorrect FIWS keyword format in the input file!");
			tempL = atof(str);
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fsat[i][j][k][0] = atof(str);
		}

		if (!fp.FileSearch("IGC")) TerM("No IGC keyword in the input file!");		//Initial global composition
		for (k = 0; k < Nz; k++)
			for (n = 0; n < Nc; n++) {
				if (!fp.ReadWord(str)) TerM("Incorrect IGC keyword format in the input file!");
				tempL = atof(str);
				for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) comp[i][j][k][n][2] = tempL;
			}

		if (!fp.FileSearch("FIGC")) TerM("No FIGC keyword in the input file!");		//Fractured Initial global composition
		for (k = 0; k < Nz; k++)
			for (n = 0; n < Nc; n++) {
				if (!fp.ReadWord(str)) TerM("Incorrect FIGC keyword format in the input file!");
				tempL = atof(str);
				for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fcomp[i][j][k][n][2] = tempL;
			}

		AllFlash();
		break;

	case 2:		//compositional grading
		if (!fp.FileSearch("REFLAYER")) TerM("No REFLAYER keyword in the input file!");
		if (!fp.ReadWord(str)) TerM("Incorrect REFLAYER keyword format in the input file!");
		refL = atoi(str);

		if (!fp.FileSearch("REFPRES")) TerM("No REFPRES keyword in the input file!");
		if (!fp.ReadWord(str)) TerM("Incorrect REFPRES keyword format in the input file!");
		tempL = atof(str);
		for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) P[i + 1][j + 1][refL + 1][1] = tempL;

		if (!fp.FileSearch("WOCHEIGHT")) TerM("No WOCHEIGHT keyword in the input file!");
		if (!fp.ReadWord(str)) TerM("Incorrect WOCHEIGHT keyword format in the input file!");
		WOCHeight = atof(str);

		if (!fp.FileSearch("REFCOMP")) TerM("No REFCOMP keyword in the input file!");
		for (n = 0; n < Nc; n++) {
			if (!fp.ReadWord(str)) TerM("Incorrect REFCOMP keyword format in the input file!");
			tempL = atof(str);
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) comp[i][j][refL][n][2] = tempL;
		}
		//GCE_NR();  //compistional grading module ->not exist yet
		break;



	case 3:		//All the same
		if (!fp.FileSearch("IPRESS")) TerM("No IPRESS keyword in the input file!");		//initial pressure
		if (!fp.ReadWord(str)) TerM("Incorrect IPRESS keyword format in the input file!");
		tempL = atof(str);
		for (k = 0; k < Nz; k++) {
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) P[i + 1][j + 1][k + 1][1] = tempL;
		}

		if (!fp.FileSearch("FIPRESS")) TerM("No FIPRESS keyword in the input file!");		//Fracture initial pressure
		if (!fp.ReadWord(str)) TerM("Incorrect FIPRESS keyword format in the input file!");
		tempL = atof(str);
		for (k = 0; k < Nz; k++) {
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fP[i + 1][j + 1][k + 1][1] = tempL;
		}

		if (!fp.FileSearch("IWS")) TerM("No IWS keyword in the input file!");		//Initial water saturation
		if (!fp.ReadWord(str)) TerM("Incorrect IWS keyword format in the input file!");
		tempL = atof(str);
		for (k = 0; k < Nz; k++) {
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) sat[i][j][k][0] = atof(str);
		}

		if (!fp.FileSearch("FIWS")) TerM("No FIWS keyword in the input file!");		//Fracture Initial water saturation
		if (!fp.ReadWord(str)) TerM("Incorrect FIWS keyword format in the input file!");
		tempL = atof(str);
		for (k = 0; k < Nz; k++) {
			for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fsat[i][j][k][0] = atof(str);
		}

		if (!fp.FileSearch("IGC")) TerM("No IGC keyword in the input file!");		//Initial global composition
		for (n = 0; n < Nc; n++) {
			if (!fp.ReadWord(str)) TerM("Incorrect IGC keyword format in the input file!");
			tempL = atof(str);
			for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) comp[i][j][k][n][2] = tempL;
		}

		if (!fp.FileSearch("FIGC")) TerM("No FIGC keyword in the input file!");		//Fractured Initial global composition
		for (n = 0; n < Nc; n++) {
			if (!fp.ReadWord(str)) TerM("Incorrect FIGC keyword format in the input file!");
			tempL = atof(str);
			for (k = 0; k < Nz; k++) for (j = 0; j < Ny; j++) for (i = 0; i < Nx; i++) fcomp[i][j][k][n][2] = tempL;
		}

		AllFlash();
		break;
	}




	/*
	Well Type
	Well X block
	Well Y block
	Well Z starting block		Well Z ending block

	Enter all well rates in positive form, the sign will be determined automatically by the well type
	*/
	if (!fp.FileSearch("WELLS")) {
		if (!MPIRank) std::cout << "Warning: No WELL keyword in the input file!\n";
	}
	fp.ReadWord(str);
	n = 0;
	for (i = 0; i < wellNO; i++) {
		for (j = 0; j < WELL_I; j++) {
			if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
			welli[i][j] = atoi(str) - 1;
		}
		(welli[i][0])++;
		switch (welli[i][0]) {
		case 0:		//Constant Pwf- Oil Production
			for (j = 0; j < 6; j++) {		//Pwf, MaxQo, Skin, Rw, T_Start, T_End, Set MaxQo=0 to lift the limit
				if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
				wellf[i][j] = atof(str);

				if (j > 3) {
					wellf[i][j] *= TIMECHFACT;
					TStepMarker[n] = wellf[i][j];
					n++;
				}
			}
			break;

		case 1:		//Constant Pwf- Gas Production
			for (j = 0; j < 6; j++) {		//Pwf, MaxQg, Skin, Rw, T_Start, T_End
				if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
				wellf[i][j] = atof(str);

				if (j > 3) {
					wellf[i][j] *= TIMECHFACT;
					TStepMarker[n] = wellf[i][j];
					n++;
				}
			}
			break;

		case 2:		//Constant Pwf- Water injection
			for (j = 0; j < 6; j++) {		//Pwf, MaxQw, Skin, Rw, T_Start, T_End
				if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
				wellf[i][j] = atof(str);

				if (j > 3) {
					wellf[i][j] *= TIMECHFACT;
					TStepMarker[n] = wellf[i][j];
					n++;
				}
			}

			break;

		case 3:		//Constant Pwf- Gas injection
			for (j = 0; j < (Nc + 6); j++) {		//Pwf, MaxQg, Skin, Rw, Nc Composition, T_Start, T_End
				if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
				wellf[i][j] = atof(str);
				if (j > (3 + Nc)) {
					wellf[i][j] *= TIMECHFACT;
					TStepMarker[n] = wellf[i][j];
					n++;
				}
			}

			break;

		case 4:		//Constant Production Oil Rate
			for (j = 0; j < 6; j++) {		//Qo, MinPwf, Skin, Rw, T_Start, T_End
				if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
				wellf[i][j] = atof(str);

				if (j > 3) {
					wellf[i][j] *= TIMECHFACT;
					TStepMarker[n] = wellf[i][j];
					n++;
				}
			}

			break;

		case 5:		//Constant Production Gas Rate
			for (j = 0; j < 6; j++) {		//Qg, MinPwf, Skin, Rw, T_Start, T_End
				if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
				wellf[i][j] = atof(str);

				if (j > 3) {
					wellf[i][j] *= TIMECHFACT;
					TStepMarker[n] = wellf[i][j];
					n++;
				}
			}

			break;

		case 6:		//Constant Production Total Rate
			for (j = 0; j < 6; j++) {		//QT, MinPwf, Skin, Rw, T_Start, T_End
				if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
				wellf[i][j] = atof(str);

				if (j > 3) {
					wellf[i][j] *= TIMECHFACT;
					TStepMarker[n] = wellf[i][j];
					n++;
				}
			}

			break;

		case 7:		//Constant Water Injection Rate
			for (j = 0; j < 6; j++) {		//Qw, MaxPwf, Skin, Rw, T_Start, T_End
				if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
				wellf[i][j] = atof(str);

				if (j > 3) {
					wellf[i][j] *= TIMECHFACT;
					TStepMarker[n] = wellf[i][j];
					n++;
				}
			}
			break;

		case 8:		//Constant Gas Injection Rate
			for (j = 0; j < (Nc + 6); j++) {		//Qg, MaxPwf, Skin, Rw, Nc Components, T_Start, T_End
				if (!fp.ReadWord(str)) TerM("Incorrect WELL keyword format in the input file!");
				wellf[i][j] = atof(str);

				if (j > (Nc + 3)) {
					wellf[i][j] *= TIMECHFACT;
					TStepMarker[n] = wellf[i][j];
					n++;
				}
			}

			break;

		default:
			TerM("Incorrect WELL keyword format in the input file!");
		}
	}


	//Time step and total time
	if (!fp.FileSearch("DT")) TerM("No DT keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect DT keyword format in the input file!");
	Dt = atof(str)*TIMECHFACT;		//Hour In
	if (!fp.FileSearch("TOTALT")) TerM("No TOTALT keyword in the input file!");
	if (!fp.ReadWord(str)) TerM("Incorrect TOTALT keyword format in the input file!");
	//totalTime = atof(str)*24*TIMECHFACT;		//Hour In
	totalTime = atof(str)*TIMECHFACT;		//Hour In
	TStepMarker[2 * wellNO] = totalTime;



	
}

void Allocation(void) {
	register int i, j, k, n;

#ifdef USE_MKL_SOLVER
	n = 2*(2 * Nc + 4) * Nx*Ny*Nz;
	if ((MKLCoefMatrix = (FType *)malloc(n*n * sizeof(FType))) == NULL) TerM("Can not allocate memory for MKLCoefMatrix");
	if ((MKLRHSVector = (FType *)malloc(n * sizeof(FType))) == NULL) TerM("Can not allocate memory for MKLRHSVector");	
#endif

	if ((OptionsDBFirst = (char **)malloc(NOOptionsDB * sizeof(char *))) == NULL) TerM("Can not allocate memory for option data base");
	for (i = 0; i < NOOptionsDB; i++) {
		if ((OptionsDBFirst[i] = (char *)malloc(MAX_STRING_LENGTH* sizeof(char))) == NULL) TerM("Can not allocate memory for option data base");
	}
	if ((OptionsDBSecond = (char **)malloc(NOOptionsDB * sizeof(char *))) == NULL) TerM("Can not allocate memory for option data base");
	for (i = 0; i < NOOptionsDB; i++) {
		if ((OptionsDBSecond[i] = (char *)malloc(MAX_STRING_LENGTH * sizeof(char))) == NULL) TerM("Can not allocate memory for option data base");
	}

	if ((SolvedXs = (FType *)malloc(2*(2*Nc+4) * Nx*Ny*Nz * sizeof(FType))) == NULL) TerM("Can not allocate memory for SolvedXs");
	//if ((ScalingVector = (FType *)malloc(Nx*Ny*Nz*(2*Nc+4) * sizeof(FType))) == NULL) TerM("Can not allocate memory for ScalingVector");
	//if ((ScalingVectorBackUp = (FType *)malloc(Nx*Ny*Nz*(2 * Nc + 4) * sizeof(FType))) == NULL) TerM("Can not allocate memory for ScalingVectorBackUp");
	if ((SchindlerList = (int *)malloc(3*Nx*Ny*Nz * sizeof(int))) == NULL) TerM("Can not allocate memory for SchindlerList");
	if ((ShiftArray = (int *)malloc(40 * Nx*Ny*Nz * sizeof(int))) == NULL) TerM("Can not allocate memory for ShiftArray");
	//if ((OffDiagonalElements = (int *)malloc(2 * Nx*Ny*Nz * sizeof(int))) == NULL) TerM("Can not allocate memory for OffDiagonalElements");
	//if ((OnDiagonalElements = (int *)malloc(2 * Nx*Ny*Nz * sizeof(int))) == NULL) TerM("Can not allocate memory for OffDiagonalElements");
	
	//if ((AllStates = (char *)malloc(2 * Nx*Ny*Nz * sizeof(char))) == NULL) TerM("Can not allocate memory for AllStates");
	//if ((StatesCopy = (char *)malloc(2 * Nx*Ny*Nz * sizeof(char))) == NULL) TerM("Can not allocate memory for StatesCopy");
	if ((PressureRow = (int *)malloc((2*Nx*Ny*Nz+1) * sizeof(int))) == NULL) TerM("Can not allocate memory for PressureRow");


	if ((WellCondition = (bool *)malloc(wellNO * sizeof(bool))) == NULL) TerM("Can not allocate memory for WellCondition");
	if ((bWellCondition = (bool *)malloc(wellNO * sizeof(bool))) == NULL) TerM("Can not allocate memory for WellCondition");
	if ((TempWellCondition = (bool *)malloc(wellNO * sizeof(bool))) == NULL) TerM("Can not allocate memory for WellCondition");
	if ((SadeqSize = (int *)malloc(4*(Nx*Ny*Nz) * sizeof(int))) == NULL) TerM("Can not allocate memory for SadeqSize");

	if ((IComp = (FType *)malloc(Nc * sizeof(FType))) == NULL) TerM("Can not allocate memory for ICopm");
	if ((IOComp = (FType *)malloc(Nc * sizeof(FType))) == NULL) TerM("Can not allocate memory for IOCopm");
	if ((IGComp = (FType *)malloc(Nc * sizeof(FType))) == NULL) TerM("Can not allocate memory for IGCopm");
	
	bCSRSize = Nx*Ny*Nz*(18 * Nc*Nc + 35 * Nc + 19) - 2 * (2 * Nc*Nc + 4 * Nc + 2)*(Nx*Ny + Ny*Nz + Nz*Nx);

	if ((gridDim = (FType *)malloc((Nx + Ny + Nz) * sizeof(FType))) == NULL) TerM("Can not allocate memory for gridDim");

	if ((porosity = (FType ***)malloc(Nx * sizeof(FType **))) == NULL) TerM("Can not allocate memory for porosity");
	for (i = 0; i<Nx; i++) {
		if ((porosity[i] = (FType **)malloc(Ny * sizeof(FType *))) == NULL) TerM("Can not allocate memory for porosity");
		for (j = 0; j<Ny; j++) if ((porosity[i][j] = (FType *)malloc(Nz * sizeof(FType))) == NULL) TerM("Can not allocate memory for porosity");
	}

	if ((fporosity = (FType ***)malloc(Nx * sizeof(FType **))) == NULL) TerM("Can not allocate memory for fporosity");
	for (i = 0; i < Nx; i++) {
		if ((fporosity[i] = (FType **)malloc(Ny * sizeof(FType *))) == NULL) TerM("Can not allocate memory for fporosity");
		for (j = 0; j < Ny; j++) if ((fporosity[i][j] = (FType *)malloc(Nz * sizeof(FType))) == NULL) TerM("Can not allocate memory for fporosity");
	}

	if ((Bift = (FType ***)malloc(Nx * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Bift");
	for (i = 0; i<Nx; i++) {
		if ((Bift[i] = (FType **)malloc(Ny * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Bift");
		for (j = 0; j<Ny; j++) if ((Bift[i][j] = (FType *)malloc(Nz * sizeof(FType))) == NULL) TerM("Can not allocate memory for Bift");
	}

	if ((fBift = (FType ***)malloc(Nx * sizeof(FType **))) == NULL) TerM("Can not allocate memory for fracture Bift");
	for (i = 0; i < Nx; i++) {
		if ((fBift[i] = (FType **)malloc(Ny * sizeof(FType *))) == NULL) TerM("Can not allocate memory for fracture Bift");
		for (j = 0; j < Ny; j++) if ((fBift[i][j] = (FType *)malloc(Nz * sizeof(FType))) == NULL) TerM("Can not allocate memory for fracture Bift");
	}

	if ((tor = (FType ***)malloc(Nx * sizeof(FType **))) == NULL) TerM("Can not allocate memory for tortuosity");
	for (i = 0; i<Nx; i++) {
		if ((tor[i] = (FType **)malloc(Ny * sizeof(FType *))) == NULL) TerM("Can not allocate memory for tortuosity");
		for (j = 0; j<Ny; j++) if ((tor[i][j] = (FType *)malloc(Nz * sizeof(FType))) == NULL) TerM("Can not allocate memory for tortuosity");
	}

	if ((perm = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for perm");
	for (i = 0; i<Nx; i++) {
		if ((perm[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for perm");
		for (j = 0; j<Ny; j++) {
			if ((perm[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for perm");
			for (k = 0; k<Nz; k++) if ((perm[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for perm");
		}
	}

	if ((fperm = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for fperm");
	for (i = 0; i < Nx; i++) {
		if ((fperm[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for fperm");
		for (j = 0; j < Ny; j++) {
			if ((fperm[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for fperm");
			for (k = 0; k < Nz; k++) if ((fperm[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for fperm");
		}
	}

	if ((fluidProp = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for components");
	for (i = 0; i<Nc; i++) {
		if ((fluidProp[i] = (FType *)malloc(COMP_PROPS * sizeof(FType))) == NULL) TerM("Can not allocate memory for components");
	}


	if ((swt = (FType **)malloc(Nswt * sizeof(FType *))) == NULL) TerM("Can not allocate memory for liquid saturation table");
	for (i = 0; i<Nswt; i++) {
		if ((swt[i] = (FType *)malloc(SAT_TABLE * sizeof(FType))) == NULL) TerM("Can not allocate memory for liquid saturation table");
	}

	if ((sgt = (FType **)malloc(Nsgt * sizeof(FType *))) == NULL) TerM("Can not allocate memory for gas saturation table");
	for (i = 0; i<Nsgt; i++) {
		if ((sgt[i] = (FType *)malloc(SAT_TABLE * sizeof(FType))) == NULL) TerM("Can not allocate memory for gas saturation table");
	}

	if ((comp = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for composition");
	for (i = 0; i<Nx; i++) {
		if ((comp[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for composition");
		for (j = 0; j<Ny; j++) {
			if ((comp[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for composition");
			for (k = 0; k<Nz; k++) {
				if ((comp[i][j][k] = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for composition");
				for (n = 0; n<Nc; n++) if ((comp[i][j][k][n] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for composition");
			}
		}
	}

	if ((fcomp = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for Fracture composition");
	for (i = 0; i < Nx; i++) {
		if ((fcomp[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture  composition");
		for (j = 0; j < Ny; j++) {
			if ((fcomp[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture composition");
			for (k = 0; k < Nz; k++) {
				if ((fcomp[i][j][k] = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture composition");
				for (n = 0; n < Nc; n++) if ((fcomp[i][j][k][n] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture composition");
			}
		}
	}
	 
	if ((EdgeComp = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for EdgeComp");
	for (i = 0; i < Nc; i++) {
		if ((EdgeComp[i] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for EdgeComp");
	}


	if ((f_m_relperm = (FType *)malloc((3) * sizeof(FType))) == NULL) TerM("Can not allocate memory for f_m_relperm");

	if ((f_m_diffusion = (FType *)malloc((3) * sizeof(FType))) == NULL) TerM("Can not allocate memory for f_m_diffusion");

	if ((f_m_drelperm = (FType **)malloc(3 * sizeof(FType *))) == NULL) TerM("Can not allocate memory for f_m_drelperm");
	for (i = 0; i < 3; i++) {
		if ((f_m_drelperm[i] = (FType *)malloc(10 * sizeof(FType))) == NULL) TerM("Can not allocate memory for f_m_drelperm");
	}


	if ((diffusion = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for diffusion");
	for (i = 0; i<Nx; i++) {
		if ((diffusion[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for diffusion");
		for (j = 0; j<Ny; j++) {
			if ((diffusion[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for diffusion");
			for (k = 0; k<Nz; k++) {
				if ((diffusion[i][j][k] = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for diffusion");
				for (n = 0; n<Nc; n++) if ((diffusion[i][j][k][n] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for diffusion");
			}
		}
	}

	if ((fdiffusion = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for fracture  diffusion");
	for (i = 0; i < Nx; i++) {
		if ((fdiffusion[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for fracture  diffusion");
		for (j = 0; j < Ny; j++) {
			if ((fdiffusion[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for fracture diffusion");
			for (k = 0; k < Nz; k++) {
				if ((fdiffusion[i][j][k] = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for fracture diffusion");
				for (n = 0; n < Nc; n++) if ((fdiffusion[i][j][k][n] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for fracture diffusion");
			}
		}
	}

	if ((P = (FType ****)malloc((Nx + 2) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for press");
	for (i = 0; i<(Nx + 2); i++) {
		if ((P[i] = (FType ***)malloc((Ny + 2) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for press");
		for (j = 0; j<(Ny + 2); j++) {
			if ((P[i][j] = (FType **)malloc((Nz + 2) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for press");
			for (k = 0; k<(Nz + 2); k++) if ((P[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for press");
		}
	}

	if ((fP = (FType ****)malloc((Nx + 2) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture press");
	for (i = 0; i < (Nx + 2); i++) {
		if ((fP[i] = (FType ***)malloc((Ny + 2) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture press");
		for (j = 0; j < (Ny + 2); j++) {
			if ((fP[i][j] = (FType **)malloc((Nz + 2) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture press");
			for (k = 0; k < (Nz + 2); k++) if ((fP[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture press");
		}
	}


	if ((sat = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for saturation");
	for (i = 0; i<Nx; i++) {
		if ((sat[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for saturation");
		for (j = 0; j<Ny; j++) {
			if ((sat[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for saturation");
			for (k = 0; k<Nz; k++) if ((sat[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for saturation");
		}
	}

	if ((fsat = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture saturation");
	for (i = 0; i < Nx; i++) {
		if ((fsat[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture saturation");
		for (j = 0; j < Ny; j++) {
			if ((fsat[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture saturation");
			for (k = 0; k < Nz; k++) if ((fsat[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture saturation");
		}
	}

	if ((bsat = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for saturation backup");
	for (i = 0; i<Nx; i++) {
		if ((bsat[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for saturation backup");
		for (j = 0; j<Ny; j++) {
			if ((bsat[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for saturation backup");
			for (k = 0; k<Nz; k++) if ((bsat[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for saturation backup");
		}
	}

	if ((fbsat = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture saturation backup");
	for (i = 0; i < Nx; i++) {
		if ((fbsat[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture saturation backup");
		for (j = 0; j < Ny; j++) {
			if ((fbsat[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture saturation backup");
			for (k = 0; k < Nz; k++) if ((fbsat[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture saturation backup");
		}
	}

	if ((bcomp = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for composition");
	for (i = 0; i<Nx; i++) {
		if ((bcomp[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for composition");
		for (j = 0; j<Ny; j++) {
			if ((bcomp[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for composition");
			for (k = 0; k<Nz; k++) {
				if ((bcomp[i][j][k] = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for composition");
				for (n = 0; n<Nc; n++) if ((bcomp[i][j][k][n] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for composition");
			}
		}
	}

	if ((fbcomp = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for Fracture composition");
	for (i = 0; i < Nx; i++) {
		if ((fbcomp[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture composition");
		for (j = 0; j < Ny; j++) {
			if ((fbcomp[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture composition");
			for (k = 0; k < Nz; k++) {
				if ((fbcomp[i][j][k] = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture composition");
				for (n = 0; n < Nc; n++) if ((fbcomp[i][j][k][n] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture composition");
			}
		}
	}

	if ((bP = (FType ****)malloc((Nx + 2) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for bpress");
	for (i = 0; i<(Nx + 2); i++) {
		if ((bP[i] = (FType ***)malloc((Ny + 2) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for bpress");
		for (j = 0; j<(Ny + 2); j++) {
			if ((bP[i][j] = (FType **)malloc((Nz + 2) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for bpress");
			for (k = 0; k<(Nz + 2); k++) if ((bP[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for bpress");
		}
	}

	if ((fbP = (FType ****)malloc((Nx + 2) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture bpress");
	for (i = 0; i < (Nx + 2); i++) {
		if ((fbP[i] = (FType ***)malloc((Ny + 2) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture bpress");
		for (j = 0; j < (Ny + 2); j++) {
			if ((fbP[i][j] = (FType **)malloc((Nz + 2) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture bpress");
			for (k = 0; k < (Nz + 2); k++) if ((fbP[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture bpress");
		}
	}

	if (wellNO) {
		if ((wellf = (FType **)malloc(wellNO * sizeof(FType *))) == NULL) TerM("Can not allocate memory for wells data");
		for (i = 0; i<wellNO; i++) {
			if ((wellf[i] = (FType *)malloc((Nc + 6) * sizeof(FType))) == NULL) TerM("Can not allocate memory for wells data");
		}
		if ((welli = (int **)malloc(wellNO * sizeof(int *))) == NULL) TerM("Can not allocate memory for wells data");
		for (i = 0; i<wellNO; i++) {
			if ((welli[i] = (int *)malloc(WELL_I * sizeof(int))) == NULL) TerM("Can not allocate memory for wells data");
		}

	}


	//////////////////////////////////////////////////////
	/////////////////////NON_INPUT////////////////////////
	//////////////////////////////////////////////////////
	if ((IFTran = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for IFTran");
	for (i = 0; i<Nx; i++) {
		if ((IFTran[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for IFTran");
		for (j = 0; j<Ny; j++) {
			if ((IFTran[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for IFTran");
			for (k = 0; k<Nz; k++) if ((IFTran[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for IFTran");
		}
	}

	if ((fIFTran = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture IFTran");
	for (i = 0; i < Nx; i++) {
		if ((fIFTran[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture IFTran");
		for (j = 0; j < Ny; j++) {
			if ((fIFTran[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture IFTran");
			for (k = 0; k < Nz; k++) if ((fIFTran[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture IFTran");
		}
	}

	if ((blockFProps = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for blockFProps");
	for (i = 0; i<Nx; i++) {
		if ((blockFProps[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for blockFProps");
		for (j = 0; j<Ny; j++) {
			if ((blockFProps[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for blockFProps");
			for (k = 0; k<Nz; k++) {
				if ((blockFProps[i][j][k] = (FType **)malloc(BLOCK_F_PROPS * sizeof(FType *))) == NULL) TerM("Can not allocate memory for blockFProps");
				for (n = 0; n<BLOCK_F_PROPS; n++) if ((blockFProps[i][j][k][n] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for blockFProps");
			}
		}
	}

	if ((fblockFProps = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for Fracture blockFProps");
	for (i = 0; i < Nx; i++) {
		if ((fblockFProps[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture blockFProps");
		for (j = 0; j < Ny; j++) {
			if ((fblockFProps[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture blockFProps");
			for (k = 0; k < Nz; k++) {
				if ((fblockFProps[i][j][k] = (FType **)malloc(BLOCK_F_PROPS * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture blockFProps");
				for (n = 0; n < BLOCK_F_PROPS; n++) if ((fblockFProps[i][j][k][n] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture blockFProps");
			}
		}
	}

	if ((trans = (FType *****)malloc((Nx + 1) * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for trans");
	for (i = 0; i<(Nx + 1); i++) {
		if ((trans[i] = (FType ****)malloc((Ny + 1) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for trans");
		for (j = 0; j<(Ny + 1); j++) {
			if ((trans[i][j] = (FType ***)malloc((Nz + 1) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for trans");
			for (k = 0; k<(Nz + 1); k++) {
				if ((trans[i][j][k] = (FType **)malloc(12 * sizeof(FType*))) == NULL) TerM("Can not allocate memory for trans");
				for (n = 0; n<12; n++) if ((trans[i][j][k][n] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for trans");
			}
		}
	}

	if ((ftrans = (FType *****)malloc((Nx + 1) * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for Fracture trans");
	for (i = 0; i < (Nx + 1); i++) {
		if ((ftrans[i] = (FType ****)malloc((Ny + 1) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture trans");
		for (j = 0; j < (Ny + 1); j++) {
			if ((ftrans[i][j] = (FType ***)malloc((Nz + 1) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture trans");
			for (k = 0; k < (Nz + 1); k++) {
				if ((ftrans[i][j][k] = (FType **)malloc(12 * sizeof(FType*))) == NULL) TerM("Can not allocate memory for Fracture trans");
				for (n = 0; n < 12; n++) if ((ftrans[i][j][k][n] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture trans");
			}
		}
	}

	if ((transS = (char *****)malloc((Nx + 1) * sizeof(char ****))) == NULL) TerM("Can not allocate memory for transS");
	for (i = 0; i<(Nx + 1); i++) {
		if ((transS[i] = (char ****)malloc((Ny + 1) * sizeof(char ***))) == NULL) TerM("Can not allocate memory for transS");
		for (j = 0; j<(Ny + 1); j++) {
			if ((transS[i][j] = (char ***)malloc((Nz + 1) * sizeof(char **))) == NULL) TerM("Can not allocate memory for transS");
			for (k = 0; k<(Nz + 1); k++) {
				if ((transS[i][j][k] = (char **)malloc(3 * sizeof(char*))) == NULL) TerM("Can not allocate memory for transS");
				for (n = 0; n<3; n++) if ((transS[i][j][k][n] = (char *)malloc(3 * sizeof(char))) == NULL) TerM("Can not allocate memory for transS");
			}
		}
	}

	if ((ftransS = (char *****)malloc((Nx + 1) * sizeof(char ****))) == NULL) TerM("Can not allocate memory for Fracture transS");
	for (i = 0; i < (Nx + 1); i++) {
		if ((ftransS[i] = (char ****)malloc((Ny + 1) * sizeof(char ***))) == NULL) TerM("Can not allocate memory for Fracture transS");
		for (j = 0; j < (Ny + 1); j++) {
			if ((ftransS[i][j] = (char ***)malloc((Nz + 1) * sizeof(char **))) == NULL) TerM("Can not allocate memory for Fracture transS");
			for (k = 0; k < (Nz + 1); k++) {
				if ((ftransS[i][j][k] = (char **)malloc(3 * sizeof(char*))) == NULL) TerM("Can not allocate memory for Fracture transS");
				for (n = 0; n < 3; n++) if ((ftransS[i][j][k][n] = (char *)malloc(3 * sizeof(char))) == NULL) TerM("Can not allocate memory for Fracture transS");
			}
		}
	}

	if ((Wtran = (FType ****)malloc((Nx + 1) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Wtran");
	for (i = 0; i<(Nx + 1); i++) {
		if ((Wtran[i] = (FType ***)malloc((Ny + 1) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Wtran");
		for (j = 0; j<(Ny + 1); j++) {
			if ((Wtran[i][j] = (FType **)malloc((Nz + 1) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Wtran");
			for (k = 0; k<(Nz + 1); k++) if ((Wtran[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Wtran");
		}
	}

	if ((fWtran = (FType ****)malloc((Nx + 1) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture Wtran");
	for (i = 0; i < (Nx + 1); i++) {
		if ((fWtran[i] = (FType ***)malloc((Ny + 1) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture Wtran");
		for (j = 0; j < (Ny + 1); j++) {
			if ((fWtran[i][j] = (FType **)malloc((Nz + 1) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture Wtran");
			for (k = 0; k < (Nz + 1); k++) if ((fWtran[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture Wtran");
		}
	}

	if ((dWtran = (FType ****)malloc((Nx + 1) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for dWtran");
	for (i = 0; i<(Nx + 1); i++) {
		if ((dWtran[i] = (FType ***)malloc((Ny + 1) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for dWtran");
		for (j = 0; j<(Ny + 1); j++) {
			if ((dWtran[i][j] = (FType **)malloc((Nz + 1) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for dWtran");
			for (k = 0; k<(Nz + 1); k++) if ((dWtran[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for dWtran");
		}
	}

	if ((fdWtran = (FType ****)malloc((Nx + 1) * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture dWtran");
	for (i = 0; i < (Nx + 1); i++) {
		if ((fdWtran[i] = (FType ***)malloc((Ny + 1) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture dWtran");
		for (j = 0; j < (Ny + 1); j++) {
			if ((fdWtran[i][j] = (FType **)malloc((Nz + 1) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture dWtran");
			for (k = 0; k < (Nz + 1); k++) if ((fdWtran[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture dWtran");
		}
	}

	if ((relPerm = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for relPerm");
	for (i = 0; i<Nx; i++) {
		if ((relPerm[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for relPerm");
		for (j = 0; j<Ny; j++) {
			if ((relPerm[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for relPerm");
			for (k = 0; k<Nz; k++) if ((relPerm[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for relPerm");
		}
	}

	if ((frelPerm = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture relPerm");
	for (i = 0; i < Nx; i++) {
		if ((frelPerm[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture relPerm");
		for (j = 0; j < Ny; j++) {
			if ((frelPerm[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture relPerm");
			for (k = 0; k < Nz; k++) if ((frelPerm[i][j][k] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture relPerm");
		}
	}

	if ((dRelPerm = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for dRelPerm");
	for (i = 0; i<Nx; i++) {
		if ((dRelPerm[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for dRelPerm");
		for (j = 0; j<Ny; j++) {
			if ((dRelPerm[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for dRelPerm");
			for (k = 0; k<Nz; k++) if ((dRelPerm[i][j][k] = (FType *)malloc(5 * sizeof(FType))) == NULL) TerM("Can not allocate memory for dRelPerm");
		}
	}

	if ((fdRelPerm = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture dRelPerm");
	for (i = 0; i < Nx; i++) {
		if ((fdRelPerm[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture dRelPerm");
		for (j = 0; j < Ny; j++) {
			if ((fdRelPerm[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture dRelPerm");
			for (k = 0; k < Nz; k++) if ((fdRelPerm[i][j][k] = (FType *)malloc(5 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture dRelPerm");
		}
	}

	/*if ((preCon=(FType **) malloc(Nx*Ny*Nz*(2*Nc+4)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for Jacobian matrix");
	for (i=0;i<(Nx*Ny*Nz*(2*Nc+4));i++) {
	if ((preCon[i]=(FType *) malloc(Nx*Ny*Nz*(2*Nc+4)*sizeof(FType)))==NULL) TerM("Can not allocate memory for Jacobian matrix");
	}*/

	if ((CSRjac = (FType *)malloc(3*bCSRSize * sizeof(FType))) == NULL) TerM("Can not allocate memory for Sparse Jacobian matrix");
	if ((CSRrow = (int *)malloc((3*Nx*Ny*Nz*(2 * Nc + 4) + 1) * sizeof(int))) == NULL) TerM("Can not allocate memory for Sparse Jacobian Row Index matrix");
	if ((CSRcol = (int *)malloc(3*bCSRSize * sizeof(int))) == NULL) TerM("Can not allocate memory for Sparse Jacobian Column Index matrix");

	if ((preCon = (FType *)malloc(3*(Nx*Ny*Nz*((2 * Nc + 4)*(2 * Nc + 3) + MRMAXNONZERO)) * sizeof(FType))) == NULL) TerM("Can not allocate memory for preconditioner matrix");
	if ((preConCSRrow = (int *)malloc(3*(Nx*Ny*Nz*(2 * Nc + 4) + 1) * sizeof(int))) == NULL) TerM("Can not allocate memory for preConCSRrow matrix");
	if ((preConCSRcol = (int *)malloc(3*(Nx*Ny*Nz*((2 * Nc + 4)*(2 * Nc + 3) + MRMAXNONZERO)) * sizeof(int))) == NULL) TerM("Can not allocate memory for preConCSRrow matrix");


	if ((preConRow = (int *)malloc(3*(Nx*Ny*Nz + 1) * sizeof(int))) == NULL) TerM("Can not allocate memory for preConRow matrix");
	if ((preConIndex = (int *)malloc(3*(Nx*Ny*Nz + 1) * sizeof(int))) == NULL) TerM("Can not allocate memory for preConIndex matrix");


	if ((ans = (FType *)malloc(3*Nx*Ny*Nz*(2 * Nc + 4) * sizeof(FType))) == NULL) TerM("Can not allocate memory for answer matrix");
	//if ((Unk=(FType *) malloc(Nx*Ny*Nz*(2*Nc+4)*sizeof(FType)))==NULL) TerM("Can not allocate memory for unknown matrix");

	if ((pJHolder = (int ***)malloc(Nx * sizeof(int **))) == NULL) TerM("Can not allocate memory for pJHolder");
	for (i = 0; i<Nx; i++) {
		if ((pJHolder[i] = (int **)malloc(Ny * sizeof(int *))) == NULL) TerM("Can not allocate memory for pJHolder");
		for (j = 0; j<Ny; j++) if ((pJHolder[i][j] = (int *)malloc(Nz * sizeof(FType))) == NULL) TerM("Can not allocate memory for pJHolder");
	}

	if ((fpJHolder = (int ***)malloc(Nx * sizeof(int **))) == NULL) TerM("Can not allocate memory for fpJHolder");
	for (i = 0; i < Nx; i++) {
		if ((fpJHolder[i] = (int **)malloc(Ny * sizeof(int *))) == NULL) TerM("Can not allocate memory for fpJHolder");
		for (j = 0; j < Ny; j++) if ((fpJHolder[i][j] = (int *)malloc(Nz * sizeof(FType))) == NULL) TerM("Can not allocate memory for fpJHolder");
	}


	if ((preProp = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for preProp");
	for (i = 0; i<Nx; i++) {
		if ((preProp[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for preProp");
		for (j = 0; j<Ny; j++) {
			if ((preProp[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for preProp");
			for (k = 0; k<Nz; k++) if ((preProp[i][j][k] = (FType *)malloc((Nc + 1) * sizeof(FType))) == NULL) TerM("Can not allocate memory for preProp");
		}
	}

	if ((fpreProp = (FType ****)malloc(Nx * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for fracture preProp");
	for (i = 0; i < Nx; i++) {
		if ((fpreProp[i] = (FType ***)malloc(Ny * sizeof(FType **))) == NULL) TerM("Can not allocate memory for fracture preProp");
		for (j = 0; j < Ny; j++) {
			if ((fpreProp[i][j] = (FType **)malloc(Nz * sizeof(FType *))) == NULL) TerM("Can not allocate memory for fracture preProp");
			for (k = 0; k < Nz; k++) if ((fpreProp[i][j][k] = (FType *)malloc((Nc + 1) * sizeof(FType))) == NULL) TerM("Can not allocate memory for fracture preProp");
		}
	}

	if ((dE = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for dE");
	for (i = 0; i<Nx; i++) {
		if ((dE[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for dE");
		for (j = 0; j<Ny; j++) {
			if ((dE[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for dE");
			for (k = 0; k<Nz; k++) {
				if ((dE[i][j][k] = (FType **)malloc((Nc + 1) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for dE");
				for (n = 0; n<(Nc + 1); n++) if ((dE[i][j][k][n] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for dE");
			}
		}
	}

	if ((fdE = (FType *****)malloc(Nx * sizeof(FType ****))) == NULL) TerM("Can not allocate memory for Fracture dE");
	for (i = 0; i < Nx; i++) {
		if ((fdE[i] = (FType ****)malloc(Ny * sizeof(FType ***))) == NULL) TerM("Can not allocate memory for Fracture dE");
		for (j = 0; j < Ny; j++) {
			if ((fdE[i][j] = (FType ***)malloc(Nz * sizeof(FType **))) == NULL) TerM("Can not allocate memory for Fracture dE");
			for (k = 0; k < Nz; k++) {
				if ((fdE[i][j][k] = (FType **)malloc((Nc + 1) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Fracture dE");
				for (n = 0; n < (Nc + 1); n++) if ((fdE[i][j][k][n] = (FType *)malloc(2 * sizeof(FType))) == NULL) TerM("Can not allocate memory for Fracture dE");
			}
		}
	}

	if ((satJac = (FType **)malloc((Nc + 1) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for Flash Jacobian matrix");
	for (i = 0; i<(Nc + 1); i++) {
		if ((satJac[i] = (FType *)malloc((Nc + 1) * sizeof(FType))) == NULL) TerM("Can not allocate memory for Flash Jacobian matrix");
	}

	if ((satAns = (FType *)malloc((Nc + 1) * sizeof(FType))) == NULL) TerM("Can not allocate memory for Flash Answer matrix");
	if ((Xm = (FType *)malloc((Nc + 1) * sizeof(FType))) == NULL) TerM("Can not allocate memory for Flash Composition matrix");
	if ((Xms = (FType *)malloc((Nc + 1) * sizeof(FType))) == NULL) TerM("Can not allocate memory for Flash Composition save matrix");

	if ((phaseStat = (signed char ***)malloc(Nx * sizeof(signed char **))) == NULL) TerM("Can not allocate memory for phaseStat");
	for (i = 0; i<Nx; i++) {
		if ((phaseStat[i] = (signed char **)malloc(Ny * sizeof(signed char *))) == NULL) TerM("Can not allocate memory for phaseStat");
		for (j = 0; j<Ny; j++) if ((phaseStat[i][j] = (signed char *)malloc(Nz * sizeof(signed char))) == NULL) TerM("Can not allocate memory for phaseStat");
	}

	if ((fphaseStat = (signed char ***)malloc(Nx * sizeof(signed char **))) == NULL) TerM("Can not allocate memory for Fracture phaseStat");
	for (i = 0; i < Nx; i++) {
		if ((fphaseStat[i] = (signed char **)malloc(Ny * sizeof(signed char *))) == NULL) TerM("Can not allocate memory for Fracture phaseStat");
		for (j = 0; j < Ny; j++) if ((fphaseStat[i][j] = (signed char *)malloc(Nz * sizeof(signed char))) == NULL) TerM("Can not allocate memory for Fracture phaseStat");
	}

	if ((bphaseStat = (signed char ***)malloc(Nx * sizeof(signed char **))) == NULL) TerM("Can not allocate memory for bphaseStat");
	for (i = 0; i<Nx; i++) {
		if ((bphaseStat[i] = (signed char **)malloc(Ny * sizeof(signed char *))) == NULL) TerM("Can not allocate memory for bphaseStat");
		for (j = 0; j<Ny; j++) if ((bphaseStat[i][j] = (signed char *)malloc(Nz * sizeof(signed char))) == NULL) TerM("Can not allocate memory for bphaseStat");
	}

	if ((fbphaseStat = (signed char ***)malloc(Nx * sizeof(signed char **))) == NULL) TerM("Can not allocate memory for Fracture bphaseStat");
	for (i = 0; i < Nx; i++) {
		if ((fbphaseStat[i] = (signed char **)malloc(Ny * sizeof(signed char *))) == NULL) TerM("Can not allocate memory for Fracture bphaseStat");
		for (j = 0; j < Ny; j++) if ((fbphaseStat[i][j] = (signed char *)malloc(Nz * sizeof(signed char))) == NULL) TerM("Can not allocate memory for Fracture bphaseStat");
	}


	if ((Ki = (FType *)malloc(2*Nc * sizeof(FType))) == NULL) TerM("Can not allocate memory for Initial Flash K matrix");

	if ((blockH = (FType ***)malloc((Nx + 2) * sizeof(FType **))) == NULL) TerM("Can not allocate memory for blockH");
	for (i = 0; i<(Nx + 2); i++) {
		if ((blockH[i] = (FType **)malloc((Ny + 2) * sizeof(FType *))) == NULL) TerM("Can not allocate memory for blockH");
		for (j = 0; j<(Ny + 2); j++) if ((blockH[i][j] = (FType *)malloc((Nz + 2) * sizeof(FType))) == NULL) TerM("Can not allocate memory for blockH");
	}	

	if ((bic = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for binary interaction coefficient matrix");
	for (i = 0; i<Nc; i++) {
		if ((bic[i] = (FType *)malloc(Nc * sizeof(FType))) == NULL) TerM("Can not allocate memory for binary interaction coefficient matrix");
	}

	if ((TStepMarker = (FType *)malloc((2 * wellNO + 1) * sizeof(FType))) == NULL) TerM("Can not allocate memory for Time Step Marker");

	if ((STcomp = (FType **)malloc(Nc * sizeof(FType *))) == NULL) TerM("Can not allocate memory for STcomp");
	for (i = 0; i<Nc; i++) {
		if ((STcomp[i] = (FType *)malloc(3 * sizeof(FType))) == NULL) TerM("Can not allocate memory for STcomp");
	}

	/*k=(2*Nc+4)*Nx*Ny*Nz;
	if ((fullPre=(FType **) malloc(k*sizeof(FType *)))==NULL) TerM("Can not allocate memory for fullPre matrix");
	for (i=0;i<k;i++) {
	if ((fullPre[i]=(FType *) malloc(k*sizeof(FType)))==NULL) TerM("Can not allocate memory for fullPre matrix");
	}*/

	if ((Jac1Row = (FType *)malloc(8 * (2 * Nc + 4) * sizeof(FType))) == NULL) TerM("Can not allocate memory for Jac1Row");
	if ((Jac1Col = (int *)malloc(8 * (2 * Nc + 4) * sizeof(int))) == NULL) TerM("Can not allocate memory for Jac1Col");


	
}

void ECLStar(char *str, int *times, FType *mvalue) {
	register int i, j;
	char Sind = 0;
	char str1[MAX_STRING_LENGTH], str2[MAX_STRING_LENGTH];

	i = 0;
	j = 0;
	while (str[i] != '\0') {
		if (str[i] == '*') {
			str1[i] = '\0';
			Sind = -1;
		}
		else {
			if (Sind) {
				str2[j] = str[i];
				j++;
			}
			else {
				str1[i] = str[i];
			}
		}
		i++;
	}

	if (Sind) {
		str2[j] = '\0';
		*times = atoi(str1);
		*mvalue = atof(str2);
	}
	else {
		str1[i] = '\0';
		*times = 1;
		*mvalue = atof(str1);
	}
}

void RestoreRST(FType *simTime) {
	register int Ix, Iy, Iz, i;
	std::ifstream fp;

	//if ((fp=fopen("test.rst","rb"))!=NULL) {
	fp.open("test.rst", std::ios::binary);
	if (fp.is_open()) {

		for (Ix = 0; Ix<(Nx + 2); Ix++)
			for (Iy = 0; Iy<(Ny + 2); Iy++)
				for (Iz = 0; Iz<(Nz + 2); Iz++)
					for (i = 0; i<3; i++) {
						//fread(&P[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
						fp.read((char *)&P[Ix][Iy][Iz][i], sizeof(FType));
						//if (ferror(fp)) puts("Error in File in!");
					}

		for (Ix = 0; Ix<Nx; Ix++)
			for (Iy = 0; Iy<Ny; Iy++)
				for (Iz = 0; Iz<Nz; Iz++) {
					for (i = 0; i<3; i++) {
						//fread(&sat[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
						fp.read((char *)&sat[Ix][Iy][Iz][i], sizeof(FType));
						//if (ferror(fp)) puts("Error in File in!");
					}

					for (i = 0; i<Nc; i++) {
						//fread(&comp[Ix][Iy][Iz][i][0], sizeof(FType), 1, fp);
						fp.read((char *)&comp[Ix][Iy][Iz][i][0], sizeof(FType));
						//if (ferror(fp)) puts("Error in File in!");
						//fread(&comp[Ix][Iy][Iz][i][1], sizeof(FType), 1, fp);
						fp.read((char *)&comp[Ix][Iy][Iz][i][1], sizeof(FType));
						//if (ferror(fp)) puts("Error in File in!");
					}
					//fread(&phaseStat[Ix][Iy][Iz], sizeof(char), 1, fp);
					fp.read((char *)&phaseStat[Ix][Iy][Iz], sizeof(char));
					//if (ferror(fp)) puts("Error in File in!");
				}
		//fread(simTime, sizeof(FType), 1, fp);
		fp.read((char *)simTime, sizeof(FType));
		//if (ferror(fp)) puts("Error in File in!");
		//fread(&Dt, sizeof(FType), 1, fp);
		fp.read((char *)&Dt, sizeof(FType));
		//if (ferror(fp)) puts("Error in File in!");
		//fread(&TSMCtrl, sizeof(int), 1, fp);
		fp.read((char *)&TSMCtrl, sizeof(int));
		//if (ferror(fp)) puts("Error in File in!");

		//fread(&SumQoProduced, sizeof(FType), 1, fp);
		//fread(&SumQgInjected, sizeof(FType), 1, fp);
		fp.read((char *)&SumQoProduced, sizeof(FType));
		fp.read((char *)&SumQgInjected, sizeof(FType));
		//fclose(fp);
		fp.close();

	}
}
