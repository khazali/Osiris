#include <iostream>
#include<string>
#include<fstream>
#include "MIfstream.h"
#include "Functions.h"
#include <time.h>
#include <iomanip>
#include <mpi.h>
#include <stdio.h>

#ifndef USE_MKL_SOLVER

int main(int argc, char *argv[]) {
	MIfstream fin;
	std::ofstream ftaout, ftaoutS;
	std::ofstream ftest, Favgf;
	FType simTime = 0;
	FType mb, runTime;
	clock_t start, stop;
	char simTStat;
	int Time_Print = 0, PrintAllCounter = 0, RSTCounter = 0, i;
	MPI_Status MPIstatus;
	//char *RatesBArray, sendBArray[2];// , bufferitoa[256], filenameitoa[1024];
	double *RatesBArray;
	bool IsRST = false;
	
	



	start = clock();
	if ((argc == 2) && (!strcmp(argv[1], "-r"))) IsRST = true;
	fin.open("test.dat", std::ios::in);
	if (!fin.is_open()) TerM("Can not open input file!");
	Data_Input(fin);
	//PJac();
	fin.close();
	SetUpPetscEnv(&argc, &argv, 2 * Nx*Ny*Nz*(2 * Nc + 4), 2 * Nc + 4);
	//if (!MPIRank) std::cin.get();
	//PetscBarrier(NULL);
	/*sprintf(bufferitoa, "%d", MPIRank);
	strcpy(filenameitoa, "All_Out");
	strcat(filenameitoa, bufferitoa);
	strcat(filenameitoa, ".dat");*/
	
	if (!MPIRank) ftaout.open("test_out_all.dat", std::ios::out);
	if (!MPIRank) if (!ftaout.is_open()) TerM("Can not open input file!");

	//ftaout.open(filenameitoa, std::ios::out);
	//if (!ftaout.is_open()) TerM("Can not open input file!");

	if (!MPIRank) ftaoutS.open("test_out_Sadeq.dat", std::ios::out);
	if (!MPIRank) if (!ftaoutS.is_open()) TerM("Can not open input file!");	
	if (IsRST) {
		if (!MPIRank) ftest.open("test_test.dat", std::ios::out | std::ios::app);
	}
	else {
		if (!MPIRank) ftest.open("test_test.dat", std::ios::out);
	}
	if (!MPIRank) if (!ftest.is_open()) TerM("Can not open output file!");
	if (!MPIRank) Favgf.open("avg_test.dat", std::ios::out);
	if (!MPIRank) if (!Favgf.is_open()) TerM("Can not open output file!");
	RatesBArray = new double[MPIsize];


	

	ManageTSMarker();
	IFK_Calc();
	Tr_Mu();
	QInitValue();
	Edge_Trans();
    Edge_fTrans();
	InitWells();

	WellCompSadeq_cum[0] = 0;
	WellCompSadeq_cum[1] = 0;
	WellCompSadeq_cum[2] = 0;
	
	SMethQ = 0;

	
	mb = MatBal(0);
	if (IsRST) simTime = ReadRST();
	UpdateDependencies(true);
	CreateBackups(simTime);


	if (!MPIRank) PrintAll(ftaout, 0.0);
	//PrintAll(ftaout, 0.0);
	//CreateBackups(simTime);

	do {
		SSimTime = simTime;		
		if (SolveNonLinear()) {
			repStat = 0;
			CreateBackups(simTime);
		}
		else {
			repStat = 1;
		}
		simTStat = TimeStepCtrl(&simTime);
		
		if (simTStat) {
			//RSTCounter++;
			//if (RSTCounter == 50) {				
				if (!MPIRank) WriteRST(simTime);
				//RSTCounter = 0;
			//}


			//sendBArray[0] = BMethQM;
			//sendBArray[1] = BMethQF;
			
			MPI_Gather(&MethQM, 1, MPI_DOUBLE, RatesBArray, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
			/*if (MPIRank) {
				if (BMethQM) MPI_Send(&MethQM, 1, MPI_DOUBLE, 0, MPITag, PETSC_COMM_WORLD);
				if (BMethQF) MPI_Send(&MethQF, 1, MPI_DOUBLE, 0, MPITag + 1, PETSC_COMM_WORLD);
			}
			else {
				for (i = 0; i < MPIsize; i++) {
					if (RatesBArray[2 * i]) {
						MSaveMPI = i;
						break;
					}
				}
				for (i = 0; i < MPIsize; i++) {
					if (RatesBArray[2 * i + 1]) {
						FSaveMPI = i;
						break;
					}
				}
				//if (MSaveMPI) MPI_Recv(&MethQM, 1, MPI_DOUBLE, MSaveMPI, MPITag, PETSC_COMM_WORLD, &MPIstatus);
				//if (FSaveMPI) MPI_Recv(&MethQF, 1, MPI_DOUBLE, FSaveMPI, MPITag+1, PETSC_COMM_WORLD, &MPIstatus);
				MPI_Recv(&MethQM, 1, MPI_DOUBLE, MSaveMPI, MPITag, PETSC_COMM_WORLD, &MPIstatus);
				MPI_Recv(&MethQF, 1, MPI_DOUBLE, FSaveMPI, MPITag+1, PETSC_COMM_WORLD, &MPIstatus);
			}
			MPITag += 2;*/
			if (!MPIRank) for (i = 0; i < MPIsize; i++) SMethQ -= RatesBArray[i] * Dt;
			//SMethQ -= (MethQM + MethQF) * Dt;
			
			if (!MPIRank) std::cout << "One Time Step!\n";
			if (!MPIRank) std::cout << "\r" << simTime * 100 / totalTime << "Completed!";

			//if ((simTime * 1000000 / totalTime) >= Time_Print) {
				Time_Print++;
				PrintAllCounter++;
				if (!MPIRank) ftaout << "\n" << mb << "\n";
				//ftaout << "\n" << mb << "\n";
				if (!MPIRank) ftest << std::setprecision(20) << simTime / TIMECHFACT << "\t" << AvgPReport() << "\t" << Dt << "\t" << SMethQ;

				if (!MPIRank) ftest << "\n";


				
				if (!MPIRank) ftest.flush();



				if (!MPIRank) PrintAll(ftaout, simTime);
				//PrintAll(ftaout, simTime);
				if (PrintAllCounter == 100) {
					if (!MPIRank) ftaout.close();
					if (!MPIRank) ftaout.open("test_out_all.dat", std::ios::out);
					if (!MPIRank) if (!ftaout.is_open()) TerM("Can not open test_out_all file!");
					PrintAllCounter = 0;
				}
				std::cout << "\n";	
				/*if (!MPIRank) {
					ftaoutS << "\n----------------------------------------------------------------------------------------\n";
					ftaoutS << std::setprecision(20) << "\n\n\ntime:" << simTime * 1000000 / 3600 << "\n\n";
					for (int ixx = 0; ixx < Nx; ixx++) {
						for (int iyy = 0; iyy < Ny; iyy++) {
							ftaoutS << std::setprecision(20) << comp[ixx][iyy][Nz / 2][0][0] << '\t';
						}
						ftaoutS << std::endl;
					}
					ftaoutS << std::endl << std::endl << std::endl;
					ftaoutS << "Line:" << std::endl;
					for (int iyy = 0; iyy < Ny; iyy++) {
						ftaoutS << std::setprecision(20) << comp[Nx/2][iyy][Nz / 2][0][0] << '\t';
					}
					ftaoutS << std::endl;
					ftaoutS << std::endl << std::endl << std::endl;
				}
				//ftaoutS.flush();
				//if (!MPIRank) PrintAll2(0, 0, 0, ftaoutS, simTime);

				FType sum = 0;
				int ttss = 0;
				for (int ixx = 0; ixx < Nx; ixx++) {
					for (int iyy = 0; iyy < Ny; iyy++) {
						for (int izz = 0; izz < Nz; izz++) {
							if (porosity[ixx][iyy][izz] > 1e-10) {
								sum += comp[ixx][iyy][izz][1][0];
								ttss++;
							}
						}
					}
				}
				sum /= ttss;
				Favgf << std::setprecision(20) << simTime * 1000000 / 3600 << '\t' << sqrt(simTime * 1000000 / 3600) << '\t' << sum << std::endl;
				//Favgf.flush();*/

			//}

			Qoprod_cum += Qoprod * Dt;
			Qginj_cum += Qginj * Dt;

			

			mb = MatBal(simTime);

			

		}


	
	} while (simTStat != 1); 
	stop = clock();
	runTime = (FType)(stop - start) / CLOCKS_PER_SEC;
	if (!MPIRank) PrintAll(ftaout, simTime);
	//PrintAll(ftaout, simTime);

	if (!MPIRank) ftaout << "\n\n\tRunTime=" << runTime;
	
	if (!MPIRank) ftaout.close();
	//ftaout.close();
	if (!MPIRank) ftest.close();
	if (!MPIRank) ftaoutS.close();
	if (!MPIRank) Favgf.close();
	delete[] RatesBArray;
	PetscFinishAll();
	TerM("\nFinished Successfully!");
	return 0;
}

#endif
