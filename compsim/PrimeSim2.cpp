#include <iostream>
#include <string>
#include <fstream>
#include "MIfstream.h"
#include "Functions.h"
#include <time.h>
#include <iomanip>

#ifdef USE_MKL_SOLVER


int main(int argc, char *argv[]) {
	MIfstream fin;
	std::ofstream ftaout, ftaoutS;
	std::ofstream ftest;
	FType simTime = 0;
	//register int NRCount;
	FType mb, runTime;
	clock_t start, stop;
	char simTStat;
	int Time_Print = 0;
	FType fssum;// , ppo, ppg;

	int NRCount;
	register int i;
	AppCtx LocalStruct;
	Vec PETScVecDummy;
	Mat PETScMatDummy;
	FType tolMain, tolMain2;

	start = clock();
	SetUpPetscEnv(&argc, &argv, 2 * Nx*Ny*Nz*(2 * Nc + 4), 2 * Nc + 4);
	fin.open("test.dat", std::ios::in);
	if (!fin.is_open()) TerM("Can not open input file!");	
	Data_Input(fin);	
	fin.close();
	MPIRank = 0;
	MPIsize = 1;

	
	if (!MPIRank) ftaout.open("test_out_all.dat", std::ios::out);
	if (!MPIRank) if (!ftaout.is_open()) TerM("Can not open input file!");

	if (!MPIRank) ftaoutS.open("test_out_Sadeq.dat", std::ios::out);
	if (!MPIRank) if (!ftaoutS.is_open()) TerM("Can not open input file!");
	
	if (!MPIRank) ftest.open("test_test.dat", std::ios::out);
	if (!MPIRank) if (!ftest.is_open()) TerM("Can not open output file!");


	ManageTSMarker();
	IFK_Calc();
	Tr_Mu();
	QInitValue();
	Edge_Trans();
	Edge_fTrans();
	InitWells();

	//TotalNXYZ = Nx * Ny*Nz;

	WellCompSadeq_cum[0] = 0;
	WellCompSadeq_cum[1] = 0;
	WellCompSadeq_cum[2] = 0;


	mb = MatBal(0);
	PrintAll(ftaout, 0.0);
	CreateBackups(simTime);
	UpdateDependencies(true);
	//ppg = blockFProps[0][0][0][RO][1];
	//ppo = fblockFProps[0][0][0][RO][0];

	
	do {
		
		//CreateBackups(simTime);
		SSimTime = simTime;
		NRCount = 0;
		
		LocalStruct.SimulationTime = simTime;
		tolMain2 = 0;
		do {			
			UpdateDependencies(false);
			tolMain2 = 0;
			
			for (i = 0; i < MatrixSize; i++) {
				LocalStruct.GlobalIndex = i;
				LocalStruct.LocalIndex = i;
				IteratorFunction(PETScVecDummy, PETScVecDummy, &LocalStruct);
				IteratorJacobian(PETScVecDummy, PETScMatDummy, &LocalStruct);
#ifdef USE_F_NORM
				tolMain2 += abs(MKLRHSVector[i]);
#endif
			}
				
					
			MKLS2(MKLCoefMatrix, MKLRHSVector, MatrixSize);
			if (repStat) break;

			tolMain = UpdateAllForMKL(MKLRHSVector, false, 1);
#ifdef USE_F_NORM
			tolMain = tolMain2 / MatrixSize;
#endif

			
			NRCount++;
			if (NRCount>MAX_NR_CYCLE) {
				//TerM("Newton-Raphson Iteration counter exceeded its limit!");		
				//std::cout<<("Newton-Raphson Iteration counter exceeded its limit!");
				if (Dt<(MINDT)) {
					repStat = -1;
					maxNR = -1;
				}
				else repStat = -1;
				break;
			}			
		} while (tolMain>TOLMAIN);

		

		simTStat = TimeStepCtrl(&simTime);
		if (simTStat) {
			CreateBackups(simTime);
			UpdateDependencies(true);
			
			std::cout << "\r" << simTime * 100 / totalTime << "Completed!";
		

			
			//if ((simTime * 10000 / totalTime) >= Time_Print) {
				Time_Print++;
				ftaout << "\n" << mb << "\n";
				ftaout << "\nNR=" << NRCount << "\n";

				for (int ixx = 0; ixx < Nx; ixx++) {
					for (int iyy = 0; iyy < Ny; iyy++) {
						ftaoutS << std::setprecision(20) << comp[ixx][iyy][Nz / 2 + 1][0][0] << '\t';
					}
					ftaoutS << std::endl;
				}
				ftaoutS << std::endl << std::endl << std::endl;
				/*
				ftest << std::setprecision(20) << simTime / 60 << "\t";// << AvgPReport() << "\t" << Dt << "\t" << -cumQo << "\t" << cumQg << "\t" << -SumQoProduced << "\t" << SumQgInjected << "\t" << wGLR << "\t" << nesbat << "\t" << MatBal(simTime) << "\t" << -Qoprod_cum << "\t" << Qginj_cum;
			

				for (unsigned int kkk = 0; kkk < Nc; kkk++) {
					fssum = 0;
					for (unsigned int jjj = 0; jjj < Nx; jjj++) fssum += blockFProps[jjj][0][0][RO][1] * 6.242796057614462e-05*comp[jjj][0][0][kkk][1];
					ftest << std::setprecision(20) << fssum / Nx << '\t';
					//ftest << std::setprecision(20) << blockFProps[0][0][0][RO][1] * 6.242796057614462e-05*comp[0][0][0][kkk][1] << '\t';
				}
				for (unsigned int kkk = 0; kkk < Nc; kkk++) {
					fssum = 0;
					for (unsigned int jjj = 0; jjj < Nx; jjj++) fssum += fblockFProps[jjj][0][0][RO][0] * 6.242796057614462e-05*fcomp[jjj][0][0][kkk][0];
					ftest << std::setprecision(20) << fssum / Nx << '\t';
					//ftest << std::setprecision(20) << fblockFProps[0][0][0][RO][0] * 6.242796057614462e-05*fcomp[0][0][0][kkk][0] << '\t';
				}
				

				ftest << "\n";
				ftest.flush();*/
				PrintAll(ftaout, simTime);
				//PrintAll2(0, 0, 0, ftaoutS, simTime);
			//}

			Qoprod_cum += Qoprod*Dt;
			Qginj_cum += Qginj*Dt;

			for (unsigned int kkkk = 0; kkkk < 3; kkkk++) {
				WellCompSadeq2[kkkk] /= (Qoprod);
			}

			WellCompSadeq_cum[0] += WellCompSadeq[0] * Dt;
			WellCompSadeq_cum[1] += WellCompSadeq[1] * Dt;
			WellCompSadeq_cum[2] += WellCompSadeq[2] * Dt;




			mb = MatBal(simTime);

			
			


		}		

	} while (simTStat != 1);

	stop = clock();
	runTime = (FType)(stop - start) / CLOCKS_PER_SEC;
	PrintAll(ftaout, simTime);

	if (!MPIRank) ftaout << "\n\n\tRunTime=" << runTime;

	///////////////////////////////////////////////////////////////
	/*ftest << std::setprecision(20) << simTime / 60 << "\t";// << AvgPReport() << "\t" << Dt << "\t" << -cumQo << "\t" << cumQg << "\t" << -SumQoProduced << "\t" << SumQgInjected << "\t" << wGLR << "\t" << nesbat << "\t" << MatBal(simTime) << "\t" << -Qoprod_cum << "\t" << Qginj_cum;


	for (unsigned int kkk = 0; kkk < Nc; kkk++) {
		fssum = 0;
		for (unsigned int jjj = 0; jjj < Nx; jjj++) fssum += blockFProps[jjj][0][0][RO][1] * 6.242796057614462e-05*comp[jjj][0][0][kkk][1];
		ftest << std::setprecision(20) << fssum / Nx << '\t';
		//ftest << std::setprecision(20) << blockFProps[0][0][0][RO][1] * 6.242796057614462e-05*comp[0][0][0][kkk][1] << '\t';
	}
	for (unsigned int kkk = 0; kkk < Nc; kkk++) {
		fssum = 0;
		for (unsigned int jjj = 0; jjj < Nx; jjj++) fssum += fblockFProps[jjj][0][0][RO][0] * 6.242796057614462e-05*fcomp[jjj][0][0][kkk][0];
		ftest << std::setprecision(20) << fssum / Nx << '\t';
		//ftest << std::setprecision(20) << fblockFProps[0][0][0][RO][0] * 6.242796057614462e-05*fcomp[0][0][0][kkk][0] << '\t';
	}


	ftest << "\n";
	ftest.flush();*/
	PrintAll(ftaout, simTime);
	PrintAll2(0, 0, 0, ftaoutS, simTime);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////







	if (!MPIRank) ftaout.close();
	if (!MPIRank) ftest.close();
	if (!MPIRank) ftaoutS.close();
	PetscFinishAll();
	TerM("\nFinished Successfully!");
	return 0;
}


FType UpdateAllForMKL(FType *XNew, bool FullTimeStep, FType Alpha) {
	int Ix, Iy, Iz, i, BLNo;
	FType normSuml, normSumg, rr;
	
	rr = 0;
	BLNo = 0;
	for (Iz = 0; Iz < Nz; Iz++)
		for (Iy = 0; Iy < Ny; Iy++)
			for (Ix = 0; Ix < Nx; Ix++) {
				if (phaseStat[Ix][Iy][Iz] == 1) {
					normSuml = 0;
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][0] += Alpha*XNew[BLNo];
						rr += fabs(XNew[BLNo]);
						if (comp[Ix][Iy][Iz][i][0] < EPS_COMP) comp[Ix][Iy][Iz][i][0] = 0;
						else if (comp[Ix][Iy][Iz][i][0] > 1) comp[Ix][Iy][Iz][i][0] = 1;

						normSuml += comp[Ix][Iy][Iz][i][0];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][0] /= normSuml;
					}

					//if ((fabs(PCoeff*Alpha*XNew[BLNo])) > (P[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) P[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF*dSgn(XNew[BLNo]);
					P[Ix + 1][Iy + 1][Iz + 1][1] += Alpha * XNew[BLNo] * PCoeff;
					//P[Ix + 1][Iy + 1][Iz + 1][1] = 11e6;
					rr += fabs(XNew[BLNo]);
					BLNo++;

					sat[Ix][Iy][Iz][0] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
					BLNo++;

					sat[Ix][Iy][Iz][1] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
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
						comp[Ix][Iy][Iz][i][1] += Alpha * XNew[BLNo];
						rr += fabs(XNew[BLNo]);
						if (comp[Ix][Iy][Iz][i][1] < EPS_COMP) comp[Ix][Iy][Iz][i][1] = 0;
						else if (comp[Ix][Iy][Iz][i][1] > 1) comp[Ix][Iy][Iz][i][1] = 1;
						normSumg += comp[Ix][Iy][Iz][i][1];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][1] /= normSumg;
					}

					//if ((fabs(PCoeff*Alpha*XNew[BLNo])) > (P[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) P[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF*dSgn(XNew[BLNo]);
					P[Ix + 1][Iy + 1][Iz + 1][1] += Alpha * XNew[BLNo] * PCoeff;
					//P[Ix + 1][Iy + 1][Iz + 1][1] = 11e6;
					rr += fabs(XNew[BLNo]);
					BLNo++;

					sat[Ix][Iy][Iz][0] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
					BLNo++;

					sat[Ix][Iy][Iz][2] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
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
						comp[Ix][Iy][Iz][i][0] += Alpha * XNew[BLNo];
						rr += fabs(XNew[BLNo]);
						if (comp[Ix][Iy][Iz][i][0] < EPS_COMP) comp[Ix][Iy][Iz][i][0] = 0;
						else if (comp[Ix][Iy][Iz][i][0] > 1) comp[Ix][Iy][Iz][i][0] = 1;
						normSuml += comp[Ix][Iy][Iz][i][0];
						BLNo++;
					}
					normSumg = 0;
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][1] += Alpha * XNew[BLNo];
						rr += fabs(XNew[BLNo]);
						if (comp[Ix][Iy][Iz][i][1] < EPS_COMP) comp[Ix][Iy][Iz][i][1] = 0;
						else if (comp[Ix][Iy][Iz][i][1] > 1) comp[Ix][Iy][Iz][i][1] = 1;
						normSumg += comp[Ix][Iy][Iz][i][1];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						comp[Ix][Iy][Iz][i][0] /= normSuml;
						comp[Ix][Iy][Iz][i][1] /= normSumg;
					}

					//if ((fabs(PCoeff*Alpha*XNew[BLNo])) > (P[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) P[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF*dSgn(XNew[BLNo]);
					P[Ix + 1][Iy + 1][Iz + 1][1] += Alpha * XNew[BLNo] * PCoeff;
					//P[Ix + 1][Iy + 1][Iz + 1][1] = 11e6;
					rr += fabs(XNew[BLNo]);
					BLNo++;

					sat[Ix][Iy][Iz][0] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
					BLNo++;

					sat[Ix][Iy][Iz][1] += Alpha * XNew[BLNo];
					BLNo++;

					sat[Ix][Iy][Iz][2] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
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
						fcomp[Ix][Iy][Iz][i][0] += Alpha * XNew[BLNo];
						rr += fabs(XNew[BLNo]);
						if (fcomp[Ix][Iy][Iz][i][0] < EPS_COMP) fcomp[Ix][Iy][Iz][i][0] = 0;
						else if (fcomp[Ix][Iy][Iz][i][0] > 1) fcomp[Ix][Iy][Iz][i][0] = 1;

						normSuml += fcomp[Ix][Iy][Iz][i][0];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][0] /= normSuml;
					}


					//if ((fabs(PCoeff*Alpha*XNew[BLNo])) > (fP[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) fP[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF*dSgn(XNew[BLNo]);
					fP[Ix + 1][Iy + 1][Iz + 1][1] += Alpha * XNew[BLNo] * PCoeff;
					//fP[Ix + 1][Iy + 1][Iz + 1][1] = 5e6;
					rr += fabs(XNew[BLNo]);
					BLNo++;

					fsat[Ix][Iy][Iz][0] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
					BLNo++;

					fsat[Ix][Iy][Iz][1] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
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
						fcomp[Ix][Iy][Iz][i][1] += Alpha * XNew[BLNo];
						rr += fabs(XNew[BLNo]);
						if (fcomp[Ix][Iy][Iz][i][1] < EPS_COMP) fcomp[Ix][Iy][Iz][i][1] = 0;
						else if (fcomp[Ix][Iy][Iz][i][1] > 1) fcomp[Ix][Iy][Iz][i][1] = 1;
						normSumg += fcomp[Ix][Iy][Iz][i][1];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][1] /= normSumg;
					}
					//if ((fabs(PCoeff*Alpha*XNew[BLNo])) > (fP[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) fP[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF*dSgn(XNew[BLNo]);
					fP[Ix + 1][Iy + 1][Iz + 1][1] += Alpha * XNew[BLNo] * PCoeff;
					//fP[Ix + 1][Iy + 1][Iz + 1][1] = 5e6;
					rr += fabs(XNew[BLNo]);
					BLNo++;

					fsat[Ix][Iy][Iz][0] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
					BLNo++;

					fsat[Ix][Iy][Iz][2] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
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
						fcomp[Ix][Iy][Iz][i][0] += Alpha * XNew[BLNo];
						rr += fabs(XNew[BLNo]);
						if (fcomp[Ix][Iy][Iz][i][0] < EPS_COMP) fcomp[Ix][Iy][Iz][i][0] = 0;
						else if (fcomp[Ix][Iy][Iz][i][0] > 1) fcomp[Ix][Iy][Iz][i][0] = 1;
						normSuml += fcomp[Ix][Iy][Iz][i][0];
						BLNo++;
					}
					normSumg = 0;
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][1] += Alpha * XNew[BLNo];
						rr += fabs(XNew[BLNo]);
						if (fcomp[Ix][Iy][Iz][i][1] < EPS_COMP) fcomp[Ix][Iy][Iz][i][1] = 0;
						else if (fcomp[Ix][Iy][Iz][i][1] > 1) fcomp[Ix][Iy][Iz][i][1] = 1;
						normSumg += fcomp[Ix][Iy][Iz][i][1];
						BLNo++;
					}
					for (i = 0; i < Nc; i++) {
						fcomp[Ix][Iy][Iz][i][0] /= normSuml;
						fcomp[Ix][Iy][Iz][i][1] /= normSumg;
					}
					//if ((fabs(PCoeff*Alpha*XNew[BLNo])) > (fP[Ix + 1][Iy + 1][Iz + 1][1] * (1 + MAX_P_COEFF))) fP[Ix + 1][Iy + 1][Iz + 1][1] += MAX_P_COEFF*dSgn(XNew[BLNo]);
					fP[Ix + 1][Iy + 1][Iz + 1][1] += Alpha * XNew[BLNo] * PCoeff;
					//fP[Ix + 1][Iy + 1][Iz + 1][1] = 5e6;
					rr += fabs(XNew[BLNo]);
					BLNo++;

					fsat[Ix][Iy][Iz][0] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
					BLNo++;

					fsat[Ix][Iy][Iz][1] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
					BLNo++;

					fsat[Ix][Iy][Iz][2] += Alpha * XNew[BLNo];
					rr += fabs(XNew[BLNo]);
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

			//comp[Nx - 1][0][0][0][0] = 1;
			//comp[Nx - 1][0][0][1][0] = 0;

	UpdateDependencies(FullTimeStep);
	return rr / BLNo;
}

#endif // USE_MKL_SOLVER