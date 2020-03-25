#pragma once
//#define USE_MKL_SOLVER
//#define USE_F_NORM
#define MAX_P_COEFF 0.01
#define WELL_MAX_P 14000000
#define	WELL_MIN_P 7000000
#define DIFF_INTERNAL_COEFF 1

#define POR_TRESHOLD 1e-2
#define MAX_EQL_ITER 20
#define PAHSE_ENV_BOUNDARY 0.95
#define EPS_SAT_FRAC 0

#define FType double
#define USE_PETSC

#include <iostream>
#include <fstream>
#include <string>
#include "MIfstream.h"

#ifdef USE_PETSC
#include <petscsnes.h>
/*extern SNES snes;         
extern KSP Petsc_ksp;          
extern PC Petsc_pc;          
extern Vec Petsc_X, Petsc_R, Petsc_F;          
extern Mat Petsc_J;*/          
extern PetscErrorCode ierr;
extern int MPIsize, MPIRank;
extern SNESConvergedReason snesreason;
struct AppCtx {
	double SimulationTime;
	int GlobalIndex;
	int LocalIndex;
}; 
typedef struct AppCtx AppCtx;
PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void *);
PetscErrorCode FormFunction(SNES, Vec, Vec, void *);
extern FType SSimTime;
PetscErrorCode SetUpPetscEnv(int *, char ***, int, int);
PetscErrorCode IteratorFunction(Vec, Vec, AppCtx *);
PetscErrorCode IteratorJacobian(Vec, Mat, AppCtx *);
PetscInt UpdateAll(Vec, bool);
void UpdateDependencies(bool);
PetscErrorCode SolveNonLinear(void);
PetscErrorCode PetscFinishAll(void);
extern AppCtx userData;
extern int MatrixSize;
PetscErrorCode InitGuess(SNES, Vec, void *);
//extern int FraqColBase;
//extern FType *ScalingVector, *ScalingVectorBackUp;

extern int *SchindlerList;
extern int SchindlerListLength;
extern FType *SolvedXs;
extern int *ShiftArray;
//extern int *OffDiagonalElements, *OnDiagonalElements;
//extern int LocalJacobianBlockSize;
extern int NewBlockSizes, SadeqIndex;
//extern int *LocalRowsAIJ;
//extern char *AllStates, *StatesCopy;
int CreateAllocationTable(void);
//void ConstructLocalSizes(int, int);
extern char SolverType[256];
extern FType KSPRelTol, KSPAbsTol, KSPDivTol;
extern PetscInt MaxIters;
extern int NOOptionsDB;
extern char **OptionsDBFirst, **OptionsDBSecond;

#ifdef USE_MKL_SOLVER
extern FType *MKLCoefMatrix, *MKLRHSVector;
#endif

void MixedCalcIFT(void);
void MixedCalcGrav(void);
void MixedCalcDiffusion(void);
void MixedCalcSatFs(void);
void MixedCalcP(void);
void MixedCalcVisco(void);
void MixedCalcTranses(void);
void MixedPreMul(void);
void MixedCalcEquilibrium(void);

PetscErrorCode InitGuess2(SNES, Vec, void *);
void MakeList(int, int);
#endif

#define SURFACE_RATES false

#define MAX_NR_CYCLE 10
#define MAX_STRING_LENGTH 80
#define WELL_I 5		//well integer properties
//#define WELL_F 3		//well FType properties

#define COMP_PROPS 15
#define SG 0
#define TB 1
#define MW 2
#define AC 3
#define PCRIT 4
#define TCRIT 5
#define VCRIT 6
#define TR 7
#define MU 8
#define EOS_A 9
#define EOS_B 10
#define COLDIA 11
#define EPDIA 12
#define PARACHOR 13
#define ZCRIT 14

#define BLOCK_F_PROPS 4
#define RO 0		//molar density
#define BLOCK_PC 1
#define BMU 2
#define BGAM 3



#define SAT_TABLE 4
#define SAT 0
#define KR 1
#define KRO 2
#define CAPILLARYPRESSURE 3
#define MAX_FLASH_N_ITR 100

#define GPU_TEMP_THRESHOLD 80
#define CPU_TEMP_THRESHOLD 80



#ifdef USE_FLOATS

#define FType float

#define G_ACC 9.80665f

#define WAT_M_RO 55.5084f
#define SI_C 9.869232667160128e-13f		// (1e-7/101325)
//Volume=Cubic meter
//Area=Square meter
//Permeability=miliDarcy
//Length=Meter
//Viscosity=centipoise
//pressure=Pascal
//Rate=Cubic meter per second
//Time=Second
#define RGAS 8.3144621f		//Universal gas constant
#define PI 3.1415926535897932384626433f
#define SQRT2 1.4142135623730950488016887f
#define EPSILON 0.001f

#define BTOL 1e-5f

#define XMTOL 0.0005f

#define EPS_SAT 0.001f		//Saturation @ bubble and dew
/*#define DKRGDSG 0;
#define DKRODSO 1;
#define DKRODSG 2;
#define DKRWDS 3;
#define DPCGDSG 4;
#define DPCWDS 5;*/

#define LN10 2.302585092994046f

#define NEWTON_TOL 1e-8f
#define FLASH_TOL 1e-8f

#define DZERO 1e-20f

#define TOLMAIN 1e-5f

#define TIMECHFACT 3600.0f		//Hour

//#define NO_ROOT_L 2000.0f

#define PZERO 1000000.0f

#define PRESSTOL 1000.0f

#define SXMTOL 1e-10f

#define MAXDP 0.01f
#define MAXPP 1e6f
#define MINDT 1000.0f

#define TIME_STEP_INC 1.05f
#define TIME_STEP_DEC 1.7f
#define XMTOLEDGE 1e-5

#else


#define FType double

#define G_ACC 9.80665
#define G_ACCG 9.80665e-3
//#define G_ACC 0
//#define G_ACCG 0
//#define BOLTZMANNK 1.38064852e-23
//#define SPLINESIZE 35


#define WAT_M_RO 55.5084
#define SI_C 9.869232667160128e-13	// (1e-7/101325)
//Volume=Cubic meter
//Area=Square meter
//Permeability=miliDarcy
//Length=Meter
//Viscosity=centipoise
//pressure=Pascal
//Rate=Cubic meter per second
//Time=Second
#define RGAS 8.3144621		//Universal gas constant
//#define BOLTZ 1.380648813e-16		//Boltzmann's Constant (erg/K)
#define PI 3.1415926535897932384626433
#define SQRT2 1.4142135623730950488016887
#define EPSILON 0.001

#define BTOL 1e-100

#define XMTOL 1e-8

#define EPS_SAT 1e-3		//Saturation @ bubble and dew
/*#define DKRGDSG 0;
#define DKRODSO 1;
#define DKRODSG 2;
#define DKRWDS 3;
#define DPCGDSG 4;
#define DPCWDS 5;*/

#define LN10 2.302585092994046

#define NEWTON_TOL 1e-8
#define FLASH_TOL 1e-14

#define DZERO 1e-20

#define TOLMAIN 1e-2

#define TIMECHFACT 3600		//Hour

//#define NO_ROOT_L 10

#define PZERO 1000000

#define PRESSTOL 1000

#define SXMTOL 1e-15

#define MAXDP 1e6
#define MAXPP 100000
#define MINDT 1e-10

#define ZERO_CONV 1
#define ZERO_DIFF 1
#define MINPPP 5000000
//#define CONSTDIFF 1.666e-4

#define TIME_STEP_INC 1.03
#define TIME_STEP_DEC 1.5
#define TSTEPTOL 1
#define TSTEPIF 1e-5
#define RELPERM0 0	//1e-11
#define MAXINITER 5000
#define NTIMESTERINC 5
#define TRIV_K 1e-6
#define MAX_MNR_ITER 20

#define MRINERITER 3
#define MRMAXNONZERO 128
#define BICGSTABTOL	1e-30

#define PSGROUP_NO 56


#define IFTBASECASE 50		//mN/m
#define IFTPOWER 7
#define MISCIBLEIFT 1		//mN/m
#define XMTOLEDGE 1e-6
#define EPS_COMP 0

//#define PCoeff 5e6


//#define PCoeff 1e7

#endif
//typedef double FType;

///////////////////GLOBALS//////////////////////////////////
extern FType PCoeff;
extern int Nx, Ny, Nz;		//Reservoir Dimension
extern FType *gridDim;	//Grid size
extern FType ***porosity;		//Porosity Distribution
extern FType ****perm;		//Permeability tensor distribution
extern FType refP, cpor, dcpor;		//Compressibility
extern int PNc, UNc, Nc;		//Number of Components
extern unsigned char PR, SRK;		//EOS Type, boolean
extern FType **fluidProp;			//Fluids components properties
extern int Nswt, Nsgt;		//Number of saturation tables inputs
extern FType **swt, **sgt;			//Saturation tables
extern int initCond;		//Initial Condition State
extern FType ****sat, ****bsat;		//Saturation and its backup
extern FType ****P;		//Pressure
extern FType *****comp;		//Phases Composition
extern int refL;		//Reference layer in initial condition
extern int wellNO;		//Number of wells
extern int **welli;	//integer well data: X-Well, Y_Well, Well_Type
extern FType **wellf;	//Double well data: Pwf or Q
extern FType ****IFTran;		//Permeability average on blocks boundary, kA/h
extern FType resTemp;		//Reservoir Temperature
extern FType *****blockFProps;		//Calculated fluids properties in each gridblock
extern FType *****trans;		//Transmissibility
extern FType ****Wtran;		//Water Transmissibility
extern FType ****relPerm;		//Relative Permeability
extern FType watRo, watMu;	//Water properties
							//FType **jac;		//Jacobian matrix//
extern FType *ans;		//Answer matrix
extern FType ****preProp;		//Multiplication of properties in previous time step
extern FType Dt;		//Timestep size, seconds
extern FType *****dE;		//Phase density derivative
extern char *****transS;		//upstream weighting condition
extern FType **satJac, *satAns, *Xm, *Xms;	//Minor NR Matrixes
extern signed char ***phaseStat;		//Number and type of the phases in GBs, 1->L, -1->G, 0->2Ph
								//int curSize;		//Jacobian Matrix current size
								//FType *Unk;		//Unknowns matrix
extern FType WOCHeight;	//WOC
extern FType totalTime;	//total simulation time, days input, seconds change
extern FType *Ki;
extern FType QoT, QgT, QwT;		//Well flow rates
extern FType cumQo, cumQg, cumQw;		//cummulative rates
extern FType ***blockH;		//Block Heights
extern int incCount;		//Time Step Control Counter
extern FType *****bcomp;		//Composition Backup
extern FType ****bP;		//Pressure Backup
extern FType Eps;		//Machine Epsilon
extern FType *CSRjac;	//Compressed Sparse Row Jacobian
extern int *CSRrow, *CSRcol;		//Compressed sparse row location holders
extern int bCSRSize;		//Compressed Sparse Row Size for a block
							//int jacIndex;		//Compressed Sparse Row Jacobian Index Holder
							//int CSRrowIndex;		//Compressed Sparse Row Index Holder
extern int ***pJHolder;		//Jacobian Position Holder: Parallel design
//extern FType PCoeff;		//Mean P for Preconditioning
extern FType *preCon;		//Preconditioner
extern int *preConRow;		//Preconditioner row holder
extern int *preConIndex;		//Preconditioner index holder
								// cl_platform_id *platformIds;
								// cl_device_id deviceId;
								// cl_context context;
								// cl_program program;
								// cl_kernel xZero, xZeroF, init, f1, MatMul, f2, f3, clBuildPrecon, dot_persist_kernel, MatMulT, Pf1, Pf2, Pf3, Pf4, dot_persist_kernelF, MatMulTF, MatMulF, ClScaling, clBinPartialSort;
								// cl_command_queue commandQueue;
								//char chWellStat=0;
extern int totalJac, totalRow;
extern FType **bic;
extern FType *TStepMarker;
extern int TSMSize;
extern char repStat;
extern FType **STcomp;
extern FType wellrate[50];
extern signed char ***bphaseStat;
extern int TSMCtrl;
extern char maxNR;
extern FType ****dRelPerm;
extern FType ****dWtran;
extern char gasInjStat;
//FType **fullPre;
extern int *preConCSRrow;
extern int *preConCSRcol;
extern double solverTimeAcc;
extern double solverTime;
//clock_t BICGStart, BICGEnd;
extern FType ***tor, *****diffusion;
extern FType Qw, Qo, Qg;
extern FType wGLR;
extern char buildPreconFlag;
//cl_mem clM, clMCSRrow, clMCSRcol;
extern FType ***Bift;
extern FType QgCumInj, QgInj;
extern FType QoCumProd, QoProd;
extern FType QoProduced, QgInjected, SumQoProduced, SumQgInjected;
extern FType nesbat;

extern long int nprocs;
extern int *CPUTempArray;
extern int *SadeqSize;
extern bool *WellCondition;
extern bool *TempWellCondition;
extern bool *bWellCondition;

extern FType Qginj, Qoprod;
extern FType Qginj_cum;
extern FType Qoprod_cum;

//////////////new Fracture//////////////////////
extern FType ***fporosity;		//Fracture Porosity Distribution
extern FType ****fperm;         //Fracture Permeability tensor distribution
extern FType ****fP;		//Pressure of fracture
extern FType *****fcomp;		//Phases Composition of fracture
extern FType *****fdE;		//fracture Phase density derivative
extern FType *****fblockFProps;     //Calculated fluids properties in each Fracture gridblock 
extern FType ****fsat, ****fbsat;		//Fractured Saturation and its backup
extern FType ****fbP;		//Fracture Pressure Backup
extern FType *****fbcomp;		//Fracture Composition Backup
extern FType ****fIFTran;		//Fracture Permeability average on blocks boundary, kA/h
extern FType *****ftrans;		//Fracture Transmissibility
extern char *****ftransS;		//Fracture upstream weighting condition
extern FType ****fWtran;		//FRacture Water Transmissibility
extern FType ****fdWtran;
extern FType ****frelPerm;		//Fracture Relative Permeability
extern FType ****fdRelPerm;
extern signed char ***fphaseStat;		//Fracture Number and type of the phases in GBs, 1->L, -1->G, 0->2Ph
extern signed char ***fbphaseStat;
extern FType ***fBift;
extern FType ****fpreProp;		//Fracture Multiplication of properties in previous time step
extern FType  *****fdiffusion;
extern int ***fpJHolder;		//Jacobian Position Holder in fractured: Parallel design
extern FType *f_m_relperm; //relative permeability between Fracture And Matrix
extern FType **f_m_drelperm; //moshtagh relative permeability between Fracture And Matrix  0 = dkrwdsw; 1 = dkrodsw; 2 = dkrodsg; 3=dkrgdsg; 4=dkrodso 
extern FType *f_m_diffusion;//Diffusion between Farcture and Matrix
extern FType **EdgeComp; // compesion in interface between matrix & Fracture
extern FType Lf;//Half of Fracture opening
extern FType Re_f;
extern bool IsConstDif;
extern FType CDiffCoef;
extern FType SigmaMF;
extern FType SMethQ, MethQM;
extern char BMethQF, BMethQM;
				//extern int FractureBase;
//extern FType Re;
extern FType *IComp, *IOComp, *IGComp;
extern FType *Jac1Row;
extern int *Jac1Col;
extern int *PressureRow;
extern FType Max_Dt;


FType ReadRST(void);
void WriteRST(FType);
void fSuccFalsh(int, int, int);
FType ShapeFactor(int, int, int);
FType ShapeFactor_1(int, int, int);
FType ShapeFactor_2(int, int, int);
FType ShapeFactor_3(int, int, int);
void fCalcSatFs(void);
void fCalc_Visco(void);
void fCalc_Transes(void);
void fCalc_Phase_P(void);
void fPhaseCtrl(void);
void fCalcPhaseGrav(void);
void fCalcIFT(void);
void fPre_Mul(void);
void fCBlock_Jac(int, int, int, double);
void fCalcDiffusion(void);
FType fDifDiffusion(int, int, int, int, int , int);
FType fDiffVisco(int, int, int, int, int);
void Equilibrium_Concentrations(int, int, int);
FType Calc_Re(int, int, int, int);
void MKLS2(FType *, FType *, int);
FType UpdateAllForMKL(FType *, bool, FType);
FType fCalcDiffIntegral(int, int, int, int, int, char);
FType CalcDiffIntegral(int, int, int, int, int, char);
int RotatingOne(int);
				
				//Parsol.cpp*******************************************
//extern void ParaSolver(FType *, int *, int *, FType *, FType *);

//pre.cpp****************************************************
//extern FType In_Product(FType *, FType *, int);
//extern void Scaling1(FType *, int *, int *, FType *);
//extern void PrebiCGStab(FType *, int *, int *, FType *);
//Data_Input.cpp***************************************************
void Allocation(void);
void ECLStar(char *, int *, FType *);
// void RestoreRST(FType *);
void Data_Input(MIfstream &);

//Parsol.hpp********************************************************
//void ParaSolver(FType *, int *, int *, FType *, FType *);

//MKL_solver.h*********************************************************
//extern void dgesv(const int*, const int*, FType*, const int*, int*, FType*, const int*, int*);
void MKLS(FType **, FType *, FType *, int);
void MKLS1(FType *, int *, int *, FType *);

/*//BiCGstab_l.h*********************************************************
void L_Mat_Mul(FType *, int *, int *, FType **, FType **, int, int, int);
void L_Mat_MulN(FType *, int *, int *, FType **, FType *, int, int);
void L_Mat_MulP(FType *, int *, int *, FType *, FType **, int, int);
void L_Mat_MulF(FType **, FType **, FType *, int, int);
void BiCGStabL(void);*/

//Do_Cycle.h************************************************************
void CalcSatFs(void);
void Calc_Visco(void);
void fCalc_Visco(void);
void Calc_Transes(void);
void fCalc_Transes(void);
void Calc_Phase_P(void);
void fCalc_Phase_P(void);
void System_Reduce(void);
void PhaseCtrl(void);
void fPhaseCtrl(void);
FType System_Update(void);
void CalcPhaseGrav(void);
void fCalcPhaseGrav(void);
//extern void DelJacRow(int, int);*
//void DelJacCol(int, int);*
void BuildPrecon(void);
double oilRateCh(int, int, int, double);
double gasRateCh(int, int, int, double);
double oilRateCh2(int, int, int, double, double);
void JacGuass(void);
void CalcIFT(void);
void fCalcIFT(void);
void CBlock_Jac(int, int, int, double);
void fCBlock_Jac(int, int, int, double);

//Do_Time_Step*******************************************************
void CreateBackups(FType);
void ZeroWellQ(void);
FType MatBal(FType);
void RestoreBackups(void);
char TimeStepCtrl(double *);
FType AvgPReport(void);
FType bAvgPReport(void);
void CreateBackups15(FType);
void Pre_Mul(void);
void fPre_Mul(void);
//Gauss_Jordan.h******************************************************
//extern void GaussJ(FType **, FType *, FType *, int);

//Globals.h*********************************************************
void TerM(char *);
FType dSgn(FType);
//extern void CLFinish(void);
//void ParaSolver(FType *, int *, int *, FType *, FType *)
FType Solve_Z(FType, FType, FType, char);

/*//GMRES.h:********************************************************
void Mat_Mul(FType *, int *, int*, FType *, FType *, int);
FType In_Product(FType *, FType *, int);
void Mat_MulGMR(FType **, FType *, FType *, int);
void Mat_MulGM2(FType **, FType *, FType *, int, int);
void Scaling(FType *, int *, int *, FType *, int);
void Mat_MulCol(FType *, int *, int *, FType **, FType *, int, int);
FType In_ProductCol(FType *, FType **, int, int);
void Scaling1(FType *, int *, int *, FType *);
void PreGMRESm(FType *, int *, int *, FType *, int );
void Mat_MulM(FType **, FType *, FType *, int );
void GMRESm(FType *, int *, int *, FType *, int );*/

//Init_calcs.h*****************************************************
void Tr_Mu(void);
void InitWells(void);
void EOS_Init(void);
//void Jac_Init(void);*
void Edge_Trans(void);
void Edge_fTrans(void);
void CPlus_Props(void);
void WaterSat(void);
void SuccFalsh(int, int, int);
void fSuccFalsh(int, int, int);
//extern void GCE(void);
//extern FType SuccFalshFug(int Ix, int Iy, int Iz, FType *Fug);
//extern char SuccFalshFugDout(int Ix, int Iy, int Iz, FType PressIn, FType *Fug, FType *dFugdP);
FType Pcgo(int, int, int);
FType Pcow(int, int, int);
void AllFlash(void);
void QInitValue(void);
void CalcBlockHeight(void);
//extern void GCE_NR(void);
//extern FType GCE_NR_eqns(FType *, FType *, FType, char);
//extern void GOC(FType *, FType *, FType, char);
void PJac(void);
void CalcPCoeff(void);
void ManageTSMarker(void);
void IFK_Calc(void);

//Mass_Tran.h*********************************************
FType DifDiffusion(int, int, int, int, int, int);
FType DiffVisco(int, int, int, int, int);
void CalcDiffusion(void);
FType fDifDiffusion(int, int, int, int, int, int);
FType fDiffVisco(int, int, int, int, int);
void fCalcDiffusion(void);

//MIfstream.h***********************************************
//Minor_NR.h**********************************************
FType Dew_NR(int, int, int);
FType Bubl_Succ(int, int, int);
FType Dew_Succ(int, int, int);
FType fBubl_Succ(int, int, int);
FType fDew_Succ(int, int, int);

FType LowerDew_Succ(int, int, int);
void critTest(int, int, int);
FType Bubl_NR(int, int, int);

//Numeric.h**********************************************
FType machine_eps(void);
FType machine_epsv(FType);
FType MachineEpsilon(void);

/*//Orthomin.h**********************************************
void Orthomin(FType **, FType *, FType *, int );
//void Mat_Mul(FType **Ls, FType *Rs, FType *As, int MSize);
void Mat_Mul2(FType **, FType **, FType *, int , int );
//FType In_Product(FType *Ls, FType *Rs, int MSize);
void GuassSiedel(FType **, FType *, FType *, int );
void Orthomin(FType **, FType *, FType *, int );*/

//Pre.h***************************************************
void NumDrop(FType *, int *, int);
void Mat_MulT(FType *, int *, int *, FType *, FType *, int);
void MRPre(void);

/*//serial_Bicgstab.h*******************************************
void PrebiCGStab(FType *, int *, int *, FType *);
void Mat_MulN(FType *, int *, int *, FType *, FType *, int);
void Scaling1(FType *, int *, int *, FType *);
void Mat_Mul(FType *, int *, int*, FType *, FType *, int);
FType In_Product(FType *, FType *, int);
void Mat_MulF(FType **, FType *, FType *, int);
void biCGStab(FType *, int *, int *, FType *);*/

/*//Temps.h********************************************************
int InitADL(void);
void GetGPUTemp(void);
void CPUInit(void);
void GetCPUTemp(void);
typedef int(*ADL_MAIN_CONTROL_CREATE)(ADL_MAIN_MALLOC_CALLBACK, int);
typedef int(*ADL_MAIN_CONTROL_DESTROY)();
typedef int(*ADL_CONSOLEMODE_FILEDESCRIPTOR_SET) (int);
typedef int(*ADL_OVERDRIVE5_TEMPERATURE_GET) (int, int, ADLTemperature*);*/

//viewAssist.h***************************************************
void PrintAll(std::ofstream &, FType);
void PrintAll2(int, int, int, std::ofstream &fp, FType);
void Print2D(FType **, FType *, int, std::ofstream &);
void Scaling(FType *, int *, int *, FType *);
void FlashInjection(FType *, FType, FType *, FType *, FType &, FType &, FType &, FType &, FType &, FType &);

extern FType WellCompSadeq[3];
extern FType WellCompSadeq_cum[3];
extern FType WellCompSadeq2[3];




#pragma once
