#ifndef MKL_SOLVER_H
#define MKL_SOLVER_H

#ifdef _WIN32
#include <mkl.h>
#include <stdio.h>
#include <string.h>
#include "Functions.h"
extern void dgesv(const int* n, const int* nrhs, FType* a, const int* lda, int* ipiv, FType* b, const int* ldb, int* info);
#define __DGESV dgesv(&n, &nrhs, AT, &lda, ipiv, b, &ldb, &info) 
#else
#include <lapacke.h>
#define __DGESV info=LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, AT, lda, ipiv, b, ldb)
#endif

#include <stdlib.h>
//#include <omp.h>


extern char repStat;

//extern void dgesv(int* n, int* nrhs, FType* a, int* lda, int* ipiv, FType* b, int* ldb, int* info);
void MKLS(FType **A, FType *b, FType *x, int MSize);

void MKLS1(FType *A, int *row, int *col, FType *b) {
	FType *AT;
	int *ipiv, MSize;
	int n, nrhs, lda, ldb, info;
	register int i, j;

	MSize = row[0];
	row[0] = 0;

	if ((AT = (FType *)malloc(MSize*MSize * sizeof(FType))) == NULL) TerM("Can not allocate memory for AT matrix");
	if ((ipiv = (int *)malloc(MSize * sizeof(int))) == NULL) TerM("Can not allocate memory for ipiv matrix");

	n = MSize;
	lda = MSize;
	ldb = MSize;
	nrhs = 1;

	for (i = 0; i<(MSize*MSize); i++) AT[i] = 0;

	for (i = 0; i<MSize; i++) {
		for (j = row[i]; j<row[i + 1]; j++) AT[MSize*col[j] + i] += A[j];
	}

	//mkl_set_dynamic(0);
	//mkl_set_num_threads(4);


	__DGESV;

	if (info) repStat = -1;
	for (i = 0; i<totalRow; i++) {
		if (ans[i] != ans[i]) {
			repStat = -1;
			break;
		}
	}

	free(AT);
	free(ipiv);
}

void MKLS(FType **A, FType *b, FType *x, int MSize) {
	FType *AT;
	int *ipiv;
	int n, nrhs, lda, ldb, info;
	register int i, j;

	if ((AT = (FType *)malloc(MSize*MSize * sizeof(FType))) == NULL) TerM("Can not allocate memory for AT matrix");
	if ((ipiv = (int *)malloc(MSize * sizeof(int))) == NULL) TerM("Can not allocate memory for ipiv matrix");

	n = MSize;
	lda = MSize;
	ldb = MSize;
	nrhs = 1;

	for (i = 0; i<MSize; i++)
		for (j = 0; j<MSize; j++) AT[MSize*i + j] = A[j][i];

	__DGESV;
	//info=LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, AT, lda, ipiv, b, ldb);

	for (i = 0; i<MSize; i++) x[i] = b[i];

	free(AT);
	free(ipiv);
}

void MKLS2(FType *AT, FType *b, int MSize) {
	int *ipiv;
	int n, nrhs, lda, ldb, info;
	register int i;

	
	if ((ipiv = (int *)malloc(MSize * sizeof(int))) == NULL) TerM("Can not allocate memory for ipiv matrix");
	
	
	n = MSize;
	lda = MSize;
	ldb = MSize;
	nrhs = 1;


	__DGESV;

	//repStat = 0;
	if (info) repStat = -1;
	for (i = 0; i < MSize; i++) {
		if (b[i] != b[i]) {
			repStat = -1;
			break;
		}
	}

	free(ipiv);	
}

#endif