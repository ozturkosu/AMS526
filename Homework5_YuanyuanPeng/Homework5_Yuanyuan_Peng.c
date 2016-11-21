#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// LAPACK function declarations
void dgesvd_(char *JOBU, char *JOBVT, int *M, int *N, double *A, int *LDA,
             double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK,
			 int *LWORK, int *INFO);
void dgemv_(char *TRANS, int *M, int *N, double *ALPHA, double *A, int *LDA,
            double *X, int *INCX, double *BETA, double *Y, int *INCY);
void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA,
			double *A,int *LDA, double *B, int *LDB, double *BETA, double *C,
			int *LDC);
void dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W,
			double *WORK, int *LWORK, int *INFO);

int mtx_get_idx(int r, int c, int m)
{
	return c*m + r;
}

void print_array(double *array, int n)
{
	printf("[");
	for (int i = 0; i < n; i++)
	{
		printf("%g, ",array[i]);
	}
	printf("]\n");
}

void compute_rsvd(double *A, int mA, int nA, double *U, double *S, double *VT)
{
	char JOBU = 'S';
	char JOBVT = 'S';
	int LDA = mA;
	int LDU = mA;
	int LDVT = nA;
	double dWORK_opt;
	int LWORK = -1;
	int INFO;
	
	// Call the function to query the WORK array.
	dgesvd_(&JOBU, &JOBVT, &mA, &nA, A, &LDA, S, U, &LDU, VT, &LDVT,
	        &dWORK_opt, &LWORK, &INFO);
	if (INFO) exit(1);
	int LWORK_opt = (int)dWORK_opt;
	double *WORK = malloc(sizeof(double)*LWORK_opt);
	dgesvd_(&JOBU, &JOBVT, &mA, &nA, A, &LDA, S, U, &LDU, VT, &LDVT,
	        WORK, &LWORK_opt, &INFO);
	free(WORK);
	if (INFO) exit(1);
}

void svd_compute_Ub(double *U, int mU, int nU, double *b, double *Ub)
{
	char TRANS = 'T';
	double ALPHA = 1.0;
	int LDA = mU;
	int INCX = 1;
	double BETA = 0.0;
	int INCY = 1;
	
	// Make the call.
	dgemv_(&TRANS, &mU, &nU, &ALPHA, U, &LDA, b, &INCX, &BETA, Ub, &INCY);
}

void compute_Vb(double *V, int mV, int nV, double *b, double *Vb)
{
	char TRANS = 'N';
	double ALPHA = 1.0;
	int LDA = mV;
	int INCX = 1;
	double BETA = 0.0;
	int INCY = 1;
	
	// Make the call.
	dgemv_(&TRANS, &mV, &nV, &ALPHA, V, &LDA, b, &INCX, &BETA, Vb, &INCY);
}

void solve_diagonal_system(double *A_diag, double nA, double *b, double *x)
{
	for (int i = 0; i < nA; i++)
	{
		x[i] = b[i] / A_diag[i];
	}
}

void solve_least_squares_svd(double *A, int mA, int nA, double *b, double *x)
{
	// Allocate some memory
	double *U = malloc(sizeof(double)*mA*nA);
	double *VT = malloc(sizeof(double)*nA*nA);
	double *S = malloc(sizeof(double)*nA);
	double *Ub = malloc(sizeof(double)*nA);
	double *w = malloc(sizeof(double)*nA);
	
	// Compute the reduced SVD of A. Assume mA > nA.
	compute_rsvd(A, mA, nA, U, S, VT);
	
	// Compute (U*)b.
	svd_compute_Ub(U, mA, nA, b, Ub);
	
	// Solve diagonal system Sigma*w = (U*)b for w.
	solve_diagonal_system(S, nA, Ub, w);
	
	// Compute x = Vw
	svd_compute_Ub(VT, nA, nA, w, x);
	
	// Memory cleanup
	free(U);
	free(VT);
	free(S);
	free(Ub);
	free(w);
}

void compute_ATA(double *A, int mA, int nA, double *ATA)
{
	char TRANSA='T';
	char TRANSB='N';
	int M = nA;
	int N = nA;
	int K = mA;
	double ALPHA = 1.0;
	int LDA = mA;
	int LDB = mA;
	double BETA = 0.0;
	int LDC = nA;
	// Call the routine.
	dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, A, &LDB, &BETA,
		   ATA, &LDC);
}

void compute_eigenvectors(double *A, int mA, double *W)
{
	char JOBZ = 'V';
	char UPLO = 'U';
	int N     = mA;
	int LDA   = mA;
	double WORK_opt;
	int LWORK = -1;
	int INFO = 20;
	
	// Calculate optimal work.
	dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, &WORK_opt, &LWORK, &INFO);
	if (INFO) exit(1);
	int LWORK_opt = (int)WORK_opt;
	double *WORK = malloc(sizeof(double)*LWORK_opt);
	
	// Call the routine.
	INFO = 20;
	dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK_opt, &INFO);
	if (INFO) exit(1);
	
	// Free memory.
	free(WORK);
}

void solve_least_squares_eigen(double *A, double mA, double nA, double *b, double *x)
{
	double *ATA = malloc(sizeof(double)*nA*nA);
	
	// Compute ATA.
	compute_ATA(A, mA, nA, ATA);
	
	// Compute eigenvalue decomp of ATA.
	double *W = malloc(sizeof(double)*nA);
	compute_eigenvectors(ATA, nA, W);
	
	// ATA is now Q.
	double *ATb = malloc(sizeof(double)*nA);
	svd_compute_Ub(A, mA, nA, b, ATb);
	double *QTATb = malloc(sizeof(double)*nA);
	svd_compute_Ub(ATA, nA, nA, ATb, QTATb);
	
	// Solve a diagonal system.
	double *w = malloc(sizeof(double)*nA);
	solve_diagonal_system(W, nA, QTATb, w);
	
	// Compute x
	compute_Vb(ATA, nA, nA, w, x);
	
	free(ATA);
	free(W);
	free(ATb);
	free(QTATb);
	free(w);
}

void generate_least_squares_problem(int m, int n, double *A, double *b)
{
	for (int i = 0; i < m; i++)
	{
		double t_i = ((double)(i))/((double)(m-1));
		double y_i = 0.0;
		for (int j = 0; j < n; j++)
		{
			double pow_ti_j = pow(t_i,(double)j); 
			y_i += pow_ti_j;
			// Build the matrix A. A is m times n.
			A[mtx_get_idx(i,j,m)] = pow_ti_j;
		}
		// Build the m-vector b.
		b[i] = y_i;
	}
}

int main()
{
	FILE *outfile = fopen("results.csv", "w");
	fprintf(outfile, "n,svd error,eigen error\n");
	for (int n = 3; n <= 15; n++)
	{
		printf("Computing n=%d\n", n);
		int m = 2*n;
		double *A = malloc(sizeof(double)*m*n);
		double *A_2 = malloc(sizeof(double)*m*n);
		double *b = malloc(sizeof(double)*m);
		double *b_2 = malloc(sizeof(double)*m);
		double *x_svd = malloc(sizeof(double)*n);
		double *x_eigen = malloc(sizeof(double)*n);
		
		generate_least_squares_problem(m,n,A,b);
		memcpy(A_2, A, sizeof(double)*m*n);
		memcpy(b_2, b, sizeof(double)*m);
		
		solve_least_squares_svd  (A, m, n, b, x_svd);
		solve_least_squares_eigen(A_2, m, n, b_2, x_eigen);
		
		// Compute the 2-norm error
		double error_2n_svd = 0.0, error_2n_eigen = 0.0;
		for (int i = 0; i < n; i++)
		{
			error_2n_svd += pow(x_svd[i] - 1.0, 2);
			error_2n_eigen += pow(x_eigen[i] - 1.0, 2);
		}
		fprintf(outfile, "%d,%e,%e\n", n, error_2n_svd, error_2n_eigen);
		
		free(A);
		free(A_2);
		free(b);
		free(b_2);
		free(x_svd);
		free(x_eigen);
	}
	fclose(outfile);
	return 0;
}