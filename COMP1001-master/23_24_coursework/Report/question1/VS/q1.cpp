/*
------------------DR VASILIOS KELEFOURAS-----------------------------------------------------
------------------COMP1001 ------------------------------------------------------------------
------------------COMPUTER SYSTEMS MODULE-------------------------------------------------
------------------UNIVERSITY OF PLYMOUTH, SCHOOL OF ENGINEERING, COMPUTING AND MATHEMATICS---
*/


#include <stdio.h>
#include <time.h>
#include <pmmintrin.h>
#include <process.h>
#include <chrono>
#include <iostream>
#include <immintrin.h>
#include <omp.h>

#define M 1024*512
#define ARITHMETIC_OPERATIONS1 3*M
#define TIMES1 1

#define N 8192
#define ARITHMETIC_OPERATIONS2 4*N*N
#define TIMES2 1


//function declaration
void initialize();
void resetY();
void resetW();
void routine1(float alpha, float beta);
void routine2(float alpha, float beta);
void routine1_vec(float alpha, float beta);
void routine2_vec(float alpha, float beta);
void checkValues(float correct[5], float checking[5]);

__declspec(align(64)) float  yReset[M], y[M], z[M] ;
__declspec(align(64)) float A[N][N], x[N], wReset[N], w[N];

int main() {

    float alpha = 0.023f, beta = 0.045f;
    double run_time, start_time;
    unsigned int t;

    initialize();

	//R1
    printf("\nRoutine1:");
    start_time = omp_get_wtime(); //start timer

    for (t = 0; t < TIMES1; t++)
        routine1(alpha, beta);

    run_time = omp_get_wtime() - start_time; //end timer
    printf("\n Time elapsed is %f secs \n %e FLOPs achieved\n", run_time, (double)(ARITHMETIC_OPERATIONS1) / ((double)run_time / TIMES1));
    
	int MLoc = M/5;
	float r1CheckValues[5] = {0, 0, 0, 0, 0};	// Was not seeming to be set to 0, manually set.
	if (MLoc == 0)	// Division does not round up. If M < 5, MLoc = 0
	{
		r1CheckValues[0] = y[MLoc];
	}
	else if (MLoc < 0)
	{
		return -1;	// Invalid M
	}
	else {		// Space inefficient, but good enough.
		 r1CheckValues[0] = y[MLoc];
		 r1CheckValues[1] = y[MLoc * 2];
		 r1CheckValues[2] = y[MLoc * 3];
		 r1CheckValues[3] = y[MLoc * 4];
		 r1CheckValues[4] = y[MLoc * 5];
	}
	 

	resetY();

	//R1_V
	printf("\nRoutine1 Vectorised:");
	start_time = omp_get_wtime(); //start timer

    for (t = 0; t < TIMES1; t++)
        routine1_vec(alpha, beta);

    run_time = omp_get_wtime() - start_time; //end timer
    printf("\n Time elapsed is %f secs \n %e FLOPs achieved\n", run_time, (double)(ARITHMETIC_OPERATIONS1) / ((double)run_time / TIMES1));

	float r1VCheckValues[5] = { 0, 0, 0, 0, 0 }; // Check values from the result of routine1_vec
												// These may be incorrect.

	if (MLoc == 0)	// Could make this a function, but I have no idea if this is allowed.
	{
		r1VCheckValues[0] = y[MLoc];
	}
	else if (MLoc < 0)
	{
		return -1;	// Invalid M
	}
	else {
		r1VCheckValues[0] = y[MLoc];
		r1VCheckValues[1] = y[MLoc * 2];
		r1VCheckValues[2] = y[MLoc * 3];
		r1VCheckValues[3] = y[MLoc * 4];
		r1VCheckValues[4] = y[MLoc * 5];
	}

	checkValues(r1CheckValues, r1VCheckValues);

	//R2
    printf("\nRoutine2:");
    start_time = omp_get_wtime(); //start timer

    for (t = 0; t < TIMES2; t++)
        routine2(alpha, beta);

    run_time = omp_get_wtime() - start_time; //end timer
    printf("\n Time elapsed is %f secs \n %e FLOPs achieved\n", run_time, (double)(ARITHMETIC_OPERATIONS2) / ((double)run_time / TIMES2));

	int NLoc = N / 5;
	float r2CheckValues[5] = { 0, 0, 0, 0, 0 };
	if (NLoc == 0)
	{
		r2CheckValues[0] = w[NLoc];
	}
	else if (NLoc < 0)
	{
		return -1;	// Invalid N
	}
	else {
		r2CheckValues[0] = w[NLoc];
		r2CheckValues[1] = w[NLoc * 2];
		r2CheckValues[2] = w[NLoc * 3];
		r2CheckValues[3] = w[NLoc * 4];
		r2CheckValues[4] = w[NLoc * 5];
	}

	resetW();

	//R2_V
	printf("\nRoutine2 Vectorised:");
	start_time = omp_get_wtime(); //start timer

	for (t = 0; t < TIMES2; t++)
		routine2_vec(alpha, beta);

	run_time = omp_get_wtime() - start_time; //end timer
	printf("\n Time elapsed is %f secs \n %e FLOPs achieved\n", run_time, (double)(ARITHMETIC_OPERATIONS2) / ((double)run_time / TIMES2));

	float r2VCheckValues[5] = { 0, 0, 0, 0, 0 };
	if (NLoc == 0)
	{
		r2VCheckValues[0] = w[NLoc];
	}
	else if (NLoc < 0)
	{
		return -1;	// Invalid N
	}
	else {
		r2VCheckValues[0] = w[NLoc];
		r2VCheckValues[1] = w[NLoc * 2];
		r2VCheckValues[2] = w[NLoc * 3];
		r2VCheckValues[3] = w[NLoc * 4];
		r2VCheckValues[4] = w[NLoc * 5];
	}

	checkValues(r2CheckValues, r2VCheckValues);		// r2Check and r2VCheck are independent.
										
    return 0;
}

void initialize() {

    unsigned int i, j;

	    //initialize routine2 arrays
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) {
			A[i][j] = (i % 99) + (j % 14) + 0.013f;
        }
	
    //initialize routine1 arrays
    for (i = 0; i < N; i++) {
        x[i] = (i % 19) - 0.01f;
        wReset[i] = (i % 5) - 0.002f;
		w[i] = wReset[i];
    }

    //initialize routine1 arrays
    for (i = 0; i < M; i++) {
        z[i] = (i % 9) - 0.08f;
        yReset[i] = (i % 19) + 0.07f;
		y[i] = yReset[i];				// Could just use std::copy, but we haven't covered it
    }
}

void resetY() {	// More efficient than initialize and safer than manual
	for (int i = 0; i < M; i++) {
		y[i] = yReset[i];
	}
}

void resetW() {
	for (int i = 0; i < N; i++) {
		w[i] = wReset[i];
	}
	
}


void routine1(float alpha, float beta) {

	unsigned int i;


	for (i = 0; i < M; i++)
		y[i] = alpha * y[i] + beta * z[i];	// Rather than editing y itself, yRes is filled.
}											// This means that both R1 and R1vec have the same array.

void routine1_vec(float alpha, float beta) {

	if (M < 8) {				// No point running this function if M is smaller than vectorisation size.
		routine1_vec(alpha, beta);
	}
    unsigned int i;

	__m256 aContain = _mm256_set1_ps(alpha);	// Sets all 8 values to alpha, rather than needing a for loop
	__m256 bContain = _mm256_set1_ps(beta);		// More performant

	__m256 finalResult;

	for (i = 0; i < M; i += 8)
	{	// Curly brackets were missing and have been added.
		__m256 aResult, bResult;		// Variables defined here so they go out of scope on loop end.
		__m256 yContain = _mm256_loadu_ps(&y[i]);
		__m256 zContain = _mm256_loadu_ps(&z[i]);
		aResult = _mm256_mul_ps(aContain, yContain);
		bResult = _mm256_mul_ps(bContain, zContain);
		finalResult = _mm256_add_ps(aResult, bResult);
		_mm256_storeu_ps(&y[i], finalResult);		// Given sets of 8 values are consistently updating, it seems to work

		//y[i] = alpha * y[i] + beta * z[i];

		if (i + 8 >= M && i != M-1)	// If next iteration puts i above M, and i is not already at its final value
		{
			i+=8;
			for (i; i < M; i++)
			{
				y[i] = alpha * y[i] + beta * z[i];
			}
		}
	}
}

void routine2(float alpha, float beta) {

    unsigned int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			w[i] = w[i] - beta + alpha * A[i][j] * x[j];
		}
			
	}
}

void routine2_vec(float alpha, float beta) {
	
	if (N < 8)
	{
		routine2(alpha, beta);
	}

	unsigned int i, j, k;
	float noVecArr[8];
	float multResultArr[8];

	__m256 alphaContain = _mm256_set1_ps(alpha);	// More performant than a for loop.

	__m256 arrAContain;		// Prevents constant redefinition.
	__m256 xContain;
	__m256 multResult;
	__m256 multResult2;

	float toStore = 0;			// Stores results if N is not a multiple of 4.

	for (i = 0; i < N; i++)
	{						// Either i or j can be vectorised (+=4), but not both
							// (Without further logic)		
		for (j = 0; j < N; j+=8)
		{
			arrAContain = _mm256_loadu_ps(&A[i][j]);
			xContain = _mm256_loadu_ps(&x[j]);

			//subResult = _mm_sub_ps(wContain, betaContain);		// This needs recalculation (w updated every j)
			multResult = _mm256_mul_ps(alphaContain, arrAContain);
			multResult2 = _mm256_mul_ps(multResult, xContain);

			_mm256_storeu_ps(&multResultArr[0], multResult2);	// The following section is not vectorised, stored back.
			for (k = 0; k < 8; k++)
			{
				w[i] = w[i] - beta + multResultArr[k];	
			}			
			//w[i] = w[i] - beta + alpha * A[i][j] * x[j];


		}
	}
	
}

void checkValues(float correct[5], float checking[5]) {
	bool isCorrect = true;
	for (int i = 0; i < 5; i++)	// Fixed i cap as the arrays are always length 5
	{
		if (abs((correct[i] - checking[i]) / checking[i]) >= 0.00001)	// FP calculation to determine if two values are equal enough
		{
			isCorrect = false;
			break;
		}
	}

	if (isCorrect) {
		printf("Values match.\n");
	}
	else {
		printf("Error: Values do not match.\n");
	}
}

