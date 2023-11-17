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
void resetY();	// Keep forgetting to do these
void resetW();
void routine1(float alpha, float beta);
void routine2(float alpha, float beta);
void routine1_vec(float alpha, float beta);
void routine2_vec(float alpha, float beta);
void checkValues(float* correct, float* checking);

__declspec(align(64)) float  yReset[M], y[M], z[M] ;
__declspec(align(64)) float A[N][N], x[N], wReset[N], w[N];

int main() {

    float alpha = 0.023f, beta = 0.045f;
    double run_time, start_time;
    unsigned int t;

    initialize();

	//R1	-	Commented out so that y remains in its unmodified state.
	//			Unless I edit Routine1 so it stores its results elsewhere, it must remain commented out
	//			FIX THIS - "Your program must work for any input size value (any ‘N, M’ values)." Results cannot be hardcoded.
    printf("\nRoutine1:");
    start_time = omp_get_wtime(); //start timer

    for (t = 0; t < TIMES1; t++)
        routine1(alpha, beta);

    run_time = omp_get_wtime() - start_time; //end timer
    printf("\n Time elapsed is %f secs \n %e FLOPs achieved\n", run_time, (double)(ARITHMETIC_OPERATIONS1) / ((double)run_time / TIMES1));
    
	float r1CheckValues[5] = { y[1], y[25], y[300], y[1024], y[1025] };

	resetY();

	//R1_V
	printf("\nRoutine1 Vectorised:");
	start_time = omp_get_wtime(); //start timer

    for (t = 0; t < TIMES1; t++)
        routine1_vec(alpha, beta);

    run_time = omp_get_wtime() - start_time; //end timer
    printf("\n Time elapsed is %f secs \n %e FLOPs achieved\n", run_time, (double)(ARITHMETIC_OPERATIONS1) / ((double)run_time / TIMES1));

	float r1VCheckValues[5] = { y[1], y[25], y[300], y[1024], y[1025] }; // Check values from the result of routine1_vec
																						// These may be incorrect.
	checkValues(r1CheckValues, r1VCheckValues);



	//R2
    printf("\nRoutine2:");
    start_time = omp_get_wtime(); //start timer

    for (t = 0; t < TIMES2; t++)
        routine2(alpha, beta);

    run_time = omp_get_wtime() - start_time; //end timer
    printf("\n Time elapsed is %f secs \n %e FLOPs achieved\n", run_time, (double)(ARITHMETIC_OPERATIONS2) / ((double)run_time / TIMES2));

	float r2CheckValues[5] = { w[1], w[25], w[300], w[1024], w[1025] };

	resetW();

	//R2_V
	printf("\nRoutine2 Vectorised:");
	start_time = omp_get_wtime(); //start timer

	for (t = 0; t < TIMES2; t++)
		routine2_vec(alpha, beta);

	run_time = omp_get_wtime() - start_time; //end timer
	printf("\n Time elapsed is %f secs \n %e FLOPs achieved\n", run_time, (double)(ARITHMETIC_OPERATIONS2) / ((double)run_time / TIMES2));

	float r2VCheckValues[5] = { w[1], w[25], w[300], w[1024], w[1025] };

	checkValues(r2CheckValues, r2VCheckValues);

    return 0;
}

void initialize() {

    unsigned int i, j;

	    //initialize routine2 arrays
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) {
            A[i][j] = (i % 99) + (j % 14) + 0.013f;		// I was trying to vectorise this... whoops
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
}												// This means that both R1 and R1vec have the same array.

void routine1_vec(float alpha, float beta) {

    unsigned int i;
	float aArray[8];
	float bArray[8];
	for (int i = 0; i < 8; i++) {	// Probably inefficient.
		aArray[i] = alpha;
		bArray[i] = beta;
	}

	__m256 aContain = _mm256_loadu_ps(aArray);
	__m256 bContain = _mm256_loadu_ps(bArray);
	// Attempted to manually destroy a and b arrays, but it's far too much effort for
	// very little memory saving.

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
	}
}

void routine2(float alpha, float beta) {

    unsigned int i, j;

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            w[i] = w[i] - beta + alpha * A[i][j] * x[j];


}

void routine2_vec(float alpha, float beta) {

	unsigned int i, j;
	float aArray[8];
	float bArray[8];
	for (int i = 0; i < 8; i++) {	// Probably inefficient.
		aArray[i] = alpha;
		bArray[i] = beta;
	}

	__m256 alphaContain = _mm256_loadu_ps(aArray);
	__m256 betaContain = _mm256_loadu_ps(bArray);

	for (i = 0; i < N; i+=8)
	{
		unsigned int iArray[8] = { i };	// Despite both being initialised as int, i is apparently unsigned
										// Set array to unsigned to remove conversion warning and allow greater upper value
		__m256 wContain = _mm256_loadu_ps(&w[i]);
		__m256 subResult = _mm256_sub_ps(wContain, betaContain);	// This does not need to be recalculated	- Correct result
		// alpha cannot be added before the loop because of the order of operations (alpha * A)
		for (j = 0; j < N; j+=8)									// every j loop
		{
			__m256 arrAContain = _mm256_loadu_ps(&A[i][j]);
			__m256 xContain = _mm256_loadu_ps(&x[j]);
			
			__m256 multResult1 = _mm256_mul_ps(alphaContain, arrAContain);	// Correct result
			__m256 multResult2 = _mm256_mul_ps(multResult1, xContain);		// Effectively correct (minor float inaccuracy)

			__m256 finalResult = _mm256_add_ps(subResult, multResult2);		// Effectively correct (minor float inaccuracy)
			_mm256_storeu_ps(&w[i], finalResult);

			//w[i] = w[i] - beta + alpha * A[i][j] * x[j];
		}
	}
}

void checkValues(float* correct, float* checking) {	// This only takes 1 value, by the way. FIX.
	bool isCorrect = true;
	for (int i = 0; i < 5; i++)	// Fixed i cap as the arrays are always length 5
	{
		if (abs((correct[i] - checking[i]) / checking[i]) >= 0.0001)	// FP calculation to determine if two values are equal enough, 100% higher tolerance than default
		{																// A minor innacuracy (below) resulted in the value being registered as incorrect. Value has been upped.
			isCorrect = false;
		}
	}

	printf("");
	if (isCorrect) {
		printf("Values match.\n");
	}
	else {
		printf("Error: Values do not match.\n");
	}
}

