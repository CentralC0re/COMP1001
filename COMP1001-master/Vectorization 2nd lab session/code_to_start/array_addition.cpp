#include "array_addition.h"

float  X1[M2], X2[M2], Y1[M2], test2[M2];

void initialization_Add() {

	float e = 0.1234f, p = 0.7264f, r = 0.11f;

	for (unsigned int j = 0; j != M2; j++) {
		Y1[j] = 0.0;
		test2[j] = 0.0;
		X1[j] = (j % 7) + r;
		X2[j] = (j % 13) + e;
	}
}


unsigned short int Add_default() {

	for (int j = 0; j < M2; j++) {
		Y1[j] = X1[j] + X2[j];
	}


	return 0;
}

unsigned short int Add_SSE() {

	//__m128 num0, num1, num2, num3, num4, num5, num6;

	__m128 container1, container2, container3;	// Forgot you could create multiple variables on one line
	int counter = 0;
	for (int i = 0; i < ((M2 / 4) * 4); i += 4) {
		container1 = _mm_loadu_ps(&X1[i]);
		container2 = _mm_loadu_ps(&X2[i]);
		container3 = _mm_add_ps(container1, container2);
		_mm_storeu_ps(&Y1[i], container3);
		counter += 4;
	}

	for (int i = counter; i < M2; i++)	// No need for if, if counter >= M2 for loop doesn't run
	{
		Y1[i] = X1[i] + X2[i];
	}
	

	return 0;
}


unsigned short int Add_AVX() {

	//__m256  ymm0, ymm1, ymm2, ymm3, ymm4;

	__m256 container1;
	__m256 container2;
	__m256 container3;
	int counter = 0;
	for (int i = 0; i < ((M2 / 8) * 8); i += 8) {		// Functionally identical, but processes 8 at a time rather than 4
		container1 = _mm256_loadu_ps(&X1[i]);
		container2 = _mm256_loadu_ps(&X2[i]);
		container3 = _mm256_add_ps(container1, container2);
		_mm256_storeu_ps(&Y1[i], container3);
		counter += 8;
	}

	for (int i = counter; i < M2; i++)	// No need for if, if counter >= M2 for loop doesn't run
	{
		Y1[i] = X1[i] + X2[i];
	}


	return 0;
}


unsigned short int Compare_Add() {


	for (int j = 0; j < M2; j++) {
		test2[j] = X1[j] + X2[j];
	}


	for (int j = 0; j < M2; j++)
		if (equal(Y1[j], test2[j]) == 1) {
			//printf("\n j=%d\n", j);
			return 1;
		}

	return 0;
}



