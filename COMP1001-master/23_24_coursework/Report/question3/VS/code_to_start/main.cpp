/*
------------------DR VASILIOS KELEFOURAS-----------------------------------------------------
------------------COMP1001 ------------------------------------------------------------------
------------------COMPUTER SYSTEMS MODULE-------------------------------------------------
------------------UNIVERSITY OF PLYMOUTH, SCHOOL OF ENGINEERING, COMPUTING AND MATHEMATICS---
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <emmintrin.h>
#include <limits.h>
#include <pmmintrin.h>
#include <immintrin.h>

//function declarations
void Gaussian_Blur(unsigned char* frame1, unsigned char* filt);
void Sobel(unsigned char* filt, unsigned char* gradient);
void SobelNoVec(unsigned char* filt, unsigned char* gradient);
//int initialize_kernel();	 Pretty sure this is for Linux, ignore.
int read_image(char* filename, unsigned char** frame1, unsigned char** filt, unsigned char** gradient);
void read_image_and_put_zeros_around(char* filename);
void write_image2(char* filename, unsigned char* output_image);
int openfile(char* filename, FILE** finput);
int getint(FILE* fp);
unsigned char* allocateArr();
char* allocatePath(unsigned int size);		// If path arrays stay static, this function becomes irrelevant, REMOVE.
void checkValues(unsigned char* grad1, unsigned char* grad2);

//CRITICAL POINT: images' paths - You need to change these paths
//#define IN "/home/wave/Desktop/comp1001/code_to_start/input_images/a15.pgm"
//#define OUT "/home/wave/Desktop/comp1001/code_to_start/output_images/blurred.pgm"
//#define OUT2 "/home/wave/Desktop/comp1001/code_to_start/output_images/edge_detection.pgm"

//IMAGE DIMENSIONS	-	Variable, should not be hardcoded.
unsigned int M = 0; //cols
unsigned int N = 0; //rows


//CRITICAL POINT:these arrays are defined statically. Consider creating these arrays dynamically instead.
//unsigned char frame1[N * M];//input image
//unsigned char filt[N * M];//output filtered image
//unsigned char gradient[N * M];//output image


const signed char Mask[5][5] = {//2d gaussian mask with integers
	{2,4,5,4,2} ,
	{4,9,12,9,4},
	{5,12,15,12,5},
	{4,9,12,9,4},
	{2,4,5,4,2}
};

const signed char GxMask[3][3] = {
	{-1,0,1} ,
	{-2,0,2},
	{-1,0,1}
};

const signed char GyMask[3][3] = {
	{-1,-2,-1} ,
	{0,0,0},
	{1,2,1}
};

char header[100];


int main(int argc, char* argv[]) {   // No order of arguments provided, assuming IN, OUT.

	if (argc == 3)
	{
		bool endOfFolder = false;
		int i = 0;
		
		unsigned int pathLength = sizeof(argv[1]);
		char iFolderPath[257];
		char oFolderPath[257];
		char oFolderPathConcat[257];

		char iFileName[16];
		char bFileName[22];	// Max int puts this at 21 chars (22 with terminator)
		char eFileName[29];	// Max int puts this at 28 chars (29 with terminator)

		int result;			// Checks if read file worked
		for (i; !endOfFolder; i++)
		{
			strcpy(iFolderPath, argv[1]);	// ifolderPath = argv[1]
			strcpy(oFolderPath, argv[2]);
			strcpy(oFolderPathConcat, oFolderPath);

			// a-.pgm	Filename from 6-16 chars
			sprintf(iFileName, "/a%d.pgm", i);
			sprintf(bFileName, "/blurred%d.pgm", i);
			sprintf(eFileName, "/edge_detection%d.pgm", i);

			
			strcat(iFolderPath, iFileName);
			strcat(oFolderPath, bFileName);
			strcat(oFolderPathConcat, eFileName);

			unsigned char* frame1;
			unsigned char* filt;
			unsigned char* gradient;
			result = read_image(iFolderPath, &frame1, &filt, &gradient);//read image from disc
			if (result == -1)
			{
				if (i == 0)
				{
					return EXIT_FAILURE;
				}
				else 
				{
					printf("All files read\n");
					return EXIT_SUCCESS;
				}
			}

			Gaussian_Blur(frame1, filt); //blur the image (reduce noise)
			
			
			unsigned char* gradientNoVec = gradient;
			Sobel(filt, gradient); //apply edge detection with vectorisation
			// filt is unchanged, but gradient changes
			// As such, it must be copied and rerun.

			SobelNoVec(filt, gradientNoVec);
			checkValues(gradientNoVec, gradient);

			write_image2(oFolderPath, filt); //store output image to the disc
			write_image2(oFolderPathConcat, gradient); //store output image to the disc

			free(frame1);
			free(filt);
			free(gradient);
		}
		return 0;
	}
	else
	{
		printf("Invalid arguments:\n");
		for (int i = 0; i < argc; i++)
			printf("%s\t", argv[i]);
		return -1;
	}
}





void Gaussian_Blur(unsigned char* frame1, unsigned char* filt) {

	int row, col, rowOffset, colOffset;
	int newPixel;
	unsigned char pix;
	//const unsigned short int size=filter_size/2;
	const unsigned short int size = 2;

	/*---------------------- Gaussian Blur ---------------------------------*/
	for (row = 0; row < N; row++) {
		for (col = 0; col < M; col++) {
			newPixel = 0;
			for (rowOffset = -size; rowOffset <= size; rowOffset++) {
				for (colOffset = -size; colOffset <= size; colOffset++) {

					if ((row + rowOffset < 0) || (row + rowOffset >= N) || (col + colOffset < 0) || (col + colOffset >= M))
						pix = 0;
					else
						pix = frame1[M * (row + rowOffset) + col + colOffset];

					newPixel += pix * Mask[size + rowOffset][size + colOffset];

				}
			}
			filt[M * row + col] = (unsigned char)(newPixel / 159);

		}
	}

}


void Sobel(unsigned char* filt, unsigned char* gradient) {

	int row, col;
	//int Gx, Gy;

	//__m128 rowOffset = _mm_set_ps(-1, 0, 1, 0);	// Not in loop because constant and performance reasons
	//__m128 colOffset = _mm_set_ps(-1, 0, 1, 0); // Need to make sure last value has no impact

	int rowOffset, colOffset, i;
	__m128i Gx, Gy, GxMaskCont, GyMaskCont, filtCont, tempStore;
	__m128i MCont = _mm_set1_epi16(M);
	__m128 sqrtVal;
	
	float tempArr[4];
	/*---------------------------- Determine edge directions and gradient strengths -------------------------------------------*/
	if (N - 1 < 8)	// If too small for vectorisation (unlikely)
	{
		SobelNoVec(filt, gradient);
	}
	else
	{
		for (row = 1; row < N - 1; row += 8) {
			for (col = 1; col < M - 1; col += 8) {
				Gx = _mm_set1_epi16(0);	// There does not appear to be a setzero.
				Gy = _mm_set1_epi16(0);

				/* Calculate the sum of the Sobel mask times the nine surrounding pixels in the x and y direction */
				rowOffset = -1;
				for (rowOffset; rowOffset <= 1; rowOffset++) {
					for (colOffset = -1; colOffset <= 1; colOffset++) {
						filtCont = _mm_loadu_epi16(&filt[M * (row + rowOffset) + col + colOffset]);	// Loads filters 0-3
						GxMaskCont = _mm_loadu_epi16(&GxMask[rowOffset + 1][colOffset + 1]);
						GyMaskCont = _mm_loadu_epi16(&GyMask[rowOffset + 1][colOffset + 1]);

						tempStore = _mm_mul_epi32(filtCont, GxMaskCont);
						Gx = _mm_add_epi16(Gx, tempStore);
						tempStore = _mm_mul_epi32(filtCont, GyMaskCont);
						Gy = _mm_add_epi16(Gy, tempStore);
						//Gx += filt[M * (row + rowOffset) + col + colOffset] * GxMask[rowOffset + 1][colOffset + 1];
						//Gy += filt[M * (row + rowOffset) + col + colOffset] * GyMask[rowOffset + 1][colOffset + 1];

					}
				}
				sqrtVal = _mm_sqrt_ps(_mm_castsi128_ps(_mm_maddubs_epi16(Gx, Gy)));
				// Maddubs muls vertically and adds horizontally
				// It produces the equivalent of Gx * Gx + Gy * Gy, I think.
				_mm_storeu_ps(&tempArr[0], sqrtVal);

				for (i = 0; i < 4; i++)
				{
					gradient[M * row + col + i] = (unsigned char)tempArr[i];
				}
				//gradient[M * row + col] = (unsigned char)sqrt(Gx * Gx + Gy * Gy); /* Calculate gradient strength		*/
				//gradient[row][col] = abs(Gx) + abs(Gy); // this is an optimized version of the above


				// Assessment:
				/* rowOffset and colOffset change, but can easily be vectorised (no complex calculation)
				* "Fully unroll the two innermost loops", only rowOffset and colOffset need vectorising
				*	Problem with the above: It seems like row and column need to be vectorised, not their offsets.
				*	Offsets change the value too much to vectorise. (-1024 to
				* filt: const Arr
				* M: const
				* row + col: int updating outside of vectorisation
				* rowOffset + colOffset: int updating inside of vectorisation (+4/8, not ++)
				* GxMask + GyMask: const Arr
				* Gx and Gy: int updating inside loop, total is the sum of the loop, they do not update
				*			 themselves, so should be simpler to implement than Q1
				*/


			}
		}
	}
}

void SobelNoVec(unsigned char* filt, unsigned char* gradient)
{
	int row, col, rowOffset, colOffset;
	int Gx, Gy;

	/*---------------------------- Determine edge directions and gradient strengths -------------------------------------------*/
	for (row = 1; row < N - 1; row++) {
		for (col = 1; col < M - 1; col++) {

			Gx = 0;
			Gy = 0;

			/* Calculate the sum of the Sobel mask times the nine surrounding pixels in the x and y direction */
			for (rowOffset = -1; rowOffset <= 1; rowOffset++) {
				for (colOffset = -1; colOffset <= 1; colOffset++) {

					Gx += filt[M * (row + rowOffset) + col + colOffset] * GxMask[rowOffset + 1][colOffset + 1];
					Gy += filt[M * (row + rowOffset) + col + colOffset] * GyMask[rowOffset + 1][colOffset + 1];
				}
			}

			gradient[M * row + col] = (unsigned char)sqrt(Gx * Gx + Gy * Gy); /* Calculate gradient strength		*/
			//gradient[row][col] = abs(Gx) + abs(Gy); // this is an optimized version of the above

		}
	}
}


int read_image(char* filename, unsigned char** frame1, unsigned char** filt, unsigned char** gradient)	// This must run for every file in the folder.
{

	int c;
	FILE* finput;
	int i, j, temp;

	printf("\nReading %s image from disk ...", filename);
	finput = NULL;
	int result = openfile(filename, &finput);
	if (result == -1)
	{
		return result;
	}

	*frame1 = allocateArr();
	*filt = allocateArr();
	*gradient = allocateArr();

	if ((header[0] == 'P') && (header[1] == '2')) {
		for (j = 0; j < N; j++) {
			for (i = 0; i < M; i++) {

				if (fscanf(finput, "%d", &temp) == EOF)	// fscanf is a security vulnerability, as is fopen.
					exit(EXIT_FAILURE);

				(*frame1)[M * j + i] = (unsigned char)temp;
			}
		}
	}
	else if ((header[0] == 'P') && (header[1] == '5')) {
		for (j = 0; j < N; j++) {
			for (i = 0; i < M; i++) {
				c = getc(finput);
				(*frame1)[M * j + i] = (unsigned char)c;	// Application did not realise what frame1 was, brackets resolve the issue.
			}
		}
	}
	else {
		printf("\n problem with reading image");
		exit(EXIT_FAILURE);
	}


	fclose(finput);
	printf("\nimage successfully read from disc\n");

	return 0;
}



void write_image2(char* filename, unsigned char* output_image)
{

	FILE* foutput;
	int i, j;

	foutput = fopen(filename, "wb");

	printf("  Writing result to disk ...\n");

	if ((foutput = fopen(filename, "wb")) == NULL) {
		fprintf(stderr, "Unable to open file %s for writing\n", filename);
		exit(-1);
	}

	fprintf(foutput, "P2\n");
	fprintf(foutput, "%d %d\n", M, N);
	fprintf(foutput, "%d\n", 255);

	for (j = 0; j < N; ++j) {
		for (i = 0; i < M; ++i) {
			fprintf(foutput, "%3d ", output_image[M * j + i]);
			if (i % 32 == 31) fprintf(foutput, "\n");
		}
		if (M % 32 != 0) fprintf(foutput, "\n");
	}
	fclose(foutput);


}




int openfile(char* filename, FILE** finput)
{
	int x0, y0, x;

	//int aa;

	if ((*finput = fopen(filename, "rb")) == NULL) {
		fprintf(stderr, "Unable to open file %s for reading\n", filename);
		return -1;
	}

	if (fscanf(*finput, "%s", header) == EOF)
		exit(EXIT_FAILURE);

	x0 = getint(*finput);//this is M
	y0 = getint(*finput);//this is N
	printf("\t header is %s, while x=%d,y=%d", header, x0, y0);

	//CRITICAL POINT: AT THIS POINT YOU CAN ASSIGN x0,y0 to M,N 
	M = x0;
	N = y0;
	printf("\n Image dimensions are M=%d,N=%d\n",M,N);


	x = getint(*finput); /* read and throw away the range info */
	//printf("\n range info is %d",x);
	return 0;
}



//CRITICAL POINT: you can define your routines here that create the arrays dynamically; now, the arrays are defined statically.
unsigned char* allocateArr()	// All 3 are the same type and length, run in a for.
{
	
	return (unsigned char*)malloc(M * N * sizeof(unsigned char));

}

char* allocatePath(unsigned int size)
{
	return (char*)malloc(size);
}


int getint(FILE* fp) /* adapted from "xv" source code */
{
	int c, i, firstchar;//, garbage;

	/* note:  if it sees a '#' character, all characters from there to end of
	   line are appended to the comment string */

	   /* skip forward to start of next number */
	c = getc(fp);
	while (1) {
		/* eat comments */
		if (c == '#') {
			/* if we're at a comment, read to end of line */
			char cmt[256], * sp;

			sp = cmt;  firstchar = 1;
			while (1) {
				c = getc(fp);
				if (firstchar && c == ' ') firstchar = 0;  /* lop off 1 sp after # */
				else {
					if (c == '\n' || c == EOF) break;
					if ((sp - cmt) < 250) *sp++ = c;
				}
			}
			*sp++ = '\n';
			*sp = '\0';
		}

		if (c == EOF) return 0;
		if (c >= '0' && c <= '9') break;   /* we've found what we were looking for */

		/* see if we are getting garbage (non-whitespace) */
	   // if (c!=' ' && c!='\t' && c!='\r' && c!='\n' && c!=',')
		//	garbage=1;

		c = getc(fp);
	}

	/* we're at the start of a number, continue until we hit a non-number */
	i = 0;
	while (1) {
		i = (i * 10) + (c - '0');
		c = getc(fp);
		if (c == EOF) return i;
		if (c < '0' || c>'9') break;
	}
	return i;
}

void checkValues(unsigned char* grad1, unsigned char* grad2) {
	bool isCorrect = true;
	for (int i = 0; i < 10; i++)	// Fixed i cap as the arrays are always length 5
	{
		if (grad1[i] != grad2[i])	// No point using the floating point calculation, as these are whole values
		{							// A character cannot be between 'a' and 'b', it is either 'a' or 'b'.
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
