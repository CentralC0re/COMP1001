/*
------------------DR VASILIOS KELEFOURAS-----------------------------------------------------
------------------COMP1001 ------------------------------------------------------------------
------------------COMPUTER SYSTEMS MODULE-------------------------------------------------
------------------UNIVERSITY OF PLYMOUTH, SCHOOL OF ENGINEERING, COMPUTING AND MATHEMATICS---
*/

// DO NOT RUN, CRASHES LAPTOP

#include <stdio.h> 
#include <stdlib.h> //for malloc
#include <windows.h>  
  
#define N 1000000000

void heap_crash();

int main() 
{ 

printf("\nHi before crashing the heap\n");

heap_crash();    

printf("\nIt looks like the heap did not crash\n");

system("pause");
return 0;
} 

void heap_crash(){

int i;

for ( i=0; i<N; i++) // do the following for N times
    {     
       int *ptr = (int *)malloc(sizeof(int));  // Allocate 8 bytes of memory without freeing it 
    }											// Causes a memory leak that continues until
												// A: Out of Memory (probable system crash)
												// B: Program crash
												// C: i == N

}

