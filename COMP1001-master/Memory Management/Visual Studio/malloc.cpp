/*
------------------DR VASILIOS KELEFOURAS-----------------------------------------------------
------------------COMP1001 ------------------------------------------------------------------
------------------COMPUTER SYSTEMS MODULE-------------------------------------------------
------------------UNIVERSITY OF PLYMOUTH, SCHOOL OF ENGINEERING, COMPUTING AND MATHEMATICS---
*/


#include <stdio.h> 
#include <stdlib.h> //this library allows for malloc to run
#include <windows.h>  
 
int main() 
{ 
  
    int* ptr;  // This pointer will hold the base address of the block created 
    int n;     //array size
    int i;
  
     
    n = 5;  // Get the number of elements for the array

    printf("\nThe number of elements is: %d\n", n); 
   
    ptr = (int*)malloc(n * sizeof(int)); // Dynamically allocate memory using malloc() 
    if (ptr == NULL) { // Check if the memory has been successfully allocated by malloc or not 
        printf("\nMemory not allocated.\n"); 
        exit(0); 
    } 

    else {  // Memory has been successfully allocated
  
         
        printf("\nMemory successfully allocated using malloc.\n"); 
  
        // initialize the array 
        for (i = 0; i < n; ++i) { 
            ptr[i] = i;				// Though watched value says it only has one element, it does update
        }							// You can see this with the memory window.
  
        // Print the array 
        printf("\nThe elements of the array are: "); 
        for (i = 0; i < n; ++i) { 
            printf("%d, ", ptr[i]);	// As array variables are pointers anyway, you can treat this as one.
        } 

        // Free the memory 
        free(ptr); 
        printf("\n\nMalloc Memory successfully freed.\n"); 
    } 
  
system("pause");
    return 0; 
} 

