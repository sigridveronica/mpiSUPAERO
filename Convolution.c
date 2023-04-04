//  gcc -Wall -std=c99 debug.c -o debug  --> Run this command to compile
//  ./debug --> Run this command to run
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define IMAGE_IN "logo.pgm"
#define IMAGE_OUT "logo2.pgm"
#define IMAGE_KERNEL_OUT "logo3.pgm"
#define M2(m,i,j,w) m[(i)*(w)+(j)]
#define KERNEL1 ((int[]) {0, 0, 0, 0, 1, 0, 0, 0, 0})//{0,0,0,0,1,0,0,0,0}//{0,-1,0,-1,5,-1,0,-1,0} // {{0,-1,0},{-1,5,-1},{0,-1,0}}
#define MINDEX(m,i,j,w) m[(i-1)*w+(j-1)] //it is found in the position i-1, j-1 from the      
#define CONVOLUTION(Kernel, MATRIX, a, index_i, index_j, w, Convolution) \
    for (int i = 0; i < (a); i++) { \
        for (int j = 0; j < (a); j++) { \
            Convolution += M2(Kernel,(i),j,(a)) * M2(MATRIX, (index_i - 1 + i),(index_j - 1 + j),w); \
        } \
    }
// a: nRows = nCols of the matrix Kernel (assuming is square)
// index_i, index_j: Center indexes for the multiplication
// w: nCols of MATRIX 
// Convolution: Result of the multiplication

int * read_image (char file_name[], int *p_h, int *p_w, int *p_levels);
void write_image (int *image, char file_name[], int h, int w, int levels);

int *allocSpace(int w, int h);
void copy_fullMatrix(int* copied_Matrix,int *new_Matrix, int w, int h);
int *create_Halo_Matrix (int * Matrix, int w, int h);
int *debugMatrix(int w, int h);
void pprintf(int *matrix, int w, int h);

int main(void){
  int w = 4;
  int h = 4;


  int *image_kernel_out = debugMatrix(w,h); //create a sample matrix
  //Call the kernel: it allocates the Minor and frees it inside
  int *kernelMatrix = KERNEL1; //ok
  int *matrix_Halo =  create_Halo_Matrix (image_kernel_out, w, h); //copies the matrix and gets a halo of 0 around
  int *hola =  allocSpace(w, h);


  int result = 0;
  for (int index_i = 1; index_i < h+1; index_i++){
      for (int index_j = 1; index_j < w+1; index_j++){
        result = 0;

          CONVOLUTION(kernelMatrix, matrix_Halo, 3, index_i, index_j, w+2, result);
          M2(hola, (index_i-1),(index_j-1),w) = result;    
    }
  }
 
  pprintf(image_kernel_out, w,h);
  free (image_kernel_out);  //ok
  free(matrix_Halo); //ok
  free(hola);

  return 0;

}

// DEFINE KERNEL 



//Matrix:
//w: n cols
//h: n rows

//space allocation:
int *allocSpace(int w, int h) {
  int *matrix;
  matrix = calloc(h*w,sizeof(int));
  if (matrix == NULL) {
    printf("Memory allocation failed\n");
    exit(1);
  }
  return matrix;
}



//this function creates a new matrix with a halo of one cell from 
//an original sample matrix
int *create_Halo_Matrix (int * Matrix, int w, int h) {
    int *matrix_withHalo =  allocSpace(w+2, h+2); //halo of 1 more cell
    // we fill the new matrix
     copy_fullMatrix(Matrix,matrix_withHalo, w, h);
  return matrix_withHalo;
}

//copy the fullMatrix exactly as it is
void copy_fullMatrix(int* copied_Matrix,int *new_Matrix, int w, int h) {   
    for (int i=0; i<h; i++){
      for (int j = 0; j < w; j++) {
        M2(new_Matrix,(i+1),(j+1),(w+2)) =  M2(copied_Matrix, i, j,w);
    } //macro defined on top
  }
}


int *debugMatrix(int w, int h) {
    int *matrix =  allocSpace(w, h); //halo of 1 more cell
    // we fill the new matrix
     for (int i=0; i<h; i++){
      for (int j=0; j<w; j++){
        M2(matrix, i,j,w) = i*w+j;
      }
     }
  return matrix;
}

void pprintf(int *matrix, int w, int h){
  for (int i=0; i<h ;i++){
    for (int j=0; j<w; j++){ 
        printf("%d  ",M2(matrix, i,j,w));
      }
      printf("\n");
    }
}


// CONVOLUTION PRODUCT:
//========================================================================================
        /* for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              result += M2(kernelMatrix,i,j,3) * M2(matrix_Halo, (index_i - 1 + i),(index_j - 1 + j),w+2);
               //if (M2(kernelMatrix,i,j,3)* M2(matrix_Halo, (index_i - 1 + i),(index_j - 1 + j),w) != 0){
                //printf(" %d %d ", M2(kernelMatrix,i,j,3), M2(matrix_Halo, (index_i - 1 + i),(index_j - 1 + j),w)) ;
               //}
               
             }
            printf("\n");
          }
          */





