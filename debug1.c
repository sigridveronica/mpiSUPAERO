#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define g2m(i,j) (i+ NX*j)
#define proc2Glob(i, j) (i+ posProcX*w_new/nProcsX + (j+posProcY*h_new/nProcsY)*w_new); //works with sqare and non square matrixes
#define prova(i,j) (j*NX +i) + NX*process_Rank
#define local(i, j) (j*accuracy +i)

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
//void resizeMatrix(int* matrix, int h, int w, int NX, int NY, int* h_new, int* w_new);
void print_matrix(int* matrix, int nrows, int ncols);
void find_close_factors(int z, int *x, int *y);
int* resizeMatrix(int* matrix, int h, int w, int NX, int NY, int* h_new, int* w_new);

#define AA 2048
#define BB 2048
unsigned char color[AA*BB];

int main(int argc, char** argv){
    


    MPI_Init(&argc, &argv);

    int i, j, h, w, levels ;
    w = 6;
    h = 5;
    int h_new, w_new; 
    int nProcsX;
    int nProcsY;
    int* matrix_Halo;
    int *N;


    //Number of processors for each cell:
    int process_Rank, nProcs;

    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

  //Number of cells in the grid
    int NX, NY;
    NX = w;//nProcsX*accuracy;
    NY = h;//nProcsY*accuracy;
  //Allocate number of processors in x and y:

  //find_close_factors(nProcs, &nProcsX, &nProcsY);
  nProcsX = 2;
  nProcsY = 2;
  nProcs = nProcsX*nProcsY;
  


//Allocation of the image in and out fields and the Kernel Matrix
//- matrix_Halo = full field, with the halo of 0 

if (process_Rank ==0) {
  
  //image_in = read_image (IMAGE_IN, &h, &w, &levels); 
  int *image_kernel_out = debugMatrix(w,h); //create a sample matrix
  //Create kernel:
  
  matrix_Halo =  create_Halo_Matrix (image_kernel_out, w, h); //copies the matrix and gets a halo of 0 around
  pprintf(matrix_Halo, w+2, h+2);
  printf("\n");
  matrix_Halo = resizeMatrix(matrix_Halo, (h+2), (w+2), nProcsX, nProcsY, &h_new, &w_new); //add zeros to matrix_HALO to make it fit with nProcsX, nProcs
  pprintf(matrix_Halo, w_new, h_new);

  // we have w_new and h_new that are divisible by nProcsX and nProcsY, respectivelly
  MPI_Bcast(&w_new, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&h_new, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  free(image_kernel_out); 
  //pprintf(matrix)
 
  //printf("%d HOLAAAA",w_new);

}

MPI_Barrier(MPI_COMM_WORLD);
int *kernelMatrix = KERNEL1; //Kernel to be applied //QUITAR: tenemos que ponerlo fuera del processor 0?
  
int w_submatrix = w_new/nProcsX;
int h_submatrix = h_new/nProcsY;
int nMessages = w_new/nProcsX*h_new/nProcsY;

int ncols =  w_new/nProcsX;
int nrows =  h_new/nProcsY;


// Allocate memory for the submatrix on each process
N = malloc(ncols * nrows * sizeof(int));
printf("ncols = %d, nrows = %d, w = %d, h= %d, wnew = %d, hnew= %d", ncols, nrows, w, h, w_new, h_new);

printf("size M %ld, size  %ld \n", sizeof(N), sizeof(matrix_Halo));
  // Calculate the send counts and displacements for MPI_Scatterv()
  int *send_counts = malloc(nProcs * sizeof(int));
  int *displs = malloc(nProcs * sizeof(int));
  for (int i = 0; i < nProcs; i++) {
      send_counts[i] = ncols * nrows;
      displs[i] = (i / nProcsY) * nrows * w + (i % nProcsX) * ncols;
  }


    // Scatter the matrix M to each process and store the submatrix in N
MPI_Scatterv(matrix_Halo, send_counts, displs, MPI_INT, N, ncols * nrows, MPI_INT, 0, MPI_COMM_WORLD);

printf("Process %d received submatrix:\n", process_Rank);
   for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            printf("%d ", N[i * ncols + j]);
        }
        printf("\n");
    }

free(send_counts);

free(displs);

MPI_Finalize();
    return 0;

  
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

//Matrix:
//w: n cols
//h: n rows
/*
* read an image from a file  and dynamic memory allocation for it
* param1 : name of the image file
* param2 : where to store the height
* param3 : where to store the weight
* param4 : where to store the number of levels
* return : address of the pixels
*/
int * read_image (char file_name[], int *p_h, int *p_w, int *p_levels){
  FILE *fin;
  int i, j, h, w, levels ;
  char buffer[80];
  int *image;

  /* open P2 image */
  fin = fopen (IMAGE_IN, "r");
  if (fin == NULL){
    printf ("file open error\n");
    exit(-1);
  } 
  fscanf (fin, "%s", buffer);
  if (strcmp (buffer, "P2")){
    printf ("the image format must be P2\n");
    exit(-1);
  }
  fgets (buffer, 80, fin);
  fgets (buffer, 80, fin);
  fscanf (fin, "%d%d", &w, &h);
  fscanf (fin, "%d", &levels);
  printf ("image reading ... h = %d w = %d\n", h, w);
  
  /* dynamic memory allocation for the pixels */
  image = malloc (h*w*sizeof(int));

  /* pixels reading */
  for (i = 0; i < h ; i++)
    for (j = 0; j < w; j++)
       fscanf (fin, "%d", image +i*w +j); 
  fclose (fin);

  *p_h = h;
  *p_w = w;
  *p_levels=levels;
  return image;
}
 /*
* write an image in a file
* param1 : address of the pixels
* param2 : name of the image file
* param3 : the height
* param4 : the weight
* param5 : the number of levels
* return : void
*/

void write_image (int *image, char file_name[], int h, int w, int levels){
  FILE *fout;
  int i, j;

  /* open the file */
  fout=fopen(IMAGE_OUT,"w");
  if (fout == NULL){
    printf ("file opening error\n");
    exit(-1);
  }
  
  /* header write */
  fprintf(fout,"P2\n# test \n%d %d\n%d\n", w, h, levels);
  /* format P2, commentary, w x h points, levels */

  /* pixels writing*/
  for (i = 0; i < h ; i++)
    for (j = 0; j < w; j++)
       fprintf (fout, "%d\n", image[i*w+j]); 

  fclose (fout);
  return;
}

//Matrix:
//w: n cols
//h: n rows

///space allocation:
int *allocSpace(int w, int h) {
  int *matrix;
  matrix = calloc(h*w,sizeof(int));
  for (int i=0;i<h; i++){
    for(int j=0; j<w; j++){
        M2(matrix, i, j,w) = 0;
    }
  }
  if (matrix == NULL) {
    printf("Memory allocation failed\n");
    exit(1);
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

int* resizeMatrix(int* matrix, int h, int w, int NX, int NY, int* h_new, int* w_new) {
    // Update new dimensions
    *h_new = h;
    *w_new = w;
    // Check if dimensions meet the requirements
    while (*w_new % NX != 0) {
        *w_new= *w_new+1;
    }
    while (*h_new % NY != 0) {
        *h_new = *h_new+1;
    }
    printf("H %d, W %d h_new %d, w_new %d \n\n",h, w, *h_new, *w_new );
    // Allocate memory for resized matrix
    int* resizedMatrix = allocSpace((*w_new), (*h_new));  
    // Copy values from original matrix to resized matrix
    for (int i = 0; i < *h_new; i++) {
        for (int j = 0; j < *w_new; j++) {
            if (i < h && j < w) {
                resizedMatrix[i * (*w_new) + j] = matrix[i * w + j];
                //printf("%d, \n ", matrix[i * w + j]);
            } else {
                resizedMatrix[i * (*w_new) + j] = 0; // Fill extra rows and columns with zeros
            }
        }
    }
    // Update matrix pointer to resized matrix
    
   return resizedMatrix;
}

//this function creates a new matrix with a halo of one cell from 
//an original sample matrix
void copy_fullMatrix(int* copied_Matrix,int *new_Matrix, int w, int h) {   
    for (int i=0; i<h; i++){
      for (int j = 0; j < w; j++) {
        M2(new_Matrix,(i+1),(j+1),(w+2)) =  M2(copied_Matrix, i, j,w);
    } //macro defined on top
  }
}

int *create_Halo_Matrix (int * Matrix, int w, int h) {
    int *matrix_withHalo =  allocSpace(w+2, h+2); //halo of 1 more cell
    // we fill the new matrix
     copy_fullMatrix(Matrix,matrix_withHalo, w, h);
  return matrix_withHalo;
}
//copy the fullMatrix exactly as it is

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


void reshape_matrix(int* i, int w, int h, int NX, int NY, int* w_new, int* h_new) {
    int new_w = w, new_h = h;
    if (w % NX != 0) {
        new_w = ((w - 1) / NX + 1) * NX;
    }
    if (h % NY != 0) {
        new_h = ((h - 1) / NY + 1) * NY;
    }
    if (new_w != w || new_h != h) {
        int* j = (int*) malloc(NX * NY * sizeof(int));
        int i_row, j_row, i_col, j_col;
        for (j_row = 0; j_row < NY; j_row++) {
            for (j_col = 0; j_col < NX; j_col++) {
                if (j_row < h && j_col < w) {
                    *(j + j_row * NX + j_col) = *(i + j_row * w + j_col);
                } else {
                    *(j + j_row * NX + j_col) = 0;
                }
            }
        }
        //print_matrix(j, NY, NX);
        free(j);
    } else {
        //print_matrix(i, h, w);
    }
}

void print_matrix(int* matrix, int nrows, int ncols) {
    int i, j;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            printf("%d ", *(matrix + i*ncols + j));
        }
        printf("\n");
    }
}

void find_close_factors(int z, int *x, int *y) {
    int i;

    for (i = 2; i < z; i++) {
        if (z % i == 0) {
            *x = i;
            *y = z / i;
            return;
        }
    }

    *x = 1;
    *y = z + 1;
}

