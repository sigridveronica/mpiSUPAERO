#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#define g2m(i,j) (i+ NX*j)
#define proc2Glob(i, j) (i+ posProcX*NX + (j+posProcY*h)*NY) //QUITAR: De posicion local a global en la matriz: REVISAR MACRO!! ncellsx, ncellsy?
//#define proc2Glob(i, j) (i+ posProcX*w + (j+posProcY*h)*NX) - Raquel: añadir un macro para no cuadradas
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


#define AA 2048
#define BB 2048
unsigned char color[AA*BB];

int main(int argc, char** argv){
    int borrar =2;
    FILE *fp;

    int i, j, h, w, levels ;
    struct timeval tdeb, tfin;
    int *image_in, *image_out;


    


      MPI_Init(&argc, &argv);

    //Number of processors for each cell:

    int process_Rank, nProcs;

    //Number of cells for each processor
    //int accuracy = 1000;



    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

//LO QUE HAGO ES SACAR LOS DIVISORES DEL NÚMERO DE PROCESADORES COGER EL DIVISOR QUE ESTÉ MÁS EN LA MITAD Y ESE ES EL NUM DE PROC EN X
    //Number of processors in X and Y to fill the grid (now it must be the same)
    int nProcsX, nProcsY, n, i, j=0;

    // create an array to store divisors
  int divisors[nProcs];

  // find divisors of n and store in the array
  for (i = 1; i <= nProcs; i++) {
    if (nProcs % i == 0) {
      divisors[j] = i;
      j++;
    }
  }
if (j%2!=0){
nProcsX = divisors[(j+1)/2];
}
else{
nProcsX = divisors[j/2];
}
nProcsY=nProcs/nProcsX;

    double *refvec;


    //Number of cells in the grid
    int NX, NY;

    NX = w;//nProcsX*accuracy;
    NY = h//nProcsY*accuracy;

  //n cells that each processor will get
    //int ncellsx = w/nProcsX;
    //int ncellsy = h/nprocsY;
    //int ncells_last_proc=

  //HE INTENTADO HACER ESTE ALGORITMO PARA DIVIDIR EL NUM DE CELDAS POR PROCESADOR HACIENDO LO QUE DIJIMOS DE DEJAR UNO LIBRE POR SI NO ES DIVISIBLE EXACTO
    int ncellsx, ncells_last_procx, ncellsy, ncells_last_procy;
    if (w%nProcsX=!0){
      ncellsx = NX/(nProcsX-1);
      ncells_last_procx = NX-ncellsx*(nProcsX-1);
    }
    else{
      ncellsx = NX/nProcs;
      ncells_last_procx = NX/nProc;
    }

    if (h%nProcsY=!0){
      ncellsy = NY/(nProcsY-1);
      ncells_last_procy = NY-ncellsy*(nProcsY-1);
    }
    else{
      ncellsx = NY/nProcs;
      ncells_last_procy = NY/nProc;
    }
    

    //Allocation of the image in and out fields and the Kernel Matrix
    //- matrix_Halo = full field, with the halo of 0 
    if (process_Rank ==0) {
    image_in = read_image (IMAGE_IN, &h, &w, &levels); 
    image_out = malloc (h*w*sizeof(int));
    int *kernelMatrix = KERNEL2; //Kernel to be applied //QUITAR: tenemos que ponerlo fuera del processor 0?
    int *matrix_Halo =  create_Halo_Matrix (image_in, w, h); //copies the matrix and gets a halo of 0 around
    
    }

  
    //Field allocations
    //GlobalField: field for the whole grid, procField: field for each of the processors
    int* procField = NULL;
    int* GlobalField;
 
    /*
    if (process_Rank ==0){
      GlobalField =  (int*)malloc(NX*NY*sizeof(int));
      refvec = GlobalField +1; //QUITAR: esto que es?
      //procField = (double*) malloc((double)accuracy*(double)accuracy*(sizeof(double))); //accuracy == number of cells per processor in X and Y axis (now is the same)
    }*/


   MPI_Barrier(MPI_COMM_WORLD); //QUITAR: para que era esto?

    int posProcX, posProcY; //Cartesian position of each processor. Ex; if 4 processors in total: (0,0) for processor 0, (0,1) for processor 1, (0,1) for processor 2, (1,0) for processor 3

        //Decide where goes each Processor //QUITAR: Cambiar esto a no cuadrados
            posProcX = process_Rank%nProcsX;
            posProcY = process_Rank/nProcsY;


       int a;


int ref;
int iter;
double *localField = NULL;


    localField  =  (int*)malloc(NX*NY*sizeof(double));

            posProcX = process_Rank%nProcsX; //determine the location of each processor in the grid
            posProcY = process_Rank/nProcsY;

      /*  for (int j=0; j<w ; j++) {
            for (int i=0; i<h; i++) {

                //ref = proc2Glob(i, j);
                iter = proc2Glob(i, j); //NEW
                //ref = mandelF(NX, 2, ref%NX, ref/NX); // this should be uncommented when it works
                iter = mandelF(NX,2, iter%NX, iter/NX); // New



                //localField[(accuracy*j+i)]= ref;
                localField[(accuracy*j+i)]= iter; //NEW



            }
        }*/

        int result = 0;
        for (int index_i = 1; index_i < h+1; index_i++){
          for (int index_j = 1; index_j < w+1; index_j++){
        result = 0;

          CONVOLUTION(kernelMatrix, matrix_Halo, 3, index_i, index_j, w+2, result);
          M2(image_out, (index_i-1),(index_j-1),w) = result; 
          proc2Glob(index_i-1, index_j-1) = result;   // QUITAR1: ESTAS aqui
         }
       }

int nMessages;
//nMessages = ncellsx*ncellsy;
nMessages = (ncellsx+ncells_last_procx)*(ncellsy+ncells_last_procy);

//MPI_Send(localField, nMessages, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
MPI_Request request;

MPI_Isend(localField, nMessages, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);

//printf("I am processor %d and indexes (%d,%d) \n", process_Rank, posProcX, posProcY);

MPI_Barrier(MPI_COMM_WORLD);



if (process_Rank==0){
 fp = fopen("color.txt","w"); //File
 for (int k=0; k<nProcs;k++) {

    MPI_Recv(localField, nMessages, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    posProcX = k%nProcsX;
    posProcY = k/nProcsY;



  for (int i=0; i<accuracy; i++) {
     for (int j=0; j<accuracy ; j++) {

                GlobalField[proc2Glob(i,j)] = localField[i+j*(accuracy)];

            }

     }
 }



     // for (int j=NY-1; j>=0; j--){
         for (int j=0; j<NX; j++){
     for(int i=0; i<NX ; i++) {



                //GlobalField[g2m(i,j)] = g2m(i,j);
                //printf("%d ", (int)GlobalField[g2m(i,j)]);

                a = (int)(GlobalField[NX*j+i]);
                //printf("%d ",a);
                fprintf(fp, "%hhu ", (unsigned char)GlobalField[j+i*NX]); // File
            }
            //printf("\n");
            fprintf(fp, "\n"); //File
        }
        fclose(fp); // File


          free(GlobalField);
    }

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