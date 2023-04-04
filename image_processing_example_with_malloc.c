/* 
  a simple example of image processing: copy an initial image to an output image
  with the help of 2 functions read_image and write_image
  begin
    read the initial image from a file and copy the pixels in an "array" of pixels (image_in)
    copy the elements of this array to a second array (image_out)
    (this step will be replaced by a more elaborated processing)
    copy the elements of the second array in a file
  end

  the format of the files is the pgm format (alias P2, ASCII file that a standard text editor car process)
  the images are in black and white
  this format is understood by a lot of image readers (eog, display, gimp, ...)
  the (graphical) display of the images is done outside this program

  the pixel "arrays" are dynamically allocated large vectors of pixels

  to obtain the executable code:
  gcc -Wall -std=c99 image_processing_example_with_malloc.c -o image_processing_with_malloc
  to execute:
  ./image_processing_with_malloc
  to display:
  eog logo2.pgm

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define IMAGE_IN "logo.pgm"
#define IMAGE_OUT "logo2.pgm"

#define IMAGE_IN "logo.pgm"
#define IMAGE_OUT "logo2.pgm"
#define IMAGE_KERNEL_OUT "logo3.pgm"
#define M2(m,i,j,w) m[(i)*(w)+(j)]
#define KERNEL1 ((int[]) {0, 0, 0, 0, 1, 0, 0, 0, 0})//{0,0,0,0,1,0,0,0,0}//{0,-1,0,-1,5,-1,0,-1,0} // {{0,-1,0},{-1,5,-1},{0,-1,0}}
#define KERNEL2 ((int[]) {0,-1,0,-1,5,-1,0,-1,0})
#define MINDEX(m,i,j,w) m[(i-1)*w+(j-1)] //it is found in the position i-1, j-1 from the      
#define CONVOLUTION(Kernel, MATRIX, a, index_i, index_j, w, Convolution) \
    for (int i = 0; i < (a); i++) { \
        for (int j = 0; j < (a); j++) { \
            Convolution += M2(Kernel,(i),j,(a)) * M2(MATRIX, (index_i - 1 + i),(index_j - 1 + j),w); \
        } \
    }



/* function declarations */
int * read_image (char file_name[], int *p_h, int *p_w, int *p_levels);
void write_image (int *image, char file_name[], int h, int w, int levels);
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

int main(void)
{

  int i, j, h, w, levels ;
  struct timeval tdeb, tfin;
  int *image_in, *image_out;

  image_in = read_image (IMAGE_IN, &h, &w, &levels); 
  image_out = malloc (h*w*sizeof(int));

  int *kernelMatrix = KERNEL2; //Kernel to be applied
  int *matrix_Halo =  create_Halo_Matrix (image_in, w, h); //copies the matrix and gets a halo of 0 around
  int *hola =  allocSpace(w, h);

  /* image processing (just a copy in this example) */
  gettimeofday(&tdeb, NULL);
/*
  for (i = 0; i < h ; i++)
    for (j = 0; j < w; j++) {
         image_out[i*w+j] = image_in[i*w+j]; 
    }
*/
  int result = 0;
  for (int index_i = 1; index_i < h+1; index_i++){
      for (int index_j = 1; index_j < w+1; index_j++){
        result = 0;

          CONVOLUTION(kernelMatrix, matrix_Halo, 3, index_i, index_j, w+2, result);
          M2(image_out, (index_i-1),(index_j-1),w) = result;    
    }
  }
  
  gettimeofday(&tfin, NULL);
  printf ("computation time (microseconds): %ld\n",  (tfin.tv_sec - tdeb.tv_sec)*1000000 + (tfin.tv_usec - tdeb.tv_usec));

  write_image (image_out, IMAGE_OUT, h, w, levels);
  free (image_in);
  free (image_out); 
  free(matrix_Halo);

  return 0;
}

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
