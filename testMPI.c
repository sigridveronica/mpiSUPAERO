#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>


#define g2m(i,j) (i+ NX*j)
#define proc2Glob(i, j) (i+ posProcX*accuracy + (j+posProcY*accuracy)*NX)
#define prova(i,j) (j*NX +i) + NX*process_Rank 
#define local(i, j) (j*accuracy +i)

int mandelF(int N, double rmax, int i, int j);

int main(int argc, char** argv){
   


      MPI_Init(&argc, &argv);

    int process_Rank, nProcs;

    //Number of cells for each processor
    int accuracy = 100;
    


    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    //Number of processors in X and Y to fill the grid (now it must be the same)
    int nProcsX, nProcsY;
    nProcsX = nProcs/2;
    nProcsY = nProcs/2;

    double *refvec;
  
    
    //Number of cells in the grid
    int NX, NY;

    NX = nProcsX*accuracy;
    NY = nProcsY*accuracy;


    //Field allocations
    //GlobalField: field for the whole grid, procField: field for each of the processors
    double* procField = NULL;
    double* GlobalField;

    if (process_Rank ==0){

    GlobalField =  (double*)malloc(NX*NY*sizeof(double));
    refvec = GlobalField +1;
    //procField = (double*) malloc((double)accuracy*(double)accuracy*(sizeof(double))); //accuracy == number of cells per processor in X and Y axis (now is the same)
    }


   MPI_Barrier(MPI_COMM_WORLD);

    int posProcX, posProcY; //Cartesian position of each processor. Ex; if 4 processors in total: (0,0) for processor 0, (0,1) for processor 1, (0,1) for processor 2, (1,0) for processor 3
       
        //Decide where goes each Processor
            posProcX = process_Rank%nProcsX;
            posProcY = process_Rank/nProcsY;


       int a;


int ref;
double *localField = NULL;
              

    localField  =  (double*)malloc(accuracy*accuracy*sizeof(double));
            
            posProcX = process_Rank%nProcsX;
            posProcY = process_Rank/nProcsY;

        for (int j=0; j<accuracy ; j++) {
            for (int i=0; i<accuracy; i++) {

                ref = proc2Glob(i, j);
                ref = mandelF(NX, 2, ref%NX, ref/NX); // this should be uncommented when it works

   
                localField[(accuracy*j+i)]= ref;


           
            }
        }

int nMessages;
nMessages =accuracy*accuracy;

//MPI_Send(localField, nMessages, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
MPI_Request request[nProcs];
            
MPI_Isend(localField, nMessages, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request[process_Rank]);

//printf("I am processor %d and indexes (%d,%d) \n", process_Rank, posProcX, posProcY);

//MPI_Barrier(MPI_COMM_WORLD);
//MPI_Waitall(nProcs, request, MPI_STATUSES_IGNORE);
  // MPI_Wait(&request[process_Rank], MPI_STATUS_IGNORE);

if (process_Rank==0){

 for (int k=0; k<nProcs;k++) {

    MPI_Recv(localField, nMessages, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   // MPI_Irecv (localField, nMessages, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, &request[k]);

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
                printf("%d ",a);
            }
            printf("\n");
        }


             
          free(GlobalField);
    }



   

    MPI_Finalize();
    return 0;   
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
int mandelF(int N, double rmax, int i, int j) {
    int result;
    double a1, a2, b1, b2, a1_2, b1_2, r2, module;

    rmax = 2; 
    r2 = rmax*rmax;


    int NX= N;
    int NY= N;

// Field of the grid to fill
double *Field;

Field =  (double*) malloc(NX*NY*sizeof(double));

//Rescale the values
double scaler = (double)10/(2*(double)N);
double scalei = (double)10/(2*(double)N);
//


//Cartesian center
double cx = (NX-1)/2*scaler;
double cy = (NY-1)/2*scalei;


 // real and imaginary part of the number on the given square of the grid
double nr;
double ni;

int iterations;

int print;
int maxIt = N;
//for  (int j=0; j<NY; j++){
   
    //for (int i=0; i<NX; i++) {

        nr = i*scaler- cx;
        ni = j*scalei- cy;

        module = r2 -1;
        a1 =0;
        b1=0;

        iterations =0;

        while (module <r2 && iterations<maxIt){
        //first z
        

        //z1
        a1_2 = a1*a1 -b1*b1;
        b1_2 = 2*a1*b1;

        //z2
        a2 = a1_2 +nr;
        b2 = b1_2 + ni;

        module = a2*a2+b2*b2;
        
        a1 = a2;
        b1 = b2;

        iterations ++;

        }// end exit

         if(iterations == maxIt) {
           // printf("*");
           result = 1;
        }
        else{
            //printf(" ");
            result = 0;
        }

    
    
   



return result;

}