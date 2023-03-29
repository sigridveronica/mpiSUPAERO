#include <stdio.h>
#include <stdlib.h>


#define g2f(i,j)= ((NX-1)*i=j)

const int N = 10;

int main(int argc, char **argv) {

    double a1, a2, b1, b2, a1_2, b1_2, rmax, r2, module;

    rmax = 2; 
    r2 = rmax*rmax;


    int NX= N;
    int NY=N;

// Field of the grid to fill
double *Field;

Field =  (double*) malloc(NX*NY*sizeof(double));

//Rescale the values
double scaler = 0.1;
double scalei = 0.1;
//


//Cartesian center
double cx = (NX-1)/2*scaler;
double cy = (NY-1)/2*scalei;


 // real and imaginary part of the number on the given square of the grid
double nr;
double ni;

int iterations;

int print;
int maxIt = 100;
for  (int j=0; j<NY; j++){
   
    for (int i=0; i<NX; i++) {

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
        b1= b2;

        iterations ++;
        //printf("module = %f\n", module);
        


        }// end exit
        //printf("(%d, %d) ", (int)nrnumber of iterations until divergence , (int)ni);double cx = 5;
        if(iterations == maxIt) {
            printf("1");
        }
        else{
            printf("0");
        }


        
    }
    printf("\n");
}


//printf("\n");

}