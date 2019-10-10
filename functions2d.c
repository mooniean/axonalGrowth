#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <ctype.h>
#include "functions2d.h"

/*Initializing a 3D matrix*/
double ** initMatrix(int Lx, int Ly){
	int i;
	double ** matrix = malloc(Lx * sizeof(double*));
	for(i = 0; i < Lx; i++){
		matrix[i] = malloc(Ly * sizeof(double));
	}
	return matrix;
}

/*Freeing a 3D matrix*/
void freeMatrix(double **matrix,int Lx, int Ly){
	int i;
	for(i = 0; i < Lx; i++){
		free(matrix[i]);
	}
	free(matrix);
}

/*Fill a matrix with zeros*/
void allZeros(double ** vector,int Lx, int Ly, double value){
	int i,j;
	for (i=0; i<Lx; i++){
		for (j=0; j<Ly; j++){
				vector[i][j]=value;				
		}
	}
}

/*Printing (region > 0 ) in a file*/
void printingWorld(FILE *f, int Lx, int Ly, double **world){
	int i,j;
	for (j=0; j<Ly; j++){
		for (i=0; i<Lx; i++){
			fprintf(f,"%lf ",world[i][j]);
		}
		fprintf(f,"\n");
	}
}

int temporaryFunction(int i, int length){
    return (i+length)%length;
}

void laplaceOperator2D(double **world, double **temporaryResult, int Lx, int Ly,double dx){
    int i, j;
    double dx2 = 1/(dx*dx);
    for (i = 0; i < Lx; i++){
        for (j = 0; j < Ly; j++){
   			temporaryResult[i][j]= dx2 * (world[temporaryFunction(i+1,Lx)][j] + 
               	           				  world[temporaryFunction(i-1,Lx)][j] + 
                   		    			  world[i][temporaryFunction(j+1,Ly)] + 
                       					  world[i][temporaryFunction(j-1,Ly)] -
                      					  4*world[i][j]);
   		}
	}
}

double dotProduct(double vectorOne[], double vectorTwo[],int n){
	double result=0.0;
	int i;
	for (i = 0; i < n; i++){
        result += vectorOne[i]*vectorTwo[i];
    }
	return result;
}

void gradient2D(double **world, double ** resX, double ** resY, int Lx, int Ly, double dx){
	int i, j;
	double dxNew = 0.5/dx;
    for (i = 0; i < Lx; i++){
        for (j = 0; j < Ly; j++){
   			resX[i][j] = dxNew * (world[temporaryFunction(i+1,Lx)][j] - 
               	           		  world[temporaryFunction(i-1,Lx)][j]) ;
   			resY[i][j] = dxNew * (world[i][temporaryFunction(j+1,Ly)] - 
                       			  world[i][temporaryFunction(j-1,Ly)]);
   		}
	}

}

void divergence2D(double **worldX, double ** worldY, double **temporaryResult, int Lx, int Ly, double dx){
	int i, j;
	double dxNew = 0.5/dx;
    for (i = 0; i < Lx; i++){
        for (j = 0; j < Ly; j++){
   			temporaryResult[i][j]= (dxNew) * (worldX[temporaryFunction(i+1,Lx)][j] - 
               	           					  worldX[temporaryFunction(i-1,Lx)][j] + 
                   		    				  worldY[i][temporaryFunction(j+1,Ly)] - 
                       						  worldY[i][temporaryFunction(j-1,Ly)]) ;
   		}
	}
}

void phigradu(double **phi, double **worldX, double ** worldY, int Lx, int Ly){
	int i, j;
    for (i = 0; i < Lx; i++){
        for (j = 0; j < Ly; j++){
   			worldX[i][j]= phi[i][j]*worldX[i][j]; 
   			worldY[i][j]= phi[i][j]*worldY[i][j];
   		}
	}
}

void phiUfunc(double **phi, double **u, double **phiu, int Lx, int Ly){
	int i, j;
    for (i = 0; i < Lx; i++){
        for (j = 0; j < Ly; j++){
   			phiu[i][j]= phi[i][j]*u[i][j];
   		}
	}
}

void totalMatrix(double **phi, double **u, double **phiu, int Lx, int Ly){
	int i, j;
    for (i = 0; i < Lx; i++){
        for (j = 0; j < Ly; j++){
   			phiu[i][j]= phi[i][j]+u[i][j];
   		}
	}
}
void alterPhi(double **phi, double **newphi, int Lx, int Ly, double max, double min){
	int i, j;// imin, jmin;
	double tempPhi;

    for (i = 0; i < Lx; i++){
        for (j = 0; j < Ly; j++){
        	tempPhi = (1+phi[i][j])*0.5;
        	if (tempPhi < 0.1){
        		newphi[i][j] = 0;
        	}
        	else{
        		newphi[i][j] = tempPhi;
        		//newphi[i][j] = (phi[i][j]-min)/(max-min);
        		//tempPhi = (phi[i][j]-min)/(max-min);// rescaling para 0 e 1
       			//newphi[i][j]= 0.5*(tempPhi+1); // otherwise, em vez de tempPhi está um phi[i][j] aqui
        	}
   		}
	}
	//newphi[imin][jmin]=pow(10,-4);
}

void assignVector(double vector[], int a, int b){
	vector[0]=a;
	vector[1]=b;
}

void iteration(double doublec[], int intc[], FILE *areaFile, FILE *volumeFile, double *** world){
	int i,j,k;
}