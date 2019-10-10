#ifndef functions2d_h
#define functions2d_h

#include <stdio.h>

double ** initMatrix(int Lx, int Ly);
void freeMatrix(double ** matrix,int Lx, int Ly);
void allZeros(double ** vector,int Lx, int Ly, double value);
void printingWorld(FILE *f, int Lx, int Ly,double **world);
int temporaryFunction(int i, int length);
void laplaceOperator2D(double **world, double **temporaryResult, int Lx, int Ly, double dx);
double dotProduct(double vectorOne[], double vectorTwo[],int n);
void assignVector(double vector[], int a, int b);
void gradient2D(double **world, double ** resX, double ** resY, int Lx, int Ly, double dx);
void divergence2D(double **worldX, double ** worldY, double **temporaryResult, int Lx, int Ly, double dx);
void phigradu(double **phi, double **worldX, double ** worldY, int Lx, int Ly);
void phiUfunc(double **phi, double **u, double **phiu, int Lx, int Ly);
void totalMatrix(double **phi, double **u, double **phiu, int Lx, int Ly);
void alterPhi(double **phi, double **newphi, int Lx, int Ly, double max, double min);
void iteration(double doublec[], int intc[], FILE *areaFile, FILE *volumeFile, double *** matrix);

#endif