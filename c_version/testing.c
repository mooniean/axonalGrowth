#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "functions2d.h"

static int row = 300, col = 300, tmax = 15, radius = 10, centreX, centreY;
static double vm=1, k=10, lambdal=0.001,lambdaf=0.005, dt=0.000005, dx = 0.5, temporaryValueFree, temporaryValueLinked, M=1, betam, gammain; //tempProteins, betam=0.01, /*dt = 0.000005,*/lambdap=0.001,  alpha=0.01,

/*
Variable definitions:
vm - 
k - Diffusion constant for free mRNA
lambdal - consumption rate of linked mRNA (rate at what it's destroyed)
lambdaf - consumption rate of free mRNA (rate at what it's destroyed)
M - Maximum Concentration of mRNA linked to microtubules
betam - binding probability
gamma - unbinding probability
*/
static char filename[100];

void mtubeFormation(double ** mtube){
  int row = 300, col = 300;

  FILE *final = fopen("mtube","r");
  FILE *coordinates = fopen("mtubeCoord", "w");

  if (final == NULL) {
    printf("Error - no matrix ");
    printf("\n");
    exit(1);
  }

  char buffer[512];
  double i,j;
  int m,n;

  while(fgets(buffer, sizeof(buffer), final) != NULL) {
    fscanf(final, "%lf %lf",&j,&i); // ficheiro guardado como y,x
    m=(int)floor(i);
    n=(int)floor(j);
    mtube[m][n] = 1;
    fprintf(coordinates,"%d %d \n",m,n);
  }

  FILE *f = fopen("shortmtube", "w");
  printingWorld(f, row,col,mtube);

  fclose(f);
  fclose(final);
  fclose(coordinates);
}


void stepVector (double ** vectorx, double ** vectory, double ** newflag, double ** phi){
  int mrow = 0, i=0, j, x,y,m,n;

  char buffer[512],c;
  double norma, tempx, tempy,dt = 0.05,D = 1, firsti, firstj;


  double ** xlap = initMatrix(row,col);
  double ** ylap = initMatrix(row,col);
  double ** checkFlag = initMatrix(row,col);
  allZeros(checkFlag,row,col,0);

  FILE *final = fopen("mtubeCoord","r");


  while( (c = fgetc(final)) != EOF ){
    if(c == '\n'){ mrow++;}  // 
  }
  
  double ** newMtube = initMatrix(mrow, 2);  // tem de estar aqui porque o mrow é definido aqui
  fclose(final);

  final = fopen("mtubeCoord","r");
  FILE *coordinates= fopen("smallTubeCoord", "w");

  while(fgets(buffer, sizeof(buffer), final) != NULL) {
    fscanf(final, "%d %d",&m,&n);
    newMtube[i][1]=m;
    newMtube[i][2]=n;

    if (i==0){
      firsti=m;
      firstj=n;
      i++;
      fprintf(coordinates,"%d %d \n",m,n);
    }
    else if( m!=newMtube[i-1][1] || n!=newMtube[i-1][2]){
      i++;
      fprintf(coordinates,"%d %d \n",m,n);
    }

  }
  fclose(final);
  fclose(coordinates);
  for (i=0; i<(mrow-5); i++){

    x=newMtube[i][1];
    y=newMtube[i][2];
    tempx=newMtube[i+5][1]-newMtube[i][1];
    tempy=newMtube[i+5][2]-newMtube[i][2];
    norma = sqrt(pow(tempx,2)+pow(tempy,2));
    vectorx[x][y]=phi[x][y]*tempx/norma;
    vectory[x][y]=phi[x][y]*tempy/norma;
    checkFlag[x][y]=1;

    if(i==0){
      continue;
    }

    m = newMtube[i-1][1];
    n = newMtube[i-1][2];
    if (fabs(vectorx[x][y]-vectorx[m][n]) > 0.3 ) {
      vectorx[x][y]= vectorx[m][n];
      vectory[x][y]= vectory[m][n];
    }
  }


  for (i=mrow-5; i<mrow; i++){
    x=newMtube[i][1];
    y=newMtube[i][2];
    if (checkFlag[x][y]==0){
      vectorx[x][y]=phi[x][y]*tempx/norma;
      vectory[x][y]=phi[x][y]*tempy/norma;
      checkFlag[x][y]=1;
    }
    // newflag[x][y]=1;
  }
  // double declive = (firstj-row/2)/(firsti-col/2), temp;
  // for (i=firsti; i<135; i++){
  // 	for (j=col/2; j<firstj; j++){
  // 		temp = (firstj-j)/(firsti-i);
  // 		if (abs(temp-declive)<0.000001){
  // 			tempx=firsti-i;
  // 			tempy=firstj-j;
  // 			norma = sqrt(pow(tempx,2)+pow(tempy,2));
  // 			vectorx[i][j]=tempx/norma;
  //  		vectory[i][j]=tempy/norma;
  //  		checkFlag[i][j]=1;

  // 		}
  // 	}
  // }
  /*Small vector diffusion*/
  double t;
  for ( t = 0; t<8; t+=dt){
    laplaceOperator2D(vectorx, xlap, row, col, 1);
    laplaceOperator2D(vectory, ylap, row, col, 1);
    for (i=1; i<col; i++){
      for (j=1; j<row; j++){
        if (checkFlag[i][j] == 1){
          continue;
        }
        if (checkFlag[i-1][j] == 1 ){
          vectorx[i][j]=vectorx[i-1][j];
          vectory[i][j]=vectory[i-1][j];
        }
        else if (checkFlag[i][j-1] == 1 ){
          vectorx[i][j]=vectorx[i][j-1];
          vectory[i][j]=vectory[i][j-1];
        }
        else{
          vectorx[i][j] = vectorx[i][j] + dt*D*xlap[i][j];
          vectory[i][j] = vectory[i][j] + dt*D*ylap[i][j];
        }
      }
    }
  }

  vectorx[0][0]=0;
  vectory[0][0]=0;

  FILE *vecx = fopen("vectorX","w");
  FILE *vecy = fopen("vectorY","w");
  FILE *flag = fopen("flag","w");
  printingWorld(flag, row,col,checkFlag);
  fclose(flag);
  printingWorld(vecx, row,col,vectorx);
  printingWorld(vecy, row,col,vectory);
  fclose(vecx);
  fclose(vecy);	

  freeMatrix(checkFlag,col,row);
  freeMatrix(newMtube,mrow, 2);
  freeMatrix(xlap,col,row);
  freeMatrix(ylap,col,row);

}

void createDensity(double ** vectorx, double ** vectory, double ** density){
  int i,j;
  for (i = 0; i<row; i++){
    for (j = 0; j<col; j++){
      density[i][j] = sqrt(pow(vectorx[i][j],2) + pow(vectory[i][j],2));
    }
  }
  char name[100];
  snprintf(name, sizeof(char) * 32, "%s_mtubedensity",filename);

  FILE *mtdens = fopen(name,"w");
  printingWorld(mtdens,row,col,density);
  fclose(mtdens);
}

void openPhiFile(double ** newphi){
  /* Opening file & getting phi matrix */
  //snprintf(name, sizeof(char) * 32, "A.%s.30000", argv[1]);
  char temporary[100];
  sprintf(temporary,"%s",filename);
  FILE *final = fopen(temporary,"r");

  int i,j;
  double ** phi = initMatrix(col,row);
  double max=0,min=0, value;
  //int imin, jmin;


  /*Calcular o Máximo/Mínimo*/
  for (j=0; j<row; j++){
    for (i=0; i<(col); i++){
      fscanf(final, "%lf",&value);
      phi[i][j]=value;
      if(phi[i][j]>1){
        phi[i][j]=1;
      }
      if (phi[i][j]>max){
        max = phi[i][j];
      }

      if (phi[i][j]<min){
        min = phi[i][j];
      }
    }
  }

  alterPhi(phi, newphi, col, row, max, min);
  char name[100];
  snprintf(name, sizeof(char) * 32, "%s_phi",filename);


  FILE *f = fopen(name, "w");
  printingWorld(f, col, row, newphi);
  fclose(f);


  freeMatrix(phi,row,col);
}

int main(int argc, char * argv[]){
  int i,j, m=0, concflag=0; char name[32]; double t, threshold = pow(10,-8), temp,concentration; 


  // Getting the variables to init the program with the nucleus in the correct place.
  FILE *initFile = fopen("initparameters.txt","r");
  fscanf(initFile, "%d",&row);
  fscanf(initFile, "%d",&col);
  fscanf(initFile, "%d",&radius);
  fscanf(initFile, "%d",&centreX);
  fscanf(initFile, "%d",&centreY);

  fclose(initFile);

  printf("Row: %d, Col: %d, radius: %d, centreX: %d, centreY: %d", row,col,radius,centreX,centreY);
  printf("\n");
  strcpy(filename, argv[1]); // This will be the name of the Phi file we're opening
  betam = atof(argv[2]); // atoi if char to int
  gammain = atof(argv[3]);

  //row = atoi(argv[2]);
  //col = atoi(argv[3]);
  //radius = atoi(argv[4]);
  //centreX = atoi(argv[5]);
  //centreY = atoi(argv[6]);
  printf("Gamma: %lf, betam: %lf \n", gammain, betam);
  printf("\n");

  double ** mtube   	  = initMatrix(row,col); // Microtubules
  double ** vectorx 	  = initMatrix(row,col); // X Component of Microtubule Vector
  double ** vectory 	  = initMatrix(row,col); // Y Component of Microtubule Vector
  double ** density 	  = initMatrix(row,col); // Microtubule Density
  double ** phi 	  	  = initMatrix(row,col); // Phi field
  double ** mfree    	  = initMatrix(row,col); // free mRNA concentration
  double ** mlinked 	  = initMatrix(row,col); // linked mRNA concentration
  double ** mfreePhi    = initMatrix(row,col); // free mRNA concentration multiplied by Phi
  double ** lapMfreePhi = initMatrix(row,col); // Laplacian of free mRNA multiplied by Phi
  double ** mlinkedVx	  = initMatrix(row,col); // X Component of linked mRNA Vector
  double ** mlinkedVy	  = initMatrix(row,col); // Y Component of linked mRNA Vector
  double ** divMlinkedV = initMatrix(row,col); // Divergence of linked mRNA Vector
  double ** newflag 	  = initMatrix(row,col); // Flag of the microtubules
  double ** oldmfreePhi = initMatrix(row,col); // Temporary phi matrix
  //double ** proteins    = initMatrix(row,col); // Proteins concentration
  //double ** protPhi     = initMatrix(row,col); // Proteins concentration multiplied by Phi
  //double ** lapprot     = initMatrix(row,col); // Laplacian of Proteins concentration multiplied by Phi
  // double ** mtotal 	  = initMatrix(row,col);

  allZeros(mtube,row,col,0);
  allZeros(vectorx,row,col,0);
  allZeros(vectory,row,col,0);
  allZeros(density,row,col,0);
  allZeros(mfree, col, row, 0);
  allZeros(mlinked, col, row, 0);
  allZeros(divMlinkedV,col,row,0);
  allZeros(mlinkedVx,col,row,0);
  allZeros(mlinkedVy,col,row,0);
  allZeros(lapMfreePhi,col,row,0);
  allZeros(mfreePhi,col,row,0);
  allZeros(phi,col,row,0);
  allZeros(newflag,col,row,0);
  allZeros(oldmfreePhi,col,row,0);
  //allZeros(proteins,col,row,0);
  //allZeros(protPhi,col,row,0);
  //allZeros(lapprot,col,row,0);

 

  openPhiFile(phi); 						// Create Phi field
  mtubeFormation(mtube); 					// Create Microtubules
  stepVector(vectorx,vectory,newflag,phi); 			// Create MT-Vectors
  createDensity(vectorx,vectory,density); // Create MT density

  printf("created stuff");
  printf("\n");

  //    double xMaxC = col*0.5 + 5; 					// Defining center
  //    double xMinC = col*0.5 - 5;
  //    double yMaxC = row*0.5 + 5;
  //    double yMinC = row*0.5 - 5;	

  //   	for (i=(xMinC); i<(xMaxC); i++){
  //   		for (j=(yMinC); j<(yMaxC); j++){
  //   			;
  //   		}
  //   	}
  FILE *new = fopen("mfreephi_final","r");
  if (new == NULL){
    for (i=(centreY-30); i<(centreY+30); i++){
      for (j=(centreX-30); j<(centreX+30); j++){
        //if (((i-row*0.5)*(i-row*0.5)+(j-col*0.5)*(j-col*0.5))<radius*radius){
        if (((i-centreY)*(i-centreY)+(j-centreX)*(j-centreX))<radius*radius){
          mfree[i][j]=1;
        }
      }
    }
  }
  else{
    for (j=0; j<row; j++){
      for (i=0; i<col; i++){
        fscanf(new, "%lf",&temp);
        oldmfreePhi[i][j]=temp;
      }
    }
    // Calcular a condição inicial do mfree
    for (j=0; j<row; j++){
      for (i=0; i<col; i++){
        if (phi[i][j] > threshold){
          mfree[i][j]=oldmfreePhi[i][j]/phi[i][j];
        }
      }
    }
    fclose(new);

    new = fopen("mlinked_final","r");
    for (j=0; j<row; j++){
      for (i=0; i<(col); i++){
        fscanf(new, "%lf",&temp);
        mlinked[i][j]=temp;
      }
    }
    fclose(new);
    /*
    new = fopen("mproteins_final","r");
    for (j=0; j<row; j++){
      for (i=0; i<(col); i++){
        fscanf(new, "%lf",&temp);
        proteins[i][j]=temp;
      }
    }
    fclose(new);
    */
  }


  int flag = 0;

  /*INITIALIZING TIME AND INTEGRATION-----------------------------------------------------------------------------------------------------------------------------*/
  for (t=0.0 ; t<(tmax); t=t+dt){
    if ((int)(t/dt)%500000==0 ){ // 500000 for dt=0.000005
      snprintf(name, sizeof(char) * 32, "%s_mfree_%d", filename,m);
      FILE *f = fopen(name, "w");
      printingWorld(f, col, row, mfree);
      fclose(f);

      snprintf(name, sizeof(char) * 32, "%s_mlinked_%d", filename,m);
      f = fopen(name, "w");
      printingWorld(f, col, row, mlinked);
      fclose(f);

      // snprintf(name, sizeof(char) * 32, "%s_mfreephi_%d", filename,m);
      // f = fopen(name, "w");
      // printingWorld(f, col, row, mfreePhi);
      // fclose(f);

      /*
      snprintf(name, sizeof(char) * 32, "%s_proteins_%d", filename,m);
      f = fopen(name, "w");
      printingWorld(f, col, row, proteins);
      fclose(f);
      */
      m++;
      // printf("m %d",m);
      // printf("\n");
    }

    //printf("t %0.1f",t);
    // if (m>20){
    //   break;
    // }

    phiUfunc(mfree, phi, mfreePhi, row, col);
    phiUfunc(mlinked, vectorx, mlinkedVx, row, col);
    phiUfunc(mlinked, vectory, mlinkedVy, row, col);
    laplaceOperator2D(mfreePhi, lapMfreePhi, row, col,dx);
    divergence2D(mlinkedVx, mlinkedVy, divMlinkedV, col, row, dx);

    //phiUfunc(proteins,phi,protPhi,row,col);
    //laplaceOperator2D(protPhi,lapprot,row,col,dx);

    for (i=0;i<row; i++){
      for (j=0;j<col;j++){
        // if (m<6){
        // if ( i>=xMinC && i<xMaxC && j>=yMinC && j<yMaxC ){
        // 	continue;
        // }
          //if (((i-row*0.5)*(i-row*0.5)+(j-col*0.5)*(j-col*0.5))<radius*radius){
        if (((i-centreY)*(i-centreY)+(j-centreX)*(j-centreX))<radius*radius){
            continue;
        }
        // }
        temporaryValueFree = 0;
        temporaryValueLinked = 0;
      //  tempProteins=0;
        temp = (M-mlinked[i][j])/M;

        if (phi[i][j] > threshold){
          temporaryValueFree = (k*lapMfreePhi[i][j] )/(phi[i][j]) - lambdaf*mfree[i][j] - density[i][j]*mfree[i][j]*temp*gammain;
          temporaryValueLinked = -vm*temp*divMlinkedV[i][j] - lambdal*mlinked[i][j] + density[i][j]*mfreePhi[i][j]*temp*gammain;
        //  tempProteins = (k*lapprot[i][j] )/(phi[i][j])-lambdap*proteins[i][j]+alpha*mfree[i][j];

          if(density[i][j] > threshold){
            temporaryValueFree+= betam*(1-density[i][j])*mlinked[i][j] / (phi[i][j]*density[i][j]);
            temporaryValueLinked-= betam*(1-density[i][j])*mlinked[i][j]/density[i][j];
          } 
        }
        mfree[i][j] = mfree[i][j] + dt*temporaryValueFree;
        mlinked[i][j] = mlinked[i][j] + dt*temporaryValueLinked;
        //proteins[i][j] = proteins[i][j] + dt*tempProteins;

        if (mlinked[i][j]<0){
          mlinked[i][j]=0;
        }

      }
    }

  }



  // freeMatrix(mtotal,col,row);
  snprintf(name, sizeof(char) * 32, "mfree_final");
  FILE *f = fopen(name, "w");
  printingWorld(f, col, row, mfree);
  fclose(f);

  snprintf(name, sizeof(char) * 32, "mlinked_final");
  f = fopen(name, "w");
  printingWorld(f, col, row, mlinked);
  fclose(f);

  snprintf(name, sizeof(char) * 32, "mfreephi_final");
  f = fopen(name, "w");
  printingWorld(f, col, row, mfreePhi);
  fclose(f);

  /*
  snprintf(name, sizeof(char) * 32, "mproteins_final");
  f = fopen(name, "w");
  printingWorld(f, col, row, proteins);
  fclose(f);
  
*/

  freeMatrix(vectorx,col,row);
  freeMatrix(vectory,col,row);
  freeMatrix(mtube,col,row);
  freeMatrix(phi,col,row);
  freeMatrix(mlinked,col,row);
  freeMatrix(mfree,col,row);
  freeMatrix(mfreePhi,col,row);
  freeMatrix(lapMfreePhi,col,row);
  freeMatrix(mlinkedVx,col,row);
  freeMatrix(mlinkedVy,col,row);
  freeMatrix(divMlinkedV,col,row);
  freeMatrix(density,col,row);
  freeMatrix(newflag,col,row);
  freeMatrix(oldmfreePhi,col,row);
  //freeMatrix(proteins,col,row);
  //freeMatrix(protPhi,col,row);
  //freeMatrix(lapprot,col,row);



  return 0;
}
