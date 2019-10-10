//imprimir ficheiro da morte

#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>

using std::cout;
using std::ofstream;
using std::ifstream;
using namespace std;

const int Lx=300, Ly=300;
static const int Arteria=40; // initial artery width
static const int centreX=Lx/2, centreY=Ly-(Arteria+5);
static const int Nprint=1000; // printing interval
static const double dt=0.001, tmax=50.001;  // integration step and final time
static const double Tang=100000; // tempo para ligar Ang2
static const double gradTmin=0.01, C=0.7; // T gradient and T value for existance of CA - Cellular Automaton

static const double Rbolas=12.0, Rbolas2=Rbolas*Rbolas; // radius, number of CA
static const int NbolasMax=200; //HUGO: Antes estava 500

static const double D1=100, epsilon=1, M=1, alpha1=D1/Rbolas2; //diffusion and decay of T, mobility, surface tension of phi

static const double D1A=70, alphaA=D1A/Rbolas2, Acut=-0.10; //Acut=0.05; Para a Angiopoietina q nao esta a funcionar!

//-------BURACOS
static const double raio_hipoxia=4.0;
static double buracos[Lx][Ly];
//-------BURACOS


static const double /*RMMP=2.0,*/ Rtip=Rbolas;
static const double D12=D1;

static const double ruido=0.0 /*D2=500,alpha2=3*/, Tfonte2=0.0, Afonte=0.0; //noise, source vegf1, source vegf2, source angiopoietin

static const double RMMP=10, Trep=0.3; // maximum velocity of CA, maximum T for reproduction
static double vmax;

static const int xinc=2, yinc=2;

static double alpha2, D2;
static int Ralimento;
static double Tfonte,Tcut;
static double vv[NbolasMax], gT[Lx][Ly][2], R[Lx][Ly], Bolas[NbolasMax][2], mu[Lx][Ly]; //maxR;
static int Nbolas, Nfontes, Fontes[10000][3], tipo[Lx][Ly];
static int xo, yo, dx, dy;
static char *fora;
static int nm=0;//numero de tipcells mortas a escrever em 'morte'
int mor[Lx*Ly][2];
double a[Lx][Ly];
//extern  int ia[Lx][Ly][2];

static double T[Lx][Ly],an[Lx][Ly],Tn[Lx][Ly],T2[Lx][Ly],T2n[Lx][Ly],dfont[Lx][Ly],ANG[Lx][Ly],ANGn[Lx][Ly];
static double **pa, **pan, **pT, **pTn,**pT2,**pT2n,**pANG,**pANGn;
static double *piT[Lx];
static double *piTn[Lx];
static double *pia[Lx];
static double *pian[Lx];
static double *piT2[Lx];
static double *piT2n[Lx];
static double *piANG[Lx];
static double *piANGn[Lx];

static double t;
static int lastt; // Flag for the opening of proteins.

static int fim;

ofstream gradientes;

static void ini(); //establece campos iniciais
static void step(); //equivalente ao evolution
static void out (int tt); //output
static int boundx (int xx); //condicoes de fronteira em X
static int boundy (int yy); //condicoes de fronteira em y
extern void flux(); //considera a passagem de sangue nos capilares

int main (int argc, char * argv[]){

  ofstream morte;
  char Mout[30];
  char outbia[100];

  int i,km;
  //command line arguments:
  fora=argv[1];
  //vmax = atof(argv[2]);
  //vmax/=100.;
    vmax=0.03;
  D2=atof(argv[2]); // chemotaxis ->250
  //alpha2=atof(argv[4]); // reproduction rate ->14, que de seguida é dividida por 10.
  //alpha2/=10.0;
  //Tfonte=atoi(argv[5]);
  //Tfonte=Tfonte/10.0;  //HUGO: Changed from 100 to 10...
    alpha2=0;
    Tfonte=1;
    Ralimento=20;
  Tcut=Tfonte*exp(-1)*0.8/(1+Ralimento/Rbolas); //expressao do artigo para o Tc; O Tfonte equivale ao Ts no artigo
  //Tcut = 0.055;  //HUGO: valor aproximado com Tfonte = 1
    Tcut=0.01;
  xo=dx=yo=dy=100;
  gradientes.open(fora);

  ini();

  t=0;
  fim=0;

  int min=2000;
  for(t=0;t<tmax;t+=dt){
    //if (((int)(t/dt))%10==0)flux();
    step();
    if(((int)(t/dt))%Nprint==0)out((int)(t/dt));
    if(((int)(t/dt))%100==0){
      for(km=0;km<Nbolas;km++){
	gradientes<<t<<"\t"<<sqrt(gT[(int)Bolas[km][0]][(int)Bolas[km][1]][0]*gT[(int)Bolas[km][0]][(int)Bolas[km][1]][0]+gT[(int)Bolas[km][0]][(int)Bolas[km][1]][1]*gT[(int)Bolas[km][0]][(int)Bolas[km][1]][1])<<"\t"<<km<<"\n";
      }
    }
    if((int)(t/dt)>=min){
      sprintf(outbia,"./newprog A.%s.%d %d %d %d %d %d",fora,(int)(t/dt), Lx, Ly, 10, centreX, centreY);
      printf("%s\n", outbia);
      system(outbia);
      min+=3000;
      lastt=(int)(t/dt);
    }
    if(fim>0)fim++;
    if(fim>50000)break;


  }

  sprintf(Mout,"M%s",fora);
  morte.open(Mout);
    morte<<nm<<"\n";
  for(i=0;i<nm;i++){
    morte<<mor[i][0]<<"\t"<<mor[i][1]<<"\n";
  }
  morte.close();

}


int notchkill(int x1,int y1, int x2,int y2){
  int i, j, n, ff, fcap, jlow, jhigh, ncount, res, jj[10][2], dj;
  double dd1, dd2, dd3, dd, ddd, tx, ty, xi, xf, yi, yf, xmed, ymed;
  double vx, vy, px, py;
  double ends[100][3], endsnew[100][3];

  ofstream fora("notas.dat");

  cout<< x1<< " "<< y1<<" "<<x2<<" "<<y2<<"\n"; //buscaerro

  dd1=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
  dd2=(x2-x1)*(x2-x1)+(y2-y1+Ly)*(y2-y1+Ly);
  dd3=(x2-x1)*(x2-x1)+(y2-y1-Ly)*(y2-y1-Ly);
    dd2=Ly;dd3=Ly;//BURACO
  if(dd1<dd2){
    if(dd1<dd3){
      dd=dd1;
    }
    else{
      dd=dd3;
      y2=y2-Ly;
    }
  }
  else{
    if(dd2<dd3){
      dd=dd2;
      y2=y2+Ly;
    }
    else{
      dd=dd3;
      y2=y2-Ly;
    }
  }

  dd=sqrt(dd); //dd é a raiz quadrada da distancia menor

  n=2;
  ends[0][0]=x1;
  ends[0][1]=y1;
  ends[1][0]=x2;
  ends[1][1]=y2;
  ends[0][2]=dd;

  ff=1;
  while(ff==1 && dd<4*Rbolas){ //enquanto estivermos dentro do capilar e a distancia for menor que 4raios

    ff=1;
    ncount=0;
    endsnew[0][0]=ends[0][0];
    endsnew[0][1]=ends[0][1];

    res=0;
    for(i=0;i<n-1;i++){
      vx=(ends[i+1][0]-ends[i][0])/ends[i][2];
      vy=(ends[i+1][1]-ends[i][1])/ends[i][2];
      fcap=1;
      jlow=0;
      for(j=1;j<ends[i][2];j++){
	tx=ends[i][0]+j*vx; //versor
	ty=ends[i][1]+j*vy;
	if(fcap==1 && pa[(int)tx][boundy((int)ty)]<0){
	  fcap=0;
	  jj[res][0]=j;
	  res++;
	}
	else if (fcap==0 && pa[(int)tx][boundy((int)ty)]>0){
	  fcap=1;
	  jj[res-1][1]=j;
	}
      }

      dj=0;
      if(res>0){
	for(j=0;j<(fcap==0?res-1:res);j++){
	  if(jj[j][1]-jj[j][0]>dj){
	    dj=jj[j][1]-jj[j][0];
	    jlow=jj[j][0];
	    jhigh=jj[j][1];
	  }
	}
	cout << jlow<< " "<<jhigh <<"\n";
	xi=ends[i][0]+jlow*vx;
	yi=ends[i][1]+jlow*vy;
	xf=ends[i][0]+jhigh*vx;
	yf=ends[i][1]+jhigh*vy;

	ddd=(jhigh+jlow)/2;
	xmed=ends[i][0]+ddd*vx;
	ymed=ends[i][1]+ddd*vy;
	px=-vy;
	py=vx;

	j=0;
	fcap=0;

	while(fcap==0){
	  j++;
	  tx=xmed+j*px;
	  ty=ymed+j*py;
	  if( pa[(int)tx][boundy((int)ty)]>0){
	    xmed=tx;
	    ymed=ty;
	    fcap=1;
	  }
	  else{
	    tx=xmed-j*px;
	    ty=ymed-j*py;
	    if( pa[(int)tx][boundy((int)ty)]>0){
	      xmed=tx;
	      ymed=ty;
	      fcap=1;
	    }
	  }
	}

	cout<< xmed<< " "<<ymed<<"\n";

	xi=ends[i][0];yi=ends[i][1];
	xf=ends[i+1][0];yf=ends[i+1][1];
	endsnew[ncount][2]=sqrt( (xmed-xi)*(xmed-xi)+(ymed-yi)*(ymed-yi) );
	ncount++;
	endsnew[ncount][0]=xmed;
	endsnew[ncount][1]=ymed;
	endsnew[ncount][2]=sqrt( (xf-xmed)*(xf-xmed)+(yf-ymed)*(yf-ymed) );
	ncount++;
	endsnew[ncount][0]=ends[i+1][0];
	endsnew[ncount][1]=ends[i+1][1];
      }
      else{
	endsnew[ncount][2]=ends[i][2];
	ncount++;
	endsnew[ncount][0]=ends[i+1][0];
	endsnew[ncount][1]=ends[i+1][1];
      }
    }

    if(ncount+1==n)ff=0;
    else{
      n=ncount+1;
      for(i=0;i<n;i++){
	ends[i][0]=endsnew[i][0];
	ends[i][1]=endsnew[i][1];
	ends[i][2]=endsnew[i][2];
	//cout<< "("<<	ends[i][0]<<", "<<	ends[i][1]<<") "<<	ends[i][2]<<"\n";
	//fora<< "("<<	ends[i][0]<<", "<<	ends[i][1]<<") "<<	ends[i][2]<<"\n";
      }
      dd=0;
      for(i=0;i<n-1;i++){
	dd+=ends[i][2];
      }
    }
  }

  res=dd<4*Rbolas?1:0; //se res=1 uma das bolas tem que morrer, porque nao pode ser uma tip cell (notch)
  return(res);
}


void step(){

  int i,j,k,m,n,x1,x2,y1,y2,ff,fm,n_criar,i_criar,j_criar,kmin,iff,id;
  double Tloc,KMM,aloc,potencial,divergencia,a_criar,dd,dd1,dd2,dd3,ddmin,D12loc,taxarep;
  double **temp;

  double nabla(double x[Lx][Ly],int i1, int i2);
  double nablabound(double **x,int i1, int i2);
  void ballmove();
  void criabola(int i10,int j10);
  void matabola(int k10);
  int notchkill(int i22,int j22,int x122,int y122);
  double prol(double trep11, double aloc11);

  //Sources disappear if are fed by capillary
  k=0;
  id=0;
  while(k<Nfontes){
    i=Fontes[k][0];j=Fontes[k][1];ff=0;iff=Fontes[k][2];
    for(m=i-Ralimento;m<i+Ralimento;m++){
      for(n=j-Ralimento;n<j+Ralimento;n++){
	if(pa[boundx(m)][boundy(n)]>0/*ia[boundx(m)][boundy(n)][0]==1 */&& (m-i)*(m-i)+(n-j)*(n-j)<Ralimento*Ralimento ){ ////NOTA NOTA: ANASTOMOSE
	  Fontes[k][2]=0;
	  ff=1;
	  goto fora;
	}
      }
    }
  // Desligar as fontes
  fora: if(ff==0)Fontes[k][2]=0;//Fontes[k][2]=1;
    //if(Fontes[k][2]!=iff){id=1;cout << "Fonte: "<< k<<"\n";}
    k++;
  }

  if(id==1){
    cout << "nova difusao\n";
    for(i=0;i<Lx;i++){
      for(j=0;j<Ly;j++){
	dfont[i][j]=Lx;
	for(k=0;k<Nfontes;k++){
	  if(Fontes[k][2]==1){
	    dd1=(i-Fontes[k][0])*(i-Fontes[k][0])+(j-Fontes[k][1])*(j-Fontes[k][1]);
	    dd2=(i-Fontes[k][0])*(i-Fontes[k][0])+(j-Fontes[k][1]+Ly)*(j-Fontes[k][1]+Ly);
	    dd3=(i-Fontes[k][0])*(i-Fontes[k][0])+(j-Fontes[k][1]-Ly)*(j-Fontes[k][1]-Ly);
          dd2=Ly;dd3=Ly;//BURACO
          dd=sqrt(dd1<dd2?(dd1<dd3?dd1:dd3):(dd2<dd3?dd2:dd3));
	    dfont[i][j]=dd<dfont[i][j]?dd:dfont[i][j];
	  }
	}
      }
    }
  }

  //Chemical potential
  for(i=0;i<Lx;i++){
    for(j=0;j<Ly;j++){
      aloc=pa[i][j];
      potencial=-aloc+aloc*aloc*aloc;
      divergencia=nablabound(pa,i,j);
      mu[i][j]=potencial-epsilon *divergencia + 3*(pa[i][j]*pa[i][j]-1)*(buracos[i][j]-2)*pow((1+buracos[i][j]),2)/16; //BURACO!!!
    }
  }

  //Gradient T
  for(j=0;j<Ly;j++){
    for(i=0;i<Lx;i++){
      aloc=pa[i][j];
      x1= boundx(i-1);
      x2= boundx(i+1);
      y1= boundy(j-1);
      y2= boundy(j+1);
      gT[i][j][0]=(pT[x2][j]+pT2[x2][j]-pT[x1][j]-pT2[x1][j])/2.;
      gT[i][j][1]=(pT[i][y2]+pT2[i][y2]-pT[i][y1]-pT2[i][y1])/2.;
    }
  }


  //CA creation
  n_criar=0; a_criar=0;
  for(j=2*Rbolas;j<Ly-2*Rbolas;j++){  //BURACO
    for(i=1;i<Lx;i++){

      if(pa[i][j]>0.9 && pANG[i][j]>Acut && pT[i][j]+pT2[i][j]-0.5*pANG[i][j] > Tcut && gT[i][j][0]*gT[i][j][0]+gT[i][j][1]*gT[i][j][1]>gradTmin*gradTmin){
	ff=1;
	for(k=0;k<Nbolas;k++){           //Notch
	  dd1=(i-Bolas[k][0])*(i-Bolas[k][0])+(j-Bolas[k][1])*(j-Bolas[k][1]);
	  dd2=(i-Bolas[k][0])*(i-Bolas[k][0])+(j-Bolas[k][1]+Ly)*(j-Bolas[k][1]+Ly);
	  dd3=(i-Bolas[k][0])*(i-Bolas[k][0])+(j-Bolas[k][1]-Ly)*(j-Bolas[k][1]-Ly);
        dd2=Ly;dd3=Ly;//BURACO
        dd=dd1<dd2?(dd1<dd3?dd1:dd3):(dd2<dd3?dd2:dd3); //dd toma o valor do menor dd1 dd2 ou dd3
	  if(dd<16*Rbolas2)ff=0; //se a distancia for menor que 4 raios vem ff=0
	}
	if(ff==1){                       //New CA is completely inside existing capillary
	  for(m=i-Rbolas; m<i+Rbolas;m++){
	    for(n=j-Rbolas;n<j+Rbolas;n++){
	      if((m-i)*(m-i)+(n-j)*(n-j)<Rbolas2 && pa[boundx(m)][boundy(n)]<0){
		ff=0;
		goto ooout;
	      }
	    }
	  }
	ooout: if(ff==1){
	    n_criar=1;
	    if(pT[i][j]>a_criar){////CHANGED!!!
	      i_criar=i;j_criar=j;a_criar=pT[i][j];////CHANGED!!!
	    }
	  }
	}
      }
    }
  }

  if(n_criar==1){
    criabola(i_criar,j_criar);  cout <<"cria\n";
  }


  ballmove();

  //CA killing
  for(k=0;k<Nbolas;k++){
    i=Bolas[k][0];j=Bolas[k][1];
    if( pT[i][j]+pT2[i][j]-0.5*pANG[i][j] < Tcut || gT[i][j][0]*gT[i][j][0]+gT[i][j][1]*gT[i][j][1]<gradTmin*gradTmin){
      matabola(k);
    }
  }

  /*
  //CA killing Notch
  m=-1;
  while(m<Nbolas-1){
    m++;
    i=Bolas[m][0];j=Bolas[m][1];
    for(k=0;k<Nbolas;k++){           //Notch
      if(k!=m){
	x1=Bolas[k][0]; y1=Bolas[k][1];
	dd1=(i-x1)*(i-x1)+(j-y1)*(j-y1);
	dd2=(i-x1)*(i-x1)+(j-y1+Ly)*(j-y1+Ly);
	dd3=(i-x1)*(i-x1)+(j-y1-Ly)*(j-y1-Ly);
	dd=dd1<dd2?(dd1<dd3?dd1:dd3):(dd2<dd3?dd2:dd3);
	if(dd<4*Rbolas2){
	  if(pT[i][j]+pT2[i][j]-0.5*pANG[i][j]<pT[x1][y1]+pT2[x1][y1]-0.5*pANG[x1][y1]) matabola(m);
	  else matabola(k);
	  m=-1;
	  break;
	}
	else if (dd<16*Rbolas2){
	  //cout<<dd<<"\n";
	  //ff=notchkill(i,j,x1,y1);
	  //cout<<ff<<"\n";
	  ff=0;
	  if(ff==1){
	    if(pT[i][j]+pT2[i][j]-0.5*pANG[i][j]<pT[x1][y1]+pT2[x1][y1]-0.5*pANG[x1][y1]) matabola(m);
	    else matabola(k);
	    m=-1;
	    break;
	  }
	}
      }
    }
  }
  */

  //New values for phi and T
  for(j=0;j<Ly;j++){
    for(i=1;i<Lx;i++){
       

      ff=0; //inside CA
      fm=0;

      ddmin=Lx*Lx;kmin=Nbolas;
      for(k=0;k<Nbolas;k++){
	dd1=(i-Bolas[k][0])*(i-Bolas[k][0])+(j-Bolas[k][1])*(j-Bolas[k][1]);
	dd2=(i-Bolas[k][0])*(i-Bolas[k][0])+(j-Bolas[k][1]+Ly)*(j-Bolas[k][1]+Ly);
	dd3=(i-Bolas[k][0])*(i-Bolas[k][0])+(j-Bolas[k][1]-Ly)*(j-Bolas[k][1]-Ly);
          dd2=Ly;dd3=Ly;//BURACO
    dd=(dd1<dd2?(dd1<dd3?dd1:dd3):(dd2<dd3?dd2:dd3));
	ddmin=dd<ddmin?dd:ddmin;
	kmin=dd<ddmin?k:kmin;
      }

      //D12loc=D12*(pow((RMMP/(RMMP+dfont[i][j])),3)+pow((Rtip/(Rtip+ddmin)),3))/2.0;
      KMM=10*exp(-dfont[i][j]/RMMP);
      D12loc=D12*(1+(KMM-C+pT2[i][j])/sqrt( pow((KMM+C-pT2[i][j]),2)+4*KMM*pT2[i][j] ))/2.0;

      if (ff==0 && tipo[i][j]!=1 && ddmin<Rbolas2){
	ff=1;
	taxarep=pT[i][j]+pT2[i][j]-0.5*pANG[i][j];
	//pan[i][j]=prol(taxarep,1.0) *3.141592*Rbolas/D2/4.0/(vv[kmin]<gradTmin?gradTmin:vv[kmin])/(1.0+ruido*R[i][j]);
          pan[i][j]=1;
	pTn[i][j]=pT[i][j]+dt*(D1*nablabound(pT,i,j));
	pT2n[i][j]=pT2[i][j]+dt*(D1*nablabound(pT2,i,j));
	pANGn[i][j]=pANG[i][j]+dt*(D1A*nablabound(pANG,i,j));
      }

      //outside CA
      if(ff==0){
	if(pa[i][j]>=0)tipo[i][j]=1;
	else tipo[i][j]=0;

	Tloc=pT[i][j];
	aloc=pa[i][j];

	pTn[i][j]=Tloc+dt*(D1*nablabound(pT,i,j)-alpha1*(aloc>0?aloc:0)*Tloc);
	pTn[i][j]=pTn[i][j]<0?0:pTn[i][j];
	taxarep=Tloc+pT2[i][j]-0.5*pANG[i][j];
	pan[i][j]=aloc+dt*(M*nabla(mu,i,j)+(pANG[i][j]>Acut? prol(taxarep,aloc)  : 0 ) );
	pT2n[i][j]=pT2[i][j]+dt*(D12loc*nablabound(pT2,i,j)-alpha1*(aloc>0?aloc:0)*pT2[i][j]);
	pT2n[i][j]=pT2n[i][j]<0?0:pT2n[i][j];
	pANGn[i][j]=pANG[i][j]+dt*(D1A*nablabound(pANG,i,j)-alphaA*(aloc>0?aloc:0)*pANG[i][j]);
	pANGn[i][j]=pANGn[i][j]<0?0:pANGn[i][j];
      }

      for(k=0;k<Nfontes;k++){
	if(i==Fontes[k][0] && j==Fontes[k][1]){
	  if(Fontes[k][2]!=0){
	    pTn[i][j]=Tfonte;pT2n[i][j]=Tfonte2;ANGn[i][j]=Afonte;
	  }
	  else{
	    if(t>Tang && i>xo && i<xo+dx && j>yo && j<yo+dy) ANGn[i][j]=Afonte;
	  }
	}
      }

        //BURACO
        /*
        if(buracos[i][j]==2) pan[i][j]=-1;
        if(buracos[i][j]==1){
            n=0;
            for(k=0;k<4;k++){
                x1=k%2==0?0:k-2;
                y1=(k+1)%2==0?0:k-1;
                pan[i][j]=0;
                if(buracos[boundx(i+x1)][boundy(j+y1)]==0){
                    pan[i][j]+=pa[boundx(i+x1)][boundy(j+y1)];
                    n++;
                }
                pan[i][j]/=n;
            }
        }*/
        //BURACO
        
        
    }
    pTn[0][j]=pT[0][j];//pTn[Lx-1][j]=pT[Lx-1][j];
    pan[0][j]=pa[0][j];//pan[Lx-1][j]=pa[Lx-1][j];
    pT2n[0][j]=pT2[0][j];
    pANGn[0][j]=pANG[0][j];
  }

  temp=pT2;
  pT2=pT2n;
  pT2n=temp;
  temp=pT;
  pT=pTn;
  pTn=temp;
  temp=pa;
  pa=pan;
  pan=temp;
  temp=pANG;
  pANG=pANGn;
  pANGn=temp;
}


void criabola(int i,int j){
  if (Nbolas<1){
  Bolas[Nbolas][0]=i;
  Bolas[Nbolas][1]=j;
  Nbolas++;
  
    if(fim>0 &&Nbolas>0)fim=0;
  }
}


void matabola(int k){
  int i,j,ve;//variavel de estado
  int r=Bolas[k][0];
  int c=Bolas[k][1];

  ve=1;
  for(i=0;i<nm;i++){
    if(nm>0 && (pow(mor[i][0]-r,2) + pow(mor[i][1]-c,2)<Rbolas2)){
      ve=0;
    }
  }
  if(ve==1){
    mor[nm][0]=r;
    mor[nm][1]=c;
    nm++;
  }
  for(j=k;j<Nbolas-1;j++){
    Bolas[j][0]=Bolas[j+1][0];
    Bolas[j][1]=Bolas[j+1][1];
  }
  Nbolas--;
  if(Nbolas==0)fim=1;
  cout <<"mata\n";
}


void ballmove(){
  int i,j,k,x1,x2,y1,y2,cflag,valy,dd1,dd2,dd3;
  double x,y,g0,gx,gy,gfx,gfy,vx,vy,vrx,vry,v,ux,uy,udotv;

    ofstream mtubout;
    mtubout.open("mtube", std::ios_base::app);

    ofstream velFile;
    velFile.open("velocityTime", std::ios_base::app);
    ofstream flagVel;
    flagVel.open("flagVel", std::ios_base::app);
    
  for(i=0;i<Nbolas;i++){
    x=Bolas[i][0];y=Bolas[i][1];

    x1=(int)(x+0.001);y1=(int)(y+0.001);

    g0=gT[x1][y1][0];
    gx=gT[boundx(x1+1)][y1][0];
    gy=gT[x1][boundy(y1+1)][0];
    gfx=g0+(x-x1)*(gx-g0)+(y-y1)*(gy-g0);

    g0=gT[x1][y1][1];
    gx=gT[boundx(x1+1)][y1][1];
    gy=gT[x1][boundy(y1+1)][1];
    gfy=g0+(x-x1)*(gx-g0)+(y-y1)*(gy-g0);

    vx=gfx; vy=gfy;

      
      //BURACO
      for(k=0;k<Nfontes;k++){
          dd1=(Fontes[k][0]-x)*(Fontes[k][0]-x)+(Fontes[k][1]-y)*(Fontes[k][1]-y);
          dd2=(Fontes[k][0]-x)*(Fontes[k][0]-x)+(Fontes[k][1]-y+Ly)*(Fontes[k][1]-y+Ly);
          dd3=(Fontes[k][0]-x)*(Fontes[k][0]-x)+(Fontes[k][1]-y-Ly)*(Fontes[k][1]-y-Ly);
          
          dd2=Ly;dd3=Ly;//BURACO
          
          
          if(dd1<dd2){
              if(dd1<dd3){
                  dd1=dd1;
                  valy=0;
              }
              else{
                  dd1=dd3;
                  valy=-1;
              }
          }
          else{
              if(dd2<dd3){
                  dd1=dd2;
                  valy=1;
              }
              else{
                  dd1=dd3;
                  valy=-1;
              }
          }
          
          if(dd1<(raio_hipoxia+Rbolas)*(raio_hipoxia+Rbolas)){
              ux=(x-Fontes[k][0])/sqrt( (x-Fontes[k][0])*(x-Fontes[k][0])+ (y+valy*Ly-Fontes[k][1])*(y+valy*Ly-Fontes[k][1]));
              uy=(y+valy*Ly-Fontes[k][1])/sqrt( (x-Fontes[k][0])*(x-Fontes[k][0])+ (y+valy*Ly-Fontes[k][1])*(y+valy*Ly-Fontes[k][1]));
              udotv=vx*ux+vy*uy;
              vx=vx-udotv*ux;
              vy=vy-udotv*uy;
              break;
          }
      }
      //BURACO
          
      

    v=sqrt(vx*vx+vy*vy);
    vv[i]=(v<gradTmin?0:v);
    
    // <Código Bia>
   /*
    char name[255];
    sprintf(name, "A.%s.%d_proteins_5",fora,lastt);
    double Dfinal=D2;
    string line;
    ifstream proteinstemp;
    proteinstemp.open(name); 
    int col=0, row=0;
    double temporary;
    //flagVel<<name<<"\n";

    if(!proteinstemp.fail()){
      flagVel<<Dfinal<<" entrou "<<t<<" "<<lastt<<endl;
      while(proteinstemp.good()){
        while(getline(proteinstemp,line)){
          istringstream streamA(line);
          col=0;
          while(streamA>>temporary){
            col++;
          }
          row++;
        }
      }
      ifstream proteins;
      proteins.open (name); 
      // Ler para uma matriz
      int px,py;
      double proteinas[col][row];
      for(px =0; px <col; px++){
        for(py = 0; py <row; py++){
          proteins>>proteinas[px][py];
        }
      }


      int proti,protj,iteration;
      double meanConc;
      meanConc = 0;
      iteration = 0;
      for (proti=x- Rbolas; proti<x+Rbolas; proti++){
        for (protj=y- Rbolas; protj<y+Rbolas; protj++){
          //if(((proti*proti)+(protj*protj))<Rbolas*Rbolas){
            meanConc += proteinas[proti][protj];
            iteration++;
          //}
        }
      }
      if (iteration > 0){
        Dfinal = meanConc/iteration;
      }
      Dfinal = Dfinal*D2;
      flagVel<<Dfinal<<" saiu "<<t<<endl;
    }*/
    // </Código Bia>
    double Dfinal=D2; // Acrescentei isto aqui para tirar a parte das proteinas deste pedaco

    Bolas[i][0]=x+dt*Dfinal*(R[x1][y1]*ruido+1)*(v<vmax?vx:vmax*vx/v);
    Bolas[i][1]=y+dt*Dfinal*(R[x1][y1]*ruido+1)*(v<vmax?vy:vmax*vy/v);

      mtubout<<Bolas[i][0]<<" "<<Bolas[i][1]<<endl;
      velFile<<Dfinal<<" "<<t<<endl;
      
      
    //Bolas[i][1]+=Bolas[i][1]<0?Ly:Bolas[i][1]>Ly?-Ly:0;
    //Bolas[i][0]=Bolas[i][0]>Lx-2?Lx-2:Bolas[i][0]<1?1:Bolas[i][0];
      //if(Bolas[i][0]>Lx-Rbolas-1)Bolas[i][0]=Lx-Rbolas-1;
      
      //BURACO
    Bolas[i][1]=Bolas[i][1]>Ly-Rbolas-1?Ly-Rbolas-1:Bolas[i][1]<Rbolas?Rbolas:Bolas[i][1];
      Bolas[i][0]=Bolas[i][0]>Lx-Rbolas-1?Lx-Rbolas-1:Bolas[i][0]<2?2:Bolas[i][0];
//BURACO
      
  }
    mtubout.close();
    velFile.close();
}


double nabla(double u[Lx][Ly], int x, int y){
  int x1,x2;
  int y1,y2;
  double lapl;

  x1= boundx(x-1);
  x2= boundx(x+1);
  y1= boundy(y-1);
  y2= boundy(y+1);

  lapl=-4*u[x][y]+u[x1][y]+u[x2][y]+u[x][y1]+u[x][y2];

  return(lapl);
}

double nablabound(double **u,int x, int y){//Laplaciano para pointer!!
  int x1,x2;
  int y1,y2;
  double lapl;
  x1= boundx(x-1);
  x2= boundx(x+1);
  y1= boundy(y-1);
  y2= boundy(y+1);

  lapl=-4*u[x][y]+u[x1][y]+u[x2][y]+u[x][y1]+u[x][y2];

  return(lapl);

}


void ini(){
  int i,j,k,m,n,dd1,dd2,dd3;
  double dd;

    ofstream outbur("buracos");
  srand(3247688);


  //initial  values for phi and T
  for(i=0;i<Lx;i++){
    for(j=0;j<Ly;j++){
      Tn[i][j]=0;
      an[i][j]=-1;
      //T[i][j]=0;
      T[i][j]=5*i/Lx;
      a[i][j]=-1;
      T2[i][j]=0;
      T2n[i][j]=0;
      ANG[i][j]=0;
      ANGn[i][j]=0;
      //ia[i][j][0]=0;
      //if(i==0)ia[i][j][0]=1;
    }
  }


  for(i=0;i<Lx;i++){
    for(j=0;j<Ly;j++){
      R[i][j]=0;
    }
  }

  //initial souces of vegf
  Nfontes=0;
  for(i=0;i<Lx;i++){
    for(j=0;j<Ly;j++){
      if( (double(rand())/RAND_MAX) /*drand48()*/<(2*(i>10?10:i)/float(150)/75)){
	Fontes[Nfontes][0]=i; Fontes[Nfontes][1]=j; Fontes[Nfontes][2]=0;
	Nfontes++;
	cout << i<<"\t"<<j<<"\t"<<Nfontes<<"\n";
      }
    }
  }


  //noise.close();

  //initial cell
  Nbolas=0;

  for(i=0;i<Lx;i++){
    for(j=0;j<Ly;j++){
      //if((i-Lx/2)*(i-Lx/2)+(j-Ly/2)*(j-Ly/2)<Arteria*Arteria){a[i][j]=1.;an[i][j]=1.;}
      if((i-centreX)*(i-centreX)+(j-centreY)*(j-centreY)<Arteria*Arteria){a[i][j]=1.;an[i][j]=1.;}
      //if(i>Ly-Arteria){a[i][j]=1.;an[i][j]=1.;}
    }
  }

  for (i=0;i<Nfontes;i++){
    cout << i<<"\t"<<Nfontes<<"\t"<< Fontes[i][0]<<"\t"<<Fontes[i][1]<<"\t"<< Fontes[i][2]<<"\n";
  }


  //sources disappear if there is food nearby
  k=0;
  while(k<Nfontes){
    i=Fontes[k][0];j=Fontes[k][1];
    for(m=i-Ralimento;m<i+Ralimento;m++){
      for(n=j-Ralimento;n<j+Ralimento;n++){
	if(a[boundx(m)][boundy(n)]==1){
        Fontes[k][0]=Fontes[Nfontes-1][0];
        Fontes[k][1]=Fontes[Nfontes-1][1];
        Fontes[k][2]=Fontes[Nfontes-1][2];
        Nfontes--;
        k--;
	  //cout << i<<"\t"<<j<<"\t"<<k<<"\t"<<Nfontes<<"\n";
	  goto foraini;
	}
      }
    }
  foraini: k++;
  }

    
    
    //-------BURACOS
    k=0;
    while(k<Nfontes){
        if(a[Fontes[k][0]][Fontes[k][1]]>-0.2){
            Fontes[k][0]=Fontes[Nfontes-1][0];
            Fontes[k][1]=Fontes[Nfontes-1][1];
            Fontes[k][2]=Fontes[Nfontes-1][2];
            Nfontes--;
            k--;
        }
        else{
        for(i=k+1;i<Nfontes;i++){
            dd1=(Fontes[k][0]-Fontes[i][0])*(Fontes[k][0]-Fontes[i][0])+(Fontes[k][1]-Fontes[i][1])*(Fontes[k][1]-Fontes[i][1]);
            dd2=(Fontes[k][0]-Fontes[i][0])*(Fontes[k][0]-Fontes[i][0])+(Fontes[k][1]-Fontes[i][1]+Ly)*(Fontes[k][1]-Fontes[i][1]+Ly);
            dd3=(Fontes[k][0]-Fontes[i][0])*(Fontes[k][0]-Fontes[i][0])+(Fontes[k][1]-Fontes[i][1]-Ly)*(Fontes[k][1]-Fontes[i][1]-Ly);
            dd2=Ly;dd3=Ly;//BURACO
            dd1=dd1<dd2?(dd1<dd3?dd1:dd3):(dd2<dd3?dd2:dd3);

            if(dd1<16*raio_hipoxia*raio_hipoxia){
                Fontes[k][0]=Fontes[Nfontes-1][0];
                Fontes[k][1]=Fontes[Nfontes-1][1];
                Fontes[k][2]=Fontes[Nfontes-1][2];
                Nfontes--;
                k--;
                break;
            }
        }
        }
        k++;
    }
    
    for(i=0;i<Lx;i++){
        for(j=0;j<Ly;j++){
            buracos[i][j]=-1;
            for(k=0;k<Nfontes;k++){
                dd1=(Fontes[k][0]-i)*(Fontes[k][0]-i)+(Fontes[k][1]-j)*(Fontes[k][1]-j);
                dd2=(Fontes[k][0]-i)*(Fontes[k][0]-i)+(Fontes[k][1]-j+Ly)*(Fontes[k][1]-j+Ly);
                dd3=(Fontes[k][0]-i)*(Fontes[k][0]-i)+(Fontes[k][1]-j-Ly)*(Fontes[k][1]-j-Ly);
                dd2=Ly;dd3=Ly;//BURACO
                dd1=dd1<dd2?(dd1<dd3?dd1:dd3):(dd2<dd3?dd2:dd3);
                if(dd1<=raio_hipoxia*raio_hipoxia){
                    buracos[i][j]=1;
                    break;
                }
            }
        }
    }
    /*
    for(i=0;i<Lx;i++){
        for(j=0;j<Ly;j++){
            if(buracos[i][j]==2 && ( buracos[boundx(i+1)][j]==0 || buracos[boundx(i-1)][j]==0 || buracos[i][boundy(j+1)]==0 || buracos[i][boundy(j-1)]==0) ) buracos[i][j]=1;
        }
    }*/
    
    for(k=0;k<15;k++){
        
        for(i=0;i<Lx;i++){
            for(j=0;j<Ly;j++){
                mu[i][j]=-buracos[i][j]+pow(buracos[i][j],3)-(buracos[boundx(i+1)][j]+buracos[boundx(i-1)][j]+buracos[i][boundy(j+1)]+buracos[i][boundy(j-1)]-4*buracos[i][j]);
            }
        }
        
        
        for(i=0;i<Lx;i++){
            for(j=0;j<Ly;j++){
                buracos[i][j]=buracos[i][j]+0.01*(mu[boundx(i+1)][j]+mu[boundx(i-1)][j]+mu[i][boundy(j+1)]+mu[i][boundy(j-1)]-4*mu[i][j]);
            }
        }
    
    }
    //-------BURACOS
    
    
    
    

  for (i=0;i<Nfontes;i++){
    cout << i<<"\t"<< Nfontes<<"\t"<<Fontes[i][0]<<"\t"<<Fontes[i][1]<<"\t"<< Fontes[i][2]<<"\n";
  }

  //T close to sources
  k=0;
  while(k<Nfontes){
    i=Fontes[k][0];j=Fontes[k][1];
    if(Fontes[k][2]!=0){
      for(m=i-10;m<i+10;m++){
	for(n=j-10;n<j+10;n++){
	  //T[boundx(m)][boundy(n)]=Tfonte*exp(-((i-m)*(i-m)+(j-n)*(j-n)));
	  T2[boundx(m)][boundy(n)]=Tfonte2*exp(-((i-m)*(i-m)+(j-n)*(j-n)));
	  ANG[boundx(m)][boundy(n)]=Afonte*exp(-((i-m)*(i-m)+(j-n)*(j-n)));
	}
      }
    }
    k++;
  }

  for(i=0;i<Lx;i++){
    for(j=0;j<Ly;j++){
      dfont[i][j]=Lx;
      for(k=0;k<Nfontes;k++){
	if(Fontes[k][2]==1){
	  dd1=(i-Fontes[k][0])*(i-Fontes[k][0])+(j-Fontes[k][1])*(j-Fontes[k][1]);
	  dd2=(i-Fontes[k][0])*(i-Fontes[k][0])+(j-Fontes[k][1]+Ly)*(j-Fontes[k][1]+Ly);
	  dd3=(i-Fontes[k][0])*(i-Fontes[k][0])+(j-Fontes[k][1]-Ly)*(j-Fontes[k][1]-Ly);
        dd2=Ly;dd3=Ly;//BURACO
      dd=sqrt(dd1<dd2?(dd1<dd3?dd1:dd3):(dd2<dd3?dd2:dd3));
	  dfont[i][j]=dd<dfont[i][j]?dd:dfont[i][j];
	}
      }
    }
  }

  for(i=0;i<Lx;i++){
    for(j=0;j<Ly;j++){
      if(a[i][j]<0)tipo[i][j]=0;
      else{ tipo[i][j]=1;}
    }
  }


    for(i=0;i<Lx;i++){
        for(j=0;j<Ly;j++){
            if(buracos[i][j]>0){
                outbur<<i<<" "<<j<<" "<<buracos[i][j]<<"\n";
            }
        }
    }

    
    
  //Inicializacao dos pointers
  for(i=0;i<Lx;i++){
    piT[i]=&T[i][0];
  }
  pT=&piT[0];

  for(i=0;i<Lx;i++){
    piTn[i]=&Tn[i][0];
  }
  pTn=&piTn[0];

  for(i=0;i<Lx;i++){
    pia[i]=&a[i][0];
  }
  pa=&pia[0];

  for(i=0;i<Lx;i++){
    pian[i]=&an[i][0];
  }
  pan=&pian[0];

  for(i=0;i<Lx;i++){
    piANG[i]=&ANG[i][0];
  }
  pANG=&piANG[0];

  for(i=0;i<Lx;i++){
    piANGn[i]=&ANGn[i][0];
  }
  pANGn=&piANGn[0];

  for(i=0;i<Lx;i++){
    piT2[i]=&T2[i][0];
  }
  pT2=&piT2[0];

  for(i=0;i<Lx;i++){
    piT2n[i]=&T2n[i][0];
  }
  pT2n=&piT2n[0];

}


void out(int tt){
  int i,j;
  char Aout[200],Tout[200],Bout[200];
  void save_ppm(int run);

  double nabla(double x[Lx][Ly], int i1, int i2);

  ofstream saia;
  ofstream sait;
  ofstream saib;

  sprintf(Aout,"A.%s.%d",fora,tt);
  sprintf(Tout,"T.%s.%d",fora,tt);
  //sprintf(Bout,"B.%s.%d",fora,tt);

  saia.open(Aout);
  sait.open(Tout);
  //saib.open(Bout);

    
  for(i=0;i<Lx;i++){
    for(j=0;j<Ly;j++){
      saia << pa[i][j] << "\t";
      sait <<pT[i][j]+pT2[i][j]<< "\t";
      //saib <<ia[i][j][0]<< "\t";
    }
    saia << "\n";
    sait << "\n";
    //saib << "\n";
  }
  save_ppm(tt);
  cout << "tempo: " << tt * dt<< "   PRINT"<< "\t" << Nbolas << "\t" << "  Posicoes:  ";
  for(i=0;i<Nbolas;i++){
    cout << "("<<Bolas[i][0]<<";"<<Bolas[i][1]<<") ";
  }
  cout <<"\n";

  saia.close();
  sait.close();
  saib.close();
}


int boundx( int xx){
  int res;
  res=xx<0?-xx:(xx>=Lx?(2*(Lx-1)-xx):xx); //volta para tras
  return(res);
}

int boundy(int yy){
  int res;
  res=yy<0?-yy:(yy>=Ly?(2*(Ly-1)-yy):yy); //volta para tras
    //res=(yy+Ly)%Ly; //entra pelo outro lado
  return(res);
}


double prol(double taxarep, double aloc){
  double z;
  //  z= alpha2* (taxarep>0? ( (aloc>0?aloc:0)*(taxarep<Trep?taxarep:(taxarep>2*Trep?(taxarep>3*Trep?0:(3*Trep-taxarep)):Trep) ) ) : ( taxarep*(aloc>-1?(aloc+1)/2:0) ) )  ;
  z= alpha2* (taxarep>0? ( (aloc>0?(aloc<1?aloc:(aloc>1.3?0:-aloc*5+6)):0)*(taxarep<Trep?taxarep:Trep) ) : 0 ) ;
  return z;
}


void save_ppm (int run)
{
  /*byte*/ unsigned char double2byte( double f );

  int i,j,l,w,xpos,ypos;
  int xsize,ysize;
  double phi;
  char filename[16] = "Pic_%s.%06i.ppm";
  char str[60],command[1000];
  /*byte*/unsigned char rbuffer[Lx*xinc][Ly*yinc],gbuffer[Lx*xinc][Ly*yinc],bbuffer[Lx*xinc][Ly*yinc];
  /*byte*/unsigned char ra,rb,rc;

  FILE *fp;

  xsize = Lx*xinc;
  ysize = Ly*yinc;

  for(w=0;w<Ly;w++){
    ypos = w*yinc;
    for(l=0;l<Lx;l++){
      xpos = l*xinc;
      phi = pa[l][w];

      rc = double2byte( ( (1.0 + phi)/2.0 ) );
      rb = double2byte((pT[l][w]+pT2[l][w])/1.1);
      ra = double2byte( ( (pANG[l][w])/1.1  )  );
      for(j=0;j<yinc;j++){
	for(i=0;i<xinc;i++){
	  rbuffer[xpos+i][ypos+j] = ra;
	  gbuffer[xpos+i][ypos+j] = rc;
	  bbuffer[xpos+i][ypos+j] = rb;
	}
      }
    }
  }

  sprintf(str,filename,fora,run);
  fp = fopen( str, "w" );
  fprintf(fp,"P6\n");
  fprintf(fp,"# A .ppm file\n");
  fprintf(fp,"%d %d\n", xsize,ysize);
  fprintf(fp,"255\n");
  for(j=ysize-1;j>=0;j--){
    for(i=0;i<xsize;i++){
      fputc( rbuffer[i][j], fp);
      fputc( gbuffer[i][j], fp);
      fputc( bbuffer[i][j], fp);
    }
  }
  fclose(fp);
}

/*byte*/unsigned char double2byte( double f )
{
  int i;
  /*byte*/ unsigned char b;

  i = (int)(f * 256.0);

  if (i <= 0) { b = (/*byte*/unsigned char)0; }
  else{	if ( i >= 255 )	{ b = (/*byte*/ unsigned char)255;}
    else						{ b = (/*byte*/ unsigned char)i;}
  }

  return b;
}
