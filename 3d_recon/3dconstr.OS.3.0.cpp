//*******************************************************************
//* 3D (Re)Construction from Radial S2(r) using Orthogonal Sampling *
//*******************************************************************


//ORTHOGONAL SAMPLING is used to generate large systems
//Both lattice-point and gray-scale method has limitations on system size
//Sometimes, small systems can not accurately reproduce the correlations...



//modified: 05/13/2012
//using a more efficient sampling method for computing S_2 for each line
// i.e., instead of actually keeping track of line segments moving through 
// the system, we focus on the black pixels as particles and directly
// compute the pair distances and then normalize them 
// This is like the orthgonal version of the lattice-point method



using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>

#define MAXX 100
#define MAXS 1000000
#define NP 446050
#define Nt MAXX/2

//double f1 = (double)NP/(double)(MAXX*MAXX*MAXX);//volume fraction of black pixels
//double f2 = 1 - f1;

int indexi; int indexj; int indexk;
int indexm; int indexn; int indexl;

int pix_position[MAXX]; //the position of black pixels
int pix_counter; //the number of pix on each line/column/height

int config[MAXX][MAXX][MAXX];

//target S2
double TargetS2[Nt];

//need to index to pin the 1D array
long int lineS2[MAXX][MAXX][Nt];
long int columeS2[MAXX][MAXX][Nt];
long int heightS2[MAXX][MAXX][Nt];

//reconstruct S2
long int SS2[Nt];
double S2[Nt];

//save for temp use
double ST2[Nt];
long int SLT[2][Nt];
long int SCT[2][Nt];
long int SHT[2][Nt];

double global_energy;

//these are for updating SS2 efficiently...
int TEMP[Nt];
int LCHV[Nt];


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//The cooling schedule
double alpha = 0.99; //cooling rate
int TN = 150; //# of temperature stage
double T = 0.0000001; //initial temperature
int Nevl = 100000; //times of evolution in each T stage 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


int abs(int a)
{
  if(a>=0) return a;
  else return -a;
}


void read_config()
{
  FILE* fp;

  if((fp = fopen("Iconfig_3D_HardSphere.txt","r"))==NULL)
    {
      printf("Can not open file Mconfig.txt! Abort!\n");
      exit(1);
    }

  int xt;
  int yt; 
  int zt;

  for(int i=0; i<NP; i++)
    {
      fscanf(fp, "%d", &xt);
      fscanf(fp, "%d", &yt);
      fscanf(fp, "%d", &zt);

      config[xt][yt][zt] = 1;
    }

  fclose(fp);
}


void sampleS2line(int index1, int index2)
{
  
  for(int r=0; r<Nt; r++)
    {
      lineS2[index1][index2][r] = 0;
    }

  //serach the line for pixel positions
  pix_counter = 0;

  for(int i=0; i<MAXX; i++)
    {
      pix_position[i] = 0;
    }

  for(int i=0; i<MAXX; i++)
    {
      if(config[i][index1][index2] == 1)
	{
	  pix_position[pix_counter] = i;
	  pix_counter++;
	}
    }

  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position[i] - pix_position[j]);
	
	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;
	
	lineS2[index1][index2][temp_dist]++;
      }
}


void sampleS2colume(int index1, int index2)
{

  for(int r=0; r<Nt; r++)
    {
      columeS2[index1][index2][r] = 0;
    }
  
  //serach the line for pixel positions
  pix_counter = 0;

  for(int i=0; i<MAXX; i++)
    {
      pix_position[i] = 0;
    }

  for(int i=0; i<MAXX; i++)
    {
      if(config[index1][i][index2] == 1)
	{
	  pix_position[pix_counter] = i;
	  pix_counter++;
	}
    }

  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position[i] - pix_position[j]);
	
	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;
	
	columeS2[index1][index2][temp_dist]++;
      }
}


void sampleS2height(int index1, int index2)
{

  for(int r=0; r<Nt; r++)
    {
      heightS2[index1][index2][r] = 0;
    }
  
  //serach the line for pixel positions
  pix_counter = 0;

  for(int i=0; i<MAXX; i++)
    {
      pix_position[i] = 0;
    }
  
  for(int i=0; i<MAXX; i++)
    {
      if(config[index1][index2][i] == 1)
	{
	  pix_position[pix_counter] = i;
	  pix_counter++;
	}
    }
  
  //now get the distance between all pixels on the line...
  int temp_dist;

  for(int i=0; i<pix_counter; i++)
    for(int j=0; j<=i; j++)
      {
	temp_dist = abs(pix_position[i] - pix_position[j]);
	
	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;
	
	heightS2[index1][index2][temp_dist]++;	
      }  
}


void calculate_S2()
{
  for(int r=0; r<Nt; r++)
    {
      SS2[r] = 0;
      S2[r] = 0;
    }

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      {
	sampleS2line(i,j);
	sampleS2colume(i,j);
	sampleS2height(i,j);

	for(int r=0; r<Nt; r++)
	  {
	    SS2[r] += lineS2[i][j][r] + columeS2[i][j][r] + heightS2[i][j][r];
	  }
      } 

  for(int r=0; r<Nt; r++)
    {
      S2[r] = (double)SS2[r]/(double)(3*MAXX*MAXX*MAXX);
    }
}


void init_config()
{
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      for(int k=0; k<MAXX; k++)
      {
        config[i][j][k] = 0;
      }

  for(int i=0; i<NP; i++)
    {
      int m = rand()%MAXX;
      int n = rand()%MAXX;
      int l = rand()%MAXX;

      while(config[m][n][l] == 1)
	{
	  m = rand()%MAXX;
	  n = rand()%MAXX;
	  l = rand()%MAXX;
	}

      config[m][n][l] = 1;
    }

  calculate_S2();
  
  //calculate the energy for the first time
  double sum = 0.0;
  
  for(int r=0; r<Nt; r++)
    {
      sum += (S2[r] - TargetS2[r])*(S2[r] - TargetS2[r]);
    } 
  
  cout<<"The initial energy is "<<sum<<endl;
  global_energy = sum;
}


void change_config()
{
  int i, j, k, m, n, l;
 
  do{   i = rand() % MAXX;
        m = rand() % MAXX;

        j = rand() % MAXX; 
        n = rand() % MAXX;

        k = rand() % MAXX; 
        l = rand() % MAXX;

     }while(config[i][j][k] == config[m][n][l]);

  int temp;
  temp = config[i][j][k];
  config[i][j][k] = config[m][n][l];
  config[m][n][l] = temp;

  indexi = i;
  indexj = j;
  indexk = k;
  indexm = m;
  indexn = n;
  indexl = l;
}


void resume_config()
{
  int temp;
  //first we resume the config
  temp = config[indexi][indexj][indexk];
  config[indexi][indexj][indexk] = config[indexm][indexn][indexl];
  config[indexm][indexn][indexl] = temp;
}


void get_temp_S2()
{
  for(int r=0; r<Nt; r++)
    {
      SLT[0][r] = lineS2[indexj][indexk][r];
      SLT[1][r] = lineS2[indexn][indexl][r];
      
      SCT[0][r] = columeS2[indexi][indexk][r];
      SCT[1][r] = columeS2[indexm][indexl][r];
      
      SHT[0][r] = heightS2[indexi][indexj][r]; 
      SHT[1][r] = heightS2[indexm][indexn][r];
    }
  
  sampleS2line(indexj, indexk);
  sampleS2line(indexn, indexl);	     
	  
  sampleS2colume(indexi, indexk);
  sampleS2colume(indexm, indexl);
	  
  sampleS2height(indexi, indexj);
  sampleS2height(indexm, indexn);
	  
	  
  //Now we compute the S2 for the new configuration...
  for(int r=0; r<Nt; r++)
    {       
      //the following method only consider the changes.. 
      //************************************************
      TEMP[r] = 0; //old value 
      LCHV[r] = 0; //current value
	      
      for(int vv=0; vv<2; vv++)
	{
	  TEMP[r] += (SLT[vv][r] + SCT[vv][r] + SHT[vv][r]);
	}
	      
      LCHV[r] = (lineS2[indexj][indexk][r] + lineS2[indexn][indexl][r]) + 
	(columeS2[indexi][indexk][r] + columeS2[indexm][indexl][r]) + 
	(heightS2[indexi][indexj][r] + heightS2[indexm][indexn][r]);  
    }
}

void update_diff_S2()
{
  for(int r=0; r<Nt; r++)
    {
      SS2[r] = SS2[r] + LCHV[r] - TEMP[r]; //new - old
      ST2[r] = (double)SS2[r]/(double)(3*MAXX*MAXX*MAXX);
    }
}

void resume_diff_S2()
{
  for(int r=0; r<Nt; r++)
    {
      SS2[r] = SS2[r] + TEMP[r] - LCHV[r];
      //now we resume lineS2, columeS2 and heightS2;
      lineS2[indexj][indexk][r] = SLT[0][r];
      lineS2[indexn][indexl][r] = SLT[1][r];
      
      columeS2[indexi][indexk][r] = SCT[0][r];
      columeS2[indexm][indexl][r] = SCT[1][r];
      
      heightS2[indexi][indexj][r] = SHT[0][r];
      heightS2[indexm][indexn][r] = SHT[1][r];		  
    }
}

double energy(double S[Nt])
{
  double E = 0.0;

  for(int i=0; i<Nt; i++)
    {
      E += (S[i] - TargetS2[i])*(S[i] - TargetS2[i]);
    }

  return E;
}

double get_dE(double S2[Nt], double ST2[Nt])
{
  double dE = 0.0;

  dE = energy(ST2) - energy(S2);

  return dE;
}


double PE(double dE, double T) // the probability that should compare ...
{
  if(dE > 0) return exp(-dE/T);
  else return 1.0;
}

void print_final_config()
{
  ofstream fout2;
  
  fout2.open("final_config.vtk");

  fout2<<"# vtk DataFile Version 2.0"<<endl;
  fout2<<"2D_to_3D example"<<endl;
  fout2<<"ASCII"<<endl;
  fout2<<"DATASET STRUCTURED_POINTS"<<endl;
  fout2<<"DIMENSIONS 100 100 100"<<endl;
  fout2<<"SPACING 1 1 1"<<endl;
  fout2<<"ORIGIN 0 0 0"<<endl;
  fout2<<"POINT_DATA 1000000"<<endl;
  fout2<<"SCALARS volume_scalars UNSIGNED_CHAR 1"<<endl;
  fout2<<"LOOKUP_TABLE default"<<endl;

  for(int k=0; k<MAXX; k++)
    for(int j=0; j<MAXX; j++)
      for(int i=0; i<MAXX; i++)
      {
	fout2<<config[i][j][k]<<endl;
      }

  fout2.close(); 
}

int main()
{
  srand(time(NULL));
  read_config();

  //now calculate S2 for the target configure
  cout<<"initial sampling S2...."<<endl;
  calculate_S2();
 
  for(int r=0; r<Nt; r++)
    {
      TargetS2[r] = S2[r];  //store the target S2
    }

  cout<<"************************************"<<endl;

  FILE* fp = fopen("TargetS2.txt", "w");
  for(int r=0; r<Nt; r++)
    {
      fprintf(fp, "%d \t  %lf \n", r, TargetS2[r]);
    }
  fclose(fp);  
  
  //Now generate a random configure and calculate the corresponding S2
  init_config();  
  int N_acc = 0; //acceptance ...
  double dE;
  
  cout<<"Staring the simulated annealing reconstruction process..."<<endl;
  for(int q=0; q<TN; q++)
    {
      T = alpha*T;
      
      N_acc = 0;
      
      cout<<"Stage "<<q+1<<" with T = "<<T<<endl;
      
      for(int ll=0; ll<Nevl; ll++)
	{
	  change_config();
	  get_temp_S2();
	  update_diff_S2();
	  dE = get_dE(S2,ST2);

	  double P = double (rand() % MAXS)/(double) MAXS;

	  if( P > PE(dE, T))
	    {
	      resume_config();
	      resume_diff_S2();
	    } 	  
	  else 
	    {	   
	      global_energy += dE;
	      N_acc++;
	    }
	  
	}
      
      printf("%d th change of temperature has finished... \n",q+1 );
      cout<<"The acceptance rate: "<<(double)N_acc/(double)Nevl<<endl;
      cout<<"The energy E = "<<global_energy<<endl;
      
      printf("*************************************************\n");
    }
  
  //*****************************************************************
  //this is the end of simulated annealing
  
  print_final_config();
  //this is the end of the codes...  
}
