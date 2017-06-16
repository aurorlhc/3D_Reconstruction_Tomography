using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>

#define MAXX 100
#define MAXY 1000000
#define NP 446050
#define Nt MAXX/2

int config[MAXX][MAXX][MAXX];

int indexi; int indexj; int indexk; 
int indexm; int indexn; int indexl;

int pix_position[MAXX];
int pix_counter; 

//target S2
int S2_x[Nt];
int S2_y[Nt];
int S2_z[Nt];
double obj[Nt];

//need to index to pin the 1D array
long int lineS2[MAXX][MAXX][Nt];
long int columeS2[MAXX][MAXX][Nt];
long int heightS2[MAXX][MAXX][Nt];

//reconstruct S2
double SS2[Nt];
double S2[Nt];
double ST2[Nt];
int SLT[2][Nt];
int SCT[2][Nt];
int SHT[2][Nt];

double global_energy;
double diff_S2[Nt];

double change[Nt];
double changeS2[Nt];

//*****************************************************
//The cooling schedule
double alpha = 0.95;
int TN = 150; // # of times lower down the temperature
double T = 0.0000001;
int Nevl = 100000; //times of evolution in each T stage 
//*****************************************************


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

void calculate_target_S2()
{
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      for(int k=0; k<MAXX; k++)
	for(int r=0; r<Nt; r++)
	  {
	    if(i+r<MAXX)
	      {
		S2_x[r] += config[i][j][k]*config[i+r][j][k];
	      }
	    else
	      {
		S2_x[r] += config[i][j][k]*config[i+r-MAXX][j][k];
	      }

	    if(j+r<MAXX)
	      {
		S2_y[r] += config[i][j][k]*config[i][j+r][k];
	      }
	    else
	      {
		S2_y[r] += config[i][j][k]*config[i][j+r-MAXX][k];
	      }

	    if(k+r<MAXX)
	      {
		S2_z[r] += config[i][j][k]*config[i][j][k+r];
	      }
	    else
	      {
		S2_z[r] += config[i][j][k]*config[i][j][k+r-MAXX];
	      }
	  }

  for(int r=0; r<Nt; r++)
    {
      obj[r] += S2_x[r] + S2_y[r] + S2_z[r];
      obj[r] = obj[r]/(double)(3*MAXX*MAXX*MAXX);
    }
}

void init_reconstr()
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
	temp_dist = abs(pix_position[i]-pix_position[j]);
	
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
      if(config[index1][i][index2]==1)
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
	temp_dist = abs(pix_position[i]-pix_position[j]);
	
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
      if(config[index1][index2][i]==1)
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
	temp_dist = abs(pix_position[i]-pix_position[j]);
	
	if(temp_dist>=MAXX/2) temp_dist = MAXX-temp_dist;
	
	heightS2[index1][index2][temp_dist]++;	
      }  
}


void calculate_temp_S2()
{
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXX; j++)
      {
	sampleS2line(i,j);
	sampleS2colume(i,j);
	sampleS2height(i,j);

	for(int r=0; r<Nt; r++)
	  {
	    SS2[r] += lineS2[i][j][r];
	    SS2[r] += columeS2[i][j][r];
	    SS2[r] += heightS2[i][j][r];
	  }
      } 

  for(int r=0; r<Nt; r++)
    {
      S2[r] = (double)SS2[r]/(double)(3*MAXX*MAXX*MAXX);
    }
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

double energy(double S[Nt])
{
  double E=0.0;

  for(int i=1; i<Nt; i++)
    {
      E = E + (S[i] - obj[i])*(S[i] - obj[i]);
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
  else return 1;
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
  calculate_target_S2();
  cout<<"Finish calculating target S2"<<endl;
  cout<<"************************************"<<endl;

  FILE* fp = fopen("TS2.txt", "w");
  for(int r=0; r<Nt; r++)
    {
      fprintf(fp, "%d \t  %lf \n", r, obj[r]);
    }
  fclose(fp);  
  
  //Now generate a random configure and calculate the corresponding S2
  init_reconstr();
  calculate_temp_S2();
  
  //calculate the energy for the first time
  double sum = 0.0;
  
  for(int r=0; r<Nt; r++)
    {
      diff_S2[r] = S2[r] - obj[r];
      sum += diff_S2[r]*diff_S2[r];
    } 
  
  cout<<"The initial energy is "<<sum<<endl;
  global_energy = sum;
  
  //**************************************************
  //simulated annealing procedure to evolve the system
  //**************************************************
  
  int N_acc = 0; //acceptance rate...
  
  int TEMP[Nt];
  int LCHV[Nt];
  
  cout<<"Staring the simulated annealing reconstruction process..."<<endl;
  for(int q=0; q<TN; q++)
    {
      T = alpha*T;
      
      N_acc = 0;
      
      cout<<"Stage "<<q+1<<" with T = "<<T<<endl;
      
      for(int i=0; i< Nevl; i++)
	{
	  change_config();

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

	      SS2[r] = (int)(SS2[r] - TEMP[r] + LCHV[r]);
	       
	      ST2[r] = (double)SS2[r]/(double)(3*MAXX*MAXX*MAXX);	      
	    }

	  double P = double (rand() % MAXY)/(double) MAXY;
	  double dE = get_dE(S2,ST2);

	  if( P > PE(dE, T))
	    {
	      resume_config();
	      
	      for(int r=0; r<Nt; r++)
		{
		  SS2[r] = (int)(SS2[r] + TEMP[r] - LCHV[r]);
		  //now we resume lineS2, columeS2 and heightS2;
		  lineS2[indexj][indexk][r] = SLT[0][r];
		  lineS2[indexn][indexl][r] = SLT[1][r];
		  
		  columeS2[indexi][indexk][r] = SCT[0][r];
		  columeS2[indexm][indexl][r] = SCT[1][r];
		  
		  heightS2[indexi][indexj][r] = SHT[0][r];
		  heightS2[indexm][indexn][r] = SHT[1][r];		  
		}
	    } 	  
	  else 
	    {	   
	      for(int r=0; r<Nt; r++)
		{
		  S2[r] = ST2[r];
		}

	      N_acc++;
	    }
	  
	  global_energy = energy(S2);
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
