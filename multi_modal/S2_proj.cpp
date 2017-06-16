#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>

using namespace std;

#define MAXX 128
#define MAXY 1000000
#define NP 199185
#define Nt MAXX/2
#define Lmax 90
#define PI 3.1415926
#define N_p 5 //number of projections, which are assumed to be evenly distributed over pi  
#define weight 1

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

double diff_S2[Nt];
double change[Nt];
double changeS2[Nt];

double S2_energy;

//projection info~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double mu_1 = 0.00004; //attenuation coefficient of phase 1
double mu_2 = 0.0; //attenuation coefficient of phase 2

double projection[N_p][2*Lmax+1][2*Lmax+1]; //the set of projections. 
double temp_projection[N_p][2*Lmax+1][2*Lmax+1];
double diff_projection[N_p][2*Lmax+1][2*Lmax+1];
int changed_location[N_p][2][2]; //the projected location of the changed pixel
double changed_projection[N_p][2]; //the associated value due to the changed pixel
double proj_energy;
double global_energy;

double x_c = (double)MAXX/2.0;
double y_c = (double)MAXX/2.0;
double z_c = (double)MAXX/2.0;

double d_theta = PI/(double)N_p;

//*****************************************************
//The cooling schedule
double alpha = 0.95;
int TN = 100; // # of times lower down the temperature
double T = 0.00000005;
int Nevl = 5000000; //times of evolution in each T stage 
//*****************************************************

int abs(int a)
{
  if(a>=0) return a;
  else return -a;
}

void read_config()
{
  FILE* fp;

  if((fp = fopen("Fconfig.txt","r"))==NULL)
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
      obj[r] = obj[r]*weight/(double)(3*MAXX*MAXX*MAXX);
    }
}

void get_projections()
{
  double temp_l; 
//double temp_bin;
  double temp_j;
  int temp_pro_x;
  int temp_pro_y;

  //loop over every angle...
  for(int i=0; i<N_p; i++)
    {
      //now check each point in the phase...
      for(int m=0; m<MAXX; m++)
	for(int n=0; n<MAXX; n++)
	  for(int l=0; l<MAXX; l++)
	    {	
	      temp_l = (m - x_c)*cos(i*d_theta) + (n - y_c)*sin(i*d_theta);
			
	      temp_j = l - z_c;

	      temp_pro_x = (int)floor(temp_l + Lmax);
	      temp_pro_y = (int)floor(temp_j + Lmax);

	      if(config[m][n][l] == 1)
		{
		  if(temp_pro_x<(2*Lmax+1) && temp_pro_x>=0 && temp_pro_y<(2*Lmax+1) && temp_pro_y>=0)
		    projection[i][temp_pro_x][temp_pro_y] += mu_1;
		}
	      else
		{
		  if(temp_pro_x<(2*Lmax+1) && temp_pro_x>=0 && temp_pro_y<(2*Lmax+1) && temp_pro_y>=0)
		    projection[i][temp_pro_x][temp_pro_y] += mu_2;
		}
	     
	    }
    }

  cout<<"finish computing the projections..."<<endl;  //get the tomography projections for the sample
}

int get_point_projection_x(int temp_x, int temp_y, int temp_z, double temp_theta)
{
  double temp_l; 
  // double temp_bin;
  int temp_pro_x;
 
  temp_l = (temp_x - x_c)*cos(temp_theta) + (temp_y - y_c)*sin(temp_theta);

  temp_pro_x = (int)floor(temp_l + Lmax);
 
  return temp_pro_x;
}

//this is actually the z-direction in the materials system
int get_point_projection_y(int temp_x, int temp_y, int temp_z, double temp_theta)
{
  double temp_i;
  int temp_pro_y;

  temp_i = temp_z - z_c;

  temp_pro_y = (int)floor(temp_i + Lmax);
 
  return temp_pro_y;
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

  double temp_theta;
  int temp_pro_x;
  int temp_pro_y;
  //now, initialize the projections..
  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<2*Lmax+1; j++)
	for(int k=0; k<2*Lmax+1; k++)
	  temp_projection[i][j][k] = 0;

      temp_theta = i*d_theta;

      //now check each point in the phase...
      for(int m=0; m<MAXX; m++)
	for(int n=0; n<MAXX; n++)
	  for(int l=0; l<MAXX; l++)
	  {
	    temp_pro_x = get_point_projection_x(m, n, l, temp_theta);
	    temp_pro_y = get_point_projection_y(m, n, l, temp_theta);
	    
	    if(config[m][n][l] == 1)
	      {
		if(temp_pro_x<(2*Lmax+1) && temp_pro_x>=0 && temp_pro_y<(2*Lmax+1) && temp_pro_y>=0)
		  temp_projection[i][temp_pro_x][temp_pro_y] += mu_1;
	      }
	    else
	      {
		if(temp_pro_x<(2*Lmax+1) && temp_pro_x>=0 && temp_pro_y<(2*Lmax+1) && temp_pro_y>=0)
		  temp_projection[i][temp_pro_x][temp_pro_y] += mu_2;
	      }
	  }
    }

  //#################################################################################
  //now get diff_projection for the first time, and then only work on this quantity
  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<2*Lmax+1; j++)
	for(int k=0; k<2*Lmax+1; k++)
	  diff_projection[i][j][k] = temp_projection[i][j][k] - projection[i][j][k];
    }

   //now compute the global energy for the first time...
  double sum = 0.0;

  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<2*Lmax+1; j++)
	for(int k=0; k<2*Lmax+1; k++)
	  sum += diff_projection[i][j][k]*diff_projection[i][j][k];
    }

  proj_energy = sum;
  cout<<"Initial proj energy is "<<proj_energy<<endl;
  //#################################################################################
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
      S2[r] = (double)SS2[r]*weight/(double)(3*MAXX*MAXX*MAXX);
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

void get_projection_difference()
{
  int temp_pro1_x, temp_pro1_y;
  int temp_pro2_x, temp_pro2_y;

  double temp_theta; 
  
  //there are the rare event that the two selected pixel have the same value...
  if(config[indexi][indexj][indexk] != config[indexm][indexn][indexl])
    {
      for(int i=0; i<N_p; i++)
	{
	  temp_theta = i*d_theta;
	  
	  temp_pro1_x = get_point_projection_x(indexi, indexj, indexk, temp_theta);
	  temp_pro1_y = get_point_projection_y(indexi, indexj, indexk, temp_theta);
	   
	  temp_pro2_x = get_point_projection_x(indexm, indexn, indexl, temp_theta);
	  temp_pro2_y = get_point_projection_y(indexm, indexn, indexl, temp_theta);
	  
	  changed_location[i][0][0] = temp_pro1_x;
	  changed_location[i][0][1] = temp_pro1_y;

	  changed_location[i][1][0] = temp_pro2_x;
	  changed_location[i][1][1] = temp_pro2_y;

	  if(changed_location[i][0][0] == changed_location[i][1][0] && changed_location[i][0][1] == changed_location[i][1][1])
	    {
	      changed_projection[i][0] = 0;
	      changed_projection[i][1] = 0;
	    }
	  
	  else if(config[indexi][indexj][indexk] == 1)
	    {
	      changed_projection[i][0] = mu_1 - mu_2;
	      changed_projection[i][1] = mu_2 - mu_1;
	    }
	  else if(config[indexi][indexj][indexk] == 0)
	    {
	      changed_projection[i][0] = mu_2 - mu_1;
	      changed_projection[i][1] = mu_1 - mu_2;
	    }
	}
    }
  else //nothing needs to be changed ...
    {
      for(int i=0; i<N_p; i++)
	for(int j=0; j<2; j++)
	  for(int k=0; k<2; k++)
	    {
	      changed_location[i][j][k] = 0;
	      changed_projection[i][j] = 0;
	    }      
    }
}

//##################################################################
void update_diff_projection()
{
  int temp_pro1_x, temp_pro1_y;
  int temp_pro2_x, temp_pro2_y;

  for(int i=0; i<N_p; i++)
    {
      temp_pro1_x = changed_location[i][0][0];
      temp_pro1_y = changed_location[i][0][1];
      temp_pro2_x = changed_location[i][1][0];
      temp_pro2_y = changed_location[i][1][1];

      if(temp_pro1_x<(2*Lmax+1) && temp_pro1_x>=0 && temp_pro1_y<(2*Lmax+1) && temp_pro1_y>=0)
	   diff_projection[i][temp_pro1_x][temp_pro1_y] += changed_projection[i][0];

      if(temp_pro2_x<(2*Lmax+1) && temp_pro2_x>=0 && temp_pro2_y<(2*Lmax+1) && temp_pro2_y>=0)
	   diff_projection[i][temp_pro2_x][temp_pro2_y] += changed_projection[i][1];
    }
}

//#################################################################
void resume_diff_projection()
{
  int temp_pro1_x, temp_pro1_y;
  int temp_pro2_x, temp_pro2_y;

  for(int i=0; i<N_p; i++)
    {
      temp_pro1_x = changed_location[i][0][0];
      temp_pro1_y = changed_location[i][0][1];

      temp_pro2_x = changed_location[i][1][0];
      temp_pro2_y = changed_location[i][1][1];

      if(temp_pro1_x<(2*Lmax+1) && temp_pro1_x>=0 && temp_pro1_y<(2*Lmax+1) && temp_pro1_y>=0)
	   diff_projection[i][temp_pro1_x][temp_pro1_y] -= changed_projection[i][0];

      if(temp_pro2_x<(2*Lmax+1) && temp_pro2_x>=0 && temp_pro2_y<(2*Lmax+1) && temp_pro2_y>=0)
	   diff_projection[i][temp_pro2_x][temp_pro2_y] -= changed_projection[i][1];
    }

}

double get_proj_dE()
{
  double dE= 0.0;

  int temp_pro1_x, temp_pro1_y;
  int temp_pro2_x, temp_pro2_y;

  for(int i=0; i<N_p; i++)
    {
      temp_pro1_x = changed_location[i][0][0];
      temp_pro1_y = changed_location[i][0][1];

      temp_pro2_x = changed_location[i][1][0];
      temp_pro2_y = changed_location[i][1][1];

      if(temp_pro1_x<(2*Lmax+1) && temp_pro1_x>=0 && temp_pro1_y<(2*Lmax+1) && temp_pro1_y>=0)
	dE += (2*diff_projection[i][changed_location[i][0][0]][changed_location[i][0][1]]*changed_projection[i][0] - changed_projection[i][0]*changed_projection[i][0]);
	  
      if(temp_pro2_x<(2*Lmax+1) && temp_pro2_x>=0 && temp_pro2_y<(2*Lmax+1) && temp_pro2_y>=0)
	dE += (2*diff_projection[i][changed_location[i][1][0]][changed_location[i][1][1]]*changed_projection[i][1] - changed_projection[i][1]*changed_projection[i][1]); 
    }

  return dE;
}

double energy(double S[Nt])
{
  double E = 0.0;

  for(int i=1; i<Nt; i++)
    {
      E = E + (S[i] - obj[i])*(S[i] - obj[i]);
    }
  
  return E;
}

double get_S2_dE(double S2[Nt], double ST2[Nt])
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

void final_proj_energy()
{

  double temp_theta;
  int temp_pro_x;
  int temp_pro_y;
  //now, initialize the projections..
  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<2*Lmax+1; j++)
	for(int k=0; k<2*Lmax+1; k++)
	  {
	    temp_projection[i][j][k] = 0;
	    diff_projection[i][j][k] = 0;
	  }
      temp_theta = i*d_theta;

      //now check each point in the phase...
      for(int m=0; m<MAXX; m++)
	for(int n=0; n<MAXX; n++)
	  for(int l=0; l<MAXX; l++)
	  {
	    temp_pro_x = get_point_projection_x(m, n, l, temp_theta);
	    temp_pro_y = get_point_projection_y(m, n, l, temp_theta);
	    
	    if(config[m][n][l] == 1)
	      {
		if(temp_pro_x<(2*Lmax+1) && temp_pro_x>=0 && temp_pro_y<(2*Lmax+1) && temp_pro_y>=0)
		  temp_projection[i][temp_pro_x][temp_pro_y] += mu_1;
	      }
	    else
	      {
		if(temp_pro_x<(2*Lmax+1) && temp_pro_x>=0 && temp_pro_y<(2*Lmax+1) && temp_pro_y>=0)
		  temp_projection[i][temp_pro_x][temp_pro_y] += mu_2;
	      }
	  }
    }

  //#################################################################################
  //now get diff_projection for the first time, and then only work on this quantity
  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<2*Lmax+1; j++)
	for(int k=0; k<2*Lmax+1; k++)
	  diff_projection[i][j][k] = temp_projection[i][j][k] - projection[i][j][k];
    }

   //now compute the global energy for the first time...
  double sum = 0.0;

  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<2*Lmax+1; j++)
	for(int k=0; k<2*Lmax+1; k++)
	  sum += diff_projection[i][j][k]*diff_projection[i][j][k];
    }

  proj_energy = sum;
  //cout<<"Initial proj energy is "<<proj_energy<<endl;
 

}

void print_final_config()
{
  ofstream fout2;
  
  fout2.open("final_config.vtk");

  fout2<<"# vtk DataFile Version 2.0"<<endl;
  fout2<<"2D_to_3D example"<<endl;
  fout2<<"ASCII"<<endl;
  fout2<<"DATASET STRUCTURED_POINTS"<<endl;
  fout2<<"DIMENSIONS 128 128 128"<<endl;
  fout2<<"SPACING 1 1 1"<<endl;
  fout2<<"ORIGIN 0 0 0"<<endl;
  fout2<<"POINT_DATA 2097152"<<endl;
  fout2<<"SCALARS volume_scalars UNSIGNED_CHAR 1"<<endl;
  fout2<<"LOOKUP_TABLE default"<<endl;

  for(int k=0; k<MAXX; k++)
    for(int j=0; j<MAXX; j++)
      for(int i=0; i<MAXX; i++)
      {
	if(config[i][j][k] == 1)
          config[i][j][k] = 255;
	fout2<<config[i][j][k]<<endl;
      }

  fout2.close(); 
}


int main()
{
  srand(time(NULL));
  read_config();
  get_projections();
  calculate_target_S2();

  /*
  cout<<"Finish calculating target S2"<<endl;
  cout<<"************************************"<<endl;

  FILE* fp = fopen("TS2.txt", "w");
  for(int r=0; r<Nt; r++)
    {
      fprintf(fp, "%d \t  %lf \n", r, obj[r]);
    }
  fclose(fp);  
  */

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
  
  cout<<"The initial S2 energy is "<<sum<<endl;
  S2_energy = sum;
  global_energy = proj_energy + S2_energy;
  cout<<"Initial global energy is "<<global_energy<<endl;

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

	  get_projection_difference();

	  update_diff_projection();

	  double P = double (rand() % MAXY)/(double) MAXY;
	  double dE = get_S2_dE(S2,ST2) + get_proj_dE();

	  if(P > PE(dE, T))
	    {
	      resume_config();
	      resume_diff_projection();

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
  double re = 0.0;
  for(int r=1; r<Nt; r++)
    {
      diff_S2[r] = S2[r] - obj[r];
      re += diff_S2[r]*diff_S2[r];
    } 
  cout<<"final_S2_energy = "<<re<<endl; 
  final_proj_energy();  
  cout<<"final_pro_ energy = "<<proj_energy<<endl;
  print_final_config();
  //this is the end of the codes...  
}
