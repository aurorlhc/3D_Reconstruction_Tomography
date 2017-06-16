//reconstructing a binary materials in 3D
//from a number of projections at different angles

//the attenuation coefficient is taken to be either 0 or 1
//the "projection" simply is proportaional to the attenuation, instead of (1 - attenuation), etc.
//we assume cone beam geomery and consider that the incoming ray is sufficiently dense, so every pixel will have a projection. The resolution is completeled determined by the detector then...



//Try to modify the code to a 3D reconstruction method using tomography technique. 
//File input method in 2D case now might not be working. It will be changed in the furture...




#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;

#define MAXS 100000 //a module for random number

#define MAXX 100 //system size
#define MAXY 100
#define MAXZ 100
#define Lmax 70 //determing the receptor resolution, the unit is one pixel
                 //the total number of bins in the recepotr 2*Lmax+1
#define L_1 200 //distance between the x-ray source and the center of the specimen.
#define L_2 100 //distance between the center of the specimen and the detector.

#define PI 3.1415926
#define N_p 30 //number of projections, which are assumed to be evenly distributed over pi  

#define NP 446050 //number of black particles

double mu_1 = 1.0; //attenuation coefficient of phase 1
double mu_2 = 0.0; //attenuation coefficient of phase 2


double config[MAXX][MAXY][MAXZ]; //a binary configuration with 0 and 1
//this is used for both storing the initial configuraion 
// and used for the recontruction after the projection is collected...
double projection[N_p][2*Lmax+1][2*Lmax+1]; //the set of projections. 
                                            //fan beam projection length equal to (2*Lmax+1) for test....

double temp_projection[N_p][2*Lmax+1][2*Lmax+1];
//this is the projection profile associated with an intermeidate configuration....
//do not explicit use this quantity, instead use the following diff_projection, to efficiently update energy



double diff_projection[N_p][2*Lmax+1][2*Lmax+1];
//the different between temp_projection and projection
//one can easily update the value of diff_projection
//one can easily compute dE using diff_projection, i.e, dE~ diff_projection*change+ change*change, in which change is only local



int changed_location[N_p][2][2]; //the projected location of the changed pixel
double changed_projection[N_p][2]; //the associated value due to the changed pixel


//need to define an energy as a global variable
double global_energy;


double x_c = (double)MAXX/2.0;
double y_c = (double)MAXY/2.0;
double z_c = (double)MAXZ/2.0;

double d_theta = PI/(double)N_p;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//for simulated annealing reconstruction...

int indexi; int indexj; int indexk; int indexm; int indexn; int indexl;

//The cooling schedule ...
int Nevl = 8000000; //maximum trials at a T-stage

double alpha = 0.98; //cooling rate

int TN = 500; //total T-stages

double T = 0.1; //initial temperature

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void read_config()
{
  FILE* fp;

  if((fp = fopen("Iconfig_3D_HardSphere.txt","r"))==NULL)
    {
      printf("Cannot open file config.txt! Abort!\n");
      exit(1);
    }
 
  int x,y,z;

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
	{
	  config[i][j][k]=0;
	}

  for(int i=0; i<NP; i++)
    {
      fscanf(fp, "%d", &x);
      fscanf(fp, "%d", &y);
      fscanf(fp, "%d", &z);

      config[x][y][z] = 1;
    }

  fclose(fp);
}

void get_projections()
{
  double temp_l; 
  double temp_bin;
  double temp_j;
  int temp_pro_x;
  int temp_pro_y;

  //loop over every angle...
  for(int i=0; i<N_p; i++)
    {
      //now check each point in the phase...
      for(int m=0; m<MAXX; m++)
	    for(int n=0; n<MAXY; n++)
	       for(int l=0; l<MAXZ; l++)
		{	
		  temp_l = (m - x_c)*cos(i*d_theta) + (n - y_c)*sin(i*d_theta);
			
			//temp_bin = dist(m, n)*sin(get_angle(m, n) + i*d_theta);

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



void print_projections()
{
  ofstream fout2;

  fout2.open("projection.txt");

  cout<<"printing projections now..."<<endl;

  for(int i=0; i<N_p; i++)
    {
      //cout<<i<<endl;

      fout2<<"### angle = "<<i*d_theta<<endl;

      for(int r=0; r<(2*Lmax+1); r++)
	for(int p=0; p<(2*Lmax+1); p++)
	{
	  fout2<<(r-Lmax)<<"     "<<(p-Lmax)<<"     "<<projection[i][r][p]<<endl;
	}

      fout2<<"#############################"<<endl;
      fout2<<endl<<endl<<endl;
    }

  fout2.close();
}


void print_config()
{
  ofstream fout2;

  fout2.open("binary_config.txt");

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
      {
	if(config[i][j][k] == 1)
	  {
	    fout2<<i<<"     "<<j<<"     "<<k<<endl;
	  }
      }

  fout2.close();
}


//given the position and angle
//get the projection of the point on the detector...
int get_point_projection_x(int temp_x, int temp_y, int temp_z, double temp_theta)
{
  double temp_l; 
  // double temp_bin;
  int temp_pro_x;
 
  temp_l = (temp_x - x_c)*cos(temp_theta) + (temp_y - y_c)*sin(temp_theta);

  //temp_bin = dist(temp_x, temp_y)*sin(get_angle(temp_x, temp_y) + temp_theta); 

  temp_pro_x = (int)floor(temp_l + Lmax);
 
  return temp_pro_x;
}

//this is actually the z-direction in the materials system
int get_point_projection_y(int temp_x, int temp_y, int temp_z, double temp_theta)
{
  //double temp_bin;
  double temp_i;
  int temp_pro_y;

  //temp_bin = dist(temp_x, temp_y)*sin(get_angle(temp_x, temp_y) + temp_theta); 

  temp_i = temp_z - z_c;

  temp_pro_y = (int)floor(temp_i + Lmax);
 
  return temp_pro_y;
}

//initilaize the initial configuration from the total number of black pixels
//it is assumed that this information can be obtained from the sum of the projection devivided by mu
//since now the only available data is considered to the projection, config is used for recosntruction
void init_reconstr()
{
  //first, clear whatever in the list...

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
      {
        config[i][j][k]=0;
      }

  srand(time(NULL));

  for(int i=0; i<NP; i++)
    {
      int m = rand()%MAXX;
      int n = rand()%MAXY;
      int l = rand()%MAXZ;

      while(config[m][n][l]==1)
	{
	  m = rand()%MAXX;
	  n = rand()%MAXY;
	  l = rand()%MAXZ;
	}

      config[m][n][l]=1;
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
	for(int n=0; n<MAXY; n++)
	  for(int l=0; l<MAXZ; l++)
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
  //#################################################################################
  //now get diff_projection for the first time, and then only work on this quantity
	for(int i=0; i<N_p; i++)
        {

	  for(int j=0; j<2*Lmax+1; j++)
	     for(int k=0; k<2*Lmax+1; k++)
	        diff_projection[i][j][k] = temp_projection[i][j][k] - projection[i][j][k];
	}


   //now compute the global energy for the first time...
	double sum = 0;

	for(int i=0; i<N_p; i++)
	  {

	    for(int j=0; j<2*Lmax+1; j++)
	      for(int k=0; k<2*Lmax+1; k++)
	        sum += diff_projection[i][j][k]*diff_projection[i][j][k];
		 //DO NOT NORMALIZE!!! 
	  }

	global_energy = sum;
  //#################################################################################
  //#################################################################################
}

void change_config()
{
  int i, j, k, m, n, l;
  //int lim = 0;

  //fist we change the lines

        i = rand()%MAXX;
        m = rand()%MAXX;
 
        j = rand()%MAXY; 
        n = rand()%MAXY; 

	k = rand()%MAXZ;
	l = rand()%MAXZ;
 
	//make sure that the two selected points are not from the same phase
	while(config[i][j][k] == config[m][n][l])
	  {
	    i = rand()%MAXX;
	    j = rand()%MAXY;
	    k = rand()%MAXZ;

	    m = rand()%MAXX;
	    n = rand()%MAXY;
	    l = rand()%MAXZ;
	  }

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

//################################################
//################################################
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

//######################################################
//get_dE, directly from diff_projection and changed_projection, changed_location
double get_dE()
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



double PE(double dE, double T) // the probability that should compare ...
{
  if(dE > 0) return exp(-dE/T);
  else return 1;
}


void print_final_config()
{
  ofstream fout2;
  
  fout2.open("final_config.vtk");

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
      {
	fout2<<config[i][j][k]<<endl;
      }

  fout2.close();
  
}


void print_final_projections()
{
  ofstream fout2;

  fout2.open("final_projection.txt");

  cout<<"printing projections now..."<<endl;

  for(int i=0; i<N_p; i++)
    {
      //cout<<i<<endl;

      fout2<<"### angle = "<<i*d_theta<<endl;

      for(int r=0; r<(2*Lmax+1); r++)
	for(int p=0; p<(2*Lmax+1); p++)
	{
	  fout2<<(r-Lmax)<<" "<<(p-Lmax)<<" "<<temp_projection[i][r][p] + diff_projection[i][r][p]<<endl;
	}

      fout2<<"#############################"<<endl;
      fout2<<endl<<endl<<endl;
    }

  fout2.close();
}


int main()
{
  //sample the projection from an input image...

  srand(time(NULL));
  
  read_config();

  cout<<"here"<<endl;

  get_projections();

  cout<<"here"<<endl;

  print_projections();

  print_config();

  //reconstruct an image from the sampled projections...
  //**********************************************************

  init_reconstr();

  int N_acc = 0; //acceptance rate...

  //double E_old, E_new;
  double E_t;
  double dE;
  double energyb = 10000000000.0;

  cout<<"Staring the simulated annealing reconstruction process..."<<endl;
  for(int q=0; q<TN; q++)
    {
      T = alpha*T;

      N_acc = 0;

      cout<<"Stage "<<q+1<<" with T = "<<T<<endl;
       
      for(int i=0; i< Nevl; i++)
	{
	  //cout<<"change "<<i<<endl;

	  change_config();
	  //sample S2 for the new configuration, using time saving methods
	   
	  get_projection_difference();
		  
	  //E_old =  get_energy();
	  
	  update_diff_projection();
	  
	  //E_new = get_energy();
	  
	  dE = get_dE();
	  //cout<<"dE = "<<dE<<endl;
	  
	  double P = (double)(rand()%MAXS)/(double)MAXS;
	  
	  if( P > PE(dE, T))
	    {
	      resume_config();
	      
	      resume_diff_projection();
	      
	      //E_t = E_old;
	      
	    } 
	  else 
	    {
	      N_acc ++;
	      
	      //E_t = E_new;
	      //###############################################
	      global_energy += dE;
	    }
	  
	  E_t = global_energy;

	  //compare and record the best energy and configuration... 
	  
	  if(E_t < energyb)
	    {
	      energyb = E_t;
	      
	    }
	  
	   //printf("%f   %d change has finished... \n", energyt, i+1 );
	  
	}
      
      printf("%d th change of temperature has finished... \n",q+1 );
      cout<<"The acceptance rate: "<<(double)N_acc/(double)Nevl<<endl;
      cout<<"The energy E = "<<energyb<<endl;
            
      printf("*************************************************\n");
      
    }
  
   //*****************************************************************
   //*****************************************************************
   //this is the end of simulated annealing

  print_final_config();

  print_final_projections();
  
}



