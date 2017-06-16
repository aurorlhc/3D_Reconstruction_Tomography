//use the sample data from Dr. Chawla's group
//reconstructing a binary materials in 2D from a number of projections at different angles

//make operations for the "projection" being proportaional to the attenuation coefficient
//use parallel beam geomery and consider that the incoming ray is sufficiently dense, so every pixel will have a projection. 
//The resolution is completeled determined by the detector


using namespace std;

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAXS 100000 //a module for random number

#define MAXX 252   //system size
#define MAXY 252
#define Lmax 126   //determing the receptor resolution, the unit is one pixel
                   //the total number of bins in the recepotr 2*Lmax+1
#define N_pr 1860  //number of projections, which are assumed to be evenly distributed over pi  
#define N_pi 2016  //number of pixels in one projection

#define N_p 30   //select 30 projections from the 1860 projections
#define N_rp 252 //reduce the resolution of pixels

#define N 0.001168   //pixel size

double mu_1 = 0.00662256; //attenuation coefficient of Al
double mu_2 = 0.0;  //attenuation coefficient of air
double mu = 5.67;

int N_tot; //the total number of black pixels...

double config[MAXX][MAXY]; //a binary configuration with 0 and 1
                           //this is used for the recontruction after the projection is collected...

double initial_projection[N_p][N_rp];
double log_projection[N_pr][N_pi];     // -ln(I/I0)=mu*L
double raw_projection[N_pr][N_pi];     //the set of raw projections
double projection[N_p][2*Lmax];        //reduce the resolution and the projection numbers

double temp_projection[N_p][2*Lmax+1];   //this is the projection profile associated with an intermeidate configuration....
int changed_location[N_p][2];          //the projected location of the changed pixel
double changed_projection[N_p][2];     //the associated value due to the changed pixel

double diff_projection[N_p][2*Lmax+1];

double global_energy = 0.0;

double x_c = (double)MAXX/2.0;
double y_c = (double)MAXY/2.0;

double d_theta = 180/(double)N_p; //pick up one projection per 62 projections to get 30 projections over 180 degree

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//for simulated annealing reconstruction...

int indexi; int indexj; int indexm; int indexn;

//The cooling schedule ...
int Nevl = 8000000; //maximum trials at a T-stage

double alpha = 0.99; //cooling rate

int TN = 400; //total T-stages

double T = 0.005; //initial temperature

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void read_data()
{
  double value;
  double number;

  ifstream fin;

  fin.open("raw_projection.txt");

  for(int i=0; i<N_pr; i++)
    for(int j=0; j<N_pi; j++)
      {
	fin>>value;
	raw_projection[i][j] = value;
	
	if(raw_projection[i][j] > 0.9)
	  raw_projection[i][j] = 1;

	log_projection[i][j] = -(double)log(raw_projection[i][j]);
      }

  fin.close();
}

//reduce resolution and projection numbers from the raw_projection
void get_projections() 
{
  for(int i=0; i<N_p; i++)
    for(int j=0; j<N_rp; j++)
      {
	for(int n=0; n<8; n++)
	  {
	    initial_projection[i][j] += 0.125*log_projection[i*N_pr/N_p][j*8+n];
	  }
      }

  for(int i=0; i<N_p; i++)
    for(int j=3; j<N_rp-3; j++)
      {
	if(initial_projection[i][j] != 0 && initial_projection[i][j-1] == 0 && initial_projection[i][j-2] == 0 && initial_projection[i][j-3] == 0 && initial_projection[i][j+1] == 0 && initial_projection[i][j+2] == 0 && initial_projection[i][j+3] == 0)
	  initial_projection[i][j] = 0;
      }

  ofstream fout;
  fout.open("initial_projection.txt");

  for(int i=0; i<N_p; i++)
    for(int j=0; j<N_rp; j++)
      {
	fout<<initial_projection[i][j]<<endl;
      }

  fout.close();

  for(int i=0; i<N_p; i++)
    for(int j=0; j<N_rp; j++)
      {
	projection[i][j] = initial_projection[i][j];
      }

}

void get_black_pixels()
{
  int temp;
  int N_temp;
  double L = 0.0;
  double n = 0.0;

  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<2*Lmax+1; j++)
	{
	  n = projection[i][j]/mu_1;
	  temp = (int)n;
	  /*
	  if(n-0.5>=temp)
	    N_temp = (int)ceil(n);
	  else
	    N_temp = (int)floor(n);
	  */
	  projection[i][j] = temp*mu_1;

	  N_tot += temp;

	  //L += projection[i][j]; //get the length of black pixels
	}
      // cout<<"The length of "<<i<<"th projection is equal to "<<L<<endl;
    }

  N_tot = (int)N_tot/N_p;
  cout<<"the average pixel number is "<<N_tot<<endl<<endl;
}

void print_projections()
{
  ofstream fout2;

  fout2.open("projection.txt");

  cout<<"printing projections now..."<<endl;

  for(int i=0; i<N_p; i++)
    {
      fout2<<"### angle = "<<i*d_theta<<endl;

      for(int r=0; r<(2*Lmax+1); r++)
	{
	  fout2<<(r-Lmax)<<" "<<projection[i][r]<<endl;
	}

      fout2<<"#############################"<<endl;
      fout2<<endl<<endl<<endl;
    }

  fout2.close();
}

int get_point_projection(int temp_x, int temp_y, double temp_theta)
{
  double temp_l; 
  int temp_bin;

  //temp_l = dist(temp_x, temp_y)*cos(get_angle(temp_x, temp_y) - temp_theta);
  temp_l = (temp_x - x_c)*cos(temp_theta) - (temp_y-y_c)*sin(temp_theta);

  temp_bin = (int)floor(temp_l+Lmax); 
  
  return temp_bin;
}

//initilaize the initial configuration from the total number of black pixels
//it is assumed that this information can be obtained from the sum of the projection devivided by mu
//since now the only available data is considered to the projection, config is used for recosntruction
void init_reconstr()
{
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      {
        config[i][j] = 0;
      }
  
  for(int i=0; i<N_tot; i++)
    {
      int m = rand()%MAXX;
      int n = rand()%MAXY;

      while(config[m][n]==1)
	{
	  m = rand()%MAXX;
	  n = rand()%MAXY;
	}

      config[m][n] = 1;
    }

  double temp_theta;
  int temp_bin;
  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<(2*Lmax+1); j++)
	temp_projection[i][j] = 0;

      temp_theta = i*d_theta;

      //now check each point in the phase...
      for(int m=0; m<MAXX; m++)
	for(int n=0; n<MAXY; n++)
	  {
	    temp_bin = get_point_projection(m, n, temp_theta);
	    
	    if(config[m][n] == 1)
	      {
		if(temp_bin<(2*Lmax+1) && temp_bin>=0)
		  temp_projection[i][temp_bin] += mu_1;
	      }
	    else
	      {
		if(temp_bin<(2*Lmax+1) && temp_bin>=0)
		  temp_projection[i][temp_bin] += mu_2;
	      }
	     
	  }
      
    }
  
//*************************************************************************
  double sum = 0.0;

  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<(2*Lmax+1); j++)
	diff_projection[i][j] = temp_projection[i][j] - projection[i][j];
    }

  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<(2*Lmax+1); j++)
	sum += diff_projection[i][j]*diff_projection[i][j];
    }

  global_energy = sum;
  cout<<"Initial_diff_energy = "<<global_energy<<endl;
 //************************************************************************

  ofstream fout;

  fout.open("initial_diff_projection.txt");

  for(int i=0; i<N_p; i++)
    {
      for(int r=0; r<(2*Lmax+1); r++)
	{
	  fout<<i<<" "<<r<<" "<<diff_projection[i][r]<<" "<<diff_projection[i][r]*diff_projection[i][r]<<endl;
	}
    }

  fout.close();

}

void change_config()
{
  int i, j, m, n;

  i = rand()%MAXX;
  j = rand()%MAXY;
 
  m = rand()%MAXX; 
  n = rand()%MAXY; 
  
  while(config[i][j] == config[m][n])
    {
      i = rand()%MAXX;
      j = rand()%MAXY;
 
      m = rand()%MAXX;
      n = rand()%MAXY;
    }

  double temp;
  
  temp = config[i][j];
  config[i][j] = config[m][n];
  config[m][n] = temp;
  
  indexi = i;
  indexj = j;
  indexm = m;
  indexn = n;

}


void resume_config()
{
  double temp;
  //first we resume the config
  temp = config[indexi][indexj];
  config[indexi][indexj] = config[indexm][indexn];
  config[indexm][indexn] = temp;

}


void get_projection_difference()
{
  int temp_bin1, temp_bin2;
  double temp_theta; 
  
  //there are the rare event that the two selected pixel have the same value...
  if(config[indexi][indexj] != config[indexm][indexn])
    {
      for(int i=0; i<N_p; i++)
	{
	  temp_theta = i*d_theta;
	  
	  temp_bin1 = get_point_projection(indexi, indexj, temp_theta);
	  
	  temp_bin2 = get_point_projection(indexm, indexn, temp_theta);
	  
	  changed_location[i][0] = temp_bin1;
	  changed_location[i][1] = temp_bin2;



	  //get rid of the effect when the two choosing points are along the same line of ray
	  //**************************************************
	  //**************************************************
	  if(changed_location[i][0] == changed_location[i][1])
	    {
	      changed_projection[i][0] = 0;
	      changed_projection[i][1] = 0;
	    }
	  //**************************************************
	  //**************************************************	  



	  //this means the old pixel is 0
	  else if(config[indexi][indexj] == 1)
	    {
	      changed_projection[i][0] = mu_1 - mu_2;
	      changed_projection[i][1] = mu_2 - mu_1;
	    }
	  else
	    {
	      changed_projection[i][0] = mu_2 - mu_1;
	      changed_projection[i][1] = mu_1 - mu_2;
	    }
	  
	}
    }

  else
    {
      for(int i=0; i<N_p; i++)
	for(int j=0; j<2; j++)
	  {
	    changed_location[i][j] = 0;
	    changed_projection[i][j] = 0;
	  }
    }

}


void update_diff_projection()
{

  int temp_bin1, temp_bin2;

  for(int i=0; i<N_p; i++)
    {
      temp_bin1 = changed_location[i][0];
      temp_bin2 = changed_location[i][1];

      if(temp_bin1<(2*Lmax+1) && temp_bin1>=0)
	diff_projection[i][temp_bin1] += changed_projection[i][0];

      if(temp_bin2<(2*Lmax+1) && temp_bin2>=0)
	diff_projection[i][temp_bin2] += changed_projection[i][1];
    }
}


void resume_diff_projection()
{
  int temp_bin1, temp_bin2;

  for(int i=0; i<N_p; i++)
    {
      temp_bin1 = changed_location[i][0];
      temp_bin2 = changed_location[i][1];

      if(temp_bin1<(2*Lmax+1) && temp_bin1>=0)
	diff_projection[i][temp_bin1] -= changed_projection[i][0];

      if(temp_bin2<(2*Lmax+1) && temp_bin2>=0)
	diff_projection[i][temp_bin2] -= changed_projection[i][1];
    }

}

double get_energy()
{
  double dE_0 = 0.0;
  double dE_1 = 0.0;

  int temp_bin1, temp_bin2;

  for(int i=0; i<N_p; i++)
    {
      temp_bin1 = changed_location[i][0];
      temp_bin2 = changed_location[i][1];

      if(temp_bin1<(2*Lmax+1) && temp_bin1>=0)
	dE_0 += 2*diff_projection[i][temp_bin1]*changed_projection[i][0] - changed_projection[i][0]*changed_projection[i][0];

      if(temp_bin2<(2*Lmax+1) && temp_bin2>=0)
	dE_1 += 2*diff_projection[i][temp_bin2]*changed_projection[i][1] - changed_projection[i][1]*changed_projection[i][1];
    }
 
  return (dE_0 + dE_1);
}


double PE(double dE, double T)
{
  if(dE > 0) return exp(-dE/T);
  else return 1;
}


void print_final_config()
{
  ofstream fout2;
  
  fout2.open("final_config.txt");

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      {
	if(config[i][j] == 1)
	  {
	    fout2<<i<<"  "<<j<<endl;
	  }
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
      fout2<<"### angle = "<<i*d_theta<<endl;

      for(int r=0; r<(2*Lmax+1); r++)
	{
	  fout2<<(r-Lmax)<<" "<<projection[i][r] + diff_projection[i][r]<<endl;
	}

      fout2<<"#############################"<<endl;
      fout2<<endl<<endl<<endl;
    }

  fout2.close();
}


void print_final_diff_projections()
{
  ofstream fout;

  fout.open("final_diff_projection.txt");

  for(int i=0; i<N_p; i++)
    {
      for(int r=0; r<(2*Lmax+1); r++)
	{
	  fout<<i<<" "<<r<<" "<<diff_projection[i][r]<<" "<<diff_projection[i][r]*diff_projection[i][r]<<endl;
	}
    }

  fout.close();
}

int main()
{
  //initiate all the matrix used in the program
  //*******************************************
  int i,j;

  for(i=0; i<N_p; i++)
    for(j=0; j<(2*Lmax+1); j++)
      {
	projection[i][j] = 0;
	temp_projection[i][j] = 0;
	diff_projection[i][j] = 0;
      }

  for(i=0; i<N_p; i++)
    for(j=0; j<2; j++)
      {
	changed_location[i][j] = 0;
	changed_projection[i][j] = 0;
      }

  for(i=0; i<MAXX; i++)
    for(j=0; j<MAXY; j++)
      {
	config[i][j] = 0;
      }
  //*******************************************

  srand(time(NULL));
  
  read_data();

  get_projections();

  get_black_pixels();

  print_projections();

  init_reconstr();

  int N_acc = 0; //acceptance rate...
  double dE;

  cout<<"The Initial Energy = "<<global_energy<<endl;

  FILE* fp = fopen("E.txt","w");
  fclose(fp);

  cout<<"Staring the simulated annealing reconstruction process..."<<endl;
  for(int q=0; q<TN; q++)
    {
      T = alpha*T;

      N_acc = 0;

      cout<<"Stage "<<q+1<<" with T = "<<T<<endl;
       
      for(int i=0; i< Nevl; i++)
	{
	  change_config();

	  get_projection_difference();
		  
	  update_diff_projection();

	  dE = get_energy();
	  
	  double P = (double)(rand()%MAXS)/(double)MAXS;

	  if( P > PE(dE, T))
	    {
	      resume_config();
	      
	      resume_diff_projection();
	      
	    } 
	  else 
	    {
	      N_acc ++;

	      global_energy += dE;
	    }
	  	  
	}

      printf("%d th change of temperature has finished... \n",q+1 );
      cout<<"The acceptance rate: "<<(double)N_acc/(double)Nevl<<endl;
      cout<<"The energy E = "<<global_energy<<endl;

      
      fp = fopen("E.txt","a");
      fprintf(fp, "%e\n", global_energy);
      fclose(fp);

      printf("*************************************************\n");
	
    }
  
  print_final_config();

  print_final_projections();

  print_final_diff_projections();
  
}



