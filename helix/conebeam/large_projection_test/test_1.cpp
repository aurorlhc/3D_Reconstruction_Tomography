//real diameter of each sphere is 0.02", which is approximately 508 micron
//construct spacial distribution of the three spheres using projections from different angle
//use the constructed structure to generate projections from different angle
//and then compare to the real projection data to see if there exists inconsistency

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
using namespace std;

#define MOD 10000 //for random number generator

#define MAXX 136
#define MAXY 136
#define MAXZ 152
#define Lmax 128
#define N_p 30   //number of projections
#define PI 3.1415926
#define L_1 1274
#define L_2 866
#define R 55    //the radius of the column for the initialization domain
#define dd 0.58  //radius of the inscribed circle for each voxel in the domain

//initial volume fraction of spheres and clay
#define IPC 0.96

//phase transition probability
#define P1 0.8 //probability of phase 1 to phase 2, then the probability of phase 1 to phase 3 is (1 - P1)
#define P2 0.6 //probability of phase 2 to phase 1, then the probability of phase 2 to phase 3 is (1 - P2)
#define P3 0.5 //probability of phase 3 to phase 1, then the probability of phase 3 to phase 2 is (1 - P3)


double config[MAXX][MAXY][MAXZ];
double initial_projection[N_p][2*Lmax][2*Lmax]; //initialization
double projection[N_p][2*Lmax][2*Lmax];         //target projections
double diff_projection[N_p][2*Lmax][2*Lmax];

double global_energy;

double mu_1 = 0.07835;
double mu_2 = 0.00223;

double changed_mu;

double x_c = (double)(MAXX-1)/2.0;
double y_c = (double)(MAXY-1)/2.0;
double z_c = (double)(MAXZ-1)/2.0;

double d_theta = PI/(double)N_p;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//for simulated annealing procedure

int indexi, indexj, indexk;
int temp;

int Nevl = 5000000;

double alpha = 0.98;

int TN = 400;

double T = 0.001;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
double abs(double x)
{
  if(x >= 0) return x;
  else return -x;
}

struct point
{
  double x, y, z;
};

struct line
{
  point p1, p2;
};


int init_point(point* a, double x, double y, double z)
{
  a->x = x;
  a->y = y;
  a->z = z;

  return 1;
}


int init_line(line* l, int x, int y, int z, int a, int b, int c)
{
  l->p1.x = x;
  l->p1.y = y;
  l->p1.z = z;

  l->p2.x = a;
  l->p2.y = b;
  l->p2.z = c;

  return 1;
}

int point2line(line l, point p, point *Q)
{
  double a,b,c;
  double A,B,C;

  a = l.p2.x - l.p1.x;
  b = l.p2.y - l.p1.y;
  c = l.p2.z - l.p1.z;
 
  A = a*p.x + b*p.y + c*p.z;
  B = b*l.p1.x - a*l.p1.y;
  C = c*l.p1.x - a*l.p1.z;
 
  if (a!=0)
    {
      Q->x = (A*a + B*b + C*c)/(a*a + b*b + c*c);
      Q->y = (b*Q->x - B)/a;
      Q->z = (c*Q->x - C)/a;
    }
  else
    {
      double D, temp;

      D = c*l.p1.y - b*l.p1.z;
      temp = b*b + c*c;
      Q->y = (A*b + D*c)/temp;
      Q->z = (A*c - D*b)/temp;
      Q->x = (B + a*Q->y)/b;
    }
  
  return 1;
}

double p2p(point a, point b)
{
  double dx, dy, dz;

  dx = a.x - b.x;
  dy = a.y - b.y;
  dz = a.z - b.z;

  return sqrt(dx*dx + dy*dy + dz*dz);
}
*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct node
{
  int projection_index;
  int detector_bin_index1;
  int detector_bin_index2;
  double length;
  node* next;
};

node* struct_array[MAXX][MAXY][MAXZ];

void initiate_struct_array()
{
  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
	struct_array[i][j][k] = NULL;
}

void update_struct_array(int voxel_index1, int voxel_index2, int voxel_index3, int projection_index, int bin_index1, int bin_index2, double length)
{
  node* pt = new node;

  pt->projection_index = projection_index;
  pt->detector_bin_index1 = bin_index1;
  pt->detector_bin_index2 = bin_index2;
  pt->length = length;

  node* pt2 = struct_array[voxel_index1][voxel_index2][voxel_index3];

  pt->next = pt2;
  struct_array[voxel_index1][voxel_index2][voxel_index3] = pt;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void initiate_config()
{
  //int range = 10000;
  double p;

  for(int i=0; i<MAXX; i++)
    for(int j=0; j<MAXY; j++)
      for(int k=0; k<MAXZ; k++)
	{
	  if((i-x_c)*(i-x_c)+(j-y_c)*(j-y_c)>R*R)
	    {
	      config[i][j][k] = 0;
	    }
	  else
	    {
	      p = (double)(rand()%MOD)/(double)MOD;
	      
	      if(p<IPC) //clay
		config[i][j][k] = 50;
	      else      //sphere
		config[i][j][k] = 255;
	    }
	}
}

void initiate_projections()
{
  for(int i=0; i<N_p; i++)
    for(int j=0; j<2*Lmax; j++)
      for(int k=0; k<2*Lmax; k++)
	{
	  initial_projection[i][j][k] = 0;
	}

  for(int i=0; i<N_p; i++)
    {
      for(int j=0; j<2*Lmax; j++)
	for(int k=0; k<2*Lmax; k++)
	  {
	    double A = j - Lmax;  //x
	    double B = L_1 + L_2; //y
	    double C = k - Lmax;  //z
	    double xu, xl, yu, yl, zu, zl;
	    double ra = 0.5*sqrt(MAXX*MAXX + MAXY*MAXY);

	    if(B*B*L_1*L_1 > (A*A + B*B)*(L_1*L_1 - ra*ra))
	      {
		yu = (B*B*L_1 + sqrt(B*B*B*B*L_1*L_1 - B*B*(A*A + B*B)*(L_1*L_1 - ra*ra)))/(A*A + B*B);
		yl = (B*B*L_1 - sqrt(B*B*B*B*L_1*L_1 - B*B*(A*A + B*B)*(L_1*L_1 - ra*ra)))/(A*A + B*B);
		xu = A*yu/B;
		xl = A*yl/B;
		zu = C*yu/B;
		zl = C*yl/B;

		xu = ceil(xu);
		xl = floor(xl);
		yu = ceil(yu);
		yl = floor(yl);
		zu = ceil(zu);
		zl = floor(zl);

		double m = A;
		double n = B;
		double p = C;

		for(int a=0; a<MAXX; a++)
		  for(int b=0; b<MAXY; b++)
		    for(int c=0; c<MAXZ; c++)
		      {
			double temp_x = (a-x_c)*cos(i*d_theta)+(b-y_c)*sin(i*d_theta);
			double temp_y = (b-y_c)*cos(i*d_theta)-(a-x_c)*sin(i*d_theta) + L_1;
			double temp_z = c-z_c;

			if(temp_x>=xl && temp_x<=xu && temp_y>=yl && temp_y<=yu && temp_z>=zl && temp_z<=zu)
			  {
			    double t = (m*temp_x + n*temp_y + p*temp_z)/(m*m + n*n + p*p);

			    double xc = m*t;
			    double yc = n*t;
			    double zc = p*t;

			    double d = sqrt((temp_x - xc)*(temp_x - xc) + (temp_y - yc)*(temp_y - yc) + (temp_z - zc)*(temp_z - zc));

			    if(d<=dd)
			      {
				double h = 2*sqrt(dd*dd - d*d); //segment length in the voxel

				update_struct_array(a, b, c, i, j, k, h);

				if(config[a][b][b] == 255)
				  initial_projection[i][j][k] += mu_1*h;
				else if(config[a][b][c] == 50)
				  initial_projection[i][j][k] += mu_2*h;
				else
				  initial_projection[i][j][k] += 0;
				//~~~~~~~~~~~~~~~~~...To be continued...~~~~~~~~~~~~~~~~~~~~~
			      }
			  }
		      }
	      }
	  }
    }      
  cout<<"finish computing the projections..."<<endl;
}

void read_projections()
{
  double value;

  ifstream fin;
  fin.open("projections_256x256x30.txt");

  for(int i=0; i<N_p; i++)
    for(int k=2*Lmax-1; k>=0; k--)
      for(int j=0; j<2*Lmax; j++)
	{
	  fin>>value;
	  projection[i][j][k] = value;
	}

  fin.close();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //calculate energy for the first time
  double sum = 0.0;

  ofstream fout;
  fout.open("initial_diff.txt");

  for(int i=0; i<N_p; i++)
    for(int j=0; j<2*Lmax; j++)
      for(int k=0; k<2*Lmax; k++)
	{
	  diff_projection[i][j][k] = initial_projection[i][j][k] - projection[i][j][k];
	  fout<<i<<" "<<j<<" "<<k<<" "<<diff_projection[i][j][k];
	  sum += diff_projection[i][j][k]*diff_projection[i][j][k];
	}
  fout.close();

  global_energy = sum;
  cout<<"Initial_energy = "<<global_energy<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}

void change_config()
{
  int i, j, k;

  i = rand()%MAXX;
  j = rand()%MAXY; 
  k = rand()%MAXZ;

  indexi = i;
  indexj = j;
  indexk = k;

  temp = config[i][j][k];
}

void resume_config()
{
  config[indexi][indexj][indexk] = temp;
}

void get_projection_difference()
{
  int i = indexi;
  int j = indexj;
  int k = indexk;

  double p = (double)(rand()%MOD)/(double)MOD;

  if(config[i][j][k] == 0)
    {
      if(p<P3)
	{
	  config[i][j][k] = 255;
	  changed_mu = mu_1 - 0;
	}
      else
	{
	  config[i][j][k] = 50;
	  changed_mu = mu_2 - 0;
	}
    }
  else if(config[i][j][k] == 50)
    {
      if(p<P2)
	{
	  config[i][j][k] = 255;
	  changed_mu = mu_1 - mu_2;
	}
      else
	{
	  config[i][j][k] = 0;
	  changed_mu = 0 - mu_2;
	}
    }
  else
    {      
      if(p<P1)
	{
	  config[i][j][k] = 50;
	  changed_mu = mu_2 - mu_1;
	}
      else
	{
	  config[i][j][k] = 0;
	  changed_mu = 0 - mu_1;
	}
    }
}

double update_diff_projection()
{
  int i = indexi;
  int j = indexj;
  int k = indexk;

  double E_o = 0.0;
  double E_n = 0.0;

  node* pt = struct_array[i][j][k];

  while(pt!=NULL)
    {
      int m,n,l;
      double length;

      m = pt->projection_index;
      n = pt->detector_bin_index1;
      l = pt->detector_bin_index2;
      length = pt->length;

      E_o += diff_projection[m][n][l]*diff_projection[m][n][l];

      diff_projection[m][n][l] += changed_mu*length;

      E_n += diff_projection[m][n][l]*diff_projection[m][n][l];

      pt = pt->next;
    }

  return (E_n - E_o);
}

void resume_diff_projection()
{
  int i = indexi;
  int j = indexj;
  int k = indexk;

  node* pt = struct_array[i][j][k];

  while(pt!=NULL)
    {
      int m,n,l;
      double length;

      m = pt->projection_index;
      n = pt->detector_bin_index1;
      l = pt->detector_bin_index2;
      length = pt->length;

      diff_projection[m][n][l] -= changed_mu*length;

      pt = pt->next;
    }
}

double p_acc(double dE, double T)
{
  if(dE<0) return 1.0;
  else return exp(-dE/T);
}


void print_final_config()
{
  ofstream fout2;
  
  fout2.open("final_config.vtk");

  fout2<<"# vtk DataFile Version 2.0"<<endl;
  fout2<<"2D_to_3D example"<<endl;
  fout2<<"ASCII"<<endl;
  fout2<<"DATASET STRUCTURED_POINTS"<<endl;
  fout2<<"DIMENSIONS 136 136 152"<<endl;
  fout2<<"SPACING 1 1 1"<<endl;
  fout2<<"ORIGIN 0 0 0"<<endl;
  fout2<<"POINT_DATA 2811392"<<endl;
  fout2<<"SCALARS volume_scalars UNSIGNED_CHAR 1"<<endl;
  fout2<<"LOOKUP_TABLE default"<<endl;


  for(int k=0; k<MAXZ; k++)
    for(int j=0; j<MAXY; j++)
      for(int i=0; i<MAXX; i++)
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
      fout2<<"### angle = "<<i*d_theta<<endl;

      for(int r=0; r<2*Lmax; r++)
	for(int p=0; p<2*Lmax; p++)
	  {
	    fout2<<(r-Lmax)<<" "<<(p-Lmax)<<" "<<projection[i][r][p] + diff_projection[i][r][p]<<endl;
	  }

      fout2<<"#############################"<<endl;
      fout2<<endl<<endl<<endl;
    }

  fout2.close();
}


int main()
{
  srand(time(NULL));
  clock_t t_1,t_2;
  initiate_struct_array();

  initiate_config();
  t_1 = clock();
  initiate_projections();
  t_2 = clock();
  read_projections();

  int N_acc = 0;

  for(int q=0; q<TN; q++)
    {
      T = alpha*T;

      N_acc = 0;

      cout<<"Stage "<<q+1<<" with T = "<<T<<endl;

      for(int i=0; i<Nevl; i++)
	{
	  double dE = 0;

	  change_config();

	  get_projection_difference();

	  dE = update_diff_projection();

	  double p_temp = (double)(rand()%MOD)/(double)MOD;

	  if(p_temp<p_acc(dE,T))
	    {
	      global_energy += dE;
	      N_acc ++;
	    }
	  else
	    {
	      resume_diff_projection();
	      resume_config();
	    }

	}

      cout<<"The acceptance rate: "<<(double)N_acc/(double)Nevl<<endl;
      cout<<"The energy E = "<<global_energy<<endl;
      cout<<"*************************************"<<endl;
    }

  print_final_config();
  print_final_projections();

  double diff = (double)(t_2 - t_1)/CLOCKS_PER_SEC;
  cout<<"Running time is "<<diff<<endl;

}
