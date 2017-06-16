#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#define MAXX 100 //system size
//#define MAXY 100
//#define MAXZ 100
#define NP 446050
#define Nt MAXX/2

int config[MAXX][MAXX][MAXX];
int S2_x[Nt];
int S2_y[Nt];
int S2_z[Nt];
double obj[Nt];

void read_config()
{
  int value;
  ifstream fin;
  fin.open("final_config.txt");
 
  for(int k=0; k<MAXX; k++)
    for(int j=0; j<MAXX; j++)
      for(int i=0; i<MAXX; i++)
	{
	  fin>>value;
	  config[i][j][k] = value;
	}

  fin.close();
}


/*
void print_config()
{
  ofstream fout2;
  fout2.open("hardsphere_config.vtk");

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

  for(int k=0; k<MAXZ; k++)
    for(int j=0; j<MAXY; j++)
      for(int i=0; i<MAXX; i++)
	{
	  fout2<<config[i][j][k]<<endl;
	}

  fout2.close();
}
*/

void calculate_S2()
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

  ofstream fout;
  fout.open("final_S2.txt");
  for(int i=0; i<Nt; i++)
    {
      fout<<i<<" "<<obj[i]<<endl;
    }
  fout.close();
}

int main()
{
  read_config();
  //print_config();
  calculate_S2();
}

