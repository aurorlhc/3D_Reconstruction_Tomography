#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#define MAX 77

int config[MAX][MAX][MAX];
int final_config[MAX][MAX][MAX];

void read_config()
{
  int value;

  ifstream fin;
  fin.open("target_config.txt");

  for(int k=0; k<MAX; k++)
    for(int j=0; j<MAX; j++)
      for(int i=0; i<MAX; i++)
	{
	  fin>>value;
	  config[i][j][k] = value;
	}

  fin.close();

  for(int k=0; k<MAX; k++)
    for(int j=0; j<MAX; j++)
      for(int i=0; i<MAX; i++)
	{
	  if(config[i][j][k] == 50)
	    config[i][j][k] = 0;
	}

  ofstream fout;
  fout.open("binary_config.vtk");

  fout<<"# vtk DataFile Version 2.0"<<endl;
  fout<<"2D_to_3D example"<<endl;
  fout<<"ASCII"<<endl;
  fout<<"DATASET STRUCTURED_POINTS"<<endl;
  fout<<"DIMENSIONS 77 77 77"<<endl;
  fout<<"SPACING 1 1 1"<<endl;
  fout<<"ORIGIN 0 0 0"<<endl;
  fout<<"POINT_DATA 456533"<<endl;
  fout<<"SCALARS volume_scalars UNSIGNED_CHAR 1"<<endl;
  fout<<"LOOKUP_TABLE default"<<endl;

  for(int k=0; k<MAX; k++)
    for(int j=0; j<MAX; j++)
      for(int i=0; i<MAX; i++)
	{
	  fout<<config[i][j][k]<<endl;
	}

  fout.close();
}

void generate_spheres()
{
  int coor[7][3] = {{17, 26, 2},
		    {23, 16, 14},
		    {34, 12, 25},
		    {44, 14, 36},
		    {53, 18, 47},
		    {59, 27, 60},
		    {61, 39, 68}};

  for(int i=0; i<MAX; i++)
    for(int j=0; j<MAX; j++)
      for(int k=0; k<MAX; k++)
	{
	  for(int m=0; m<7; m++)
	    {
	      if((i-coor[m][0])*(i-coor[m][0]) + (j-coor[m][1])*(j-coor[m][1]) + (k-coor[m][2])*(k-coor[m][2]) <= 25)
		{
		  final_config[i][j][k] = 255;
		}
	    }
	}

  config[61][44][74] = 255;

  for(int i=0; i<MAX; i++)
    for(int j=0; j<MAX; j++)
      {
	if((i-61)*(i-61) + (j-44)*(j-44) <= 9)
	  final_config[i][j][MAX-1] = 255;
      }

  ofstream fout1;
  fout1.open("final_config.vtk");

  fout1<<"# vtk DataFile Version 2.0"<<endl;
  fout1<<"2D_to_3D example"<<endl;
  fout1<<"ASCII"<<endl;
  fout1<<"DATASET STRUCTURED_POINTS"<<endl;
  fout1<<"DIMENSIONS 77 77 77"<<endl;
  fout1<<"SPACING 1 1 1"<<endl;
  fout1<<"ORIGIN 0 0 0"<<endl;
  fout1<<"POINT_DATA 456533"<<endl;
  fout1<<"SCALARS volume_scalars UNSIGNED_CHAR 1"<<endl;
  fout1<<"LOOKUP_TABLE default"<<endl;

  for(int k=0; k<MAX; k++)
    for(int j=0; j<MAX; j++)
      for(int i=0; i<MAX; i++)
	{
	  fout1<<final_config[i][j][k]<<endl;
	}

  fout1.close();
}


int main()
{
  read_config();
  generate_spheres();
}

