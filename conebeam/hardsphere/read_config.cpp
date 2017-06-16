#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#define MAXX 100 //system size
#define MAXY 100
#define MAXZ 100
#define NP 446050

int config[MAXX][MAXY][MAXZ];

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
	  config[i][j][k] = 0;
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

int main()
{
  read_config();
  print_config();
}

