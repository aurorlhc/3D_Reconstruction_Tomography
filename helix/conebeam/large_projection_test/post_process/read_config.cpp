#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
#define MAX 10000
#define MAXX 136
#define MAXY 136
#define MAXZ 152
#define r 10
#define R 62

double x_c = MAXX/2;
double y_c = MAXY/2;
double z_c = MAXZ/2;

int config[MAXX][MAXY][MAXZ];

void read_config()
{
  ifstream fin;
  fin.open("final_config.txt");

  int value;

  for(int k=0; k<MAXZ; k++)
    for(int j=0; j<MAXY; j++)
      for(int i=0; i<MAXX; i++)
	{
	  fin>>value;
	  config[i][j][k] = value;
	}

  fin.close();
}

void processing_image()
{
  for(int i=1; i<MAXX-1; i++)
    for(int j=1; j<MAXY-1; j++)
      for(int k=1; k<MAXZ-1; k++)
	{
	  int counter = 0;

	  for(int m=-1; m<=1; m++)
	    for(int n=-1; n<=1; n++)
	      for(int l=-1; l<=1; l++)
		{
		  if(config[i+m][j+n][k+l] == 255)
		    {
		      counter++;
		    }
		}

	  if(counter <= 6)
	    {
	      /*
	      double p = (double)(rand()%MAX)/(double)MAX;

	      if(p<0.5)
		config[i][j][k] = 50;
	      else
		config[i][j][k] = 0;
	      */

	      if((i - x_c)*(i - x_c) + (j - y_c)*(j - y_c) <= R*R)
		config[i][j][k] = 50;
	      else
		config[i][j][k] = 0;

	    }
	}

    for(int k=0; k<MAXZ; k++)
      for(int j=0; j<MAXY; j++)
	for(int i=0; i<MAXX; i++)
	  {
	    if((i - x_c)*(i - x_c) + (j - y_c)*(j - y_c) <= r*r)
	      {
		config[i][j][k] = 50;
	      }
	  }

    for(int k=74, k<=77; k++)
      {


}

void output()
{
  ofstream fout2;
  fout2.open("after_processing.vtk");

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

int main()
{
  srand(time(NULL));
  read_config();
  processing_image();
  output();
}
