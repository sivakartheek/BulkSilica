#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "angular.h"
#include "ewald.h"
#include<vector>


double mindis(double dx,double dy,double dz,double a,double b,double c)
{
  return norm(dx - a * round(dx/a),dy - b * round(dy/b),dz - c * round(dz/c));
}

int Rand_INT(int min,int max) //[min,max)
{
  return rand()%(max-min)+min;
}

double Rand_DOUBLE() //(0,1)
{
  return (double) rand() / (double) RAND_MAX;
}

double min(double a,double b)
{
  if(a>b)
    {
      return b;
    }
  else
    {
      return a;
    }
}


double angle_btwn_3points(vector<double> const &r,int i,int j1,int j2,double a,double b,double c)
{
  double x1,x2;
  double y1,y2;
  double z1,z2;

  x1 = r[j1+0] - r[i];
  x2 = r[i+0] - r[j2];
  y1 = r[j1+1] - r[i+1];
  y2 = r[i+1] - r[j2+1];
  z1 = r[j1+2] - r[i+2];
  z2 = r[i+2] - r[j2+2]; 

  x1 = x1 - a*round(x1/a);
  x2 = x2 - a*round(x2/a);
  y1 = y1 - b*round(y1/b);
  y2 = y2 - b*round(y2/b);
  z1 = z1 - c*round(z1/c);
  z2 = z2 - c*round(z2/c);

  double top =  (x1*x2+y1*y2+z1*z2) - 0.0001 ;
  double bot =  norm(x1,y1,z1) * norm(x2,y2,z2);
  if( top == bot )
    {
      return 180; 
    }
  else
    {
      return  180 - (acos(top/bot) * 57.296);
    }
}


void angular_distribution_SiOSi(double L, double **siosi, int a_size, vector<double> const &r, double q1, double q2)
{
     const int end  = r.size();
     for(int i = 0 ; i < end;i = i + 4)
      {
        if(r[i+3] == q2)
        {
          int j1 = -1,j2 = -1;
          double min = 1000;
          for(int j = 0 ; j < end;j = j + 4)
          {
            if(i == j){}
            else if(r[j+3] == q1)
            {
              double rij = mindis(r[i]-r[j],r[i+1]-r[j+1],r[i+2]-r[j+2],L,L,L);
              if(rij < min && rij < 2.3)
              {
                min = rij;
                j1 = j;
              }
            }
          }

          //second oxygen
          min = 1000;
          for(int j = 0 ; j < end;j = j + 4)
          {
            if(i == j || j == j1){}
            else if(r[j+3] == q1)
            {
              double rij = mindis(r[i]-r[j],r[i+1]-r[j+1],r[i+2]-r[j+2],L,L,L);
              if(rij < min && rij < 2.3 )
              {
                min = rij;
                j2 = j;
              }
            }
          }

          double angle  = angle_btwn_3points(r,i,j1,j2,L,L,L);

          for(int kk = 0;kk < a_size;kk++)
          {
            if(angle <= siosi[0][kk] + 0.25 && angle > siosi[0][kk] - 0.25)
            {
              siosi[1][kk] += 1;
            }
          } 
        }
      }
}

void angular_distribution_OSiO( double L, double **osio, int a_size, vector<double> const &r, double q1, double q2)
{
    const int end            = r.size();
    for(int i = 0 ; i < end;i = i + 4)
    {
      if(r[i+3] == q1 )
      {
        int j1 = -1,j2 = -1,j3 = -1,j4 = -1;
        double min = 1000;
        for(int j = 0 ; j < end;j = j + 4)
        {
          if(i == j){}
          else if(r[j+3] == q2)
          {
            double rij = mindis(r[i]-r[j],r[i+1]-r[j+1],r[i+2]-r[j+2],L,L,L);
            if(rij < min && rij < 2.5)
            {
              min = rij;
              j1 = j;
            }
          }
        }

        //second oxygen
        min = 1000;
        for(int j = 0 ; j < end;j = j + 4)
        {
          if(i == j || j == j1){}
          else if(r[j+3] == q2)
          {
            double rij = mindis(r[i]-r[j],r[i+1]-r[j+1],r[i+2]-r[j+2],L,L,L);
            if(rij < min && rij < 2.5 )
            {
              min = rij;
              j2 = j;
            }
          }
        }

        //third oxygen
        min = 1000;
        for(int j = 0 ; j < end;j = j + 4)
        {
          if(i == j || j == j1 || j == j2){}
          else if(r[j+3] == q2)
          {
            double rij = mindis(r[i]-r[j],r[i+1]-r[j+1],r[i+2]-r[j+2],L,L,L);
            if(rij < min && rij < 2.5 )
            {
              min = rij;
              j3 = j;
            }
          }
        }

        // Fourth oxygen
        min = 1000;
        for(int j = 0 ; j < end;j = j + 4)
        {
          if(i == j || j == j1 || j == j2 || j == j3){}
          else if(r[j+3] == q2)
          {
            double rij = mindis(r[i]-r[j],r[i+1]-r[j+1],r[i+2]-r[j+2],L,L,L);
            if(rij < min && rij < 2.5 )
            {
              min = rij;
              j4 = j;
            }
          }
        }

  double angle1  = angle_btwn_3points(r,i,j1,j2,L,L,L);
  double angle2  = angle_btwn_3points(r,i,j1,j3,L,L,L);
  double angle3  = angle_btwn_3points(r,i,j1,j4,L,L,L);
  double angle4  = angle_btwn_3points(r,i,j2,j3,L,L,L);
  double angle5  = angle_btwn_3points(r,i,j2,j4,L,L,L);
  double angle6  = angle_btwn_3points(r,i,j3,j4,L,L,L);

  for(int kk = 0;kk < a_size;kk++)
  {
    if(angle1 < osio[0][kk] + 0.25 && angle1 > osio[0][kk] - 0.25)
    {
      osio[1][kk] += 1;
    }
    if(angle2 < osio[0][kk] + 0.25 && angle2 > osio[0][kk] - 0.25)
    {
      osio[1][kk] += 1;
    }
    if(angle3 < osio[0][kk] + 0.25 && angle3 > osio[0][kk] - 0.25)
    {
      osio[1][kk] += 1;
    }
    if(angle4 < osio[0][kk] + 0.25 && angle4 > osio[0][kk] - 0.25)
    {
      osio[1][kk] += 1;
    }
    if(angle5 < osio[0][kk] + 0.25 && angle5 > osio[0][kk] - 0.25)
    {
      osio[1][kk] += 1;
    }
    if(angle6 < osio[0][kk] + 0.25 && angle6 > osio[0][kk] - 0.25)
    {
      osio[1][kk] += 1;
    }
  } 

  }
  }
}
