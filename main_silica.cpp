// @main.cpp determines the structural and thermodynamic properties 
// of silica glass material for given temperature and pressure
// conditions. The PBC is used to make the limited number of atoms
// as a bulk. The interatomic interactions were modeled using the 
// Coulomb and Buckingham potential.

// @author - Naveen Kumar Kaliannan

#include<iostream>
#include<cmath>
#include<time.h>
#include<fstream>
#include<stdlib.h>
#include<vector>
#include"ewald.h"
#include<mpi.h>

using namespace std;

const double Pi = 3.141593;
const double Unit_Coulomb = 1389.55;
const double Unit_eV = 96.49930;
const double R        = 0.0083144598; // Gas constant
const double Temp     = 2750;         //Temperature of the system
const double Pr       = 0.6023;      //Pressure of the system


// generates a random  integer between [min,max)
int Rand_INT(int min,int max) //[min,max)
{
  return rand()%(max-min)+min;
}

// generates the random number between 0 and 1
double Rand_DOUBLE() //(0,1)
{
  return (double) rand() / (double) RAND_MAX;
}

// finds the minimum (a,b)
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

// Maximum of 3 elements
template <class T>
T max3(T a, T b, T c)
{
 return (a < c && b < c) ? c : (a < b) ? b : a;
}

// Short-range interactions - BKS
double BKS(vector<double> const &r, double L, double q1, double q2, int size, int rank)
{
  double u = 0, rij  = 0;  
  // MPI variables
  int from = 0, to = 0;
  from = rank * r.size() / size ;
  to = (rank + 1) * r.size() / size ;
  if(rank + 1 == size)
  {
    to = to - 4;
  }

  for(unsigned int i = from;i < to ;i += 4)
    {
      for(unsigned int j = i+4;j < r.size();j += 4)
        {
          rij = sqrt( pow((r[i+0]-r[j+0]) - L * round((r[i+0]-r[j+0])/L),2) 
                    + pow((r[i+1]-r[j+1]) - L * round((r[i+1]-r[j+1])/L),2) 
                    + pow((r[i+2]-r[j+2]) - L * round((r[i+2]-r[j+2])/L),2) );
          // Potential truncation
          if(rij < 5)
            {
              if(r[i+3] == r[j+3] && r[i+3] == q2)//oxygen - oxygen
                {
                  u += (134015.6224*exp(-rij*2.7600)) - (16887.3775/pow(rij,6));
                }
              if(r[i+3] != r[j+3]) //silicon - oxygen
                { 
                  u += (1737349.967*exp(-rij*4.87318)) - (12599.51795/pow(rij,6));
                }
            }
        }
    }
  return u;
}

// Short-range interactions - TCD potential
double TCD(vector<double> const &r, double L, double q1, double q2, int size, int rank)
{
  double u = 0, rij  = 0;  
  // MPI variables
  int from = 0, to = 0;
  from = rank * r.size() / size ;
  to = (rank + 1) * r.size() / size ;
  if(rank + 1 == size)
  {
    to = to - 4;
  }

  for(unsigned int i = from;i < to ;i += 4)
    {
      for(unsigned int j = i+4;j < r.size();j += 4)
        {
          rij = sqrt( pow((r[i+0]-r[j+0]) - L * round((r[i+0]-r[j+0])/L),2) 
                    + pow((r[i+1]-r[j+1]) - L * round((r[i+1]-r[j+1])/L),2) 
                    + pow((r[i+2]-r[j+2]) - L * round((r[i+2]-r[j+2])/L),2) );
          // Potential truncation
          if(rij < 5)
            {
              if(r[i+3] == r[j+3] && r[i+3] == q2)//oxygen - oxygen
                {
                  u += (Unit_eV * 2029.2204*exp(-rij/0.343645)) - (Unit_eV * 192.58/pow(rij,6));
                }
              if(r[i+3] != r[j+3]) //silicon - oxygen
                { 
                  u += (Unit_eV * 13702.905*exp(-rij/0.193817)) - (Unit_eV * 54.681/pow(rij,6));
                }
            }
        }
    }
  return u;
}



// Main implementation
int main(int argc, char** argv)
{

  int size, rank, root = 0; 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Variables
  double L = 0, V = 0, V_new = 0;
  unsigned int nx = 0, ny = 0, nz = 0, mean = 0 ;
  double r0 = 0;
  unsigned int N = 0, N_ = 0;
  unsigned int per = 0;
  double q1 = 0 , q2 = 0;
  vector<double> r;
  vector<double> r_new;

  // Properties variables
  //Radial Distribution Function(RDF)
  unsigned int RDF_size = 600;
  double RDF_h = 0.01;
  double* rad = new double[RDF_size];
  double* RDF_all = new double[RDF_size];
  double* RDF_Si = new double[RDF_size];
  double* RDF_O  = new double[RDF_size];
  double* RDF_SiO =  new double[RDF_size];
  for(unsigned int i = 1;i < RDF_size ;i++)
    {
      rad[i]          = i * 0.01;
      RDF_all[i]    = 0;
      RDF_Si[i]     = 0;
      RDF_O [i]     = 0;
      RDF_SiO[i]    = 0;
    }
  ofstream outfile1;
  ofstream outfile2;

  // Processor 0 - input/output processing
  if(rank == root) 
    {
      double Lx = 0, Ly = 0, Lz = 0;
      ifstream infile(argv[1]);
      infile >> nx >> ny >> nz ;
      infile >> r0 ;
      infile >> N ;
      infile >> per ;
      per *= nx * ny * nz;
      infile >> q1 >> q2 ;
      infile >> Lx >> Ly >> Lz;
      N_ = N * 4 * nx * ny * nz ;
      r.resize(N*4*nx*ny*nz);
      r_new.resize(N*4*nx*ny*nz);

      for(unsigned int i = 0; i < (N*4); i += 4)
        {
          infile >> r[i+0] >> r[i+1] >> r[i+2] >> r[i+3];
          r[i+0] *= Lx;  r[i+1] *= Ly; r[i+2] *= Lz; 
        }
      infile.close();
      infile.clear();

      unsigned int i = 0;
      for(unsigned int x = 0; x < nx; ++x)
        {
          for(unsigned int y = 0; y < ny; ++y)
            {
              for(unsigned int z = 0; z < nz; ++z)
                {
                  for(unsigned int j = 0;j < (N*4); j += 4)
                    {
                      r[i+0] = Lx * x + r[j+0];
                      r[i+1] = Ly * y + r[j+1];
                      r[i+2] = Lz * z + r[j+2];
                      r[i+3] = r[j+3];
                      //cout << i << " " <<  r[i+0] << " " << r[i+1] << " " << r[i+2] <<" " << r[i+3] << endl;
                      i += 4;
                    }
                }
            }  
        }
      r_new = r;
      Lx *= nx;
      Ly *= ny;
      Lz *= nz;
      L = max3(Lx,Ly,Lz);
      V = pow(L,3);
      V_new = V;
      outfile1.open(argv[2]);
      outfile2.open(argv[3]);
    }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&N_, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  r.resize(N_);
  r_new.resize(N_);
  MPI_Bcast(&q1, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&q2, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&V, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&L, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&r[0], r.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);

  double U_local = 0, U_global = 0, kappa = 0.2150, u = 0, u_new = 0;
  U_local = ( RealandReciprocalSpace(r, L, L, L, kappa, 1, size, rank) 
             + kappa * PointEnergy(r, size, rank) / sqrt(Pi) ) * Unit_Coulomb + TCD(r, L, q1, q2, size, rank);
  MPI_Reduce(&U_local, &u, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
  u += (Unit_Coulomb *  2 * Pi * pow(Dipole(r),2) / (3 * V)) ;

  srand(time(NULL));
  for(int k = 0;k < 5000000; k++)
    {
      if(rank == root)
        {
          //1. randomly select an atom
          unsigned int P = Rand_INT(0, r.size()/4) * 4;

          //2. move the selected atom randomly
          r_new[P+0] += 0.1*(2*Rand_DOUBLE()-1); 
          r_new[P+1] += 0.1*(2*Rand_DOUBLE()-1);
          r_new[P+2] += 0.1*(2*Rand_DOUBLE()-1);
          
          //3. Volume change
          V_new += 10 * (2*Rand_DOUBLE()-1);
          L = pow(V_new,0.333333);
          double ratio = pow(V_new,0.3333333)/pow(V,0.3333333);
          for(unsigned int i = 0 ; i < r.size();i = i + 4)
            {
              r_new[i+0] = r_new[i+0]*ratio;
              r_new[i+1] = r_new[i+1]*ratio;
              r_new[i+2] = r_new[i+2]*ratio;
            }
          U_global = 0;
        }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&L, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&r_new[0], r_new.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&V_new, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

      //4. Calculate the potential energy of the configurations
      U_global = 0;
      U_local = ( RealandReciprocalSpace(r_new, L, L, L, kappa, 1, size, rank) 
             + kappa * PointEnergy(r_new, size, rank) / sqrt(Pi) ) * Unit_Coulomb + TCD(r_new, L, q1, q2, size, rank);

      MPI_Reduce(&U_local, &U_global, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
     if (rank == root)
       {
         u_new = U_global + (Unit_Coulomb *  2.0 * Pi * pow(Dipole(r_new),2) / (3.0 * V_new));
        
        //5. Minimum Metropolis condition
         double weight = exp((r.size()/4.0)*log(V_new/V)- (1/(R*Temp))*((u_new-u)+Pr*(V_new-V)));
         if(min(1,weight) >= Rand_DOUBLE())
           {
             r = r_new;
             V = V_new;
             u = u_new;
             L = pow(V,0.33333);
           }
         else
           {
             r_new = r;
             V_new = V;
             L = pow(V,0.33333);
           }
         
         // Averages of the properties
         if(k % 100 == 1)
           {
             outfile1 << k <<" "<< ((r.size()/4) * 33.252/V) << "  " <<  u/(96.49930*(r.size()/4)/3) << "\n";
           }      
         if(k > 3000000 && k%10000 == 1)
           {
             mean += 1;
             // 1. Radial Distribution Function
             for(unsigned int i = 1;i < RDF_size-1;++i)
               {
                 unsigned int N_All = 0,N_Si = 0,N_O = 0,N_SiO = 0;
                 double rij = 0;
                 for(unsigned int j = 0;j < r.size() - 4 ;j = j + 4)
                   {
                     for(unsigned int k_ = j + 4;k_ < r.size() ;k_ = k_ + 4)
                       {
                         rij = sqrt( pow((r[k_+0]-r[j+0]) - L * round((r[k_+0]-r[j+0])/L),2) 
                                   + pow((r[k_+1]-r[j+1]) - L * round((r[k_+1]-r[j+1])/L),2) 
                                   + pow((r[k_+2]-r[j+2]) - L * round((r[k_+2]-r[j+2])/L),2) );

                         if(rij <= rad[i] + (RDF_h/2) && rij > rad[i] - (RDF_h/2) )
                           {
                             N_All += 1;
                           }
                         if(rij <= rad[i] + (RDF_h/2) && rij > rad[i] - (RDF_h/2) && r[k_+3] == r[j+3] && r[k_+3] == q2)
                           {
                             N_O += 1;
                           }
                         if(rij <= rad[i] + (RDF_h/2) && rij > rad[i] - (RDF_h/2) && r[k_+3] == r[j+3] && r[k_+3] == q1 )
                           {
                             N_Si += 1;
                           }
                         if(rij <= rad[i] + (RDF_h/2) && rij > rad[i] - (RDF_h/2) && r[k_+3] != r[j+3] && r[k_+3] == q2)
                           {
                             N_SiO += 1;
                           }
                       }
                   }

                 RDF_all[i] += 2 * N_All  /( ((r.size()/4) / V) * 4 * Pi * rad[i] * rad[i] *  RDF_h);
                 RDF_Si[i]  += 2 * N_Si  / ( (((r.size()/4)*1/3) / V) * 4 * Pi * rad[i] * rad[i] *  RDF_h);
                 RDF_O [i]  += 2 * N_O   / ( (((r.size()/4)*2/3) / V) * 4 * Pi * rad[i] * rad[i] *  RDF_h);
                 RDF_SiO[i] += 2 * N_SiO / ( (((r.size()/4)*2/3) / V) * 4 * Pi * rad[i] * rad[i] *  RDF_h);
               }
           }
      }
    }

  if(rank == root)
    {
      for(unsigned int i = 0; i < RDF_size; ++i)
        {
          outfile2 << rad[i] << " " << RDF_all[i] / (mean * r.size()/4) << " " << RDF_Si[i] / (mean * (r.size()/4)/3) 
                   << " " << RDF_O [i] / (mean * (r.size()/4) * 2/3) << " " << RDF_SiO[i] / (mean * (r.size()/4) * 2/3) << "\n";
        }
    }

  MPI_Finalize();  
  return 0;
}
