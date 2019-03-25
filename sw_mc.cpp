#include "sw.hpp"


int main(int argc, char* argv[])
{
  if(argc!=3)
  {
    std::cerr << "Enter (1) lattice size, (2) number of SW sweeps. " << endl;
    exit(1);
  }

  size = std::atoi(argv[1]);
  int no_sweeps = std::atoi(argv[2]);
  int no_therm = no_sweeps/2;
  int no_msmt = no_sweeps-no_therm;
  long idum = time(NULL);

  sw** lattice = new sw* [size];
  for(int it=0; it<size; it++) lattice[it] = new sw [size];

  //initialize spins randomly
  for(int xc=0; xc < size; xc++)
  {
    for(int yc=0; yc<size; yc++)
    {
      int rand_spin = (ran0(&idum)<=0.5)?-1:1;
      lattice[xc][yc].set_S(rand_spin);
    }
  }

  // print_lattice(lattice);

  std::ofstream outfile_data("results.dat");
  for(double temperature = 6.0; temperature >= 0.1; temperature -= 0.04)
  {
    for(int sweep=1; sweep < no_therm; sweep++)
    {
      swendsen_flip(lattice, temperature);
    } 
    
    double total_M = 0.0;
    double total_E = 0.0;
    double total_E2 = 0.0;
    double total_M2 = 0.0;

    for(int sweep= no_therm; sweep < no_sweeps; sweep++)
    {
      swendsen_flip(lattice, temperature);
      total_M += magnetization(lattice);
      total_E += energy(lattice);
      total_E2 += pow(energy(lattice),2);
      total_M2 += pow(magnetization(lattice),2);
    }

    double avg_M = total_M/no_msmt;
    double avg_E = total_E/no_msmt;
    double avg_M2 = total_M2/no_msmt;
    double avg_E2 = total_E2/no_msmt;

    double err_E = sqrt(avg_E2 - pow(avg_E,2))/sqrt(no_msmt);
    double err_M = sqrt(avg_M2 - pow(avg_M,2))/sqrt(no_msmt);

    double cv = (avg_E2 - pow(avg_E,2)) /pow(temperature,2); 
    double xi = (avg_M2 - pow(avg_M,2)) /temperature; 

    outfile_data << temperature << " " << avg_M <<  " " << err_M << " " << avg_E << " " << err_E << " " << cv << " " << xi << endl;
  }

  // print_lattice(lattice);

  return 0;
}



