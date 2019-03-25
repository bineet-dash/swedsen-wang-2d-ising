#ifndef _SW_HPP_DEFINED_
#define _SW_HPP_DEFINED_

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <clocale>
#include <iomanip>

using std::cout;
using std::exp;
using std::endl;

extern int size;
extern double J;

int size;
double J=1;
long idum = time(NULL);
inline int pbc(int a, int b){return (a+b+size)%size;}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
double ran0(long *idum)
{
	long  k;
	double ans;

	*idum ^= MASK;
	k = (*idum)/IQ;
	*idum = IA*(*idum - k*IQ) - IR*k;
	if(*idum < 0) *idum += IM;
	ans=AM*(*idum);
	*idum ^= MASK;
	return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK


void print_arrow(int s)
{
  const char  *ptr = NULL;
  const wchar_t uparrow[] = L"\u2191";
  const wchar_t downarrow[] = L"\u2193";
  freopen(ptr, "w", stdout);
  setlocale(LC_ALL, "");
  if(s==1)   std::wcout << uparrow;
  else if(s==-1) std::wcout << downarrow;
  else  std::wcout << "-" << " ";
  freopen(ptr, "w", stdout);
}

class sw
{
private:
  int spin; bool in_cluster; int cluster;
public:
  sw(){spin=1; in_cluster=false; cluster=-1;}
  sw(int s){spin=s; in_cluster=false;}
  void set_cl(int c){in_cluster=true; cluster=c;}
  bool bonded(void){return in_cluster;}
  void reset_bond(void){in_cluster=false; cluster=-1;}
  int S(void){return spin;}
  void set_S(int s) {spin=s;}
  int cl(void){return cluster;}
  void flip(void) {spin *= -1;}
  double bond_prob(sw s2, double beta){return (1-exp(-1*beta*J))*spin*s2.S(); }
  void print(void){ cout << ( (cluster<10)? "(0" : "(") << cluster << ") "; print_arrow(spin); }
};

void merge_clusters(sw** lattice, int cluster_id1, int cluster_id2 )
{
  int merged_cluster_id = std::min(cluster_id1, cluster_id2);
  int max_cluster_id = std::max(cluster_id1, cluster_id2); 
  for(int xc = 0; xc < size; xc++)
  {
    for(int yc=0; yc<size; yc++)
    {
      if(lattice[xc][yc].cl()== max_cluster_id) lattice[xc][yc].set_cl(merged_cluster_id); 
    }
  }
}

void make_bonds(sw** lattice, int xc, int yc, double temperature)
{
  double beta = 1/temperature;
  sw central = lattice[xc][yc];

  if(ran0(&idum) <= central.bond_prob(lattice[pbc(xc,1)][yc], beta)) 
  {
    if(!lattice[pbc(xc,1)][yc].bonded()) 
    {
      lattice[pbc(xc,1)][yc].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[pbc(xc,1)][yc].cl(), central.cl()) ;
    }
  }
  
  if(ran0(&idum) <= central.bond_prob(lattice[pbc(xc,-1)][yc], beta)) 
  {
    if(!lattice[pbc(xc,-1)][yc].bonded()) 
    {
      lattice[pbc(xc,-1)][yc].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[pbc(xc,-1)][yc].cl(), central.cl()) ;
    }
  }
   
  if(ran0(&idum) <= central.bond_prob(lattice[xc][pbc(yc,1)], beta)) 
  {
    if(!lattice[xc][pbc(yc,1)].bonded()) 
    {
      lattice[xc][pbc(yc,1)].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[xc][pbc(yc,1)].cl(), central.cl()) ;
    }
  }
    
  if(ran0(&idum) <= central.bond_prob(lattice[xc][pbc(yc,-1)], beta)) 
  {
    if(!lattice[xc][pbc(yc,-1)].bonded()) 
    {
      lattice[xc][pbc(yc,-1)].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[xc][pbc(yc,-1)].cl(), central.cl()) ;
    }
  }
}

int make_cluster(sw** lattice, double temperature)
{
  int cluster_id = 0;
  for(int xc=0; xc<size; xc++)
  {
    for(int yc=0; yc<size; yc++)
    {
      if(!lattice[xc][yc].bonded()) 
      {
        lattice[xc][yc].set_cl(cluster_id);
      }

      make_bonds(lattice, xc, yc, temperature);
      
      if(lattice[xc][yc].cl()==cluster_id)
      {
        cluster_id++;
      }
    }
  }
  return cluster_id;
}

void flip_cluster(sw** lattice, bool* decision_list)
{
  for(int xc=0; xc<size; xc++)
  {
    for(int yc=0; yc<size; yc++)
    {
      if(decision_list[lattice[xc][yc].cl()]) lattice[xc][yc].flip();
    }
  }
}

void swendsen_flip(sw** lattice, double temperature)
{
  int no_clusters = make_cluster(lattice, temperature);
  bool* flip_decision = new bool [no_clusters];
  for(int i=0; i<no_clusters; i++) flip_decision[i] = (ran0(&idum)<=0.5);
  flip_cluster(lattice, flip_decision);
  delete[] flip_decision;
  for(int xc=0; xc<size; xc++)
  {
    for(int yc=0; yc<size; yc++)
    {
      lattice[xc][yc].reset_bond();
    }
  }
}

double magnetization(sw** lattice)
{
  double mag = 0.0;
  for(int xc=0; xc<size; xc++)
  {
    for(int yc=0; yc<size; yc++)
    {
      mag += double(lattice[xc][yc].S());
    }
  }
  return abs(mag)/(double(size*size));
}

double energy(sw** lattice)
{
  double E = 0.0;
  for(int xc=0; xc < size; xc++)
  {
    for(int yc=0; yc<size; yc++)
    {
      E += -J*lattice[xc][yc].S()*(lattice[pbc(xc,1)][yc].S() + lattice[pbc(xc,-1)][yc].S() + lattice[xc][pbc(yc,1)].S() + lattice[xc][pbc(yc,-1)].S());
    }
  }
  return E/(size*size);
}

void print_lattice(sw** lattice)
{    
  for(int xc=0; xc < size; xc++)
  {
    for(int yc=0; yc<size; yc++)
    {
      lattice[xc][yc].print(); cout << " ";
      // cout << std::setw(3)  << lattice[xc][yc].cl() << " ";
    }
    cout << endl;
  }
}


#endif