#include <iostream>
#include <fstream> 
#include <string> 
#include <vector> 
#include <cmath> 
#include <complex> // for complex numbers
#include <cstdlib> 

using namespace std;

// declaration de triangular_solve
// voir en dessous de main() pour la definition
// resout  un system A * solution = rhs si la matrice A est tridiagonal
template <class T> 
void triangular_solve(const vector<T>& diag,
		      const vector<T>& lower,
		      const vector<T>& upper,
		      const vector<T>& rhs,
		      vector<T>& solution);


// this is specialication function for the simplified output
// of std::valarray's without the need of writing all the
// elements oneself:
template <class T> 
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) 
{
  int len = v.size();
  
  for(int i = 0; i < len - 1; ++i)
    o << v[i] << " ";
  
  if(len > 0)
    o << v[len-1];
  
  return o;
}


// declaration de probability: calculer la probabilite de trouver la particule entre les points nl et nr du maillage avec alpha x = dx
// voir en dessous de main() pour la definition
template <class T> 
double probability(const T& psi, int nl, int nr, double dx);
double getenergy(const vector<complex<double> >& psi, const vector<complex<double> >& diagH, const vector<complex<double> >& upperH, const vector<complex<double> >& lowerH, const double &);
double getmeanx(const vector<complex<double> >& psi, const vector<double>& x, const double &);
double getmeanx2(const vector<complex<double> >& psi, const vector<double>& x, const double &);
double getmeanp(const vector<complex<double> >& psi, const vector<double>& x, const double & dx);
double getmeanp2(const vector<complex<double> >& psi, const vector<double>& x, const double & dx);

// fonction qui calcule le potentiel
double getpotential(double alpha, double V0, double Rnucleus, double x)
{
  return (x<-Rnucleus)? -alpha/x : V0;
}


int main(int argc,char **argv) 
{
    double hbar = 1.;
  double mass = 0.5;
  double xl=-256;
  double xr=0.0;
  double L= xr-xl;

  // the imaginary i
  complex<double> complex_I = complex<double> (0,1);

  
  double n=32;
  cerr << "n " << endl;
  cin >> n;
  
  double alpha;
  cerr << "alpha " << endl;
  cin >> alpha;
  double zk = 2.0 * M_PI * n/L;

  // energy cinetique approximative de l'onde
  double E0 = hbar * hbar * zk * zk/(2. * mass);
  cerr << "E0 = " << E0 << endl;
  
  
  double Rnucleus=50;
  cerr << "Rnucleus " << endl;
  cin >> Rnucleus;
  
  double V0=-1;
  cerr << "V0 " << endl;
  cin >> V0;
  
  double x0=-128;
  cerr << "x0 " << endl;
  cin >> x0;
  
  double sigmanorm = 0.05;
  cerr << "sigmanorm " << endl;
  cin >> sigmanorm;  
  double sigma0=sigmanorm*L;
  
  // temps
  double dt = 1.0;
  cerr << "pas de temps: dt= " << endl;
  cin >> dt;

  double tfinal = 2500.;
  cerr << "t final " << endl;
  cin >> tfinal;

  int ndx = 1024;   // nombre d'intervalles
  cerr << "nombre d'intervales en x: ndx= " << endl;
  cin >> ndx;

  int nx = ndx + 1;  // nombre de points 
  double dx = (xr-xl)/ndx;
  vector<double> x(nx);
  
  // les points du maillage
  for(int i = 0; i < nx; ++i)
    x[i] = xl + i * dx;


  // les fonctions d'onde
  vector<complex<double> > psi_next(nx), psi_now(nx);
  
  // the diagonal and the next-to-diagonal elements of the Hamiltonian matrix:
  // Hamiltonian matrix is needed for energy calculation
  vector<complex<double> > dH(nx), aH(ndx, -(hbar*hbar)/(2*mass*dx*dx)), cH(ndx, -(hbar*hbar)/(2*mass*dx*dx));
  
  // les 3 diagonales des matrices A et B dans l'equation (4.80)
  // d = diagonal
  // a = lower diagonal
  // c = upper diagonal
  vector<complex<double> > dA(nx), dB(nx), aA(ndx, -(hbar*dt*1i)/(4*mass*dx*dx)), aB(ndx, (hbar*dt*1i)/(4*mass*dx*dx)), cA(ndx, -(hbar*dt*1i)/(4*mass*dx*dx)), cB(ndx, (hbar*dt*1i)/(4*mass*dx*dx));

  for(int i = 0; i < nx; i++)
  {
      if(x[i] < -Rnucleus)
      {
        dH(i) = (hbar*hbar)/(2*mass*dx*dx)+alpha/x[i];
      }
      else
      {
        dH(i) = (hbar*hbar)/(2*mass*dx*dx)+V0;
      }
      dA(i) = 1+1i*dt*dH(i)/(2*hbar);
      dB(i) = 1-1i*dt*dH(i)/(2*hbar);
  }

  // prepare la fonction d'onde initiale;
  for(int i = 0; i < nx; ++i) 
    {
      psi_next[i] = 0.0;
      psi_now[i] = exp(complex_I * zk * x[i]) * exp( -(x[i] - x0) *(x[i] - x0)/(2.0 * sigma0 * sigma0) );
    }
  
  // normaliser la fonction d'onde
  double pre_norm = sqrt(probability(psi_now,0,nx-1,dx));
  if(pre_norm==0){pre_norm=1;}
  
  for(int i = 0; i < nx; ++i)
    psi_now[i] /= pre_norm;
  
  cerr << "pre_norm = " << pre_norm << endl;
  // verifier que la norme vaut un
  cerr << "Norm = " << probability(psi_now,0,nx-1,dx) << endl;
  // ecrire l'energie moyenne initiale
  cerr << "<H> = " << getenergy(psi_now, dH, aH, cH, dx) << endl;
  // ecrire la position moyenne initiale
  cerr << "<x> = " << getmeanx(psi_now, x, dx) << endl;
  cerr << "V(0) = " << getpotential(alpha, V0, Rnucleus , 0.0) << endl;
  cerr << "<p> = "<< getmeanp(psi_now,x,dx) << endl;
  cerr << "<pp> = "<< getmeanp2(psi_now,x,dx) << endl;
  
  cout << 0. << " "
	   << probability(psi_next,xl,-Rnucleus,dx) << " "  // proba "`a gauche"
	   << probability(psi_next,-Rnucleus,xr,dx) << " " // proba "\`a droite"
       << getmeanx(psi_now,x,dx) << " "             // mean position 
       << getenergy(psi_now,dH,aH,cH,dx) << " "   // mean energy
       << getmeanp(psi_now,x,dx) << " "
       << getmeanp2(psi_now,x,dx)<< " " // proba "\`a droite"
       << getmeanx2(psi_now,x,dx) << endl;
       
  // ecrire la fonction d'onde initiale dans le fichier "psi.dat"
  ostringstream oss;
  oss << "psi.dat";

  ofstream ofs(oss.str().c_str()); // asci
    for(int i = 0; i < nx; ++i)
      ofs << 0 << " " << x[i] << " " << abs(psi_now[i]) * abs(psi_now[i]) << endl;
  
    for(double time = 0; time < tfinal - dt/2.; time += dt)
    {
        psi_next[0] = dB[0]*psi_now[0] + cB[0]*psi_now[1];
        for(int i = 1; i < ndx; i++)
        {
            psi_next[i] = dB[i]*psi_now[i] + cB[i]*psi_now[i+1] + aB[i-1]*psi_now[i-1];
        }
        psi_next[ndx] = dB[ndx]*psi_now[ndx] + aB[ndx-1]*psi_now[ndx-1];

        triangular_solve(dA, aA, cA, psi_next, psi_next);
      
      // output the probabilities "left" and "right", mean position and mean energy
      cout << time+dt << " "
	   << probability(psi_next,xl,-Rnucleus,dx) << " "  // proba "`a gauche"
	   << probability(psi_next,-Rnucleus,xr,dx) << " " // proba "\`a droite"
	   << getmeanx(psi_next,x,dx) << " "            // mean position 
           << getmeanx2(psi_next,x,dx) << " "
	   << getenergy(psi_next,dH,aH,cH,dx) << " "  // mean energy
           << getmeanp(psi_next,x,dx) << " "
           << getmeanp2(psi_next,x,dx)<< endl; // proba "\`a droite"
           
      for(int i = 0; i < nx; ++i)
	 {
	  // asci mode
	  ofs << time << " " << x[i] << " " << abs(psi_next[i]) * abs(psi_next[i]) << endl;
	 }
      
      psi_now = psi_next;

    } // end of time evolution loop

  ofs.close();

} // end of main(...)


//
// Calcule la probabilite de trouver la particule dans
// un intervalle, entre les points de maillage nl et nr
//
template <class V> 
double probability(const V& psi, int nl, int nr, double dx) 
{
  double retval = 0;
  if(nr >= psi.size())
  {
      nr = psi.size()-1;
  }
  for(int i = nl; i < nr; i++)
  {
      retval += norm(psi[i])+norm(psi[i+1]);
  }

  return retval*dx/2;
}

double getenergy(const vector<complex<double> >& psi, const vector<complex<double> >& diagH, const vector<complex<double> >& upperH, const vector<complex<double> >& lowerH, const double & dx)
{
  vector<complex<double> > psi_tmp(psi.size());
  double energy=0.;
  psi_tmp[0] = diagH[0]*psi[0] + upperH[0]*psi[1];
  for(int i = 1; i < upperH.size(); i++)
  {
      psi_tmp[i] = diagH[i]*psi[i] + upperH[i]*psi[i+1] + lowerH[i-1]*psi[i-1];
  }
  psi_tmp[upperH.size()] = diagH[upperH.size()]*psi[upperH.size()] + lowerH[upperH.size()-1]*psi[upperH.size()-1];

  for(int i = 0; i < upperH.size(); i++)
  {
      energy += conj(psi[i])*psi_tmp[i]+conj(psi[i+1])*psi_tmp[i+1];
  }

  return energy*dx/2;
}

double getmeanx(const vector<complex<double> >& psi, const vector<double>& x, const double & dx)
{
  double meanx=0.;
  //TODO
  return meanx;
}

double getmeanx2(const vector<complex<double> >& psi, const vector<double>& x, const double & dx)
{
  double meanx2=0.;
  //TODO
  return meanx2;
}

double getmeanp(const vector<complex<double> >& psi, const vector<double>& x, const double & dx)
{
  double meanp= 0;
  //TODO
  return meanp;
}

double getmeanp2(const vector<complex<double> >& psi, const vector<double>& x, const double & dx)
{
  double meanp2=0;
  //TODO
  return meanp2;
}


//
// definition of triangular_solve
//
template <class T> 
void triangular_solve(const vector<T>& diag,
		      const vector<T>& lower,
		      const vector<T>& upper,
		      const vector<T>& rhs,
		      vector<T>& solution) 
{
  vector<T> new_diag = diag;
  vector<T> new_rhs = rhs;
  
  // forward elimination
  for(int i = 1; i < diag.size(); ++i) 
    {
      T pivot = lower[i-1]/new_diag[i-1];
      new_diag[i] -= pivot * upper[i-1];
      new_rhs[i] -= pivot * new_rhs[i-1];
    }
  
  solution.resize(diag.size());
  
  // solve last equation
  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];
  
  // back substitution
  for(int i = diag.size() - 2; i >= 0; --i) 
    {
      solution[i] = (new_rhs[i] - upper[i] * solution[i+1]) / new_diag[i];
    }
}

