/*
 * =====================================================================================
 *
 *
 *    Description:  Fourier transform.
 *                  Be careful with indices: 
 *                  C++ and C alike languages start at 0.
 *
 *        Version:  1.0
 *        Created:  07/21/2016 08:32:14 PM
 *       Revision:  none
 *
 *         Author:  Daniel (dmc), danielmail7@gmail.com
 *   Organization:  University of Oviedo
 *
 * =====================================================================================
 */

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

// Fourier to Real
vector< complex<double> > fourier_to_real(vector< complex<double> >& c)
{
  const unsigned int M = c.size();
  const unsigned int N = (M-1)/2;
  vector< complex<double> > u(M, 0.0); // Initialize with zeros. 

  int k;
  for (unsigned int l=0; l<M; l++)
  {
    for (unsigned int k1=0; k1<M; k1++)
    {
      k = k1 - N;
      u.at(l) += c.at(k1) * polar<double>(1.0, 2.0 * M_PI * k * l / M);
    }
  }
  return u;
}
// Real to Fourier
vector< complex<double> > real_to_fourier(vector< complex<double> >& u)
{
  const unsigned int M = u.size();
  const unsigned int N = (M-1)/2;
  vector< complex<double> > c(M, 0.0); // Initialize with zeros. 

  int k;
  for (unsigned int k1=0; k1<M; k1++)
  {
    for (unsigned int l=0; l<M; l++)
    {
      k = k1 - N;
      c.at(k1) += (1.0/M) * u.at(l) * polar<double>(1.0, - 2.0 * M_PI * k * l / M);
    }
  }
  return c;
}


int main(void)
{
  //vector< complex<double> > c(21, complex<double>(1.0,0.0));
  const unsigned int N = 10;
  const unsigned int M = 2*N+1;

  vector< double > x(M, 0.0);
  vector< complex<double> > u(M, 1.0);
  vector< complex<double> > c(M, 0.0);
  vector< complex<double> > unew(M, 0.0);

  for (unsigned int l=0; l<M; ++l)
  {
    x.at(l) = 2 * M_PI * l / M;
    //Solution: u = exp(cos(x));
    u.at(l) = exp(cos(x.at(l)));
  }
  
  // Try fourier transform
  c = real_to_fourier(u);
  unew = fourier_to_real(c);
  cout << "#x    " << " u_calc   " << " u   " << endl;
  for (unsigned int l=0; l<M; ++l)
  {
    cout << x.at(l) << "  " 
         << real(unew.at(l)) << "  " 
         << real(u.at(l)) << endl;
  }

  return EXIT_SUCCESS;

}

  //u = for_each(x.begin(), x.end(), [](double &n){ exp(cos(n)); });
  //for (auto& uel: u)
  //  uel = exp(cos(x));
  //
  //for (const auto i: unew)    // C++11
  //  cout << real(i) << endl;

