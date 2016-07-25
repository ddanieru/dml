/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Fourier transform
 *
 *        Version:  1.0
 *        Created:  06/23/2016 01:26:33 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Daniel (dmc), danielmail7@gmail.com
 *   Organization:  University of Oviedo
 *
 * =====================================================================================
 */

#include <iostream>
#include <armadillo>
#include <fftw3.h>

using namespace std;
using namespace arma;

/*******************************************************************************
 * Fourier to Real
 * *****************************************************************************/
cx_vec fourier_to_real(cx_vec c)
{ // Zero frequency shifted to the origin
  const unsigned long M = c.size();
  const unsigned long N = (M-1)/2;
  //cx_vec u = zeros<cx_vec>(M);
  cx_vec u(M);

  long k;
  for (unsigned long l=0; l<M; l++)
  {
    for (unsigned long kpositive=0; kpositive<M; kpositive++)
    {
      k = kpositive - N;
      u(l) += c(kpositive) * polar<double>(1.0, 2.0 * M_PI * k * l / M);
    }
  }
  return u;
}
/*******************************************************************************
 * Real to Fourier
 * *****************************************************************************/
cx_vec real_to_fourier(cx_vec u)
{ // Zero frequency shifted to the origin
  const unsigned long M = u.size();
  const unsigned long N = (M-1)/2;
  //cx_vec c = zeros<cx_vec>(M);
  cx_vec c(M);

  long k;
  for (unsigned long kpositive=0; kpositive<M; kpositive++)
  {
    for (unsigned long l=0; l<M; l++)
    {
      k = kpositive - N;
      c(kpositive) += (1.0/M) * u(l) * polar<double>(1.0, - 2.0 * M_PI * k * l / M);
    }
  }
  return c;
}
/*******************************************************************************
 * Compute Lower triangular matrix
 * *****************************************************************************/
mat computeL(unsigned long M)
{
  /* <*,*> is the inner product
   * e_k(x) = e^(ikx) is the Fourier basis 
   * L_{ij} = < e_i, (-e_j'' + e_j) > = 2 pi delta_{ij}*(j^2+1)
   * or simply L_{ii} = 2 pi (j^2 + 1)
   * */
  const unsigned long N = (M-1)/2;
  //mat L = zeros<mat>(M,M);
  mat L(M,M);

  /*
  long j, i;
  for (unsigned i1=0; i1<M; i1++)
  {
    i = i1 - N;
    for (unsigned j1=0; j1<M; j1++)
    {
      j = j1 - N;
      if (i==j) 
      {
        L(i1,j1) = 2.0 * M_PI * (pow(j,2) + 1);
      }
    }
  }
  */
  long j;
  for (unsigned j1=0; j1<M; j1++)
  {
    j = j1 - N;
    L(j1,j1) = 2.0 * M_PI * (pow(j,2) + 1);
  }
  return L;
}

/*******************************************************************************
 * Shift the zero-frequency to the center before fft
 * *****************************************************************************/
cx_vec ifftshift(cx_vec anyvector)
{
  long n;
  size_t M = anyvector.size();

  if (M % 2 ==1) 
  { /*  M is odd */ 
    n = M/2 + 1;
  } else { 
    n = M/2;
  }
  anyvector = shift(anyvector,n); // M=21 => +11

  return anyvector;

}

/*******************************************************************************
 * Shift the zero-frequency to the center after fft
 * *****************************************************************************/
cx_vec fftshift(cx_vec anyvector)
{
  long n;
  size_t M = anyvector.size();

  if (M % 2 ==1) 
  { /*  M is odd */ 
    n = M/2;
  } else {
    n = M/2;
  }
  anyvector = shift(anyvector,n); // M=21 => +10

  return anyvector;

}
/*******************************************************************************
 * Solve linear system with constant coefficients.
 * *****************************************************************************/
cx_vec solve_constant_coeff(cx_vec f)
{
  /* Solve c in:
   * Lc = f
   * */
  const unsigned long M = f.size();
  //mat L     = zeros<mat>(M,M);
  //cx_vec c  = zeros<cx_vec>(M);
  //vec rec   = zeros<vec>(M);
  mat L(M,M);
  cx_vec c(M);
  vec rec(M);

  L = computeL(M);

  // FFT
  fftw_complex* in = reinterpret_cast<fftw_complex*> (f.memptr());
  fftw_plan plan = fftw_plan_dft_1d(M, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
  // Put data in samples (other than ESTIMATE)
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  f /= M;                  // Normalization
  f = fftshift(f); 
  f *= - 2.0 * M_PI ;

  vec rfi   = real(f);
  rec = solve(trimatl(L), -rfi);
  c.set_real(rec);

  return c;
}


/*******************************************************************************
 * main()
 * *****************************************************************************/
int main()
{
  cout.precision(8);
  cout.setf(ios::fixed);
  const unsigned long N = 5000;
  const unsigned long M = 2*N+1;

  //cx_vec x    = zeros<cx_vec>(M);
  //cx_vec u    = zeros<cx_vec>(M);
  //cx_vec c    = zeros<cx_vec>(M);
  //cx_vec unum = zeros<cx_vec>(M);
  //cx_vec f    = zeros<cx_vec>(M);

  cx_vec x(M);
  cx_vec u(M);
  cx_vec c(M);
  cx_vec unum(M);
  cx_vec f(M);


  for (unsigned long l=0; l<M; ++l)
  {
    x(l) = 2 * M_PI * l / M;
  }
  /********************************* 
   * Solve:
   *         -u''(x) + u'(x) = f(x)
   ********************************* 
   * Knowing that u = exp(cos(x)) then
   * f = u (-sin(x) sin(x) + cos(x) + 1)
   * and now we solve u from this f
   * with the Galerkin method.
   ************************************* 
   * */
  u = exp(cos(x)); // Analytical solution.
  f=u % (-sin(x) % sin(x)+cos(x) + 1);

  c = solve_constant_coeff(f);

  // inverse FFT
  c = ifftshift(c);
  fftw_complex* in = reinterpret_cast<fftw_complex*> (c.memptr());
  fftw_complex* out = reinterpret_cast<fftw_complex*> (unum.memptr());
  fftw_plan plan = fftw_plan_dft_1d(M, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  // Put data in samples. 
  fftw_execute(plan);
  /* Free memory */
  fftw_destroy_plan(plan);
  //fftw_free(in);
  //fftw_free(out);


  // Pretty printing
  mat A(M,3);
  A.col(0) = real(x);
  A.col(1) = real(unum);
  A.col(2) = real(u);
  A.print("#       x     unum        u");

  cout << "#Armadillo version: " << arma_version::as_string() << endl;

  return 0;
}
