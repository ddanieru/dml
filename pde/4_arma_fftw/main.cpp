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
  mat L(M,M);

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
template <class Tvector>
Tvector ifftshift(Tvector anyvector)
{
  long n;
  size_t M = anyvector.size();

  if (M % 2 ==1) 
  { /*  M is odd */ 
    n = M/2 + 1;
  } else { 
    n = M/2+1;
  }
  anyvector = shift(anyvector,n); // M=21 => +11

  return anyvector;

}

/*******************************************************************************
 * Shift the zero-frequency to the center after fft
 * *****************************************************************************/
template <class Tvector>
Tvector fftshift(Tvector anyvector)
{
  long n;
  size_t M = anyvector.size();

  if (M % 2 ==1) 
  { /*  M is odd */ 
    n = M/2;
  } else {
    n = M/2-1;
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
  f = fftshift<cx_vec>(f); 
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
  const unsigned long M = pow(2,10); //2*N+1;

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
  c = ifftshift<cx_vec>(c);
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
