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

using namespace std;
using namespace arma;

/*******************************************************************************
 * Fourier to Real
 * *****************************************************************************/
cx_vec fourier_to_real(cx_vec c)
{ // Zero frequency shifted to the origin
  const unsigned int M = c.size();
  const unsigned int N = (M-1)/2;
  cx_vec u = zeros<cx_vec>(M);

  int k;
  for (unsigned int l=0; l<M; l++)
  {
    for (unsigned int kpositive=0; kpositive<M; kpositive++)
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
  const unsigned int M = u.size();
  const unsigned int N = (M-1)/2;
  cx_vec c = zeros<cx_vec>(M);

  int k;
  for (unsigned int kpositive=0; kpositive<M; kpositive++)
  {
    for (unsigned int l=0; l<M; l++)
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
mat computeL(unsigned int M)
{
  /* <*,*> is the inner product
   * e_k(x) = e^(ikx) is the Fourier basis 
   * L_{ij} = < e_i, (-e_j'' + e_j) > = 2 pi delta_{ij}*(j^2+1)
   * or simply L_{ii} = 2 pi (j^2 + 1)
   * */
  const unsigned int N = (M-1)/2;
  mat L = zeros<mat>(M,M);

  int j, i;
  for (unsigned int i1=0; i1<M; i1++)
  {
    i = i1 - N;
    for (unsigned int j1=0; j1<M; j1++)
    {
      j = j1 - N;
      if (i==j) 
      {
        L(i1,j1) = 2.0 * M_PI * (pow(j,2) + 1);
      }
    }
  }
  return L;
}
/*******************************************************************************
 * Solve linear system with constant coefficients.
 * *****************************************************************************/
cx_vec solve_constant_coeff(cx_vec f)
{
  /* Solve c in:
   * Lc = f
   * */
  const unsigned int M = f.size();
  mat L     = zeros<mat>(M,M);
  cx_vec fi = zeros<cx_vec>(M);
  cx_vec c  = zeros<cx_vec>(M);
  vec rec  = zeros<vec>(M);

  L = computeL(M);
  fi = - 2.0 * M_PI * real_to_fourier(f);
  vec rfi = real(fi);
  //rfi.print("rfi: ");
  rec = solve(trimatl(L), -rfi);
  c.set_real(rec);
  //c = solve(trimatl(L),-fi);


  return c;
}


/*******************************************************************************
 * main()
 * *****************************************************************************/
int main()
{
  cout.precision(8);
  cout.setf(ios::fixed);
  const unsigned int N = 100;
  const unsigned int M = 2*N+1;

  cx_vec x    = zeros<cx_vec>(M);
  cx_vec u    = zeros<cx_vec>(M);
  cx_vec c    = zeros<cx_vec>(M);
  cx_vec unum = zeros<cx_vec>(M);
  cx_vec f    = zeros<cx_vec>(M);

  for (unsigned int l=0; l<M; ++l)
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

  //c = real_to_fourier(u);
  c = solve_constant_coeff(f);
  unum = fourier_to_real(c); // Numerical solution

  // Pretty printing
  mat A(M,3);
  A.col(0) = real(x);
  A.col(1) = real(unum);
  A.col(2) = real(u);
  A.print("#       x     unum        u");
 

  
  return 0;
}
