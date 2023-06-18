#pragma once

#include <iostream>
#include <vector>
#include <assert.h>
#include <chrono>
#include <algorithm>
#include <memory>
#include <numeric>
// #include <boost/python.hpp>
// #include"cnpy.h"
// #include <Eigen/Dense>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
// #include <boost/math/special_functions/factorials.hpp>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
// #include "gnuplot-iostream.h"

// #include <matplot/matplot.h>

typedef double real_t;
using std::size_t;

std::vector<real_t> linspace(const real_t start, const real_t stop, const size_t num);
std::vector<real_t> arange(const real_t start, const real_t stop, const real_t step);
size_t unravel_index(std::vector<size_t> indices, std::vector<size_t> shape);
size_t _unravel_index(std::vector<size_t> indices, std::vector<size_t> shape, size_t rek);
// void::draw_line
std::vector<std::vector<real_t>> hyperslab2d(std::vector<real_t> input, std::vector<size_t> shape, std::vector<size_t> start, std::vector<size_t> end);


template<typename T>
void print_vector(std::vector<T> input){
  std::cout << "{";
  for(size_t i = 0; i < input.size(); i++){
    std::cout << std::to_string(input[i]);
    if(i < input.size()-1){std::cout << ", ";};
  }
  std::cout << "}" ;
}

template<typename T, size_t D>
void print_vector(std::array<T, D> input){
  std::cout << "{";
  for(size_t i = 0; i < input.size(); i++){
    std::cout << std::to_string(input[i]);
    if(i < input.size()-1){std::cout << ", ";};
  }
  std::cout << "}" ;
}

template<size_t D>
void print_vector(std::array<std::string, D> input){
  std::cout << "{";
  for(size_t i = 0; i < input.size(); i++){
    std::cout << input[i];
    if(i < input.size()-1){std::cout << ", ";};
  }
  std::cout << "}" ;
}




template<typename T>
Eigen::VectorXd coefficients(std::vector<T> stencil, size_t derivative){
  size_t n = stencil.size();
  Eigen::MatrixXd matrix (n,n);
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < n; j++)
    {
      matrix(i,j) = pow(stencil[j], i);
    }
  }
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
  // rhs(derivative) = boost::math::factorial<real_t>(derivative);
  rhs(derivative) = std::tgamma(n+1);
  Eigen::VectorXd coeffs = matrix.colPivHouseholderQr().solve(rhs);
  // Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
  // real_t cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
  // std::cout << "cond: " << cond << std::endl;
  // Eigen::VectorXd coeffs = matrix.fullPivHouseholderQr().solve(rhs);
  return coeffs;
}

template<typename T>
Eigen::SparseMatrix<T> finite_difference_uneven_spaced_1d(std::vector<T> h, size_t acc=2, bool norm=false){
  Eigen::VectorXd normed_h = Eigen::VectorXd::Map(h.data(), h.size());
  if(norm){
    normed_h.normalize();
    normed_h *= h.size();
  }
  long n = h.size();
  long stencil_size = (acc*2+1);
  size_t nnz = n * stencil_size;
  // std::vector<Eigen::Triplet<T>> (nnz); // row, column, value
  Eigen::SparseMatrix<T> coefficent_sparse_matrix(n,n);
  coefficent_sparse_matrix.reserve(nnz);
  // row = np.zeros((nnz))
  // col = np.zeros((nnz))
  // data = np.zeros((nnz))
  size_t counter = 0;
  for (long i = 0; i < n; i++)
  {
    std::vector<T> stencil (stencil_size);
    long k = 0;
    long left = i-acc;
    long right = i+acc;
    while (left < 0)
    {
      k++;
      left++;
    }
    while (right >= n)
    {
      k--;
      right--;
    }
    long center_index = acc-k;
    // std::cout << "stencil for i: " << i << " and center index: " << center_index << std::endl;
    for (long j = 0; j < stencil_size; j++)
    {
      // std::cout << "j: " << j << " ";
      if (j == center_index)
      {
        stencil[j] = 0.;
        // std::cout << 0;
      }
      else if(j < center_index){
        stencil[j] = -std::abs(normed_h[i+j-center_index] - normed_h[i]);
        // std::cout << "-abs(" << h[i+j-center_index] << " - " << h[i] << ")";
      }
      else if(j > center_index){
        stencil[j] = std::abs(normed_h[i+j-center_index] - normed_h[i]);
        // std::cout << "abs(" << h[i+j-center_index] << " - " << h[i] << ")";
      }
      // else{
      //   stencil[j] = h[i+j-center_index] - h[i];
      // }
      // std::cout << std::endl;
    }
    // std::cout << "stencil for i: " << i << std::endl;
    // print_vector(stencil); std::cout << std::endl;
    Eigen::VectorXd coeffs = coefficients(stencil, 1);
    for (long j = 0; j < stencil_size; j++)
    {
      // Eigen::Triplet<T> triplet = {i, i+j-center_index, coeffs[j]};
      // coefficients[counter] = triplet;
      if(k==0 && j==center_index)
        coeffs[j] = 0.;
      coefficent_sparse_matrix.insert(i, i+j-center_index) = coeffs[j];
      counter ++;
    }
  }
  return coefficent_sparse_matrix;
}

template<typename T>
Eigen::SparseMatrix<T> finite_difference_uneven_spaced_2d_dx(size_t ny, std::vector<T> x, size_t acc=2, bool norm=false){
  Eigen::SparseMatrix<T> one_dim_x = finite_difference_uneven_spaced_1d(x, acc, norm);
  Eigen::SparseMatrix<T> I(ny, ny);
  I.reserve(ny);
  for (size_t i = 0; i < ny; i++)
  {
    I.insert(i, i) = 1.;
  }
  Eigen::SparseMatrix<T> Dx = Eigen::kroneckerProduct(I, one_dim_x);
  return Dx;
}

template<typename T>
Eigen::SparseMatrix<T> finite_difference_uneven_spaced_2d_dy(size_t nx, std::vector<T> y, size_t acc=2, bool norm=false){
  Eigen::SparseMatrix<T> one_dim_y = finite_difference_uneven_spaced_1d(y, acc, norm);
  Eigen::SparseMatrix<T> I(nx, nx);
  I.reserve(nx);
  for (size_t i = 0; i < nx; i++)
  {
    I.insert(i, i) = 1.;
  }
  Eigen::SparseMatrix<T> Dy = Eigen::kroneckerProduct(one_dim_y, I);
  return Dy;
}

template<size_t D>
std::array<size_t, D> ravel_index(const size_t &index, const std::array<size_t, D> &subspaces){
  std::array<size_t, D> indices;
  size_t remainder = index;
  for(size_t rek = 0; rek < D; rek++){
    indices[rek] = remainder/subspaces[rek];
    remainder = remainder%subspaces[rek];
  }
  return indices;
}

template<size_t D>
inline std::array<size_t, D> ravel_index(const size_t* subspaces, const size_t &index){
  std::array<size_t, D> indices;
  size_t remainder = index;
  for(size_t rek = 0; rek < D; rek++){
    indices[rek] = remainder/subspaces[rek];
    remainder = remainder%subspaces[rek];
  }
  return indices;
}

template<size_t D>
inline size_t unravel_index(const std::array<size_t, D> &indices, const std::array<size_t, D> &subspaces){
  size_t index = 0;
  for(size_t rek = 0; rek != D; rek++){

    index += (indices[rek] * subspaces->at(rek));
  }

  return index;
}

template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
    const std::vector<T>& vec,
    const Compare& compare)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
    return p;
}

template <typename T>
std::vector<T> apply_permutation(
    const std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](std::size_t i){ return vec[i]; });
    return sorted_vec;
}




// Given three collinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(real_t p_x, real_t p_y, real_t q_x, real_t q_y, real_t r_x, real_t r_y);
// // To find orientation of ordered triplet (p, q, r).
// // The function returns following values
// // 0 --> p, q and r are collinear
// // 1 --> Clockwise
// // 2 --> Counterclockwise
int orientation(real_t p_x, real_t p_y, real_t q_x, real_t q_y, real_t r_x, real_t r_y);

// // The main function that returns true if line segment 'p1q1'
// // and 'p2q2' intersect.
bool doIntersect(real_t p1_x, real_t p1_y, real_t q1_x, real_t q1_y, real_t p2_x, real_t p2_y, real_t q2_x, real_t q2_y);

/// @brief 
/// @param tau Lagearameter/Erwarungswert
/// @param sigma sigma squared is the variance
/// @return 
real_t normal_distribution(const real_t tau, const real_t sigma);
