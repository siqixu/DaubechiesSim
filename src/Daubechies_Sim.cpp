#include <RcppArmadillo.h>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

double wv_value(double x, List wv, int filter_number)
{
    NumericVector x_vec = as<NumericVector>(wv["x"]);
	NumericVector y_vec = as<NumericVector>(wv["y"]);
	if (x > (2 * filter_number - 1) || x < 0)
	{
        return(0);
    }
	else
	{
        return(y_vec[which_min(abs(x_vec - x))]);
    }
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

NumericVector wv_values(NumericVector X, List wv, int filter_number)
{
    NumericVector x_vec = as<NumericVector>(wv["x"]);
	NumericVector y_vec = as<NumericVector>(wv["y"]);
	NumericVector wv_x(X.size());
	for (unsigned i = 0; i != X.size(); ++i)
	{
		if (X[i] > (2 * filter_number - 1) || X[i] < 0)
		{
			wv_x[i] = 0;
		}
		else
		{
			wv_x[i] = y_vec[which_min(abs(x_vec - X[i]))];
		}
	}
	return wv_x;
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

NumericVector wvs_dilation_shift(NumericVector X, List wv, int filter_number, int j, int k)
{
	return pow(2, j / 2.0) * wv_values(pow(2, j) * X - k, wv, filter_number);
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

double wv_dilation_shift(double x, List wv, int filter_number, int j, int k)
{
	return pow(2, j / 2.0) * wv_value(pow(2, j) * x - k, wv, filter_number);
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

NumericVector wv_trans_fast(NumericVector X, List wv, int filter_number, double R, NumericVector gammas, double beta)
{
	NumericVector b(X.size());
	for (unsigned i = 0; i != X.size(); ++i)
	{
		for (unsigned g = 0; g!= gammas.size(); ++g)
		{
			double gamma = gammas[g];
			double temp = 0;
			int kk = floor(pow(2, gamma) * X[i]);
			for (unsigned k = kk - 2 * filter_number + 1; k != kk + 2; ++k)
			{
				temp = temp + wv_dilation_shift(X[i], wv, filter_number, gamma, k);
			}
			b[i] = b[i] + temp * R * pow(2, - gamma * (beta + 0.5));
		}
	}
	return b;
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

NumericVector wv_trans(NumericVector X, List wv, int filter_number, double R, NumericVector gammas, double beta)
{
    NumericVector b(X.size());
	for (unsigned i = 0; i != gammas.size(); ++i)
	{
		double gamma = gammas[i];
		NumericVector temp(X.size());
		for (unsigned k = - 2 * filter_number + 1; k != pow(2, gamma) + 1; ++k)
		{
			temp = temp + wvs_dilation_shift(X, wv, filter_number, gamma, k);
		}
		b = b + temp * R * pow(2, - gamma * (beta + 0.5));
	}
	return b;
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

NumericMatrix wv_analysis(NumericVector X, List wv, int filter_number, int level, NumericVector ks)
{
	NumericMatrix wv_mat(X.size(), ks.size());
	for (unsigned k = 0; k != ks.size(); ++k)
	{
		wv_mat(_, k) = wvs_dilation_shift(X, wv, filter_number, level, k + 1 - 5);
	}
	return wv_mat;
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat outer_prod(NumericVector X, arma::vec X_space_g_density, List wv, int filter_number, int level, NumericVector ks, int res, double norm_const)
{
  arma::mat cov = as<arma::mat>(wv_analysis(X, wv, filter_number, level, ks));
  return cov.t() * cov * X_space_g_density[X[X.size() - 1] * pow(10, res)] / norm_const;
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat generate_cov(NumericVector X, arma::vec X_space_g_density, List wv, int filter_number, int level, NumericVector ks, int res, double norm_const)
{
	const int dim = ks.size();
	arma::mat Cov(dim, dim);
	Cov.fill(0);
	NumericVector v(1);
	for (unsigned i = 0; i != X.size(); ++i)
	{
		v = X[i];
		Cov = Cov + outer_prod(v, X_space_g_density, wv, filter_number, level, ks, res, norm_const);
	}
	return Cov;
}
