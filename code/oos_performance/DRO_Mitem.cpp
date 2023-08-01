#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>    
#include <iostream>
#include <random>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "solveNVSAA.h"
#include "solveNVDROW.h"
#include "solveNVDROM.h"
#include "solveNVKDE.h"
#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define ABS(x) (((x) > 0 ) ? (x) : -(x))	
#define EPSILON	 0.000000000001    //Added to the fractional solution obtained at the node
#define ENOUGH ((CHAR_BIT * sizeof(int) - 1) / 3 + 42)
FILE   *Open_File(const char *name, const char *mode);
double *prob_;
void compute_prob(double *fea, double **fea_train, int M, int n, int mk);
using namespace std;
double *h;
namespace Eigen {
	namespace internal {
		template<typename Scalar>
		struct scalar_normal_dist_op
		{
			static boost::mt19937 rng;    // The uniform pseudo-random algorithm
			mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

			EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

				template<typename Index>
			inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
		};

		template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

		template<typename Scalar>
		struct functor_traits<scalar_normal_dist_op<Scalar> >
		{
			enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false };
		};
	} // end namespace internal
} // end namespace Eigen


void main(int  argc, char *argv[]) {
	int CVK = 4;
	int testN = 1;
	int Ms = 3;// number of group scenarios
	int size = 3; // Dimensionality (rows)
	int *M_ins;
	M_ins = create_int_vector(Ms);//number of scenarios
	M_ins[0] = 8;
	M_ins[1] = 52;
	M_ins[2] = 100;
	h = create_double_vector(Ms);
	int num_Method = 7;
	int expk, Expk = 1000;
	double ***all_results2;
	all_results2 = create_double_tmatrix(Expk, num_Method, Ms);
	double ***all_results2_;
	all_results2_ = create_double_tmatrix(Expk, num_Method, Ms);
	double *oos_profit;
	oos_profit = create_double_vector(num_Method);

	int n;
	int num_G;//number of gammas
	int M_oos;//number of oos sample
	double *Gammas;
	double *fea_fix;
	double **fea_train_all;
	char stri[ENOUGH];
	sprintf(stri, "./data_oos/input.dat");
	FILE		*out;
	out = Open_File(stri, "a+");
	fscanf(out, "%d", &n);
	fscanf(out, "%d", &num_G);
	fscanf(out, "%d", &M_oos);
	fea_fix = create_double_vector(n);
	Gammas = create_double_vector(num_G);
	for (int i = 0; i < num_G; i++) {
		fscanf(out, "%lf", &Gammas[i]);
		//printf("gamma is %.2lf\t", Gammas[i]);
	}
	for (int i = 0; i < Ms; i++) {
		fscanf(out, "%lf", &h[i]);		
	}
	fclose(out);


	double *ubar, *bs, *cs;
	ubar = create_double_vector(n);
	double *uu;
	uu = create_double_vector(n);
	bs = create_double_vector(n);
	cs = create_double_vector(n);
	int max_Ms = 0;
	for (int i = 0; i < Ms; i++) {
		if (max_Ms < M_ins[i]) {
			max_Ms = M_ins[i];
		}
	}
	double **oos_xis;
	oos_xis = create_double_matrix(n, M_oos);
	double d;
	double **ins_xis_all;
	ins_xis_all = create_double_matrix(n, max_Ms);
	double ***all_fvals;
	all_fvals = create_double_tmatrix(testN, num_Method, Ms);
	double ***all_fvals_;
	all_fvals_ = create_double_tmatrix(testN, num_Method, Ms);
	double ***fvals_valid;
	fvals_valid = create_double_tmatrix(4, num_G, CVK);
	double *fvals_valid2;
	fvals_valid2 = create_double_vector(num_G);
	fea_train_all = create_double_matrix(n, max_Ms);
	for (int k = 0; k < n; k++) {
		for (int s = 0; s < max_Ms; s++) {
			fea_train_all[k][s] = 0.;
		}
	}
	char strM[ENOUGH];
	sprintf(strM, "%d_results_M.dat", n);
	out = Open_File(strM, "a+");
	char strW[ENOUGH];
	sprintf(strW, "%d_results_W.dat", n);
	FILE		*out1;
	out1 = Open_File(strW, "a+");
	for (expk = 0; expk < Expk; expk++) {
		printf("exp is %d\t", expk);
		while (1)
		{
			srand(time(0));
			typedef boost::mt19937 RNGType;
			RNGType rng(time(0));
			Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
			Eigen::internal::scalar_normal_dist_op<double>::rng.seed(time(0)); // Seed the rng

																			   // Define mean and covariance of the distribution
			Eigen::VectorXd mean(n);
			Eigen::MatrixXd covar(n, n);
			mean << 0., 0.;
			covar << 1., 0.,
				0., 1.;
			
			
			Eigen::MatrixXd normTransform(n, n);

			Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

			// We can only use the cholesky decomposition if 
			// the covariance matrix is symmetric, pos-definite.
			// But a covariance matrix might be pos-semi-definite.
			// In that case, we'll go to an EigenSolver
			if (cholSolver.info() == Eigen::Success) {
				// Use cholesky solver
				normTransform = cholSolver.matrixL();
			}
			else {
				// Use eigen solver
				Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
				normTransform = eigenSolver.eigenvectors()
					* eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
			}
			for (int k = 0; k < n; k++) {
				Eigen::MatrixXd samples = (normTransform
					* Eigen::MatrixXd::NullaryExpr(n, 1, randN)).colwise()
					+ mean;
				for (int i = 0; i < n; i++) {
					fea_fix[i] = samples(i, 0);
				}
			}

			boost::uniform_real<> zero_to_five(0, 5);
			boost::variate_generator< RNGType, boost::uniform_real<> >
				dice(rng, zero_to_five);
			for (int i = 0; i < n; i++) {
				bs[i] = dice();
			}

			for (int i = 0; i < n; i++)
				cs[i] = 1.;
			d = 5. * n;


			for (int i = 0; i < n; i++) {
				double sum_fea = fea_fix[i];

				double phi = exp(sum_fea);
				std::default_random_engine generator;
				std::normal_distribution<double> distribution(3., phi);
				for (int s = 0; s < M_oos; s++) {
					do {
						oos_xis[i][s] = distribution(generator);
					} while (oos_xis[i][s] < 0.);
					//printf("oos is %.2lf\t", oos_xis[i][s]);
				}

				//printf("\n");
			}
			//double Gamma = 0.;
			//oos_profit[3] = solveNVDROM(n, M_oos, ubar, Gamma, bs, oos_xis, cs, d, M_oos, oos_xis,prob);
			oos_profit[4] = solveNVSAA(n, M_oos, bs, oos_xis, cs, d, M_oos, oos_xis);

			if (oos_profit[4]>10.) {
				break;
			}
		}

		srand(time(0));
		typedef boost::mt19937 RNGType;
		RNGType rng(time(0));

		Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
		Eigen::internal::scalar_normal_dist_op<double>::rng.seed(time(0)); // Seed the rng

																		   // Define mean and covariance of the distribution
		Eigen::VectorXd mean(n);
		Eigen::MatrixXd covar(n, n);
		mean << 0., 0.;
		covar << 1., 0.,
			0., 1.;
		
		
		Eigen::MatrixXd normTransform(n, n);

		Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

		// We can only use the cholesky decomposition if 
		// the covariance matrix is symmetric, pos-definite.
		// But a covariance matrix might be pos-semi-definite.
		// In that case, we'll go to an EigenSolver
		if (cholSolver.info() == Eigen::Success) {
			// Use cholesky solver
			normTransform = cholSolver.matrixL();
		}
		else {
			// Use eigen solver
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
			normTransform = eigenSolver.eigenvectors()
				* eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
		}
		for (int j = 0; j < max_Ms; j++) {

			Eigen::MatrixXd samples = (normTransform
				* Eigen::MatrixXd::NullaryExpr(n, 1, randN)).colwise()
				+ mean;

			for (int i = 0; i < n; i++) {
				fea_train_all[i][j] = samples(i, 0);
			}

		}


		
		double *upper_u, *lower_u;
		upper_u = create_double_vector(n);
		lower_u = create_double_vector(n);
		for (int i = 0; i < n; i++) {
			upper_u[i] = 0.;
			lower_u[i] = 10000.;

		}

		for (int i = 0; i < n; i++) {
			for (int s = 0; s < max_Ms; s++) {
				double sum_fea = fea_train_all[i][s];

				double phi = exp(sum_fea);
				std::default_random_engine generator;
				std::normal_distribution<double> distribution(3., phi);

				do {
					ins_xis_all[i][s] = distribution(generator);
				} while (ins_xis_all[i][s] < 0.);

			}

			//printf("\n");
		}
		for (int i = 0; i < n; i++) {
			for (int s = 0; s < max_Ms; s++) {
				if (ins_xis_all[i][s] > upper_u[i]) {
					upper_u[i] = ins_xis_all[i][s];
				}
				if (ins_xis_all[i][s] < lower_u[i]) {
					lower_u[i] = ins_xis_all[i][s];
				}
			}
		}
		for (int i = 0; i < n; i++) {
			ubar[i] = upper_u[i];
			uu[i] = lower_u[i];
		}

		free_and_null((char **)&upper_u);
		free_and_null((char **)&lower_u);

		
		for (int mk = 0; mk < Ms; mk++) {
			printf("mk is %d\t", mk);
			int M = M_ins[mk];
			double **ins_xis_;
			ins_xis_ = create_double_matrix(n, M);
			for (int i = 0; i < n; i++) {
				for (int s = 0; s < M; s++) {
					ins_xis_[i][s] = ins_xis_all[i][s];
				}
			}
			double **fea_train_;
			fea_train_ = create_double_matrix(n, M);
			for (int k = 0; k < n; k++) {
				for (int s = 0; s < M; s++) {
					fea_train_[k][s] = fea_train_all[k][s];
				}
			}
			prob_ = create_double_vector(M);
			compute_prob(fea_fix, fea_train_, M, n, mk);
			double *prob;
			prob = create_double_vector(M);
			int ss = 0;
			double sum_prob = 0.;
			for (int s = 0; s < M; s++) {
				if (prob_[s] > 0.00000001) {
					prob[ss] = prob_[s];
					sum_prob += prob_[s];
					ss++;
				}
			}
			for (int s = 0; s < ss; s++) {
				prob[s] = prob[s] / sum_prob;				
			}
			int M_old = M;
			M = ss;
			ss = 0;
			double **ins_xis;
			ins_xis = create_double_matrix(n, M);
			
			for (int s = 0; s < M_old; s++) {
				if (prob_[s] > 0.00000001) {
					for (int i = 0; i < n; i++) {
						ins_xis[i][ss] = ins_xis_[i][s];
					}
					ss++;
				}
			}
			if (M / CVK <= 0) {
				CVK = M;
			}
			double **xis_train, **xis_valid;
			xis_train = create_double_matrix(n, M / CVK);
			xis_valid = create_double_matrix(n, M - M / CVK);
			double *prob_train, *prob_valid;
			prob_train = create_double_vector(M / CVK);
			prob_valid = create_double_vector(M - M / CVK);
			double *prob_train_, *prob_valid_;
			prob_train_ = create_double_vector(M / CVK);
			prob_valid_ = create_double_vector(M - M / CVK);
			int *index_s;
			index_s = create_int_vector(M / CVK);
			int *index_s_valid;
			index_s_valid = create_int_vector(M - M / CVK);
			for (int j = 0; j < num_G; j++) {
				double Gamma = Gammas[j];
				for (int kk = 0; kk < CVK; kk++) {
					int num_s = 0;
					for (int s = 0; s < M / CVK; s++) {
						index_s[s] = s + kk * M / CVK;
					}
					int ii = 0;
					for (int s = 0; s < M - M / CVK; s++) {
						if (kk == 0) {
							index_s_valid[s] = M / CVK + s;

						}
						else {
							if (ii >= M / CVK * kk) {
								index_s_valid[ii] = M / CVK + s;
								ii++;
							}
							else {
								index_s_valid[ii] = s;
								ii++;
							}

						}

					}
					
					for (int i = 0; i < n; i++) {
						for (int s = 0; s < M / CVK; s++) {
							xis_train[i][s] = ins_xis[i][index_s[s]];
						}
					}
					double sum_prob = 0.;
					for (int s = 0; s < M / CVK; s++) {
						prob_train[s] = prob[index_s[s]];
						sum_prob += prob_train[s];
					}
					if (sum_prob < 0.001) {
						break;
					}
					for (int s = 0; s < M / CVK; s++) {
						prob_train[s] = prob_train[s] / sum_prob;
						prob_train_[s] = CVK * 1. / M;
						//printf("prob is %.6lf\t", prob_train[s]);
					}
					for (int i = 0; i < n; i++) {
						for (int s = 0; s < M - M / CVK; s++) {
							xis_valid[i][s] = ins_xis[i][index_s_valid[s]];
						}
					}
					sum_prob = 0.;
					for (int s = 0; s < M - M / CVK; s++) {
						prob_valid[s] = prob[index_s_valid[s]];
						sum_prob += prob_valid[s];
					}
					if (sum_prob < 0.001) {
						break;
					}
					for (int s = 0; s < M - M / CVK; s++) {
						prob_valid[s] = prob_valid[s] / sum_prob;
						prob_valid_[s] = 1. / (M - M / CVK);
					}

					double valid_profit = solveNVDROW(n, M / CVK, ubar, Gamma, bs, xis_train, cs, d, M - M / CVK, xis_valid, prob_train, uu, prob_valid);
					fvals_valid[0][j][kk] = valid_profit;
					//printf("w is %.2lf\t", valid_profit);

					//printf("solv M");
					valid_profit = solveNVDROM(n, M / CVK, ubar, Gamma, bs, xis_train, cs, d, M - M / CVK, xis_valid, prob_train, uu, prob_valid);
					fvals_valid[1][j][kk] = valid_profit;
					//printf("m is %.2lf\t", valid_profit);

					valid_profit = solveNVDROW(n, M / CVK, ubar, Gamma, bs, xis_train, cs, d, M - M / CVK, xis_valid, prob_train_, uu, prob_valid_);
					fvals_valid[2][j][kk] = valid_profit;
					//printf("w is %.2lf\t", valid_profit);

					//printf("solv M");
					valid_profit = solveNVDROM(n, M / CVK, ubar, Gamma, bs, xis_train, cs, d, M - M / CVK, xis_valid, prob_train_, uu, prob_valid_);
					fvals_valid[3][j][kk] = valid_profit;

				}
			}
			double max_fvals = 0.;
			int jstar=0;
			for (int j = 0; j < num_G; j++) {
				fvals_valid2[j] = 0.;
				for (int kk = 0; kk < CVK; kk++) {
					fvals_valid2[j] += fvals_valid[0][j][kk];
				}
				fvals_valid2[j] = fvals_valid2[j] / CVK;
				if (fvals_valid2[j] > max_fvals) {
					max_fvals = fvals_valid2[j];
					jstar = j;
				}
			}
			int jstarM = jstar;
			/**/
			max_fvals = 0.;
			for (int j = 0; j < num_G; j++) {
				fvals_valid2[j] = 0.;
				for (int kk = 0; kk < CVK; kk++) {
					fvals_valid2[j] += fvals_valid[1][j][kk];
				}
				fvals_valid2[j] = fvals_valid2[j] / CVK;
			//	printf("fval is %.2lf\t", fvals_valid2[j]);
				if (fvals_valid2[j] > max_fvals) {
					max_fvals = fvals_valid2[j];
					jstarM = j;
				}
			}
			
			int jstar_=0;
			/**/
			max_fvals = 0.;
			for (int j = 0; j < num_G; j++) {
				fvals_valid2[j] = 0.;
				for (int kk = 0; kk < CVK; kk++) {
					fvals_valid2[j] += fvals_valid[2][j][kk];
				}
				fvals_valid2[j] = fvals_valid2[j] / CVK;
				//	printf("fval is %.2lf\t", fvals_valid2[j]);
				if (fvals_valid2[j] > max_fvals) {
					max_fvals = fvals_valid2[j];
					jstar_ = j;
				}
			}

			int jstarM_ = 0;
			/**/
			max_fvals = 0.;
			for (int j = 0; j < num_G; j++) {
				fvals_valid2[j] = 0.;
				for (int kk = 0; kk < CVK; kk++) {
					fvals_valid2[j] += fvals_valid[3][j][kk];
				}
				fvals_valid2[j] = fvals_valid2[j] / CVK;
				//	printf("fval is %.2lf\t", fvals_valid2[j]);
				if (fvals_valid2[j] > max_fvals) {
					max_fvals = fvals_valid2[j];
					jstarM_ = j;
				}
			}

			printf("j star is %.5lf\t%.5lf\n", Gammas[jstar], Gammas[jstarM]);
			fprintf(out, "%.6lf\t", Gammas[jstarM]);
			fprintf(out1, "%.6lf\t", Gammas[jstar]);
			double saa_gamma = 0.;
			double *prob_oos;
			prob_oos = create_double_vector(M_oos);
			for (int s = 0; s < M_oos; s++) {
				prob_oos[s] = 1. / M_oos;
			}
			double gamma_kde = 0.;
			double *prob_;
			prob_ = create_double_vector(M);
			for (int s = 0; s < M; s++) {
				prob_[s] = 1. / M;
			}
			//printf("n is %d\t", n);
			//oos_profit[0] = solveNVDROM(n, M, ubar, gamma_kde, bs, ins_xis, cs, d, M_oos, oos_xis, prob_, uu, prob_oos);//solveNVKDE(n, M, bs, ins_xis, cs, d, M_oos, oos_xis, prob);
			oos_profit[0] = solveNVSAA(n, M, bs, ins_xis, cs, d, M_oos, oos_xis);//solveNVDROM(n, M, ubar, saa_gamma, bs, ins_xis, cs, d, M_oos, oos_xis, prob,uu);
			oos_profit[2] = solveNVDROW(n, M, ubar, Gammas[jstar], bs, ins_xis, cs, d, M_oos, oos_xis, prob, uu, prob_oos);
			oos_profit[3] = solveNVDROM(n, M, ubar, Gammas[jstarM], bs, ins_xis, cs, d, M_oos, oos_xis, prob, uu, prob_oos);
			oos_profit[1] = solveNVDROW(n, M, ubar, gamma_kde, bs, ins_xis, cs, d, M_oos, oos_xis, prob, uu, prob_oos);//solveNVKDE(n, M, bs, ins_xis, cs, d, M_oos, oos_xis, prob);
			oos_profit[5] = solveNVDROW(n, M, ubar, Gammas[jstar_], bs, ins_xis, cs, d, M_oos, oos_xis, prob_, uu, prob_oos);
			oos_profit[6] = solveNVDROM(n, M, ubar, Gammas[jstarM_], bs, ins_xis, cs, d, M_oos, oos_xis, prob_, uu, prob_oos);
			//oos_profit[1] = solveNVKDE(n, M, bs, ins_xis, cs, d, M_oos, oos_xis, prob);
			printf("saa_ is %.2lf\t", oos_profit[0]);
			printf("kde is %.2lf\t", oos_profit[1]);
			printf("w is %.2lf\t", oos_profit[2]);
			printf("M is %.2lf\t", oos_profit[3]);
			printf("w_ is %.2lf\t", oos_profit[5]);
			printf("M_ is %.2lf\t", oos_profit[6]);
			all_fvals[0][0][mk] = oos_profit[0];
			all_fvals[0][1][mk] = oos_profit[1];
			all_fvals[0][2][mk] = oos_profit[2];
			all_fvals[0][3][mk] = oos_profit[3];
			all_fvals[0][4][mk] = oos_profit[4];
			all_fvals[0][5][mk] = oos_profit[5];
			all_fvals[0][6][mk] = oos_profit[6];
			for (int i = 0; i < n; i++) {
				free_and_null((char **)&ins_xis[i]);
				free_and_null((char **)&xis_train[i]);
				free_and_null((char **)&xis_valid[i]);
				free_and_null((char **)&fea_train_[i]);
			}
			free_and_null((char **)&ins_xis);
			free_and_null((char **)&xis_train);
			free_and_null((char **)&xis_valid);
			free_and_null((char **)&prob_train);
			free_and_null((char **)&prob_valid);
			free_and_null((char **)&prob_train_);
			free_and_null((char **)&prob_valid_);
			free_and_null((char **)&index_s);
			free_and_null((char **)&index_s_valid);
			free_and_null((char **)&fea_train_);
		}


		for (int it = 0; it < testN; it++) {
			for (int i = 0; i < num_Method; i++) {
				for (int s = 0; s < Ms; s++) {
					all_fvals_[it][i][s] = all_fvals[it][i][s];
					all_fvals[it][i][s] = all_fvals[it][i][s] / oos_profit[4];
				}
			}
		}
		for (int i = 0; i < num_Method; i++) {
			for (int s = 0; s < Ms; s++) {
				all_results2_[expk][i][s] = all_fvals_[0][i][s];
				all_results2[expk][i][s] = all_fvals[0][i][s];
			}
		}

	}
	fclose(out);
	fclose(out1);
	char str1[ENOUGH];
	sprintf(str1, "%d_results_oos1000_106d.dat", n);
	FILE		*prog_;
	prog_ = Open_File(str1, "a+");
	for (int i = 0; i < Expk; i++) {
		for (int j = 0; j < num_Method; j++) {
			for (int k = 0; k < Ms; k++) {
				fprintf(prog_, "%.8lf\t", all_results2[i][j][k]);
			}
			fprintf(prog_, "\n");
		}
	}
	fclose(prog_);
	char str2[ENOUGH];
	sprintf(str2, "%d_results_oos10_4_.dat", n);
	prog_ = Open_File(str2, "a+");
	for (int i = 0; i < Expk; i++) {
		for (int j = 0; j < num_Method; j++) {
			for (int k = 0; k < Ms; k++) {
				fprintf(prog_, "%.8lf\t", all_results2_[i][j][k]);
			}
			fprintf(prog_, "\n");
		}
	}
	fclose(prog_);
	//}
	free_and_null((char **)&ubar);
	free_and_null((char **)&bs);
	free_and_null((char **)&cs);
	free_and_null((char **)&Gammas);
	for (int i = 0; i < n; i++) {
		free_and_null((char **)&oos_xis[i]);
		free_and_null((char **)&ins_xis_all[i]);
	}
	free_and_null((char **)&oos_xis);
	free_and_null((char **)&ins_xis_all);
	for (int i = 0; i < testN; i++) {
		for (int j = 0; j < 5; j++) {
			free_and_null((char **)&all_fvals[i][j]);
			free_and_null((char **)&all_fvals_[i][j]);
		}
		free_and_null((char **)&all_fvals[i]);
		free_and_null((char **)&all_fvals_[i]);
	}
	free_and_null((char **)&all_fvals);
	free_and_null((char **)&all_fvals_);
	free_and_null((char **)&M_ins);
	for (int i = 0; i < Expk; i++) {
		for (int j = 0; j < 5; j++) {
			free_and_null((char **)&all_results2[i][j]);
			free_and_null((char **)&all_results2_[i][j]);
		}
		free_and_null((char **)&all_results2[i]);
		free_and_null((char **)&all_results2_[i]);
	}
	free_and_null((char **)&all_results2);
	free_and_null((char **)&all_results2_);
	free_and_null((char **)&oos_profit);
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < num_G; j++) {
			free_and_null((char **)&fvals_valid[2][j]);
		}
		free_and_null((char **)&fvals_valid[2]);
	}
	free_and_null((char **)&fvals_valid);
	free_and_null((char **)&fvals_valid2);
}
void compute_prob(double *fea, double **fea_train, int M, int n, int mk) {
	double *sort_d;
	sort_d = create_double_vector(M);
	for (int nn = 0; nn < M; nn++) {
		sort_d[nn] = 0.;
		for (int k = 0; k < n; k++) {
			sort_d[nn] += pow((fea[k] - fea_train[k][nn]) / h[mk], 2);
		}
	}


	double val_sum = 0.;
	for (int nn = 0; nn < M; nn++) {
		val_sum += exp(-sort_d[nn] / 2) / (pow(2 * 3.1416, 0.5)*h[mk]);
	}

	//printf("sum w is %.2lf\t", val_sum);

	double sum_w = 0.;
	for (int nn = 0; nn < M; nn++) {
		prob_[nn] = 0.;
		double weight = exp(-sort_d[nn] / 2) / (pow(2 * 3.1416, 0.5)*h[mk]);
		prob_[nn] = weight / val_sum;
		//printf("prob is %.4lf\t", prob[nn]);
	}
	free_and_null((char **)&sort_d);
}
FILE *Open_File(const char *name, const char *mode)  // This function opens the file "name" under mode "mode" (reading, writing, etc)
{
	FILE *file;
	int OK;
	if ((file = fopen(name, mode)) == NULL) {
		OK = 1;
		printf("\nError: File cannot be opened %d\n", OK);

	}
	return file;
}
