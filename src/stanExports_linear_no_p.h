// Generated by rstantools.  Do not edit by hand.

/*
    treeRingSplinesStan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    treeRingSplinesStan is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with treeRingSplinesStan.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_linear_no_p_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_linear_no_p");
    reader.add_event(78, 76, "end", "model_linear_no_p");
    return reader;
}
#include <stan_meta_header.hpp>
class model_linear_no_p
  : public stan::model::model_base_crtp<model_linear_no_p> {
private:
        int K;
        int n;
        matrix_d X;
        vector_d y;
        int n_tree;
        std::vector<int> tree_idx;
        int n_pred;
        matrix_d X_pred;
        std::vector<int> tree_idx_pred;
public:
    model_linear_no_p(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_linear_no_p(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_linear_no_p_namespace::model_linear_no_p";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            check_greater_or_equal(function__, "K", K, 0);
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "n", "int", context__.to_vec());
            n = int(0);
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            n = vals_i__[pos__++];
            check_greater_or_equal(function__, "n", n, 0);
            current_statement_begin__ = 6;
            validate_non_negative_index("X", "n", n);
            validate_non_negative_index("X", "K", K);
            context__.validate_dims("data initialization", "X", "matrix_d", context__.to_vec(n,K));
            X = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, K);
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_2_max__ = K;
            size_t X_j_1_max__ = n;
            for (size_t j_2__ = 0; j_2__ < X_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                    X(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("y", "n", n);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(n));
            y = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 9;
            context__.validate_dims("data initialization", "n_tree", "int", context__.to_vec());
            n_tree = int(0);
            vals_i__ = context__.vals_i("n_tree");
            pos__ = 0;
            n_tree = vals_i__[pos__++];
            check_greater_or_equal(function__, "n_tree", n_tree, 0);
            current_statement_begin__ = 10;
            validate_non_negative_index("tree_idx", "n", n);
            context__.validate_dims("data initialization", "tree_idx", "int", context__.to_vec(n));
            tree_idx = std::vector<int>(n, int(0));
            vals_i__ = context__.vals_i("tree_idx");
            pos__ = 0;
            size_t tree_idx_k_0_max__ = n;
            for (size_t k_0__ = 0; k_0__ < tree_idx_k_0_max__; ++k_0__) {
                tree_idx[k_0__] = vals_i__[pos__++];
            }
            size_t tree_idx_i_0_max__ = n;
            for (size_t i_0__ = 0; i_0__ < tree_idx_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "tree_idx[i_0__]", tree_idx[i_0__], 1);
            }
            current_statement_begin__ = 12;
            context__.validate_dims("data initialization", "n_pred", "int", context__.to_vec());
            n_pred = int(0);
            vals_i__ = context__.vals_i("n_pred");
            pos__ = 0;
            n_pred = vals_i__[pos__++];
            check_greater_or_equal(function__, "n_pred", n_pred, 0);
            current_statement_begin__ = 13;
            validate_non_negative_index("X_pred", "n_pred", n_pred);
            validate_non_negative_index("X_pred", "K", K);
            context__.validate_dims("data initialization", "X_pred", "matrix_d", context__.to_vec(n_pred,K));
            X_pred = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n_pred, K);
            vals_r__ = context__.vals_r("X_pred");
            pos__ = 0;
            size_t X_pred_j_2_max__ = K;
            size_t X_pred_j_1_max__ = n_pred;
            for (size_t j_2__ = 0; j_2__ < X_pred_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_pred_j_1_max__; ++j_1__) {
                    X_pred(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 14;
            validate_non_negative_index("tree_idx_pred", "n_pred", n_pred);
            context__.validate_dims("data initialization", "tree_idx_pred", "int", context__.to_vec(n_pred));
            tree_idx_pred = std::vector<int>(n_pred, int(0));
            vals_i__ = context__.vals_i("tree_idx_pred");
            pos__ = 0;
            size_t tree_idx_pred_k_0_max__ = n_pred;
            for (size_t k_0__ = 0; k_0__ < tree_idx_pred_k_0_max__; ++k_0__) {
                tree_idx_pred[k_0__] = vals_i__[pos__++];
            }
            size_t tree_idx_pred_i_0_max__ = n_pred;
            for (size_t i_0__ = 0; i_0__ < tree_idx_pred_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "tree_idx_pred[i_0__]", tree_idx_pred[i_0__], 1);
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 19;
            num_params_r__ += 1;
            current_statement_begin__ = 20;
            validate_non_negative_index("beta0_t_tilde", "n_tree", n_tree);
            num_params_r__ += (1 * n_tree);
            current_statement_begin__ = 21;
            num_params_r__ += 1;
            current_statement_begin__ = 24;
            validate_non_negative_index("beta", "K", K);
            num_params_r__ += K;
            current_statement_begin__ = 27;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_linear_no_p() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 19;
        if (!(context__.contains_r("mu_beta0")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu_beta0 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu_beta0");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu_beta0", "double", context__.to_vec());
        double mu_beta0(0);
        mu_beta0 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(mu_beta0);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu_beta0: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 20;
        if (!(context__.contains_r("beta0_t_tilde")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta0_t_tilde missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta0_t_tilde");
        pos__ = 0U;
        validate_non_negative_index("beta0_t_tilde", "n_tree", n_tree);
        context__.validate_dims("parameter initialization", "beta0_t_tilde", "double", context__.to_vec(n_tree));
        std::vector<double> beta0_t_tilde(n_tree, double(0));
        size_t beta0_t_tilde_k_0_max__ = n_tree;
        for (size_t k_0__ = 0; k_0__ < beta0_t_tilde_k_0_max__; ++k_0__) {
            beta0_t_tilde[k_0__] = vals_r__[pos__++];
        }
        size_t beta0_t_tilde_i_0_max__ = n_tree;
        for (size_t i_0__ = 0; i_0__ < beta0_t_tilde_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_unconstrain(beta0_t_tilde[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta0_t_tilde: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 21;
        if (!(context__.contains_r("s_beta0_t")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable s_beta0_t missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("s_beta0_t");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "s_beta0_t", "double", context__.to_vec());
        double s_beta0_t(0);
        s_beta0_t = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, s_beta0_t);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable s_beta0_t: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 24;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "K", K);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(K));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(K);
        size_t beta_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 27;
        if (!(context__.contains_r("sigma_y")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma_y missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma_y");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma_y", "double", context__.to_vec());
        double sigma_y(0);
        sigma_y = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma_y);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma_y: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 19;
            local_scalar_t__ mu_beta0;
            (void) mu_beta0;  // dummy to suppress unused var warning
            if (jacobian__)
                mu_beta0 = in__.scalar_constrain(lp__);
            else
                mu_beta0 = in__.scalar_constrain();
            current_statement_begin__ = 20;
            std::vector<local_scalar_t__> beta0_t_tilde;
            size_t beta0_t_tilde_d_0_max__ = n_tree;
            beta0_t_tilde.reserve(beta0_t_tilde_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < beta0_t_tilde_d_0_max__; ++d_0__) {
                if (jacobian__)
                    beta0_t_tilde.push_back(in__.scalar_constrain(lp__));
                else
                    beta0_t_tilde.push_back(in__.scalar_constrain());
            }
            current_statement_begin__ = 21;
            local_scalar_t__ s_beta0_t;
            (void) s_beta0_t;  // dummy to suppress unused var warning
            if (jacobian__)
                s_beta0_t = in__.scalar_lb_constrain(0, lp__);
            else
                s_beta0_t = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 24;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(K, lp__);
            else
                beta = in__.vector_constrain(K);
            current_statement_begin__ = 27;
            local_scalar_t__ sigma_y;
            (void) sigma_y;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_y = in__.scalar_lb_constrain(0, lp__);
            else
                sigma_y = in__.scalar_lb_constrain(0);
            // transformed parameters
            current_statement_begin__ = 32;
            validate_non_negative_index("beta0_t", "n_tree", n_tree);
            std::vector<local_scalar_t__> beta0_t(n_tree, local_scalar_t__(0));
            stan::math::initialize(beta0_t, DUMMY_VAR__);
            stan::math::fill(beta0_t, DUMMY_VAR__);
            current_statement_begin__ = 33;
            validate_non_negative_index("mu", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu(n);
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 35;
            for (int t = 1; t <= n_tree; ++t) {
                current_statement_begin__ = 36;
                stan::model::assign(beta0_t, 
                            stan::model::cons_list(stan::model::index_uni(t), stan::model::nil_index_list()), 
                            (mu_beta0 + (s_beta0_t * get_base1(beta0_t_tilde, t, "beta0_t_tilde", 1))), 
                            "assigning variable beta0_t");
            }
            current_statement_begin__ = 39;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 40;
                stan::model::assign(mu, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (get_base1(beta0_t, get_base1(tree_idx, i, "tree_idx", 1), "beta0_t", 1) + multiply(get_base1(X, i, "X", 1), beta)), 
                            "assigning variable mu");
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 32;
            size_t beta0_t_k_0_max__ = n_tree;
            for (size_t k_0__ = 0; k_0__ < beta0_t_k_0_max__; ++k_0__) {
                if (stan::math::is_uninitialized(beta0_t[k_0__])) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: beta0_t" << "[" << k_0__ << "]";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable beta0_t: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 33;
            size_t mu_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(mu(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: mu" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable mu: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 48;
            lp_accum__.add(normal_log<propto__>(mu_beta0, 0, 100));
            current_statement_begin__ = 49;
            lp_accum__.add(cauchy_log<propto__>(s_beta0_t, 0, 2.5));
            current_statement_begin__ = 52;
            lp_accum__.add(normal_log<propto__>(beta, 0, 100));
            current_statement_begin__ = 55;
            lp_accum__.add(cauchy_log<propto__>(sigma_y, 0, 2.5));
            current_statement_begin__ = 58;
            lp_accum__.add(normal_log<propto__>(y, mu, sigma_y));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("mu_beta0");
        names__.push_back("beta0_t_tilde");
        names__.push_back("s_beta0_t");
        names__.push_back("beta");
        names__.push_back("sigma_y");
        names__.push_back("beta0_t");
        names__.push_back("mu");
        names__.push_back("y_rep");
        names__.push_back("y_pred");
        names__.push_back("log_lik");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_tree);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_tree);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_pred);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_linear_no_p_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double mu_beta0 = in__.scalar_constrain();
        vars__.push_back(mu_beta0);
        std::vector<double> beta0_t_tilde;
        size_t beta0_t_tilde_d_0_max__ = n_tree;
        beta0_t_tilde.reserve(beta0_t_tilde_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < beta0_t_tilde_d_0_max__; ++d_0__) {
            beta0_t_tilde.push_back(in__.scalar_constrain());
        }
        size_t beta0_t_tilde_k_0_max__ = n_tree;
        for (size_t k_0__ = 0; k_0__ < beta0_t_tilde_k_0_max__; ++k_0__) {
            vars__.push_back(beta0_t_tilde[k_0__]);
        }
        double s_beta0_t = in__.scalar_lb_constrain(0);
        vars__.push_back(s_beta0_t);
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_constrain(K);
        size_t beta_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        double sigma_y = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma_y);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 32;
            validate_non_negative_index("beta0_t", "n_tree", n_tree);
            std::vector<double> beta0_t(n_tree, double(0));
            stan::math::initialize(beta0_t, DUMMY_VAR__);
            stan::math::fill(beta0_t, DUMMY_VAR__);
            current_statement_begin__ = 33;
            validate_non_negative_index("mu", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> mu(n);
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 35;
            for (int t = 1; t <= n_tree; ++t) {
                current_statement_begin__ = 36;
                stan::model::assign(beta0_t, 
                            stan::model::cons_list(stan::model::index_uni(t), stan::model::nil_index_list()), 
                            (mu_beta0 + (s_beta0_t * get_base1(beta0_t_tilde, t, "beta0_t_tilde", 1))), 
                            "assigning variable beta0_t");
            }
            current_statement_begin__ = 39;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 40;
                stan::model::assign(mu, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (get_base1(beta0_t, get_base1(tree_idx, i, "tree_idx", 1), "beta0_t", 1) + multiply(get_base1(X, i, "X", 1), beta)), 
                            "assigning variable mu");
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t beta0_t_k_0_max__ = n_tree;
                for (size_t k_0__ = 0; k_0__ < beta0_t_k_0_max__; ++k_0__) {
                    vars__.push_back(beta0_t[k_0__]);
                }
                size_t mu_j_1_max__ = n;
                for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                    vars__.push_back(mu(j_1__));
                }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 63;
            validate_non_negative_index("y_rep", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> y_rep(n);
            stan::math::initialize(y_rep, DUMMY_VAR__);
            stan::math::fill(y_rep, DUMMY_VAR__);
            current_statement_begin__ = 64;
            validate_non_negative_index("y_pred", "n_pred", n_pred);
            Eigen::Matrix<double, Eigen::Dynamic, 1> y_pred(n_pred);
            stan::math::initialize(y_pred, DUMMY_VAR__);
            stan::math::fill(y_pred, DUMMY_VAR__);
            current_statement_begin__ = 65;
            validate_non_negative_index("log_lik", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> log_lik(n);
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 67;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 68;
                stan::model::assign(y_rep, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            normal_rng((get_base1(beta0_t, get_base1(tree_idx, i, "tree_idx", 1), "beta0_t", 1) + multiply(get_base1(X, i, "X", 1), beta)), sigma_y, base_rng__), 
                            "assigning variable y_rep");
                current_statement_begin__ = 69;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            normal_log(get_base1(y, i, "y", 1), (get_base1(beta0_t, get_base1(tree_idx, i, "tree_idx", 1), "beta0_t", 1) + multiply(get_base1(X, i, "X", 1), beta)), sigma_y), 
                            "assigning variable log_lik");
            }
            current_statement_begin__ = 71;
            for (int i = 1; i <= n_pred; ++i) {
                current_statement_begin__ = 72;
                stan::model::assign(y_pred, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            normal_rng((get_base1(beta0_t, get_base1(tree_idx_pred, i, "tree_idx_pred", 1), "beta0_t", 1) + multiply(get_base1(X_pred, i, "X_pred", 1), beta)), sigma_y, base_rng__), 
                            "assigning variable y_pred");
            }
            // validate, write generated quantities
            current_statement_begin__ = 63;
            size_t y_rep_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < y_rep_j_1_max__; ++j_1__) {
                vars__.push_back(y_rep(j_1__));
            }
            current_statement_begin__ = 64;
            size_t y_pred_j_1_max__ = n_pred;
            for (size_t j_1__ = 0; j_1__ < y_pred_j_1_max__; ++j_1__) {
                vars__.push_back(y_pred(j_1__));
            }
            current_statement_begin__ = 65;
            size_t log_lik_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
                vars__.push_back(log_lik(j_1__));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_linear_no_p";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu_beta0";
        param_names__.push_back(param_name_stream__.str());
        size_t beta0_t_tilde_k_0_max__ = n_tree;
        for (size_t k_0__ = 0; k_0__ < beta0_t_tilde_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta0_t_tilde" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "s_beta0_t";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_y";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t beta0_t_k_0_max__ = n_tree;
            for (size_t k_0__ = 0; k_0__ < beta0_t_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta0_t" << '.' << k_0__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t mu_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        size_t y_rep_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < y_rep_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_rep" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t y_pred_j_1_max__ = n_pred;
        for (size_t j_1__ = 0; j_1__ < y_pred_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_pred" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t log_lik_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu_beta0";
        param_names__.push_back(param_name_stream__.str());
        size_t beta0_t_tilde_k_0_max__ = n_tree;
        for (size_t k_0__ = 0; k_0__ < beta0_t_tilde_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta0_t_tilde" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "s_beta0_t";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_y";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t beta0_t_k_0_max__ = n_tree;
            for (size_t k_0__ = 0; k_0__ < beta0_t_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta0_t" << '.' << k_0__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t mu_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        size_t y_rep_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < y_rep_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_rep" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t y_pred_j_1_max__ = n_pred;
        for (size_t j_1__ = 0; j_1__ < y_pred_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_pred" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t log_lik_j_1_max__ = n;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_linear_no_p_namespace::model_linear_no_p stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
