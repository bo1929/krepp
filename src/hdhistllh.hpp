#ifndef _HDHISTLLH_H
#define _HDHISTLLH_H

#include <cstdint>
#include <lbfgsb.h>

namespace optimize {
  const Vector d_init{{0.01}};             // Initial guess
  const Vector d_lb{{0.0000000000000001}}; // Lower bounds on d_llh
  const Vector d_ub{{0.5}};                // Upper bounds on d_llh

  class HDistHistLLH : public Function
  {
  private:
    uint32_t hdist_th;
    uint32_t nmers;
    uint32_t k;
    uint32_t h;
    double rho = 1;
    uint32_t* mc_ptr = nullptr;
    std::vector<uint64_t> binom_coef_k;

  public:
    HDistHistLLH(uint32_t h, uint32_t k, uint32_t hdist_th, uint32_t nmers)
      : h(h)
      , k(k)
      , hdist_th(hdist_th)
      , nmers(nmers)
    {
      binom_coef_k.resize(k + 1);
      binom_coef_k[0] = 1;
      for (uint32_t i = 1; i <= k; ++i) {
        uint64_t vc = k - i + 1;
        uint64_t vp;
        for (uint64_t j = 1; j < i; ++j) {
          vc = (vp = vc) * (k - i + 1 + j) / (j + 1);
        }
        binom_coef_k[i] = vc;
      }
    }

    double prob_collide(const uint32_t x)
    {
      return 1.0 - (1.0 - pow(1.0 - static_cast<float>(x) / static_cast<float>(k), h));
    }

    double prob_mutate(const double d, const uint32_t x)
    {
      return pow((1.0 - d), (k - x)) * pow(d, x) * binom_coef_k[x];
    }

    double prob_miss(const double d)
    {
      double p = 0;
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        p += (1.0 - prob_collide(x)) * prob_mutate(d, x);
      }
      for (uint32_t x = hdist_th + 1; x <= k; ++x) {
        p += prob_mutate(d, x);
      }
      return rho * p + 1.0 - rho;
    }

    double prob_hit(const double d, const uint32_t x)
    {
      return rho * prob_collide(x) * prob_mutate(d, x);
    }

    Scalar computeValue(const Vector& d) override
    {
      uint32_t uc = nmers;
      double sum = 0;
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        uc -= *(mc_ptr + x);
        sum -= log(prob_hit(d[0], x)) * (*(mc_ptr + x));
      }
      return sum - log(prob_miss(d[0])) * uc;
    }

    Vector computeGradient(const Vector& d) override
    {
      uint32_t uc = nmers;
      double grad = 0;
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        uc -= *(mc_ptr + x);
        grad -= (x / d[0] - (k - x) / (1.0 - d[0])) * (*(mc_ptr + x));
      }
      double mg_n = 0;
      for (uint32_t x = 0; x <= k; ++x) {
        if (x <= hdist_th) {
          mg_n += ((1.0 - prob_collide(x)) * prob_mutate(d[0], x) * (static_cast<double>(x) - k) /
                   (1.0 - d[0]));
          mg_n += ((1.0 - prob_collide(x)) * prob_mutate(d[0], x) * (x) / d[0]);
        } else {
          mg_n += (prob_mutate(d[0], x) * (static_cast<double>(x) - k) / (1.0 - d[0]));
          mg_n += (prob_mutate(d[0], x) * (x) / d[0]);
        }
      }
      mg_n = rho * uc * mg_n;
      double mg_d = prob_miss(d[0]);
      Vector grad_v(d.size());
      grad_v[0] = grad - mg_n / mg_d;
      /* #pragma omp critical */
      /*       { */
      /*         std::cout << "d : " << d[0] << std::endl; */
      /*         std::cout << "hllhg : " << grad << std::endl; */
      /*         std::cout << "grad : " << grad - mg_n / mg_d << std::endl; */
      /*         std::cout << "mg_n : " << mg_n << std::endl; */
      /*         std::cout << "mg_d : " << mg_d << std::endl; */
      /*       } */
      return grad_v;
    }
    void set_mc(uint32_t* curr_mc_ptr) { mc_ptr = curr_mc_ptr; }
    void set_ro(double curr_ro) { rho = curr_ro; }
  };
}
/* Vector grad(d.size()); */
/* uint32_t uc = nmers; */
/* grad[0] = 0; */
/* for (uint32_t i = 0; i <= hdist_th; ++i) { */
/*   uc -= *(mc_ptr + i); */
/*   grad[0] -= (i / d[0] - (k - i) / (1 - d[0])) * (*(mc_ptr + i)); */
/* } */
/* // Better modeling of Hamming distance threshold: */
/* double tmp = 0; */
/* double f_bth = 0, g_bth = 0; */
/* for (uint32_t i = 0; i <= hdist_th; ++i) { */
/*   tmp = pow(d[0], i) * pow(1 - d[0], k - i) * binom_coef_k[i] * (1 - pow(1 - d[0], h)); */
/*   f_bth += tmp; */
/*   g_bth += (tmp * i) / d[0] + (tmp * i - tmp * k) / (1 - d[0]); */
/* } */
/* double f_ath = 0, g_ath = 0; */
/* for (uint32_t i = hdist_th + 1; i <= k; ++i) { */
/*   tmp = pow(d[0], i) * pow(1 - d[0], k - i) * binom_coef_k[i]; */
/*   f_ath += tmp; */
/*   g_ath += (tmp * i) / d[0] + (tmp * i - tmp * k) / (1 - d[0]); */
/* } */
/* double numer = uc * rho * (g_ath + g_bth); */
/* double denom = (rho * (f_ath + f_bth) + 1 - rho); */
/* grad[0] -= numer / denom; */
/* // Simple likelihood for miss: */
/* /1* double ps = pow(1 - d[0], h - 1); *1/ */
/* /1* grad[0] -= (h * uc * rho * ps) / (1 - rho * ps * (1 - d[0])); *1/ */
/* return grad; */

#endif
