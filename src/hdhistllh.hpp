#ifndef _HDHISTLLH_H
#define _HDHISTLLH_H

#include <cmath>

// TODO: Make this templated.
namespace optimize {
  class HDistHistLLH
  {
  private:
    uint32_t hdist_th;
    uint32_t k;
    uint32_t h;
    double rho = 1;
    double uc;
    double* mc_ptr = nullptr;
    std::vector<uint64_t> binom_coef_k;
    std::vector<uint64_t> binom_coef_hnk;

  public:
    double prob_elude(const uint32_t x)
    {
      /* return (1.0 - pow(1.0 - static_cast<double>(x) / static_cast<double>(k), h)); */
      return 1.0 - static_cast<double>(binom_coef_hnk[x]) / static_cast<double>(binom_coef_k[x]);
    }

    double prob_collide(const uint32_t x)
    {
      /* return pow(1.0 - static_cast<double>(x) / static_cast<double>(k), h); */
      return static_cast<double>(binom_coef_hnk[x]) / static_cast<double>(binom_coef_k[x]);
    }

    double prob_mutate(const double d, const uint32_t x) { return pow((1.0 - d), (k - x)) * pow(d, x) * binom_coef_k[x]; }

    double prob_miss(const double d)
    {
      double p = 0;
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        p += prob_elude(x) * prob_mutate(d, x);
      }
      for (uint32_t x = hdist_th + 1; x <= k; ++x) {
        p += prob_mutate(d, x);
      }
      return rho * p + 1.0 - rho;
    }

    double prob_hit(const double d, const uint32_t x) { return rho * prob_collide(x) * prob_mutate(d, x); }

    HDistHistLLH() {}

    HDistHistLLH(uint32_t h, uint32_t k, uint32_t hdist_th)
      : h(h)
      , k(k)
      , hdist_th(hdist_th)
    {
      uint64_t vc = 1;
      uint32_t nh = k - h;
      binom_coef_k.resize(k + 1);
      binom_coef_hnk.resize(hdist_th + 1);
      binom_coef_k[0] = 1;
      binom_coef_hnk[0] = 0;
      for (int32_t i = 0; i < k; ++i) {
        binom_coef_k[i + 1] = (binom_coef_k[i] * (k - i)) / (i + 1);
      }
      for (uint32_t i = 1; i <= hdist_th; ++i) {
        vc = (vc * (nh - i + 1)) / i;
        binom_coef_hnk[i] = binom_coef_k[i] - vc;
      }
    }

    double operator()(double const& d)
    {
      double sum = 0.0, lv_m = 0.0;
      double powdc = pow((1.0 - d), k);
      double logdn = log(1.0 - d);
      double logdp = log(d) - logdn;
      logdn *= k;
      double dratio = d / (1.0 - d);
      for (uint32_t x = 0; x <= k; ++x) {
        if (x <= hdist_th) {
          sum -= (logdn + x * logdp) * (*(mc_ptr + x));
          lv_m += binom_coef_hnk[x] * powdc;
        } else {
          lv_m += powdc * binom_coef_k[x];
        }
        powdc *= dratio;
      }
      return sum - log(rho * lv_m + 1.0 - rho) * uc;
    }

    void set_parameters(double* curr_mc_ptr, double curr_uc, double curr_rho)
    {
      mc_ptr = curr_mc_ptr;
      uc = curr_uc;
      rho = curr_rho;
    }
  };
}

#endif
