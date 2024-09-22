#ifndef _HDHISTLLH_H
#define _HDHISTLLH_H

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
      return 1.0 - (1.0 - pow(1.0 - static_cast<double>(x) / static_cast<double>(k), h));
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
      Vector grad(d.size());
      uint32_t uc = nmers;
      grad[0] = 0;
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        uc -= *(mc_ptr + x);
        grad[0] -= (x / d[0] - (k - x) / (1.0 - d[0])) * (*(mc_ptr + x));
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
      grad[0] = grad[0] - mg_n / mg_d;
      return grad;
    }
    void set_mc(uint32_t* curr_mc_ptr) { mc_ptr = curr_mc_ptr; }
    void set_ro(double curr_ro) { rho = curr_ro; }
  };
}

#endif
