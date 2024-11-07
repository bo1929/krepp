#ifndef _HDHISTLLH_H
#define _HDHISTLLH_H

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
    std::vector<double> prob_elude_v;

  public:
    double prob_collide(const uint32_t x)
    {
      return 1.0 - (1.0 - pow(1.0 - static_cast<double>(x) / static_cast<double>(k), h));
    }

    double prob_elude(const uint32_t x)
    {
      return 1.0 - pow(1.0 - static_cast<double>(x) / static_cast<double>(k), h);
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

    HDistHistLLH() {}

    HDistHistLLH(uint32_t h, uint32_t k, uint32_t hdist_th)
      : h(h)
      , k(k)
      , hdist_th(hdist_th)
    {
      binom_coef_k.resize(k + 1);
      prob_elude_v.resize(hdist_th + 1);
      binom_coef_k[0] = 1;
      for (uint32_t i = 1; i <= k; ++i) {
        uint64_t vc = k - i + 1;
        uint64_t vp;
        for (uint64_t j = 1; j < i; ++j) {
          vc = (vp = vc) * (k - i + 1 + j) / (j + 1);
        }
        binom_coef_k[i] = vc;
      }
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        prob_elude_v[x] = prob_elude(x);
      }
    }

    /* Scalar computeValue(const Vector& d) override */
    /* { */
    /*   double sum = 0; */
    /*   for (uint32_t x = 0; x <= hdist_th; ++x) { */
    /*     // sum -= (log(rho) + log(prob_collide(x)) + log(prob_mutate(d, x))) * (*(mc_ptr + x)); */
    /*     // sum -= log(prob_hit(d, x)) * (*(mc_ptr + x)); */
    /*     sum -= ((k - x) * log(1.0 - d[0]) + x * log(d[0])) * (*(mc_ptr + x)); */
    /*   } */
    /*   return sum - log(prob_miss(d[0])) * uc; */
    /* } */

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
          lv_m += prob_elude_v[x] * (powdc * binom_coef_k[x]);
        } else {
          lv_m += powdc * binom_coef_k[x];
        }
        powdc *= dratio;
      }
      return sum - log(rho * lv_m + 1.0 - rho) * uc;
    }

    /* Vector computeGradient(const Vector& d) override */
    /* { */
    /*   Vector grad(d.size()); */
    /*   grad[0] = 0; */
    /*   double mg_n = 0, mg_d = 0; */
    /*   double prob_dxc, prob_dxm, dhitx; */
    /*   for (uint32_t x = 0; x <= k; ++x) { */
    /*     dhitx = (x / d[0] - (k - x) / (1.0 - d[0])); */
    /*     prob_dxm = prob_mutate(d[0], x); */
    /*     prob_dxc = prob_dxm * dhitx; */
    /*     if (x <= hdist_th) { */
    /*       grad[0] -= dhitx * (*(mc_ptr + x)); */
    /*       mg_n += prob_elude_v[x] * prob_dxc; */
    /*       mg_d += prob_elude_v[x] * prob_dxm; */
    /*     } else { */
    /*       mg_n += prob_dxc; */
    /*       mg_d += prob_dxm; */
    /*     } */
    /*   } */
    /*   mg_n = rho * uc * mg_n; */
    /*   mg_d = rho * mg_d + 1.0 - rho; */
    /*   grad[0] = grad[0] - mg_n / mg_d; */
    /*   return grad; */
    /* } */

    /* double operator()(double const& d) */
    /* { */
    /*   double grad = 0; */
    /*   double mg_n = 0, mg_d = 0; */
    /*   double prob_dxc, prob_dxm, dhitx; */
    /*   for (uint32_t x = 0; x <= k; ++x) { */
    /*     dhitx = (x / d - (k - x) / (1.0 - d)); */
    /*     prob_dxm = prob_mutate(d, x); */
    /*     prob_dxc = prob_dxm * dhitx; */
    /*     if (x <= hdist_th) { */
    /*       grad -= dhitx * (*(mc_ptr + x)); */
    /*       mg_n += prob_elude_v[x] * prob_dxc; */
    /*       mg_d += prob_elude_v[x] * prob_dxm; */
    /*     } else { */
    /*       mg_n += prob_dxc; */
    /*       mg_d += prob_dxm; */
    /*     } */
    /*   } */
    /*   mg_n = rho * uc * mg_n; */
    /*   mg_d = rho * mg_d + 1.0 - rho; */
    /*   grad = grad - mg_n / mg_d; */
    /*   return grad; */
    /* } */

    void set_parameters(double* curr_mc_ptr, double curr_uc, double curr_rho)
    {
      mc_ptr = curr_mc_ptr;
      uc = curr_uc;
      rho = curr_rho;
    }
  };
}

#endif
