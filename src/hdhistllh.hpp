#ifndef _HDHISTLLH_H
#define _HDHISTLLH_H

#include <lbfgsb.h>

namespace optimize {
  const Vector d_init{{0.01}}; // Initial guess
  const Vector d_lb{{0.0}};    // Lower bounds on d_llh
  const Vector d_ub{{0.5}};    // Upper bounds on d_llh

  class HDistHistLLH : public Function
  {
  private:
    double ro = 1;
    uint32_t* mc_ptr = nullptr;
    uint32_t hdist_th;
    uint32_t nmers;
    uint32_t k;
    uint32_t h;

  public:
    HDistHistLLH(uint32_t h, uint32_t k, uint32_t hdist_th, uint32_t nmers)
      : h(h)
      , k(k)
      , hdist_th(hdist_th)
      , nmers(nmers)
    {}
    double operator()(const Vector& d, Vector& grad)
    {
      double fx = 0.0;
      uint32_t uc = nmers;
      double ps = pow(1 - d[0], h - 1);
      grad[0] = 0;
      for (uint32_t i = 0; i <= hdist_th; ++i) {
        uc -= *(mc_ptr + i);
        fx -= (log(1 - d[0]) * (k - i) + (log(d[0]) * i)) * (*(mc_ptr + i));
        grad[0] -= (i / d[0] - (k - i) / (1 - d[0])) * (*(mc_ptr + i));
      }
      fx -= log(ro * (1 - ps * (1 - d[0])) + (1 - ro)) * uc;
      grad[0] -= (h * uc * ro * ps) / (1 - ro * ps * (1 - d[0]));
      return fx;
    }

    Scalar computeValue(const Vector& d) override
    {
      double fx = 0.0;
      uint32_t uc = nmers;
      for (uint32_t i = 0; i <= hdist_th; ++i) {
        uc -= *(mc_ptr + i);
        fx -= (log(1 - d[0]) * (k - i) + (log(d[0]) * i)) * (*(mc_ptr + i));
      }
      fx -= log(ro * (1 - pow(1 - d[0], h - 1) * (1 - d[0])) + (1 - ro)) * uc;
      return fx;
    }

    Vector computeGradient(const Vector& d) override
    {
      Vector grad(d.size());
      uint32_t uc = nmers;
      double ps = pow(1 - d[0], h - 1);
      grad[0] = 0;
      for (uint32_t i = 0; i <= hdist_th; ++i) {
        uc -= *(mc_ptr + i);
        grad[0] -= (i / d[0] - (k - i) / (1 - d[0])) * (*(mc_ptr + i));
      }
      grad[0] -= (h * uc * ro * ps) / (1 - ro * ps * (1 - d[0]));
      return grad;
    }
    void set_mc(uint32_t* curr_mc_ptr) { mc_ptr = curr_mc_ptr; }
    void set_ro(double curr_ro) { ro = curr_ro; }
  };
}

#endif
