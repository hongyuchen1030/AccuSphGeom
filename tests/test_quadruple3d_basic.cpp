#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

#include "spip/predicates/quadruple3d.hpp"

using spip::predicates::Sign;

static const char* to_cstr(Sign s) {
  switch(s) {
    case Sign::Positive: return "Positive";
    case Sign::Negative: return "Negative";
    default:             return "Zero";
  }
}

// A simple helper to call the generated filter directly (for debug)
static int filter_only(const double* a, const double* b, const double* c, const double* d) {
  return quadruple_3d_filter(a,b,c,d);
}

int main() {
  // --- Test 1: simple, non-degenerate, known sign ---
  // Choose orthonormal basis:
  // a=e1, b=e2, c=e1, d=e2
  // (a×b)·(c×d) = (e3)·(e3) = +1
  std::array<double,3> e1{1,0,0};
  std::array<double,3> e2{0,1,0};

  {
    Sign s = spip::predicates::quadruple3d(e1, e2, e1, e2);
    if(s != Sign::Positive) {
      std::cerr << "Test1 failed: expected Positive, got " << to_cstr(s) << "\n";
      return 1;
    }
  }

  // --- Test 2: sign flip by swapping c and d ---
  // (a×b)·(d×c) = (e3)·(-e3) = -1
  {
    Sign s = spip::predicates::quadruple3d(e1, e2, e2, e1);
    if(s != Sign::Negative) {
      std::cerr << "Test2 failed: expected Negative, got " << to_cstr(s) << "\n";
      return 1;
    }
  }

  // --- Test 3: try to trigger filter uncertainty (fallback) ---
  // Make (a·c)(b·d) - (a·d)(b·c) extremely small and positive.
  //
  // Let a = e1, b = e2.
  // Let c = e1 + t e2, d = t e1 + e2.
  //
  // Then:
  // a·c = 1
  // b·d = 1
  // a·d = t
  // b·c = t
  // Delta = 1*1 - t*t = 1 - t^2  (positive but very close to 0 if t≈1)
  //
  // We set t = 1 - eps with eps tiny so Delta ≈ 2*eps.
  //
  // This often makes the filter return FPG_UNCERTAIN_VALUE, forcing exact fallback.
  {
    const double eps = std::ldexp(1.0, -52);          // ~2.22e-16
    const double t   = 1.0 - eps;

    std::array<double,3> a{1,0,0};
    std::array<double,3> b{0,1,0};
    std::array<double,3> c{1,t,0};
    std::array<double,3> d{t,1,0};

    int f = filter_only(a.data(), b.data(), c.data(), d.data());
    Sign s = spip::predicates::quadruple3d(a,b,c,d);

    // True sign is Positive (Delta = 1 - t^2 > 0)
    if(s != Sign::Positive) {
      std::cerr << "Test3 failed: expected Positive, got " << to_cstr(s) << "\n";
      std::cerr << "Filter returned: " << f << " (uncertain=" << FPG_UNCERTAIN_VALUE << ")\n";
      return 1;
    }

    // Not a hard assertion, but useful info: did we likely hit fallback?
    if(f == FPG_UNCERTAIN_VALUE) {
      std::cout << "Test3: filter uncertain -> fallback exercised (good).\n";
    } else {
      std::cout << "Test3: filter decided sign directly (still ok). filter=" << f << "\n";
    }
  }

  std::cout << "All quadruple3d tests passed.\n";
  return 0;
}