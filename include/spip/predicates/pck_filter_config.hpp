#pragma once
#include <cmath>   // for std::fabs if you want, but generated code calls fabs()

#include <cstdlib> // optional

// Use an out-of-band sentinel; do NOT use 0.
#ifndef SPIP_FPG_UNCERTAIN_VALUE
#define SPIP_FPG_UNCERTAIN_VALUE 2147483647
#endif

#ifndef FPG_UNCERTAIN_VALUE
#define FPG_UNCERTAIN_VALUE SPIP_FPG_UNCERTAIN_VALUE
#endif