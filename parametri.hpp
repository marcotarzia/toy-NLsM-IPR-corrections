#ifndef PARAMETRI_HPP
#define PARAMETRI_HPP

// Simulation parameters

#include <vector>

// Finite-H equation
#ifndef PAR_M
#define PAR_M       2                     // power m in the equation (only m = 2 is supported)
#endif
#ifndef PAR_J
#define PAR_J       0.026540159736761493  // J at the critical point
#endif
#ifndef PAR_LM_MAX
#define PAR_LM_MAX  10.0                  // right end of the t grid
#endif
#ifndef PAR_ORDINE
#define PAR_ORDINE  6                     // order of the Gregory quadrature (2, 4 or 6)
#endif
#ifndef PAR_TOLL
#define PAR_TOLL    1e-13                 // convergence threshold on max|psi_new - psi_old|
#endif
#ifndef PAR_MAXITER
#define PAR_MAXITER 200000                // maximum number of iterations
#endif
#ifndef PAR_CHOP
#define PAR_CHOP    1e-10                 // values below this are set to 0 in the initial guess
#endif
#ifndef PAR_LM_TEST
#define PAR_LM_TEST 14.0                  // default Lm (overridden in main.cpp)
#endif
#ifndef PAR_NN_BASE
#define PAR_NN_BASE 175                   // default NN for salva_psi_tre_dt (dt = 14/175 = 0.08). Optional, not used by main.cpp
#endif

// H = 0 reference solution (left tail fixed to 1 - c e^{b t}; Lm0 must be >= every Lm in LM_LIST)
#ifndef PAR_LM0
#define PAR_LM0     100.0                 // |t_min| for H = 0
#endif
#ifndef PAR_LM0_MAX
#define PAR_LM0_MAX 10.0                  // t_max for H = 0
#endif
#ifndef PAR_NN0
#define PAR_NN0     625                   // number of intervals for H = 0 (dt0 = 100/625 = 0.16)
#endif
#ifndef PAR_B
#define PAR_B       0.5                   // exponent b (1/2 at the critical point)
#endif
#ifndef PAR_C
#define PAR_C       0.61                  // tail coefficient c
#endif
#ifndef PAR_TOLL0
#define PAR_TOLL0   1e-7                  // convergence threshold for H = 0
#endif

namespace par {

    // Finite-H equation
    constexpr int    m        = PAR_M;
    constexpr double J        = PAR_J;
    constexpr double LM       = PAR_LM_MAX;
    constexpr int    ORDINE   = PAR_ORDINE;
    constexpr double TOLL     = PAR_TOLL;
    constexpr int    MAX_ITER = PAR_MAXITER;
    constexpr double CHOP     = PAR_CHOP;
    constexpr double LM_TEST  = PAR_LM_TEST;
    constexpr int    NN_BASE  = PAR_NN_BASE;    // optional, not used by main.cpp

    // Values of Lm (H = e^{-Lm}) to solve for
    inline const std::vector<double> LM_LIST = { 50, 54 };

    // H = 0 reference solution
    constexpr double Lm0   = PAR_LM0;
    constexpr double LM0   = PAR_LM0_MAX;
    constexpr int    NN0   = PAR_NN0;
    constexpr double b     = PAR_B;
    constexpr double c     = PAR_C;
    constexpr double TOLL0 = PAR_TOLL0;

    constexpr const char* PREFISSO_FILE = "psi";   // output file prefix. Optional, not used by main.cpp

}

#endif  // PARAMETRI_HPP
