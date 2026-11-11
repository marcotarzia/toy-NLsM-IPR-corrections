#ifndef GRUZBERG_HPP
#define GRUZBERG_HPP

// Iterative solver for the Gruzberg-Mirlin integral equation.
//   - Gruzberg:   psi(t,H) at finite H = e^{-Lm}, for a single step dt = Lm/NN.
//                 For t <= 0 the substitution psi(t) = 1 - deltapsi(t) e^{b t} is used (m = 2 only)
//   - GruzbergH0: reference solution P(t,0) at H = 0, with the left tail fixed to 1 - c e^{b t}
// All parameters are set in parametri.hpp

#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include "parametri.hpp"

// Modified Bessel function I1: std::cyl_bessel_i if available (C++17 special math) otherwise Boost
#if defined(__cpp_lib_math_special_functions)
  #define GM_BESSEL_I1(x) std::cyl_bessel_i(1.0, (x))
#else
  #include <boost/math/special_functions/bessel.hpp>
  #define GM_BESSEL_I1(x) boost::math::cyl_bessel_i(1.0, (x))
#endif

namespace gm {

// Parameters of the finite-H solver (defaults from parametri.hpp)
struct GMParams {
    int    m       = par::m;
    double J       = par::J;
    double Lm      = par::LM_TEST;
    double LM      = par::LM;
    int    ordine  = par::ORDINE;
    double toll    = par::TOLL;
    int    maxiter = par::MAX_ITER;
    double chop    = par::CHOP;
    double b       = par::b;   // exponent of the substitution psi = 1 - deltapsi e^{b t} for t <= 0
};

// Number of grid points in [start, end] with step dt
inline int numero_punti(double start, double end, double dt) {
    return static_cast<int>((end - start) / dt + 1e-9) + 1;
}

// Gregory end-point correction weights for the given order (2, 4 or 6)
inline std::vector<double> pesi_gregory(int ordine) {
    switch (ordine) {
        case 2: return { 1.0/2.0 };
        case 4: return { 3.0/8.0, 7.0/6.0, 23.0/24.0 };
        case 6: return { 95.0/288.0, 317.0/240.0, 23.0/30.0, 793.0/720.0, 157.0/160.0 };
        default: return { 1.0/2.0 };
    }
}

// log I1(x) using the asymptotic expansion for x >= 700 to avoid overflow
inline double log_I1(double x) {
    constexpr double DUE_PI = 6.283185307179586476925287;
    if (x < 700.0) return std::log(GM_BESSEL_I1(x));
    double ix = 1.0/x;
    double S  = 1.0 - 0.375*ix - (15.0/128.0)*ix*ix - (105.0/1024.0)*ix*ix*ix;
    return x - 0.5*std::log(DUE_PI*x) + std::log(S);
}

// sqrt(2 pi x) e^{-x} I1(x) with its asymptotic form 1 - 3/(8x) for large x
inline double calcola_fac(double x) {
    constexpr double DUE_PI = 6.283185307179586476925287;
    if (std::log(x) < 15.0)
        return std::exp(log_I1(x) - x + 0.5*std::log(DUE_PI) + 0.5*std::log(x));
    return 1.0 - 0.375/x;
}

// 1 - e^{-x} computed as 2 e^{-x/2} sinh(x/2) to avoid cancellation for small x
// (for x >= 700 the result is 1 to double precision)
inline double uno_meno_exp_meno(double x) {
    if (x < 700.0) return 2.0*std::exp(-x*0.5)*std::sinh(x*0.5);
    return 1.0 - std::exp(-x);
}

// Solution psi(t,H) on the grid t in [-Lm, LM] with step dt = Lm/NN
class Gruzberg {
public:
    GMParams prm;
    int    NN, N, iterazioni = 0;
    long double dt, H, dmax = 1.0, pref, eHm, expJ;
    std::vector<double> t, psi, npsi, stt, q34, q14, em, psim, wfull, w;

    long double deltapsi0;                   // deltapsi(-Lm), exact boundary value
    std::vector<double> deltapsi, ndeltapsi; // deltapsi at current / next iteration
    std::vector<double> bj;                  // integrand factor replacing psi^2 e^{-e^t} for t <= 0
    std::vector<double> unmenoem;            // 1 - e^{-e^t}
    std::vector<double> bordoneg;            // analytic boundary term of deltapsi for t <= 0

    // Builds the grid, precomputes all t-dependent factors and sets the initial guess
    // psi_0(t) = e^{H/(m-1)} e^{-e^t}
    Gruzberg(const GMParams& p, int NN_) : prm(p), NN(NN_) {
        const double DUE_PI = 6.283185307179586476925287;
        if (NN < 1) NN = 1;
        // The expansion of psi^m in deltapsi is implemented for m = 2 only.
        if (prm.m != 2) {
            throw std::runtime_error("Gruzberg: the deltapsi substitution for t <= 0 is implemented only for m = 2");
        }
        dt = prm.Lm / NN;
        N  = numero_punti(-prm.Lm, prm.LM, dt);
        H    = std::exp(-prm.Lm);
        expJ = std::exp(prm.J);
        eHm  = std::exp(H / (prm.m - 1));
        pref = std::sqrt(prm.J) * expJ / std::sqrt(DUE_PI);

        // 1 - e^{H/(m-1)} = 2 e^{H/(2(m-1))} sinh(-H/(2(m-1))), written without cancellation.
        long double Hm1 = H / (prm.m - 1);
        long double due_sinh = 2.0L*std::exp(Hm1/2.0L)*std::sinh(-Hm1/2.0L);
        deltapsi0 = due_sinh * std::pow(H, -prm.b);

        t.assign(N,0); psi.assign(N,0); npsi.assign(N,0);
        stt.assign(N,0); q34.assign(N,0); q14.assign(N,0);
        em.assign(N,0); psim.assign(N,0); wfull.assign(N,1.0);
        deltapsi.assign(N,0); ndeltapsi.assign(N,0); bj.assign(N,0);
        unmenoem.assign(N,0); bordoneg.assign(N,0);

        const double H2 = H*H;
        for (int k=0;k<N;++k){
            t[k] = -prm.Lm + k*dt;
            double u = 1.0 - H2*std::exp(-2.0*t[k]);
            if (u<0.0) u=0.0;
            stt[k]=std::sqrt(u);
            q34[k]=(u>0.0)?std::pow(u,-0.75):0.0;
            q14[k]=std::pow(u,0.25);
            em[k]=std::exp(-std::exp(t[k]));
            unmenoem[k] = uno_meno_exp_meno(std::exp(t[k]));
            bordoneg[k] = std::exp(prm.J - prm.J*std::exp(t[k]+prm.Lm))
                          * std::exp(-prm.b*t[k]) * (double)due_sinh;
        }
        w = pesi_gregory(prm.ordine);
        int nw=(int)w.size();
        for (int k=1;k<nw && k<N;++k) wfull[k]=w[k];

        // Initial guess for psi
        for (int k=0;k<N;++k){
            double v = em[k];
            if (std::fabs(v) < prm.chop) v=0.0;
            npsi[k]=eHm*v; psi[k]=npsi[k];
        }
        // Initial guess for deltapsi, derived from psi
        // the boundary value is set exactly
        for (int k=0;k<N;++k){
            ndeltapsi[k] = (1.0 - psi[k]) * std::exp(-prm.b * t[k]);
            deltapsi[k]  = ndeltapsi[k];
        }
        deltapsi[0] = (double)deltapsi0; ndeltapsi[0] = (double)deltapsi0;
    }

    // Fixed-point iteration until max|psi_new - psi_old| < toll or maxiter is reached
    // For t <= 0 the equation is solved for deltapsi, for t > 0 directly for psi
    void itera_solve() {
        const double H2 = H*H;
        while (dmax > prm.toll && iterazioni < prm.maxiter) {
            ++iterazioni;

            // bj = 1 - psi^2 e^{-e^t} expressed through deltapsi (exact for m = 2)
            for (int j=1;j<N;++j){
                double tj = t[j];
                double dpj = deltapsi[j];
                bj[j] = unmenoem[j]
                        + 2.0*std::exp(-std::exp(tj)+prm.b*tj)*dpj
                        - std::exp(-std::exp(tj)+2.0*prm.b*tj)*dpj*dpj;
            }

            for (int k=0;k<N;++k) psim[k]=std::pow(psi[k],prm.m);
            npsi[0]=eHm;
            ndeltapsi[0]=(double)deltapsi0;

            // t <= 0 (indices 1..NN): update deltapsi, then psi = 1 - deltapsi e^{b t}
            for (int i=1;i<=NN;++i)
            {
                double ti=t[i]; double integ=0.0;
                double base_end = expJ*deltapsi0*(stt[i]*stt[i])
                                  *(prm.J*prm.J*std::exp(2.0*(ti+prm.Lm)))
                                  *std::exp(-prm.J*std::exp(ti+prm.Lm));
                integ += w[0]*0.5*base_end;
                for (int j=1;j<N;++j){
                    double tj=t[j];
                    double x = prm.J*std::exp(ti+tj)/H2 * stt[j]*stt[i];
                    double y = -prm.J/(1.0+stt[j]*stt[i])
                               *( std::exp(ti-tj)+std::exp(tj-ti)-H2*std::exp(-ti-tj) );
                    double fac = calcola_fac(x);
                    double g = pref*std::exp(-prm.b*ti)*bj[j]*q34[j]*q14[i]
                               *std::exp((ti-tj)*0.5+y)*fac;
                    integ += wfull[j]*g;
                }
                ndeltapsi[i] = bordoneg[i] + integ*dt;
                npsi[i] = 1.0 - std::exp(prm.b*ti)*ndeltapsi[i];
            }

            // t > 0 (indices NN+1..N-1): update psi directly
            // deltapsi is kept for the next bj
            for (int i=NN+1;i<N;++i){
                double ti=t[i]; double integ=0.0;
                double base_end = expJ*eHm*(stt[i]*stt[i])
                                  *(prm.J*prm.J*std::exp(2.0*(ti+prm.Lm)))
                                  *std::exp(-prm.J*std::exp(ti+prm.Lm));
                integ += w[0]*0.5*base_end;
                for (int j=1;j<N;++j){
                    double tj=t[j];
                    double x = prm.J*std::exp(ti+tj)/H2 * stt[j]*stt[i];
                    double y = -prm.J/(1.0+stt[j]*stt[i])
                               *( std::exp(ti-tj)+std::exp(tj-ti)-H2*std::exp(-ti-tj) );
                    double fac = calcola_fac(x);
                    double g = pref*psim[j]*em[j]*q34[j]*q14[i]
                               *std::exp((ti-tj)*0.5+y)*fac;
                    integ += wfull[j]*g;
                }
                npsi[i] = std::exp(prm.J - prm.J*std::exp(ti+prm.Lm))*eHm + integ*dt;
                ndeltapsi[i] = (1.0 - npsi[i]) * std::exp(-prm.b*ti);
            }

            dmax=0.0;
            for (int k=0;k<N;++k){ double d=std::fabs(psi[k]-npsi[k]); if(d>dmax)dmax=d; }
            for (int k=0;k<N;++k) psi[k]=npsi[k];
            for (int k=0;k<N;++k) deltapsi[k]=ndeltapsi[k];
        }
    }
};

// Optional: not used by main.cpp
// Extrapolation to dt -> 0 from three runs at dt, dt/2, dt/4
// Usage: gm::salva_psi_tre_dt(prm, NN, "file.dat")
// Extrapolates psi to dt = 0 by a linear least-squares fit of psi vs dt^6
// on the three values at dt, dt/2, dt/4
inline double estrapola_dt6(double y1, double y2, double y4, double dt) {
    double x[3] = { std::pow(dt,6), std::pow(dt/2.0,6), std::pow(dt/4.0,6) };
    long double y[3] = { y1, y2, y4 };
    double Sx=0,Sy=0,Sxx=0,Sxy=0;
    for (int k=0;k<3;++k){ Sx+=x[k]; Sy+=y[k]; Sxx+=x[k]*x[k]; Sxy+=x[k]*y[k]; }
    double b = (3*Sxy - Sx*Sy) / (3*Sxx - Sx*Sx);
    return (Sy - b*Sx) / 3.0;
}

// Solves at dt, dt/2, dt/4 and writes, on the common grid:
// t  psi_dt  psi_dt/2  psi_dt/4  psi_extrapolated  dt  dt/2  dt/4
inline void salva_psi_tre_dt(const GMParams& prm, int NN, const std::string& nome) {
    Gruzberg s1(prm, NN);     s1.itera_solve();
    Gruzberg s2(prm, 2*NN);   s2.itera_solve();
    Gruzberg s4(prm, 4*NN);   s4.itera_solve();
    std::ofstream out(nome); out << std::setprecision(15);
    out << "# t   psi_dt   psi_dt/2   psi_dt/4   psi_estrap   dt   dt/2   dt/4\n";
    for (int q=1; q<=s1.N; ++q) {
        int i1=q-1, i2=2*q-2, i4=4*q-4;
        if (i2>=s2.N || i4>=s4.N) break;
        long double pe = estrapola_dt6(s1.psi[i1], s2.psi[i2], s4.psi[i4], s1.dt);
        out << s1.t[i1] << "  "
            << s1.psi[i1] << "  " << s2.psi[i2] << "  " << s4.psi[i4] << "  "
            << pe << "  "
            << s1.dt << "  " << s2.dt << "  " << s4.dt << "\n";
    }
    out.close();
}
// End of optional part 

// Reference solution P(t,0) at H = 0 with fixed step dt = Lm0/NN0
// The left tail (first 10 units of t) is fixed analytically to 1 - c e^{b t}
class GruzbergH0 {
public:
    int    m, NN0, N, iterazioni = 0, maxiter;
    long double J, Lm0, LM0, b, c, toll0, chop, dt, dmax = 1.0, pref;
    std::vector<long double> t, psi, npsi, em, psim;

    // Builds the grid and sets the initial guess psi_0(t) = e^{-e^t}.
    GruzbergH0(int m_, double J_, double Lm0_, double LM0_, int NN0_,
               double b_, double c_, double toll0_, int maxiter_, double chop_)
        : m(m_), NN0(NN0_), maxiter(maxiter_), J(J_), Lm0(Lm0_), LM0(LM0_),
          b(b_), c(c_), toll0(toll0_), chop(chop_) {
        const double DUE_PI = 6.283185307179586476925287;
        dt = Lm0 / NN0;
        N  = numero_punti(-Lm0, LM0, dt);
        pref = std::sqrt(J)*std::exp(J)/std::sqrt(DUE_PI);
        t.assign(N,0); psi.assign(N,0); npsi.assign(N,0); em.assign(N,0); psim.assign(N,0);
        for (int k=0;k<N;++k){ t[k]=-Lm0+k*dt; em[k]=std::exp(-std::exp(t[k])); }
        for (int k=0;k<N;++k){ double v=em[k]; if(std::fabs(v)<chop)v=0.0; npsi[k]=v; psi[k]=v; }
    }

    // Fixed-point iteration until max|psi_new - psi_old| < toll0 or maxiter is reached
    void itera_solve() {
        int n_coda = static_cast<int>(10.0/dt + 1e-9) + 1;
        while (dmax > toll0 && iterazioni < maxiter) {
            ++iterazioni;
            for (int k=0;k<N;++k) psim[k]=std::pow(psi[k],m);
            npsi[0]=1.0;
            for (int i=1;i<N;++i){
                double ti=t[i]; double integ=0.0;
                for (int j=1;j<N;++j){
                    double tj=t[j];
                    double y = -J*(std::exp(ti-tj)+std::exp(tj-ti))*0.5;   // -J cosh(ti-tj)
                    integ += pref*psim[j]*em[j]*std::exp((ti-tj)*0.5+y);
                }
                npsi[i]=integ*dt;
            }
            for (int k=0;k<n_coda && k<N;++k) npsi[k]=1.0 - c*std::exp(b*t[k]);
            dmax=0.0;
            for (int k=0;k<N;++k){ double d=std::fabs(psi[k]-npsi[k]); if(d>dmax)dmax=d; }
            for (int k=0;k<N;++k) psi[k]=npsi[k];
        }
    }
};

// Writes the H = 0 solution with columns: t  P(t,0)  1-P(t,0)
inline void salva_psiH0(const GruzbergH0& ref, const std::string& nome) {
    std::ofstream out(nome); out << std::setprecision(15);
    out << "# t   P(t,0)   1-P(t,0)\n";
    for (int k=0;k<ref.N;++k)
        out << ref.t[k] << "  " << ref.psi[k] << "  " << (1.0-ref.psi[k]) << "\n";
    out.close();
}

}  // namespace gm
#endif
