// Solves the finite-H equation at fixed dt = 0.01 for each Lm in par::LM_LIST
// and writes psi_Lm<Lm>.dat with columns: t  psi  1-psi  dt
// Also writes the H = 0 reference solution to psi_H0.dat

#include "gruzberg.hpp"
#include <iostream>
#include <cmath>

int main() {
    const long double DT = 0.01L;   // fixed time step

    for (double Lm : par::LM_LIST) {
        gm::GMParams prm;
        prm.Lm = Lm;
        int NN = static_cast<int>(std::llround(Lm / DT));   // NN such that dt = Lm/NN = 0.01
        if (NN < 1) NN = 1;
        gm::Gruzberg s(prm, NN);
        s.itera_solve();

        std::string nome = "psi_Lm" + std::to_string(static_cast<int>(Lm)) + ".dat";
        std::ofstream out(nome); out << std::setprecision(19);
        out << "# t   psi   1-psi   dt\n";
        for (int k = 0; k < s.N; ++k) {
            // 1-psi is computed as deltapsi e^{b t} to avoid cancellation when psi is close to 1
            long double uno_meno_psi = (long double)s.deltapsi[k] * std::exp(prm.b * s.t[k]);
            out << s.t[k] << "  " << s.psi[k] << "  " << uno_meno_psi << "  " << s.dt << "\n";
        }
        out.close();
        std::cout << "wrote " << nome << "  (dt=" << (double)s.dt << ", N=" << s.N << ")\n";
    }

    // H = 0 reference solution
    gm::GruzbergH0 ref0(par::m, par::J, par::Lm0, par::LM0, par::NN0,
                        par::b, par::c, par::TOLL0, par::MAX_ITER, par::CHOP);
    ref0.itera_solve();
    gm::salva_psiH0(ref0, "psi_H0.dat");
    std::cout << "wrote psi_H0.dat\n";
    return 0;
}
