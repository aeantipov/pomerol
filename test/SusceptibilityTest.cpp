//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2012 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2012 Igor Krivenko <igor@shg.ru>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.

#include <pomerol/Misc.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/StatesClassification.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/DensityMatrix.hpp>
#include <pomerol/MonomialOperator.hpp>
#include <pomerol/EnsembleAverage.hpp>
#include <pomerol/Susceptibility.hpp>

#include "catch2/catch-pomerol.hpp"

#include <cmath>

using namespace Pomerol;

TEST_CASE("Susceptibilities of a single Hubbard atom", "[Susceptibility]") {
    RealType U = 1.0;
    RealType mu = 0.4;
    RealType h_field = 0.01;
    RealType beta = 10.0;
    int n_iw = 20;

    auto omega = [beta](int n) { return M_PI * (2*n) / beta; };

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu) + Magnetization("A", -h_field);
    INFO("Hamiltonian\n" << HExpr);

    auto IndexInfo = MakeIndexClassification(HExpr);
    INFO("Indices\n" << IndexInfo);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);
    INFO("Energy levels " << H.getEigenValues());
    INFO("The value of ground energy is " << H.getGroundEnergy());

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    ParticleIndex up_index = IndexInfo.getIndex("A",0,up);
    ParticleIndex dn_index = IndexInfo.getIndex("A",0,down);

    // Quadratic operators of form c^+ c
    QuadraticOperator s_plus(IndexInfo, HS, S, H, up_index, dn_index);
    QuadraticOperator s_minus(IndexInfo, HS, S, H, dn_index, up_index);
    QuadraticOperator n_up(IndexInfo, HS, S, H, up_index, up_index);
    QuadraticOperator n_dn(IndexInfo, HS, S, H, dn_index, dn_index);

    for(auto * op : {&s_plus, &s_minus, &n_up, &n_dn}) {
        op->prepare(HS);
        op->compute();
    }

    // Reference statistical weights of states
    RealVectorType weights_ref(4); // {0, up, down, 2}
    weights_ref << 1.0,
                   exp(-beta*(-mu-h_field)),
                   exp(-beta*(-mu+h_field)),
                   exp(-beta*(-2*mu+U));
    weights_ref /= weights_ref.sum();
    RealType wu = weights_ref[1];
    RealType wd = weights_ref[2];
    RealType w2 = weights_ref[3];

    SECTION("Ensemble averages") {
        EnsembleAverage s_plus_aver(S, H, s_plus, rho);
        s_plus_aver.compute();
        REQUIRE_THAT(s_plus_aver.getResult(), IsCloseTo(0, 1e-14));

        EnsembleAverage s_minus_aver(S, H, s_minus, rho);
        s_minus_aver.compute();
        REQUIRE_THAT(s_minus_aver.getResult(), IsCloseTo(0, 1e-14));

        EnsembleAverage n_up_aver(S, H, n_up, rho);
        n_up_aver.compute();
        RealType n_up_ref = weights_ref[1] + weights_ref[3];
        REQUIRE_THAT(n_up_aver.getResult(), IsCloseTo(n_up_ref, 1e-14));

        EnsembleAverage n_dn_aver(S, H, n_dn, rho);
        n_dn_aver.compute();
        RealType n_dn_ref = weights_ref[2] + weights_ref[3];
        REQUIRE_THAT(n_dn_aver.getResult(), IsCloseTo(n_dn_ref, 1e-14));
    }

    SECTION("<S_+; S_->") {
        Susceptibility Chi(S, H, s_plus, s_minus, rho);
        Chi.prepare();
        Chi.compute();
        Chi.subtractDisconnected();

        auto ref = [&](int n) {
            ComplexType g = 0;
            if(std::abs(wu - wd) < 1e-8 && n == 0) // E_up == E_down
                g += wu * beta;
            else
                g += -(wu - wd) / (I*omega(n) - 2*h_field);
            return g;
        };
        for(int n = 0; n < n_iw; ++n)
            REQUIRE_THAT(Chi(n), IsCloseTo(ref(n), 1e-14));
    }

    SECTION("<n_up; n_up>") {
        Susceptibility Chi(S, H, n_up, n_up, rho);
        Chi.prepare();
        Chi.compute();
        Chi.subtractDisconnected();

        auto ref = [&](int n) {
            return n == 0 ? (wu + w2) * (1 - wu - w2) * beta : 0;
        };
        for(int n = 0; n < n_iw; ++n)
            REQUIRE_THAT(Chi(n), IsCloseTo(ref(n), 1e-14));
    }

    SECTION("<n_up; n_dn>") {
        Susceptibility Chi(S, H, n_up, n_dn, rho);
        Chi.prepare();
        Chi.compute();
        Chi.subtractDisconnected();

        auto ref = [&](int n) {
            return n == 0 ? (w2 - (wu + w2) * (wd + w2)) * beta : 0;
        };
        for(int n = 0; n < n_iw; ++n)
            REQUIRE_THAT(Chi(n), IsCloseTo(ref(n), 1e-14));
    }
}
