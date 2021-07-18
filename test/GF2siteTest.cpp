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

/** \file tests/green.cpp
** \brief Test of a Green's function calculation (1 s-orbital).
**
** \author Igor Krivenko (igor@shg.ru)
*/

#include <pomerol/Misc.hpp>
#include <pomerol/Operators.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/StatesClassification.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/DensityMatrix.hpp>
#include <pomerol/FieldOperatorContainer.hpp>
#include <pomerol/GreensFunction.hpp>

#include "catch2/catch-pomerol.hpp"

using namespace Pomerol;

TEST_CASE("Green's function of a Hubbard dimer", "[GF2site]") {
    RealType U = 1.0;
    RealType mu = 0.5;
    RealType beta = 10.0;

    // Reference Green's function
    ComplexVectorType G_ref(10);
    G_ref << -2.53021005e-01*I,
             -4.62090702e-01*I,
             -4.32482782e-01*I,
             -3.65598615e-01*I,
             -3.07785174e-01*I,
             -2.62894141e-01*I,
             -2.28274316e-01*I,
             -2.01170772e-01*I,
             -1.79539602e-01*I,
             -1.61950993e-01*I;

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu) + CoulombS("B", U, -mu);
    HExpr += Hopping("A","B", -1.0);
    INFO("Hamiltonian\n" << HExpr);

    auto IndexInfo = MakeIndexClassification(HExpr);
    INFO("Indices\n" << IndexInfo);

    auto A_dn = IndexInfo.getInfo(0);
    auto A_up = IndexInfo.getInfo(1);
    auto B_dn = IndexInfo.getInfo(2);
    auto B_up = IndexInfo.getInfo(3);

    auto N = Operators::N(std::vector<decltype(A_up)>{A_up, A_dn, B_up, B_dn});
    INFO("N = " << N);

    auto Sz = Operators::Sz(std::vector<decltype(A_up)>{A_up, B_up},
                            std::vector<decltype(A_dn)>{A_dn, B_dn});
    INFO("Sz = " << Sz);

    REQUIRE(Sz * N == N * Sz);
    REQUIRE(HExpr * N == N * HExpr);
    REQUIRE(HExpr * Sz == Sz * HExpr);

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
    for(QuantumState i = 0; i < S.getNumberOfStates(); ++i)
        INFO("Weight " << i << " = " << rho.getWeight(i));

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex A_down_index = IndexInfo.getIndex("A",0,down);

    auto const& c_map = Operators.getCreationOperator(A_down_index).getBlockMapping();
    for(auto c_map_it = c_map.right.begin(); c_map_it != c_map.right.end(); c_map_it++)
    {
        INFO(c_map_it->first << "->" << c_map_it->second);
    }

    GreensFunction GF(S,
                      H,
                      Operators.getAnnihilationOperator(A_down_index),
                      Operators.getCreationOperator(A_down_index),
                      rho);

    GF.prepare();
    GF.compute();

    for(int n = 0; n < G_ref.size(); ++n) {
        auto result = GF(n);
        auto ref = G_ref[n];
        // The tolerance has to be fairly large as Pomerol discards
        // some contributions to the GF.
        REQUIRE_THAT(result, IsCloseTo(ref, 1e-8));
    }
}
