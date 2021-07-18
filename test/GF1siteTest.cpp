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
#include <pomerol/LatticePresets.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/StatesClassification.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/DensityMatrix.hpp>
#include <pomerol/FieldOperatorContainer.hpp>
#include <pomerol/GreensFunction.hpp>
#include <pomerol/GFContainer.hpp>

#include "catch2/catch-pomerol.hpp"

using namespace Pomerol;

TEST_CASE("Green's function of a Hubbard atom", "[GF1site]") {
    RealType U = 1.0;
    RealType mu = 0.4;
    RealType beta = 10.0;

    // Reference Green's function
    auto G_ref = [U,mu,beta](int n)
    {
        RealType omega = M_PI*(2*n+1)/beta;

        RealType w0 = 1.0;
        RealType w1 = exp(beta*mu);
        RealType w2 = exp(-beta*(-2*mu+U));
        RealType Z = w0 + 2*w1 +w2;
        w0 /= Z; w1 /= Z; w2 /= Z;

        return (w0+w1)/(I*omega+mu) + (w1+w2)/(I*omega+mu-U);
    };

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", U, -mu);
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
    for(QuantumState i = 0; i < S.getNumberOfStates(); ++i)
        INFO("Weight " << i << " = " << rho.getWeight(i));

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex down_index = IndexInfo.getIndex("A",0,down);

    auto const& c_map = Operators.getCreationOperator(down_index).getBlockMapping();
    for (auto c_map_it = c_map.right.begin(); c_map_it != c_map.right.end(); c_map_it++)
    {
        INFO(c_map_it->first << "->" << c_map_it->second);
    }

    GreensFunction GF(S,
                      H,
                      Operators.getAnnihilationOperator(down_index),
                      Operators.getCreationOperator(down_index),
                      rho);

    GF.prepare();
    GF.compute();

    for(int n = 0; n < 100; ++n) {
        auto result = GF(n);
        auto ref = G_ref(n);
        REQUIRE_THAT(result, IsCloseTo(ref, 1e-14));
    }

    SECTION("GFContainer") {
        GFContainer G(IndexInfo,S,H,rho,Operators);

        std::set<IndexCombination2> indices;
        indices.insert(IndexCombination2(0,0));
        indices.insert(IndexCombination2(0,1));
        indices.insert(IndexCombination2(1,0));
        indices.insert(IndexCombination2(1,1));

        G.prepareAll(indices);
        G.computeAll();

        for(int n = -100; n < 100; ++n) {
            auto ref = G_ref(n);
            REQUIRE_THAT(G(0,0)(n), IsCloseTo(ref, 1e-14));
            REQUIRE_THAT(G(0,1)(n), IsCloseTo(0, 1e-14));
            REQUIRE_THAT(G(1,0)(n), IsCloseTo(0, 1e-14));
            REQUIRE_THAT(G(1,1)(n), IsCloseTo(ref, 1e-14));
        }
    }
}
