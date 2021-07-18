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

/** \file tests/IndexTerm.cpp
** \brief Test of the Symmetrizer::IndexPermutation.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include <pomerol/LatticePresets.hpp>
#include <pomerol/Index.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/Operators.hpp>
#include <pomerol/StatesClassification.hpp>
#include <pomerol/HamiltonianPart.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/MonomialOperator.hpp>

#include "catch2/catch-pomerol.hpp"

using namespace Pomerol;

TEST_CASE("Simple Hamiltonian test", "[hamiltonian]") {
    using namespace LatticePresets;

    auto HExpr = CoulombS("A", 1.0, -0.5);
    HExpr += CoulombS("B", 2.0, -1.0);
    HExpr += Hopping("A", "B", -1.0);

    auto IndexInfo = MakeIndexClassification(HExpr);
    INFO(IndexInfo);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);

    INFO(H.getPart(BlockNumber(4)));
    H.compute(MPI_COMM_WORLD);
    INFO(H.getPart(BlockNumber(4)));

    SECTION("Ground state energy") {
        RealType E_ref = -2.8860009;
        RealType E = H.getGroundEnergy();

        REQUIRE_THAT(E, IsCloseTo(E_ref, 1e-7));
    }

    SECTION("Monomial operators") {
        ParticleIndex op_index = IndexInfo.getIndex("B", 0, up);
        BlockNumber test_block = 4;

        CreationOperator op(IndexInfo, HS, S, H, op_index);
        op.prepare(HS);
        op.compute();
        BlockNumber result_block = op.getLeftIndex(test_block);

        INFO("Acting with rotated cdag_" << op_index << " on block " << test_block << " and receiving " << result_block);

        using LOperatorType = libcommute::loperator<RealType, libcommute::fermion>;

        HamiltonianPart HpartRHS(LOperatorType(HExpr, HS.getFullHilbertSpace()), S, test_block);
        HpartRHS.prepare();
        HpartRHS.compute();
        INFO(HpartRHS);

        HamiltonianPart HpartLHS(LOperatorType(HExpr, HS.getFullHilbertSpace()), S, result_block);
        HpartLHS.prepare();
        HpartLHS.compute();
        INFO(HpartLHS);

        auto cdag1op = LOperatorType(Operators::c_dag("B", (unsigned short)0, up), HS.getFullHilbertSpace());
        MonomialOperatorPart Cdag1(cdag1op, S, HpartRHS, HpartLHS);
        Cdag1.compute();
        INFO(Cdag1);

        auto c1op = LOperatorType(Operators::c("B", (unsigned short)0, up), HS.getFullHilbertSpace());
        MonomialOperatorPart C1(c1op, S, HpartLHS, HpartRHS);
        C1.compute();
        INFO(C1);

        // Check transposition
        auto diff1 = (Cdag1.getRowMajorValue<false>() - C1.getColMajorValue<false>().transpose()).eval();
        diff1.prune(1e-12);
        REQUIRE(diff1.nonZeros() == 0);

        auto diff2 = (Cdag1.getColMajorValue<false>() - C1.getRowMajorValue<false>().transpose()).eval();
        diff2.prune(1e-12);
        REQUIRE(diff2.nonZeros() == 0);
    }
}
