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

#include "Misc.hpp"
#include "LatticePresets.hpp"
#include "Operators.hpp"
#include "IndexClassification.hpp"
#include "HilbertSpace.hpp"
#include "StatesClassification.hpp"
#include "HamiltonianPart.hpp"
#include "Hamiltonian.hpp"
#include "DensityMatrix.hpp"
#include "FieldOperatorContainer.hpp"
#include "GFContainer.hpp"

#include "./Utility.hpp"

#include <cstdlib>

using namespace Pomerol;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    using namespace LatticePresets;

    auto HExpr = CoulombS("A", 1.0, -0.5);
    HExpr += CoulombS("B", 2.0, -1.1);
    HExpr += CoulombS("C", 3.0, -0.7);
    HExpr += CoulombS("D", 4.0, -1.1);

    HExpr += Hopping("A", "B", -1.3);
    HExpr += Hopping("B", "C", -0.45);
    HExpr += Hopping("C", "D", -0.127);
    HExpr += Hopping("A", "D", -0.255);

    auto IndexInfo = MakeIndexClassification(HExpr);
    print_section("Indices");
    std::cout << IndexInfo << std::endl;

    INFO("Hamiltonian");
    INFO(HExpr);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();

    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);
    if(pMPI::rank(MPI_COMM_WORLD) == 0) {
        INFO("Energy levels " << H.getEigenValues());
        INFO("The value of ground energy is " << H.getGroundEnergy());
    }

    RealType beta = 10.0;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();
    for (QuantumState i=0; i<S.getNumberOfStates(); ++i) INFO(rho.getWeight(i));


    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    auto c_map = Operators.getCreationOperator(0).getBlockMapping();
    for (auto c_map_it = c_map.right.begin(); c_map_it != c_map.right.end(); c_map_it++)
    {
        INFO(c_map_it->first << "->" << c_map_it->second);
    }

    ParticleIndex down_index = IndexInfo.getIndex("A", 0, down);
    GreensFunction GF(S,
                      H,
                      Operators.getAnnihilationOperator(down_index),
                      Operators.getCreationOperator(down_index),
                      rho);

    GF.prepare();
    GF.compute();

    ComplexVectorType GF_ref(10);

    // Exact reference results
    GF_ref <<  0.00515461461  -0.191132319*I,
              -0.0129218293   -0.35749415*I,
              -0.0063208255   -0.364571553*I,
              -0.00244599255  -0.326995909*I,
              -0.000938220077 -0.285235829*I,
              -0.000360621591 -0.248974505*I,
              -0.000129046261 -0.219206946*I,
              -3.20102701e-05 -0.194983212*I,
               9.51503858e-06 -0.175149329*I,
               2.68929175e-05 -0.158732731*I;

    bool result = true;
    auto compare = [](ComplexType a, ComplexType b) {
        // The tolerance has to be fairly large as Pomerol discards
        // some contributions to the GF.
        return std::abs(a - b) < 1e-6;
    };

    for(int n = 0; n<10; ++n) {
        DEBUG(GF(n) << " " << GF_ref(n));
        result = result && compare(GF(n), GF_ref(n));
    }

    MPI_Finalize();

    return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
