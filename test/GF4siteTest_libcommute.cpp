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

#include "Misc.h"
//#include "Lattice.h"
//#include "LatticePresets.h"
#include "Operators.h"
#include "IndexClassification.h"
//#include "OperatorPresets.h"
#include "HilbertSpace.h"
#include "StatesClassification.h"
#include "HamiltonianPart.h"
#include "Hamiltonian.h"
#include "DensityMatrix.h"
#include "FieldOperatorContainer.h"
#include "GFContainer.h"

#include<cstdlib>

using namespace Pomerol;

bool compare(ComplexType a, ComplexType b)
{
    return abs(a-b) < 1e-5;
}

void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
  std::cout << std::string(str.size(),'=') << std::endl;
}

int main(int argc, char* argv[])
{

    MPI_Init(&argc, &argv);

    // TODO: Remove Lattice (?) and rework LatticePresets
    /*
    Lattice L;
    L.addSite(new Lattice::Site("A",1,2));
    LatticePresets::addCoulombS(&L, "A", 1.0, -0.5);
    L.addSite(new Lattice::Site("B",1,2));
    LatticePresets::addCoulombS(&L, "B", 2.0, -1.1);
    L.addSite(new Lattice::Site("C",1,2));
    LatticePresets::addCoulombS(&L, "C", 3.0, -0.7);
    L.addSite(new Lattice::Site("D",1,2));
    LatticePresets::addCoulombS(&L, "D", 4.0, -1.1);

    LatticePresets::addHopping(&L, "A","B", -1.3);
    LatticePresets::addHopping(&L, "B","C", -0.45);
    LatticePresets::addHopping(&L, "C","D", -0.127);
    LatticePresets::addHopping(&L, "A","D", -0.255);
    INFO("Sites");
    L.printSites();
    INFO("Terms with 2 operators");
    L.printTerms(2);
    INFO("Terms with 4 operators");
    L.printTerms(4);
    */

    using namespace Operators;

    auto HExpr = -0.5 * (n("A", up) + n("A", down)) + 1.0 * n("A", up) * n("A", down);
    HExpr += -1.1 * (n("B", up) + n("B", down)) + 2.0 * n("B", up) * n("B", down);
    HExpr += -0.7 * (n("C", up) + n("C", down)) + 3.0 * n("C", up) * n("C", down);
    HExpr += -1.1 * (n("D", up) + n("D", down)) + 4.0 * n("D", up) * n("D", down);
    for(spin s : {up, down}) {
        HExpr += -1.3 * c_dag("A", s) * c("B", s) + hc;
        HExpr += -0.45 * c_dag("B", s) * c("C", s) + hc;
        HExpr += -0.127 * c_dag("C", s) * c("D", s) + hc;
        HExpr += -0.255 * c_dag("A", s) * c("D", s) + hc;
    }

    auto IndexInfo = MakeIndexClassification(HExpr);
    print_section("Indices");
    IndexInfo.printIndices();

    INFO("Hamiltonian");
    INFO(HExpr);

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();

    StatesClassification S;
    S.compute(HS);

    // FIXME: Remove
    for(int i = 0; i < S.getNumberOfStates(); i++) {
      std::cout << i << ": " << S.getBlockNumber(i) << "/" << S.getInnerState(i) << std::endl;
    }
    for(int b = 0; b < S.getNumberOfBlocks(); b++) {
      std::cout << "size(block("  << b << ")) = " << S.getBlockSize(b);
      std::cout << "; states = ";
      for(auto const& s : S.getFockStates(b)) std::cout << s << ", ";
      std::cout << std::endl;
    }

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

    ParticleIndex down_index = IndexInfo.getIndex("A", down);
    GreensFunction GF(S,
                      H,
                      Operators.getAnnihilationOperator(down_index),
                      Operators.getCreationOperator(down_index),
                      rho);

    GF.prepare();
    GF.compute();

    ComplexVectorType GF_ref(10);

    GF_ref << 0.00515461461  -0.191132319*I,
              -0.0129218293  -0.35749415*I,
              -0.0063208255  -0.364571553*I,
              -0.00244599255  -0.326995909*I,
              -0.000938220077  -0.285235829*I,
              -0.000360621591  -0.248974505*I,
              -0.000129046261  -0.219206946*I,
              -3.20102701e-05  -0.194983212*I,
              9.51503858e-06  -0.175149329*I,
              2.68929175e-05  -0.158732731*I;

    bool result = true;
    for(int n = 0; n<10; ++n) {
        DEBUG(GF(n) << " " << GF_ref(n));
        result = (result && compare(GF(n),GF_ref(n)));
    }

    MPI_Finalize();

    return result ? EXIT_SUCCESS : EXIT_FAILURE;
}