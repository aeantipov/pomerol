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

/** \file tests/Susceptibility1siteTest.cpp
** \brief Test of a calculation of dynamical susceptibility (1 s-orbital).
**
** \author Igor Krivenko (igor@shg.ru)
** \author Junya Otsuki (j.otsuki@okayama-u.ac.jp)
*/

#include "Misc.h"
#include "LatticePresets.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Operators.h"
#include "StatesClassification.h"
#include "Hamiltonian.h"
#include "FieldOperatorContainer.h"
#include "Susceptibility.h"

#include "./Utility.h"

#include <cmath>
#include <cstdlib>

using namespace Pomerol;

RealType U = 1.0;
RealType mu = 0.4;
RealType h_field = 0.01;
//RealType h_field = 0;

bool compare(ComplexType a, ComplexType b)
{
    return std::abs(a-b) < 1e-14;
}

// Reference Green's function

struct Weights{
    RealType omega;
    RealType w0, wu, wd, w2;

    Weights(int n, RealType beta){
        omega = M_PI*(2*n)/beta;
        w0 = 1.0;
        wu = exp(beta*(mu+h_field));
        wd = exp(beta*(mu-h_field));
        w2 = exp(-beta*(-2*mu+U));
        RealType Z = w0 + wu + wd + w2;
        w0 /= Z; w2 /= Z; wu /= Z; wd /= Z;
    }
};

// <S_+; S_->
ComplexType Gref_pm(int n, RealType beta)
{
    Weights W(n, beta);

    ComplexType g = 0;
    if(std::abs(W.wu - W.wd) < 1e-8 ){  // E_up == E_down
        if(n==0)  g += W.wu * beta;
    }
    else{
        g += -(W.wu - W.wd) / (I*W.omega - 2*h_field);
    }
    return g;
}

// <n_u; n_u>
ComplexType Gref_uu(int n, RealType beta)
{
    Weights W(n, beta);

    ComplexType g = 0;
//    if(n==0)  g += (W.wu + W.w2) * beta;
    if(n==0)  g += (W.wu + W.w2) * (1 - W.wu - W.w2) * beta;
    return g;
}

// <n_u; n_d>
ComplexType Gref_ud(int n, RealType beta)
{
    Weights W(n, beta);

    ComplexType g = 0;
//    if(n==0)  g += W.w2 * beta;
    if(n==0)  g += (W.w2 - (W.wu + W.w2) * (W.wd + W.w2)) * beta;
    return g;
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    using namespace LatticePresets;

    // h_field (n_down - n_up)
    // addHopping double-counts diagonal term, and divide h_field by 2.
    auto HExpr = CoulombS("A", U, -mu) + Magnetization("A", -h_field);

    INFO("Hamiltonian");
    INFO(HExpr);
    auto IndexInfo = MakeIndexClassification(HExpr);
    print_section("Indices");
    std::cout << IndexInfo << std::endl;

    auto HS = MakeHilbertSpace(IndexInfo, HExpr);
    HS.compute();
    StatesClassification S;
    S.compute(HS);

    Hamiltonian H(S);
    H.prepare(HExpr, HS, MPI_COMM_WORLD);
    H.compute(MPI_COMM_WORLD);

    RealType beta = 10.0;

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    ParticleIndex dn_index = IndexInfo.getIndex("A",0,down);
    ParticleIndex up_index = IndexInfo.getIndex("A",0,up);

    // quadratic operators, c^+ c
    QuadraticOperator s_plus(IndexInfo, HS, S, H, up_index, dn_index);
    QuadraticOperator s_minus(IndexInfo, HS, S, H, dn_index, up_index);
    QuadraticOperator n_up(IndexInfo, HS, S, H, up_index, up_index);
    QuadraticOperator n_dn(IndexInfo, HS, S, H, dn_index, dn_index);

    std::vector<QuadraticOperator*> quad_ops{&s_plus, &s_minus, &n_up, &n_dn};
    for(auto op : quad_ops){
        op->prepare(HS);
        op->compute();
    }

    // compute 3 susceptibilities
    //  < S_+ ; S_- >
    //  < n_up ; n_up >
    //  < n_up ; n_dn >
    std::vector< std::pair<QuadraticOperator*, QuadraticOperator*> > op_pairs;
    op_pairs.emplace_back(std::make_pair(&s_plus, &s_minus));
    op_pairs.emplace_back(std::make_pair(&n_up, &n_up));
    op_pairs.emplace_back(std::make_pair(&n_up, &n_dn));

    // for print
    std::vector<std::string> names = {"< S_+ ; S_- >", "< n_up ; n_up >", "< n_up ; n_dn >"};

    // reference data
    ComplexType (*Grefs[3])(int n, RealType beta) = {Gref_pm, Gref_uu, Gref_ud};

    // compute susceptibilities, and compare them with reference data
    bool result = true;
    for(int i=0; i<op_pairs.size(); i++){
        print_section(names[i]);
        Susceptibility Chi(S,H, *(op_pairs[i].first), *(op_pairs[i].second), rho);
        Chi.prepare();
        Chi.compute();
        Chi.subtractDisconnected();

        // check if results are correct
        for(int n = 0; n<20; ++n) {
            INFO(Chi(n) << " == " << Grefs[i](n, beta));
            result = (result && compare(Chi(n), Grefs[i](n,beta)));
        }
    }

    MPI_Finalize();
    return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
