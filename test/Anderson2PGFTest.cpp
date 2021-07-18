#include <pomerol/Misc.hpp>
#include <pomerol/Operators.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/StatesClassification.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/DensityMatrix.hpp>
#include <pomerol/FieldOperatorContainer.hpp>
#include <pomerol/TwoParticleGFContainer.hpp>

#include "catch2/catch-pomerol.hpp"

#include <cmath>
#include <string>
#include <tuple>
#include <vector>

using namespace Pomerol;

TEST_CASE("Two-particle GF of the Anderson model", "[Anderson2PGF]") {
    RealType U = 0.5;
    RealType mu = 0.25;
    std::vector<RealType> levels = {1.02036910873357, -1.02036910873357};
    std::vector<RealType> hoppings = {0.296439333614347, 0.296439333614347};
    RealType beta = 26;

    RealType reduce_tol = 1e-5;
    RealType coeff_tol = 1e-8;

    RealVectorType chi_ref(10);
    chi_ref << -2.342841271771e+01,
                0.000000000000e+00,
                6.932231165814e-03,
                2.037522082872e-03,
               -2.150424835716e-03,
               -4.384848776411e-03,
               -5.253420668000e-03,
               -5.370700986029e-03,
               -5.126175681822e-03,
               -4.732777836189e-03;

    using namespace LatticePresets;

    auto HExpr = CoulombS("C", U, -mu);
    for(int i = 0; i < levels.size(); ++i)
    {
        auto bath_name = "b" + std::to_string(i);
        HExpr += Level(bath_name, levels[i]);
        HExpr += Hopping("C", bath_name, hoppings[i]);
    }
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

    DensityMatrix rho(S, H, beta);
    rho.prepare();
    rho.compute();

    ParticleIndex d0 = IndexInfo.getIndex("C", 0, down);
    ParticleIndex u0 = IndexInfo.getIndex("C", 0, up);

    std::set<ParticleIndex> f = {u0, d0};
    FieldOperatorContainer Operators(IndexInfo, HS, S, H, f);
    Operators.prepareAll(HS);
    Operators.computeAll();

    std::set<IndexCombination4> indices4 = {
        IndexCombination4(u0,u0,u0,u0),
        IndexCombination4(u0,d0,u0,d0),
        IndexCombination4(d0,d0,d0,d0)
    };
    TwoParticleGFContainer Chi4(IndexInfo,S,H,rho,Operators);
    Chi4.ReduceResonanceTolerance = reduce_tol;
    Chi4.CoefficientTolerance = coeff_tol;
    Chi4.MultiTermCoefficientTolerance = 1e-6;
    Chi4.prepareAll(indices4);
    MPI_Barrier(MPI_COMM_WORLD);

    ComplexType Omega = I*2.*M_PI/beta;
    ComplexType omega = I*M_PI/beta;

    std::vector<std::tuple<ComplexType, ComplexType, ComplexType>> freqs;

    SECTION("Chi4.computeAll() for arbitrary frequencies")
    {
        Chi4.computeAll(false, freqs, MPI_COMM_WORLD, true);

        const TwoParticleGF& chi_uuuu = Chi4(IndexCombination4(u0,u0,u0,u0));
        const TwoParticleGF& chi_dddd = Chi4(IndexCombination4(d0,d0,d0,d0));

        for(int i = 0; i < chi_ref.size(); ++i)
        {
            ComplexType w_p = I*(2.*i+1.)*M_PI/beta;
            INFO("i = " << i << ", w_p = " << w_p);
            auto ref = chi_ref[i];
            ComplexType chi_uuuu_val =  chi_uuuu(omega+Omega, w_p, omega);
            ComplexType chi_dddd_val =  chi_dddd(omega+Omega, w_p, omega);
            REQUIRE_THAT(chi_uuuu_val, IsCloseTo(ref, 1e-6));
            REQUIRE_THAT(chi_dddd_val, IsCloseTo(ref, 1e-6));
        }
    }

    SECTION("Chi4.computeAll() with precomputation for specific frequencies")
    {
        freqs.resize(chi_ref.size());
        for(int i = 0; i < chi_ref.size(); ++i) {
            ComplexType w_p = I*(2.*i+1.)*M_PI/beta;
            freqs[i] = std::make_tuple(omega+Omega, w_p, omega);
        }

        auto computed_data = Chi4.computeAll(true, freqs, MPI_COMM_WORLD, true);
        auto chi_uuuu = computed_data[IndexCombination4(u0,u0,u0,u0)];
        auto chi_dddd = computed_data[IndexCombination4(d0,d0,d0,d0)];

        for(int i = 0; i < chi_ref.size(); ++i)
        {
            ComplexType w_p = I*(2.*i+1.)*M_PI/beta;
            INFO("i = " << i << ", w_p = " << w_p);
            auto ref = chi_ref[i];
            ComplexType chi_uuuu_val =  chi_uuuu[i];
            ComplexType chi_dddd_val =  chi_dddd[i];
            REQUIRE_THAT(chi_uuuu_val, IsCloseTo(ref, 1e-6));
            REQUIRE_THAT(chi_dddd_val, IsCloseTo(ref, 1e-6));
        }
    }
}
