#include <pomerol/Misc.hpp>
#include <pomerol/LatticePresets.hpp>
#include <pomerol/IndexClassification.hpp>
#include <pomerol/HilbertSpace.hpp>
#include <pomerol/StatesClassification.hpp>
#include <pomerol/Hamiltonian.hpp>
#include <pomerol/DensityMatrix.hpp>
#include <pomerol/FieldOperatorContainer.hpp>
#include <pomerol/GreensFunction.hpp>

#include "catch2/catch-pomerol.hpp"

#include <cmath>
#include <vector>

using namespace Pomerol;

// Parameters
double mu = 1.2;
double U = 2.0;
double beta = 10;

// Reference Green's function
struct GF_ref {
    RealType phi;
    std::vector<RealType> E;
    std::vector<RealType> w;
    RealType Z = 0;
    GF_ref(ComplexType J) :
        phi(std::arg(J)),
        E{0, -mu-std::abs(J), -mu+std::abs(J), -2*mu+U} {
        for(int i = 0; i < E.size(); ++i) {
            w.push_back(std::exp(-beta*E[i]));
            Z += w.back();
        }
        for(int i = 0; i < E.size(); ++i) w[i] /= Z;
    }

    ComplexType operator()(spin s1, spin s2, int n) const {
        ComplexType iw = ComplexType(0,M_PI*(2*n+1)/beta);

        ComplexType res = 0;
        if(s1 == up && s2 == up) {
            res += 0.5*(w[1] + w[0])/(iw - (E[1] - E[0]));
            res += 0.5*(w[2] + w[0])/(iw - (E[2] - E[0]));
            res += 0.5*(w[3] + w[1])/(iw - (E[3] - E[1]));
            res += 0.5*(w[3] + w[2])/(iw - (E[3] - E[2]));
        } else if(s1 == up && s2 == down) {
            ComplexType u = std::exp(ComplexType(0,phi));
            res += -u*0.5*(w[1] + w[0])/(iw - (E[1] - E[0]));
            res += u*0.5*(w[2] + w[0])/(iw - (E[2] - E[0]));
            res += u*0.5*(w[3] + w[1])/(iw - (E[3] - E[1]));
            res += -u*0.5*(w[3] + w[2])/(iw - (E[3] - E[2]));
        } else if(s1 == down && s2 == up) {
            ComplexType u = std::exp(ComplexType(0,-phi));
            res += -u*0.5*(w[1] + w[0])/(iw - (E[1] - E[0]));
            res += u*0.5*(w[2] + w[0])/(iw - (E[2] - E[0]));
            res += u*0.5*(w[3] + w[1])/(iw - (E[3] - E[1]));
            res += -u*0.5*(w[3] + w[2])/(iw - (E[3] - E[2]));
        } else {
            res += 0.5*(w[1] + w[0])/(iw - (E[1] - E[0]));
            res += 0.5*(w[2] + w[0])/(iw - (E[2] - E[0]));
            res += 0.5*(w[3] + w[1])/(iw - (E[3] - E[1]));
            res += 0.5*(w[3] + w[2])/(iw - (E[3] - E[2]));
        }
        return res;
    }
};

TEST_CASE("Bare Anderson atom with complex spin mixing", "[AndersonComplex]") {

    // Execute this test case for a few values of 'J'
    auto J = GENERATE(ComplexType(0.1, 0),
                      ComplexType(-0.1, 0),
                      ComplexType(0, 0.1),
                      ComplexType(0,-0.1),
                      ComplexType(0.1, 0.1),
                      ComplexType(0.1,-0.1),
                      ComplexType(-0.1,0.1),
                      ComplexType(-0.1,-0.1)
                      );

    using namespace LatticePresets;

    ComplexExpr HExpr = CoulombS("C", U, -mu);
    HExpr += Hopping("C", "C", J, 0, 0, up, down);
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

    FieldOperatorContainer Operators(IndexInfo, HS, S, H);
    Operators.prepareAll(HS);
    Operators.computeAll();

    // Reference
    GF_ref G_ref(J);

    for(spin s1 : {down, up}){
        for(spin s2 : {down, up}){
            ParticleIndex index1 = IndexInfo.getIndex("C",0,s1);
            ParticleIndex index2 = IndexInfo.getIndex("C",0,s2);
            GreensFunction GF(S,H,
                Operators.getAnnihilationOperator(index1),
                Operators.getCreationOperator(index2),
            rho);
            GF.prepare();
            GF.compute();

            for(int n = 0; n < 10; ++n) {
                auto result = GF(n);
                auto ref = G_ref(s1, s2, n);
                REQUIRE_THAT(result, IsCloseTo(ref, 1e-12));
            }
        }
    }
}
