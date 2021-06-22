//
// This file is a part of pomerol - a scientific ED code for obtaining
// properties of a Hubbard model on a finite-size lattice
//
// Copyright (C) 2010-2012 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2012 Igor Krivenko <Igor.S.Krivenko@gmail.com>
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
#ifndef POMEROL_INCLUDE_TWOPARTICLEGFCONTAINER_H
#define POMEROL_INCLUDE_TWOPARTICLEGFCONTAINER_H

#include"Misc.hpp"
#include"Index.hpp"
#include"TwoParticleGF.hpp"
#include"FieldOperatorContainer.hpp"
#include"IndexContainer4.hpp"

#include <mpi.h>

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace Pomerol {

class TwoParticleGFContainer: public IndexContainer4<TwoParticleGF,TwoParticleGFContainer>, public Thermal
{
public:
    /** A difference in energies with magnitude less than this value is treated as zero. default = 1e-8. */
    RealType ReduceResonanceTolerance = 1e-8;
    /** Minimal magnitude of the coefficient of a term to take it into account. default = 1e-16. */
    RealType CoefficientTolerance = 1e-16;
    /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. default = 1e-5. */
    RealType MultiTermCoefficientTolerance = 1e-5;

    template<typename... IndexTypes>
    TwoParticleGFContainer(const IndexClassification<IndexTypes...>& IndexInfo,
                           const StatesClassification &S,
                           const Hamiltonian &H,
                           const DensityMatrix &DM,
                           const FieldOperatorContainer& Operators) :
        IndexContainer4<TwoParticleGF, TwoParticleGFContainer>(this, IndexInfo), Thermal(DM),
        S(S), H(H), DM(DM), Operators(Operators)
    {}

    void prepareAll(const std::set<IndexCombination4>& InitialIndices = {});
    std::map<IndexCombination4,std::vector<ComplexType>> computeAll(
        bool clearTerms = false,
        FreqVec const& freqs = {},
        const MPI_Comm& comm = MPI_COMM_WORLD,
        bool split = true
    );
    std::map<IndexCombination4,std::vector<ComplexType>> computeAll_nosplit(
        bool clearTerms,
        FreqVec const& freqs = {},
        const MPI_Comm& comm = MPI_COMM_WORLD
    );
    std::map<IndexCombination4,std::vector<ComplexType>> computeAll_split(
        bool clearTerms,
        FreqVec const& freqs = {},
        const MPI_Comm& comm = MPI_COMM_WORLD
    );

protected:

    friend class IndexContainer4<TwoParticleGF,TwoParticleGFContainer>;
    std::shared_ptr<TwoParticleGF> createElement(const IndexCombination4& Indices) const;

    const StatesClassification &S;

    const Hamiltonian &H;
    const DensityMatrix &DM;
    const FieldOperatorContainer &Operators;
};

} // namespace Pomerol

#endif // #ifndef POMEROL_INCLUDE_TWOPARTICLEGFCONTAINER_H