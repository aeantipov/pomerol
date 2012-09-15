//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2011 Andrey Antipov <antipov@ct-qmc.org>
// Copyright (C) 2010-2011 Igor Krivenko <igor@shg.ru>
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


/** \file OperatorPresets.h
**  \brief Declarations of the OperatorPresets class. This is some workaround to generate easy classes
** 
**  \author    Andrey Antipov (antipov@ct-qmc.org)
*/

#ifndef __INCLUDE_OPERATOR_PRESETS_H__
#define __INCLUDE_OPERATOR_PRESETS_H__

#include "Misc.h"
#include "Operator.h"

namespace Pomerol { 
namespace OperatorPresets {

class N : public Operator {
private:
    const ParticleIndex Nmodes;
public:
    N(ParticleIndex Nmodes);
    std::map <FockState,MelemType> actRight(const FockState &ket) const;
    MelemType getMatrixElement(const FockState &bra, const FockState &ket) const;
    MelemType getMatrixElement(const FockState &ket) const;
};

class Sz : public Operator {
private:
    const int Nmodes;
    bool SimpleOrdering; // If true, spins up are in the beginning, spins down are in the end.
    std::vector<ParticleIndex> SpinUpIndices; 
    std::vector<ParticleIndex> SpinDownIndices; 
public:
    Sz(ParticleIndex Nmodes);
    Sz(const std::vector<ParticleIndex> & SpinUpIndices, const std::vector<ParticleIndex> & SpinDownIndices);
    std::map <FockState,MelemType> actRight(const FockState &ket) const;
    MelemType getMatrixElement(const FockState &bra, const FockState &ket) const;
    MelemType getMatrixElement(const FockState &ket) const;
};

class Cdag : public Operator {
private:
    ParticleIndex index;
public:
    Cdag(ParticleIndex index);
    friend std::ostream& operator<< (std::ostream& output, const Cdag& out) { output << "C^+_" << out.index; return output; };
};

class C : public Operator {
private:
    ParticleIndex index;
public:
    C(ParticleIndex index);
    friend std::ostream& operator<< (std::ostream& output, const C& out) { output << "C_" << out.index; return output; };
};


} // end of namespace OperatorPresets
} // end of namespace Pomerol

#endif // endif :: #ifndef __INCLUDE_OPERATOR_PRESETS_H__
