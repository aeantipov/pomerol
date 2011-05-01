// This file is part of pomerol ED code
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


/** \file src/Thermal.h
** \brief Thermal object (an object which has sense only for a finite temperature).
**
** \author Igor Krivenko (igor@shg.ru)
*/
#ifndef __INCLUDE_THERMAL_H
#define __INCLUDE_THERMAL_H

#include "Misc.h"
#include "HDF5Storage.h"

struct Thermal {
    const RealType beta;

    Thermal(RealType beta);
};

#endif // endif :: #ifndef __INCLUDE_THERMAL_H
