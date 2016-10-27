//
// Copyright (c) 2015, Roland Kaminski
//
// This file is part of Potassco. See http://potassco.sourceforge.net/
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef LIBFOUNDED_ASPIFC_H_INCLUDED
#define LIBFOUNDED_ASPIFC_H_INCLUDED
#include <potassco/aspif.h>
#include <vector>

using ConditionVec = std::vector<std::vector<Potassco::Lit_t>>;

class AspifCInput : public Potassco::AspifInput {
public:
    AspifCInput(Potassco::LpElement& out, ConditionVec &conditions, Potassco::TheoryData& theory)
    : Potassco::AspifInput(out, &theory)
    , conditions_(conditions) { }
    virtual ~AspifCInput() noexcept = default;
protected:
    virtual unsigned newTheoryCondition(const Potassco::LitSpan& lits) {
        if (lits.size == 0) { return 0; }
        else {
            conditions_.emplace_back(Potassco::begin(lits), Potassco::end(lits));
            return conditions_.size();
        }
    }
private:
    ConditionVec &conditions_;
};

#endif
