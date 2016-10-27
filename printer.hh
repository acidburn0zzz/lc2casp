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
#ifndef LIBFOUNDED_PRINTER_H_INCLUDED
#define LIBFOUNDED_PRINTER_H_INCLUDED
#include <potassco/basic_types.h>
#include <potassco/theory_data.h>
#include <gringo/output/theory.hh>

using ConditionVec = std::vector<std::vector<Potassco::Lit_t>>;

class Printer : public Potassco::LpElement {
public:
    class Impl;
    Printer(std::ostream& out, ConditionVec &conditions, Potassco::TheoryData &data);
    Printer(const Printer&) = delete;
    Printer& operator=(const Printer&) = delete;
    virtual void initProgram(bool);
    virtual void beginStep();
    virtual void rule(const Potassco::HeadView& head, const Potassco::BodyView& body);
    virtual void minimize(Potassco::Weight_t prio, const Potassco::WeightLitSpan& lits);
    virtual void output(const Potassco::StringSpan& str, const Potassco::LitSpan& cond);
    virtual void assume(const Potassco::LitSpan& lits);
    virtual void external(Potassco::Atom_t a, Potassco::Value_t v);
    virtual void project(const Potassco::AtomSpan& atoms);
    virtual void acycEdge(int s, int t, const Potassco::LitSpan& condition);
    virtual void heuristic(Potassco::Atom_t a, Potassco::Heuristic_t t, int bias, unsigned prio, const Potassco::LitSpan& condition);
    virtual void endStep();
    virtual ~Printer() noexcept;
private:
    std::unique_ptr<Impl> impl_;
};

#endif
