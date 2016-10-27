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

class Printer : public Potassco::AbstractProgram {
public:
    class Impl;
    Printer(std::ostream& out);
    Printer(const Printer&) = delete;
    Printer& operator=(const Printer&) = delete;
    void initProgram(bool) override;
    void beginStep() override;
    void rule(Potassco::Head_t ht, const Potassco::AtomSpan& head, const Potassco::LitSpan& body) override;
    void rule(Potassco::Head_t ht, const Potassco::AtomSpan& head, Potassco::Weight_t bound, const Potassco::WeightLitSpan& body) override;
    void minimize(Potassco::Weight_t prio, const Potassco::WeightLitSpan& lits) override;
    void output(const Potassco::StringSpan& str, const Potassco::LitSpan& cond) override;
    void assume(const Potassco::LitSpan& lits) override;
    void external(Potassco::Atom_t a, Potassco::Value_t v) override;
    void project(const Potassco::AtomSpan& atoms) override;
    void acycEdge(int s, int t, const Potassco::LitSpan& condition) override;
    void heuristic(Potassco::Atom_t a, Potassco::Heuristic_t t, int bias, unsigned prio, const Potassco::LitSpan& condition) override;
    void endStep() override;

    void theoryTerm(Potassco::Id_t termId, int number) override;
    void theoryTerm(Potassco::Id_t termId, const Potassco::StringSpan& name) override;
    void theoryTerm(Potassco::Id_t termId, int cId, const Potassco::IdSpan& args) override;
    void theoryElement(Potassco::Id_t elementId, const Potassco::IdSpan& terms, const Potassco::LitSpan& cond) override;
    void theoryAtom(Potassco::Id_t atomOrZero, Potassco::Id_t termId, const Potassco::IdSpan& elements) override;
    void theoryAtom(Potassco::Id_t atomOrZero, Potassco::Id_t termId, const Potassco::IdSpan& elements, Potassco::Id_t op, Potassco::Id_t rhs) override;

    ~Printer() noexcept;
private:
    std::unique_ptr<Impl> impl_;
};

#endif
