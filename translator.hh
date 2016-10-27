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
#ifndef LIBFOUNDED_FOUNDED_H_INCLUDED
#define LIBFOUNDED_FOUNDED_H_INCLUDED

#include <potassco/basic_types.h>
#include <potassco/theory_data.h>
#include <gringo/output/theory.hh>

using ConditionVec = std::vector<std::vector<Potassco::Lit_t>>;


class FoundedOutput : public Potassco::LpElement {
    enum class Op { Add, Sub, Mul };
    struct Define;
    struct SimpleDefine;
    struct GeneralDefine;
    struct LinearTerm;
    struct Disjunction;
    struct Assignment;
    struct Variable {
        using Domain = std::vector<std::pair<int, int>>;

        Variable(Potassco::Atom_t atom);

        Variable(Variable &&);
        Variable &operator=(Variable &&) noexcept;
        ~Variable() noexcept;

        static bool bounded(int left, int right);
        bool bounded() const;
        void extend(int left, int right);
        void unbind();

        // TODO: consider making this a function that fails if the variable is defined
        Potassco::Atom_t atom;
        Domain domain;
        bool defined = false;
    };
    using VariableSet = std::set<Potassco::Id_t>;
    using Disjunctions = std::vector<Disjunction>;
    // TODO: better put the variable id into the variable class and make this a set...
    using VariableMap = std::unordered_map<Potassco::Id_t, Variable>;
    using ShowTable = std::set<std::pair<char const *, int>>;
    using Facts = std::unordered_set<Potassco::Atom_t>;
public:
    FoundedOutput(std::ostream& out, ConditionVec &conditions, Potassco::TheoryData &data, int min, int max);
    FoundedOutput(const FoundedOutput&) = delete;
    FoundedOutput& operator=(const FoundedOutput&) = delete;
    virtual ~FoundedOutput() noexcept;
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
private:
    void rewriteDom(Potassco::TheoryAtom const &atom);
    void rewriteConstraint(Gringo::Output::TheoryData &data, Potassco::TheoryAtom const &atom);
    void rewriteShow(Potassco::TheoryAtom const &atom);
    void rewriteMinimize(Gringo::Output::TheoryData &data, Potassco::TheoryAtom const &atom);
    Potassco::Id_t rewriteTerm(Gringo::Output::TheoryData &data, Potassco::Id_t term);
    Potassco::Id_t rewriteTerm(Gringo::Output::TheoryData &data, LinearTerm const &term);
    template <class ElemFilter>
    Potassco::Atom_t rewriteAtom(Gringo::Output::TheoryData &data, Potassco::TheoryAtom const &atom, bool reMap, ElemFilter f);
    Potassco::Atom_t rewriteAtom(Gringo::Output::TheoryData &data, Potassco::TheoryAtom const &atom, bool reMap);
    VariableSet collectVariables(Potassco::TheoryAtom const &atom) const;
    void collectVariables(VariableSet &variables, Potassco::Id_t termId) const;
    void collectVariablesWeightPrio(VariableSet &vars, Potassco::Id_t termId) const;
    bool isOp(char const *op) const;
    void require(bool exp, char const *message) const;
    Potassco::TheoryElement const &requireEmptyCondition(Potassco::Id_t elemId) const;
    Variable &mapVar(Potassco::Id_t var);
    Potassco::Atom_t addSum(Gringo::Output::TheoryData &data, Potassco::Id_t term, char const *rel, Potassco::Id_t rhs);
    Potassco::Id_t addSum(Gringo::Output::TheoryData &data, LinearTerm const &term, char const *rel, Potassco::Id_t rhs);
    void addDom(Gringo::Output::TheoryData &data, Potassco::Id_t var, std::vector<std::pair<int, int>> const &dom);
    void showVariable(Gringo::Output::TheoryData &data, Potassco::Id_t varId, Variable &var, std::vector<Potassco::Id_t> &elems);
    Potassco::Id_t requireNotOperator(Potassco::Id_t termId) const;
    Potassco::Id_t requireVariable(Potassco::Id_t termId) const;
    Potassco::Id_t requireWeight(Potassco::Id_t termId) const;
    LinearTerm combine(LinearTerm &&a, LinearTerm &&b, Op op);
    LinearTerm parseLinearTerm_(Potassco::Id_t ti);
    LinearTerm parseLinearTerm(Potassco::Id_t ti);
    void printAssign(Gringo::Output::TheoryData &data, Disjunction const &assign);
    bool isFact() const;

    std::ostream& out_;
    Potassco::TheoryData &data_;
    ConditionVec &conditions_;
    Potassco::Atom_t atoms_;
    VariableMap varMap_;
    ShowTable showTable_;
    Disjunctions assign_;
    Facts facts_;
    int min_;
    int max_;
};

#endif
