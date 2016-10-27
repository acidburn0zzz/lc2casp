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
#include "printer.hh"
#include <cstring>

using namespace Potassco;
using ConditionVec = std::vector<std::vector<Potassco::Lit_t>>;

namespace {

struct Rule {
    std::vector<Atom_t> head;
    std::vector<WeightLit_t> body;
    Head_t headType;
    Body_t bodyType;
    Weight_t bound;
};

struct Output {
    std::string str;
    std::vector<Lit_t> body;
};

struct Minimize {
    Weight_t prio;
    std::vector<WeightLit_t> body;
};

struct Assume {
    std::vector<Lit_t> body;
};

struct External {
    Atom_t a;
    Value_t v;
};

struct Project {
    std::vector<Atom_t> body;
};

struct Acyc {
    int s;
    int t;
    std::vector<Lit_t> body;
};

struct Heuristic {
    Atom_t a;
    Heuristic_t t;
    int bias;
    unsigned prio;
    std::vector<Lit_t> body;
};

} // namespace

class Printer::Impl {
public:
    Impl(std::ostream& out)
    : out_(out) { }

    void add(Rule &&rule) {
        rules_.emplace_back(std::move(rule));
    }
    void add(Minimize &&m) {
        minimize_.emplace_back(std::move(m));
    }
    void add(Output &&output) {
        if (output.body.size() == 1 && output.body.front() > 0) {
            named_.emplace(output.body.front(), output.str);
        }
        else {
            output_.emplace_back(std::move(output));
        }
    }
    void add(Assume &&a) {
        assume_.emplace_back(std::move(a));
    }
    void add(Project &&a) {
        project_.emplace_back(std::move(a));
    }
    void add(External &&a) {
        external_.emplace_back(std::move(a));
    }
    void add(Heuristic &&a) {
        heuristic_.emplace_back(std::move(a));
    }
    void add(Acyc &&a) {
        acyc_.emplace_back(std::move(a));
    }
    void print() {
        for (auto &a : data_) {
            if (a->atom()) {
                atoms_[a->atom()] = a;
            }
            else {
                print(*a);
                out_ << ".\n";
            }
        }
        for (auto &r : rules_) { print(r); }
        for (auto &m : minimize_) { print(m); }
        for (auto &a : assume_) { print(a); }
        for (auto &h : heuristic_) { print(h); }
        for (auto &p : project_) { print(p); }
        for (auto &e : external_) { print(e); }
        for (auto &a : acyc_) { print(a); }
    }
    Potassco::TheoryData &data() { return data_; }
    Potassco::Id_t addCondition(const Potassco::LitSpan& cond) {
        if (cond.size > 0) {
            conditions_.emplace_back(cond.first, cond.first + cond.size);
            return conditions_.size();
        }
        return 0;
    }
private:
    void print(Rule const &r) {
        if (r.headType == Head_t::Choice) {
            out_ << "{ ";
        }
        bool comma = false;
        for (auto &h : r.head) {
            if (comma) { out_ << "; "; }
            else       { comma = true; }
            print(h);
        }
        if (r.headType == Head_t::Choice) {
            out_ << " }";
        }
        if (r.bodyType == Body_t::Sum || !r.body.empty() || (r.head.empty() && r.headType == Head_t::Disjunctive)) {
            out_ << " :- ";
            comma = false;
            if (r.bodyType == Body_t::Sum) {
                out_ << "[ ";
            }
            for (auto &b : r.body) {
                if (comma) { out_ << ", "; }
                else       { comma = true; }
                if (r.bodyType == Body_t::Sum) { print(b); }
                else { print(b.lit); }

            }
            if (r.bodyType == Body_t::Sum) {
                out_ << "] >= " << r.bound;
            }
        }
        out_ << ".\n";
    }
    void print(Output const &o) {
        out_ << "#show " << o.str;
        print(o.body);
        out_ << ".\n";
    }
    void print(Minimize const &o) {
        bool comma = false;
        out_ << "#minimize [ ";
        for (auto &b : o.body) {
            if (comma) { out_ << ", "; }
            else       { comma = true; }
            print(b);
        }
        out_ << " ] @ " << o.prio << ".\n";
    }
    void print(Assume const &a) {
        out_ << "#assume { ";
        print(a.body, false);
        out_ << " }.\n";
    }
    void print(Project const &p) {
        out_ << "#project { ";
        bool comma = false;
        for (auto &a : p.body) {
            if (comma) { out_ << ", "; }
            else       { comma = true; }
            print(a);
        }
        out_ << " }.\n";
    }
    void print(External const &e) {
        out_ << "#external ";
        print(e.a);
        out_ << "=";
        switch (e.v) {
            case Value_t::True:    { out_ << "true"; break; }
            case Value_t::False:   { out_ << "false"; break; }
            case Value_t::Free:    { out_ << "free"; break; }
            case Value_t::Release: { out_ << "release"; break; }
        }
        out_ << ".\n";
    }
    void print(Acyc const &a) {
        out_ << "#acyc ";
        out_ << "(" << a.s << "," << a.t << ")";
        print(a.body);
        out_ << ".\n";
    }
    void print(Heuristic const &a) {
        out_ << "#heuristic ";
        print(a.a);
        out_ << " ";
        print(a.body);
        out_ << "[ ";
        switch (a.t) {
            case Heuristic_t::Level:  { out_ << "level"; break; }
            case Heuristic_t::Sign:   { out_ << "sign"; break; }
            case Heuristic_t::Factor: { out_ << "factor"; break; }
            case Heuristic_t::Init:   { out_ << "init"; break; }
            case Heuristic_t::True:   { out_ << "true"; break; }
            case Heuristic_t::False:  { out_ << "false"; break; }
        }
        out_ << "," << a.bias << "@" << a.prio << " ].\n";
    }
    void print(std::vector<Lit_t> const &lits, bool cond = true) {
        if (!lits.empty()) {
            if (cond) { out_ << ": "; }
            bool comma = false;
            for (auto &l : lits) {
                if (comma) { out_ << ", "; }
                else       { comma = true; }
                print(l);
            }
        }
    }
    void print(WeightLit_t l) {
        print(l.lit);
        out_ << "=" << l.weight;
    }
    void print(Lit_t l) {
        if (l < 0) { out_ << "not "; }
        print(atom(l));
    }
    void print(Atom_t a) {
        auto it = named_.find(a);
        if (it != named_.end()) {
           out_ << it->second;
        }
        else {
            auto jt = atoms_.find(a);
            if (jt != atoms_.end()) {
                print(*jt->second);
            }
            else {
                out_ << "x_" << a;
            }
        }
    }
    void print(TheoryAtom const &atom) {
        out_ << "&";
        print(data_.getTerm(atom.term()));
        out_ << " { ";
        bool comma = false;
        for (auto &elemId : atom) {
            if (comma) { out_ << "; ";}
            else { comma = true; }
            print(data_.getElement(elemId));
        }
        out_ << " }";
        if (atom.guard()) {
            out_ << " ";
            print(data_.getTerm(*atom.guard()));
            out_ << " ";
            print(data_.getTerm(*atom.rhs()));
        }
    }

    void print(TheoryElement const &elem) {
        Gringo::print_comma(out_, elem, ", ", [this](std::ostream &, Id_t termId){ print(data_.getTerm(termId)); });
        if (elem.size() == 0 && !elem.condition()) {
            out_ << ": ";
        }
        else if (elem.condition()) {
            print(conditions_[elem.condition() - 1]);
        }
    }

    void print(TheoryTerm const &term) const {
        switch (term.type()) {
            case Theory_t::Number: {
                if (term.number() < 0) { out_ << "("; }
                out_ << term.number();
                if (term.number() < 0) { out_ << ")"; }
                break;
            }
            case Theory_t::Symbol: { out_ << term.symbol(); break; }
            case Theory_t::Compound: {
                Tuple_t x;
                auto parens = toString(term.isTuple() ? term.tuple() : Tuple_t::Paren);
                bool isOp = false;
                if (term.isFunction()) {
                    auto &name = data_.getTerm(term.function());
                    char buf[2] = { *name.symbol(), 0 };
                    isOp = term.size() <= 2 && std::strpbrk(buf, "/!<=>+-*\\?&@|:;~^.");
                    if (!isOp) { print(data_.getTerm(term.function())); }
                }
                out_ << parens[0];
                if (isOp && term.size() <= 1) {
                    print(data_.getTerm(term.function()));
                }
                Gringo::print_comma(out_, term, isOp ? data_.getTerm(term.function()).symbol() : ",", [this](std::ostream &, Id_t termId){ print(data_.getTerm(termId)); });
                if (term.isTuple() && term.tuple() == Tuple_t::Paren && term.size() == 1) { out_ << ","; }
                out_ << parens[1];
                break;
            }
        }
    }

private:
    std::ostream& out_;
    TheoryData data_;
    ConditionVec conditions_;
    std::vector<Rule> rules_;
    std::vector<Output> output_;
    std::vector<Minimize> minimize_;
    std::vector<Assume> assume_;
    std::vector<External> external_;
    std::vector<Heuristic> heuristic_;
    std::vector<Acyc> acyc_;
    std::vector<Project> project_;
    std::unordered_map<Atom_t, std::string> named_;
    std::unordered_map<Atom_t, TheoryAtom const*> atoms_;
};

Printer::Printer(std::ostream& out)
: impl_(new Printer::Impl(out)) { }

void Printer::initProgram(bool) {
}

void Printer::beginStep() {
}

void Printer::rule(Potassco::Head_t ht, const Potassco::AtomSpan& head, const Potassco::LitSpan& body) {
    std::vector<WeightLit_t> wlits;
    for (auto &x : body) {
        wlits.emplace_back(Potassco::WeightLit_t{x, 1});
    }
    impl_->add({
        {head.first, head.first + head.size},
        wlits,
        ht, Potassco::Body_t::Normal,
        0});
}

void Printer::rule(Potassco::Head_t ht, const Potassco::AtomSpan& head, Potassco::Weight_t bound, const Potassco::WeightLitSpan& body) {
    impl_->add({
        {head.first, head.first + head.size},
        {body.first, body.first + body.size},
        ht, Potassco::Body_t::Sum,
        bound});
}

void Printer::minimize(Weight_t prio, const WeightLitSpan& lits) {
    impl_->add({
        prio,
        {lits.first, lits.first + lits.size}});
}

void Printer::output(const StringSpan& str, const LitSpan& lits) {
    impl_->add(Output{
        {str.first, str.first + str.size},
        {lits.first, lits.first + lits.size}});
}

void Printer::assume(const LitSpan& lits) {
    impl_->add(Assume{{lits.first, lits.first + lits.size}});
}

void Printer::external(Atom_t a, Value_t v) {
    impl_->add(External{a, v});
}

void Printer::project(const AtomSpan& atoms) {
    impl_->add(Project{{atoms.first, atoms.first + atoms.size}});
}

void Printer::acycEdge(int s, int t, const LitSpan& body) {
    impl_->add(Acyc{s, t, {body.first, body.first + body.size}});
}

void Printer::heuristic(Atom_t a, Heuristic_t t, int bias, unsigned prio, const LitSpan& body) {
    impl_->add(Heuristic{a, t, bias, prio, {body.first, body.first + body.size}});
}

void Printer::endStep() {
    impl_->print();
}

void Printer::theoryTerm(Potassco::Id_t termId, int number) {
    impl_->data().addTerm(termId, number);
}

void Printer::theoryTerm(Potassco::Id_t termId, const Potassco::StringSpan& name) {
    impl_->data().addTerm(termId, name);
}

void Printer::theoryTerm(Potassco::Id_t termId, int cId, const Potassco::IdSpan& args) {
    impl_->data().addTerm(termId, cId, args);
}

void Printer::theoryElement(Potassco::Id_t elementId, const Potassco::IdSpan& terms, const Potassco::LitSpan& cond) {
    impl_->data().addElement(elementId, terms, impl_->addCondition(cond));
}

void Printer::theoryAtom(Potassco::Id_t atomOrZero, Potassco::Id_t termId, const Potassco::IdSpan& elements) {
    impl_->data().addAtom(atomOrZero, termId, elements);
}

void Printer::theoryAtom(Potassco::Id_t atomOrZero, Potassco::Id_t termId, const Potassco::IdSpan& elements, Potassco::Id_t op, Potassco::Id_t rhs) {
    impl_->data().addAtom(atomOrZero, termId, elements, op, rhs);
}

Printer::~Printer() noexcept = default;

