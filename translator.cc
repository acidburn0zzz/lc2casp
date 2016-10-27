#include "translator.hh"
#include <potassco/theory_data.h>
#include <gringo/output/literals.hh>
#include <ostream>
#include <iostream>
#include <cstring>

#define ASSIGN "assign"

using namespace Potassco;

namespace {

// {{{1 Printing

template <class T>
struct PrintWrapper {
    PrintWrapper(T const &value, Atom_t &atoms)
    : value(value), atoms(atoms) { }
    T const &value;
    Atom_t &atoms;
};

template <class T>
PrintWrapper<T> p(T const &value, Atom_t &atoms) {
    return {value, atoms};
}

std::ostream &operator<<(std::ostream &out, PrintWrapper<LitSpan> const &p) {
    out << p.value.size;
    for (auto &&x : p.value) {
        assert(x != 0);
        out << " " << x;
        p.atoms = std::max<Atom_t>(p.atoms, std::abs(x) + 1);
    }
    return out;
}

std::ostream &operator<<(std::ostream &out, PrintWrapper<AtomSpan> const &p) {
    out << p.value.size;
    for (auto &&x : p.value) {
        assert(x != 0);
        out << " " << x;
        p.atoms = std::max(p.atoms, x + 1);
    }
    return out;
}

std::ostream &operator<<(std::ostream &out, PrintWrapper<WeightLitSpan> const &p) {
    out << p.value.size;
    for (auto &&x : p.value) {
        assert(x.lit != 0);
        out << " " << x.lit << " " << x.weight;
        p.atoms = std::max<Atom_t>(p.atoms, std::abs(x.lit) + 1);
    }
    return out;
}

std::ostream &operator<<(std::ostream &out, PrintWrapper<HeadView> const &p) {
    return out << p.value.type << " " << ::p(p.value.atoms, p.atoms);
}

std::ostream &operator<<(std::ostream &out, PrintWrapper<BodyView> const &p) {
    out << p.value.type;
    if (p.value.type != Body_t::Normal) { out << " " << p.value.bound; }
    out << " " << p.value.lits.size;
    for (auto &&x : p.value.lits) {
        assert(x.lit != 0);
        out << " " << x.lit;
        p.atoms = std::max<Atom_t>(p.atoms, std::abs(x.lit) + 1);
        if (p.value.type == Body_t::Sum) {
            out << " " << x.weight;
        }
    }
    return out;
}

class TheoryPrinter {
public:
    TheoryPrinter(Gringo::Output::TheoryData const &data, std::ostream &out)
    : data_(data)
    , out_(out) { }

    void printTerm(Potassco::Id_t termId) {
        if (seenTerms_.size() <= termId) { seenTerms_.resize(termId + 1, false); }
        if (!seenTerms_[termId]) {
            seenTerms_[termId] = true;
            auto &term = data_.data().getTerm(termId);
            switch (term.type()) {
                case Potassco::Theory_t::Number: {
                    out_ << Potassco::Directive_t::Theory << " " << Potassco::Theory_t::Number << " " << termId << " " << term.number() << "\n";
                    break;
                }
                case Potassco::Theory_t::Symbol: {
                    out_ << Potassco::Directive_t::Theory<< " " << Potassco::Theory_t::Symbol << " " << termId << " " << std::strlen(term.symbol()) << " " << term.symbol() << "\n";
                    break;
                }
                case Potassco::Theory_t::Compound: {
                    for (auto &termId : term) { printTerm(termId); }
                    if (term.isFunction()) { printTerm(term.function()); }
                    out_ << Potassco::Directive_t::Theory << " " << Potassco::Theory_t::Compound << " " << termId << " " << term.compound() << " " << term.size();
                    for (auto &termId : term) { out_ << " " << termId; }
                    out_ << "\n";
                    break;
                }
            }
        }
    }

    void debugTheoryAtom(std::ostream &out, Potassco::TheoryAtom const &atom) {
        out << "% &";
        data_.printTerm(out, atom.term());
        out << "{";
        bool comma = false;
        for (auto &elemId : atom) {
            if (comma) { out << "; ";}
            else { comma = true; }
            data_.printElem(out, elemId, [](std::ostream &out, Gringo::Output::LiteralId lit){ out << lit.offset(); });
        }
        out << "}";
        if (atom.guard()) {
            out << " ";
            data_.printTerm(out, *atom.guard());
            out << " ";
            data_.printTerm(out, *atom.rhs());
        }
        out << " <=> " << atom.atom() << std::endl;
    }
    void debugTheoryTerm(std::ostream &out, Id_t term) {
        data_.printTerm(out, term);
    }
    void printTheoryAtom(Potassco::TheoryAtom const &atom) {
        //debugTheoryAtom(std::cerr, atom);
        printTerm(atom.term());
        for (auto &elemId : atom) {
            if (seenElems_.size() <= elemId) { seenElems_.resize(elemId + 1, false); }
            if (!seenElems_[elemId]) {
                seenElems_[elemId] = true;
                auto &elem = data_.data().getElement(elemId);
                for (auto &termId : elem) { printTerm(termId); }
                std::vector<Lit_t> cond;
                for (auto &&lit : data_.getCondition(elemId)) { cond.emplace_back(lit.offset()); }
                out_ << Potassco::Directive_t::Theory << " " << Potassco::Theory_t::Element << " " << elemId << " " << elem.size();
                for (auto &termId : elem) { out_ << " " << termId; }
                out_ << " " << cond.size();
                for (auto &lit : cond) { out_ << " " << lit; }
                out_ << "\n";
            }
        }
        if (atom.guard()) {
            printTerm(*atom.rhs());
            printTerm(*atom.guard());
        }
        out_ << Potassco::Directive_t::Theory << " " << (atom.guard() ? Potassco::Theory_t::AtomWithGuard : Potassco::Theory_t::Atom) << " " << atom.atom() << " " << atom.term() << " " << atom.size();
        for (auto &elemId : atom) { out_ << " " << elemId; }
        if (atom.guard()) { out_ << " " << *atom.guard() << " " << *atom.rhs(); }
        out_ << "\n";
    }

private:
    Gringo::Output::TheoryData const &data_;
    std::ostream &out_;
    std::vector<bool> seenTerms_;
    std::vector<bool> seenElems_;
};

// {{{1 Helpers

template <class C, class F>
void groupBy(C &c, F f) {
    auto ib = std::begin(c);
    auto ie  = std::end(c);
    if (ib != ie) {
        std::sort(ib, ie);
        auto it = ib;
        while (++ib != ie) {
            if (!f(*it, *ib) && ++it != ib) {
                *it = std::move(*ib);
            }
        }
        c.erase(++it, ie);
    }
}

// }}}1

} // namespace

// {{{1 FoundedOutput::Variable

FoundedOutput::Variable::Variable(Atom_t atom)
: atom(atom) { }

FoundedOutput::Variable::Variable(Variable &&) = default;
FoundedOutput::Variable &FoundedOutput::Variable::operator=(Variable &&) noexcept = default;
FoundedOutput::Variable::~Variable() noexcept = default;

bool FoundedOutput::Variable::bounded(int left, int right) {
    return left != std::numeric_limits<int>::min() && right != std::numeric_limits<int>::max();
}

bool FoundedOutput::Variable::bounded() const {
    return domain.size() != 1 || bounded(domain.front().first, domain.front().second);
}

void FoundedOutput::Variable::extend(int left, int right) {
    if (bounded()) {
        domain.emplace_back(left, right);
    }
}

void FoundedOutput::Variable::unbind() {
    domain.clear();
    domain.emplace_back(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
}


// {{{1 FoundedOutput::LinearTerm

struct FoundedOutput::LinearTerm {
    LinearTerm(int fixed)
    : fixed(fixed) { }
    LinearTerm(int fixed, std::initializer_list<std::pair<Potassco::Id_t, int>> terms)
    : terms(terms)
    , fixed(fixed) { }

    LinearTerm(LinearTerm &&) = default;
    LinearTerm &operator=(LinearTerm &&) noexcept = default;
    ~LinearTerm() noexcept = default;

    void simplify() {
        using E = std::pair<Id_t, int>;
        groupBy(terms, [](E &a, E &b){
            if (a.first == b.first) {
                a.second+= b.second;
                return true;
            }
            return false;
        });
    }
    bool constant() const { return terms.empty(); }
    void collect(VariableSet &vars) const {
        for (auto &&term : terms) {
            vars.emplace(term.first);
        }
    }
    bool operator<(FoundedOutput::LinearTerm const &b) const {
        if (this->fixed != b.fixed) { return this->fixed < b.fixed; }
        return this->terms < b.terms;
    }

    bool operator==(FoundedOutput::LinearTerm const &b) const {
        return this->terms == b.terms && this->fixed == b.fixed;
    }

    bool operator!=(FoundedOutput::LinearTerm const &b) const {
        return !(*this == b);
    }

    std::vector<std::pair<Potassco::Id_t, int>> terms;
    int fixed;
};

// {{{1 FoundedOutput::Assignment

struct FoundedOutput::Assignment {
    Assignment(Potassco::Id_t var, LinearTerm &&left, LinearTerm &&right)
    : var(var)
    , left(std::move(left))
    , right(std::move(right)) { }

    Assignment(Assignment &&) = default;
    Assignment &operator=(Assignment &&) noexcept = default;
    ~Assignment() noexcept = default;

    bool operator==(FoundedOutput::Assignment const &b) const {
        return var == b.var && left == b.left && right == b.right;
    }
    bool operator!=(FoundedOutput::Assignment const &b) const {
        return !(*this == b);
    }
    bool operator<(FoundedOutput::Assignment const &b) const {
        if (var != b.var) { return var < b.var; }
        if (left != b.left) { return left < b.left; }
        return right < b.right;
    }

    Potassco::Id_t var;
    LinearTerm left;
    LinearTerm right;
};

// {{{1 FoundedOutput::Disjunction

struct FoundedOutput::Disjunction {
    using Elements = std::vector<Assignment>;
    Disjunction(Atom_t atom)
    : atom(atom) { }

    Disjunction(Disjunction &&) = default;
    Disjunction &operator=(Disjunction &&) noexcept = default;
    ~Disjunction() noexcept = default;

    void add(Potassco::Id_t var, LinearTerm &&left, LinearTerm &&right) {
        left.simplify();
        right.simplify();
        elems.emplace_back(var, std::move(left), std::move(right));
    }
    bool defines(FoundedOutput &out, Id_t &var) {
        if (out.facts_.find(atom) == out.facts_.end()) { return false; }
        if (elems.empty()) { return false; }
        var = elems.front().var;
        for (auto &elem : elems) {
            if (elem.var != var) { return false; }
        }
        return true;
    }
    bool operator==(FoundedOutput::Disjunction const &b) const {
        return this->atom == b.atom && elems == b.elems;
    }

    Atom_t atom;
    Elements elems;
};

// {{{1 FoundedOutput::Define

struct FoundedOutput::Define {
    virtual void encode(Gringo::Output::TheoryData &data, FoundedOutput &out, Id_t var, Variable &variable, Atom_t c) = 0;
};

// {{{1 FoundedOutput::SimpleDefine

struct FoundedOutput::SimpleDefine : FoundedOutput::Define {
    SimpleDefine(int left, int right)
    : left(left)
    , right(right) { }
    void encode(Gringo::Output::TheoryData &data, FoundedOutput &out, Id_t var, Variable &variable, Atom_t c) override {
        // Note: this translation introduces a lot of loops
        //       which are most likely unnecessary but cannot be removed by equivalence preprocessing
        // c <=> v && d
        // % is equivalent to:
        // d :- &sum { v } >= l, &sum { v } <= r.
        // v :- c.
        // d :- c.
        // c :- v, d.
        // % is equivalent to:
        // v :- c.
        //   :- c, not &sum { v } >= r.
        //   :- c, not &sum { v } <= r.
        // c :- v,     &sum { v } >= l, &sum { v } <= r.
        Atom_t l = out.addSum(data, var, ">=", data.addTerm(left));
        Atom_t r = out.addSum(data, var, "<=", data.addTerm(right));
        WeightLit_t body[3] = {{lit(c), 1}, {-lit(l), 1}, {-lit(r), 1}};
        if (!variable.defined) {
            out.rule({Head_t::Disjunctive, {&variable.atom, 1}}, {Body_t::Normal, 1, {body, 1}});
        }
        out.rule({Head_t::Disjunctive, {nullptr, 0}}, {Body_t::Normal, 1, {body, 2}});
        std::swap(body[1], body[2]);
        out.rule({Head_t::Disjunctive, {nullptr, 0}}, {Body_t::Normal, 1, {body, 2}});
        body[0] = {lit(variable.atom), 1};
        body[1].lit*= -1;
        body[2].lit*= -1;
        out.rule({Head_t::Disjunctive, {&c, 1}}, {Body_t::Normal, 1, {!variable.defined ? body : body + 1, size_t(!variable.defined ? 3 : 2)}});
    }
    int left;
    int right;
};

// {{{1 FoundedOutput::GeneralDefine

struct FoundedOutput::GeneralDefine : FoundedOutput::Define {
    GeneralDefine(LinearTerm const &left, LinearTerm const &right)
    : left(left)
    , right(right) { }
    void rule(FoundedOutput &out, std::initializer_list<Atom_t> head, std::initializer_list<Lit_t> body) {
        Atom_t h[head.size()];
        WeightLit_t b[body.size()];
        int i = 0;
        for (auto &&a : head) { h[i++] = a; }
        i = 0;
        for (auto &&l : body) { b[i++] = {l, 1}; }
        out.rule({Head_t::Disjunctive, {h, head.size()}}, {Body_t::Normal, 1, {b, body.size()}});
    }
    void rule(FoundedOutput &out, std::initializer_list<Atom_t> head, std::vector<WeightLit_t> const &body) {
        Atom_t h[head.size()];
        int i = 0;
        for (auto &&a : head) { h[i++] = a; }
        out.rule({Head_t::Disjunctive, {h, head.size()}}, {Body_t::Normal, 1, toSpan(body)});
    }
    void encode(Gringo::Output::TheoryData &data, FoundedOutput &out, Id_t var, Variable &variable, Atom_t c) override {
        // c <=> ~~a & (a => v & d)
        // % is equivalent to:
        // c => ~~a & (a => v & d)
        // ~~a & (a => v & d) => c
        // % is equivalent to:
        // c => ~~a
        // c => (a => v & d)
        // ~~a & (a => v & d) => c
        // % is equivalent to:
        // c & ~a => false
        // c & a => v
        // c & a & ~d => false
        // (a => v & d) => c | ~a
        // % is equivalent to:
        // c & ~a => false
        // c & a => v
        // c & a & ~d => false
        // ((v & d) | ~a) => (c | ~a)
        // c | ~a | a | ~v | ~d
        // % is equivalent to:
        // c & ~a => false
        // c & a => v
        // c & a & ~d => false
        // v & d & ~~a => c
        // ~~a & ~~v & ~~d => c | a
        // % is equivalent to:
        //       :- c, ~a.
        //       :- c,  a, ~d.
        //     v :- c,  a.
        //     c :-   v, d, ~~a.
        // a | c :- ~~v, d, ~~a.
        // TODO: if v is defined, then the element of the disjunction can be shifted into the rule body
        // this should be done before calling encode
        // ... :- B, (~a; ~d).
        VariableSet vars;
        left.collect(vars);
        right.collect(vars);
        Atom_t na = out.atoms_++;
        std::vector<Atom_t> a;
        for (auto &v : vars) {
            auto &&var = out.mapVar(v);
            if (!var.defined) {
                a.emplace_back(var.atom);
                rule(out, {na}, {-lit(var.atom)});
            }
        }
        // :- c, ~a.
        rule(out, {}, {lit(c), lit(na)});
        // :- c,  a, ~d.
        Atom_t l = out.addSum(data, left, "<=", var);
        Atom_t r = out.addSum(data, right, ">=", var);
        std::vector<WeightLit_t> body;
        body.push_back({lit(c), 1});
        for (auto &&l : a) {
            body.push_back({lit(l), 1});
        }
        body.push_back({-lit(l), 1});
        rule(out, {}, body);
        body.back().lit = -lit(r);
        rule(out, {}, body);
        // v :- c,  a.
        body.pop_back();
        if (!variable.defined) {
            rule(out, {variable.atom}, body);
        }
        // c :- v, d, ~~a.
        body.clear();
        if (!variable.defined) {
            body.push_back({lit(variable.atom), 1});
        }
        body.push_back({lit(l), 1});
        body.push_back({lit(r), 1});
        body.push_back({-lit(na), 1});
        rule(out, {c}, body);
        // a | c :- ~~v, d, ~~a.
        body.clear();
        if (!variable.defined) {
            Atom_t nv = out.atoms_++;
            rule(out, {nv}, {-lit(variable.atom)});
            body.push_back({-lit(nv), 1});
        }
        body.push_back({lit(l), 1});
        body.push_back({lit(r), 1});
        body.push_back({-lit(na), 1});
        for (auto &&l : a) {
            rule(out, {l, c}, body);
        }
    }
    LinearTerm const &left;
    LinearTerm const &right;
};

// {{{1 FoundedOutput

FoundedOutput::FoundedOutput(std::ostream& out, ConditionVec &conditions, TheoryData &data, int min, int max)
: out_(out)
, data_(data)
, conditions_(conditions)
, atoms_(0)
, min_(min)
, max_(max) { }
FoundedOutput::~FoundedOutput() noexcept = default;

void FoundedOutput::initProgram(bool incremental) {
    if (incremental) { throw std::runtime_error("incremental programs are not supported at the moment"); }
    out_ << "asp 1 0 0" << (incremental ? " incremental" : "") << "\n";
}

void FoundedOutput::beginStep() {
}

void FoundedOutput::rule(const HeadView& head, const BodyView& body) {
    if (head.type == Head_t::Disjunctive && head.atoms.size == 1 && body.type == Body_t::Normal && body.lits.size == 0) {
        facts_.emplace(*head.atoms.first);
    }
    out_ << Directive_t::Rule << " " << p(head, atoms_) << " " << p(body, atoms_) << "\n";
}

void FoundedOutput::minimize(Weight_t prio, const WeightLitSpan& lits) {
    out_ << Directive_t::Minimize << " " << prio << " " << p(lits, atoms_) << "\n";
}

void FoundedOutput::output(const StringSpan& str, const LitSpan& cond) {
    out_ << Directive_t::Output << " " << str.size << " ";
    out_.write(str.first, str.size) << " " << p(cond, atoms_) << "\n";
}

void FoundedOutput::assume(const LitSpan& lits) {
    out_ << Directive_t::Assume << " " << p(lits, atoms_) << "\n";
}

void FoundedOutput::external(Atom_t a, Value_t v) {
    atoms_ = std::max(atoms_, a + 1);
    out_ << Directive_t::External << " " << a << " " << v << "\n";
}

void FoundedOutput::project(const AtomSpan& atoms) {
    out_ << Directive_t::Project << " " << p(atoms, atoms_) << "\n";
}

void FoundedOutput::acycEdge(int s, int t, const LitSpan& condition) {
    out_ << Directive_t::Edge << " " << s << " " << t << " " << p(condition, atoms_) << "\n";
}

void FoundedOutput::heuristic(Atom_t a, Heuristic_t t, int bias, unsigned prio, const LitSpan& condition) {
    atoms_ = std::max(atoms_, a + 1);
    out_ << Directive_t::Heuristic << " " << t << " " << a << " " << bias << " " << prio << " " << p(condition, atoms_) << "\n";
}

void FoundedOutput::require(bool exp, char const *message) const {
    if (!exp) { throw std::runtime_error(message); }
}

// Note: for now theory elements with conditions are not supported
//       because we do not have a semantics for such language constructs yet
TheoryElement const &FoundedOutput::requireEmptyCondition(Id_t elemId) const {
    auto &&elem = data_.getElement(elemId);
    require(elem.condition() == 0, "non empty conditions are not supported");
    return elem;
}

bool FoundedOutput::isOp(char const *op) const {
    // is there an algorithm for this???
    return strcmp(op, "+") == 0 || strcmp(op, "-") == 0 || strcmp(op, "*") == 0;
}

Id_t FoundedOutput::requireNotOperator(Id_t termId) const {
    auto &&term = data_.getTerm(termId);
    if (term.type() == Theory_t::Compound) {
        require(term.isFunction() && !isOp(data_.getTerm(term.function()).symbol()), "unexpected operator");
        for (auto &&t : term) { requireNotOperator(t); }
        return termId;
    }
    return termId;
}

Id_t FoundedOutput::requireVariable(Id_t termId) const {
    auto &&term = data_.getTerm(termId);
    switch (term.type()) {
        case Theory_t::Number: { require(false, "variable expected"); }
        case Theory_t::Symbol: { break; }
        case Theory_t::Compound: {
            require(term.isFunction() && !isOp(data_.getTerm(term.function()).symbol()), "variable expected");
            for (auto &&t : term) { requireNotOperator(t); }
            return termId;
        }
    }
    return termId;
}

Id_t FoundedOutput::requireWeight(Id_t termId) const {
    auto &&term = data_.getTerm(termId);
    if (term.type() == Theory_t::Compound) {
        require(term.isFunction(), "weight expected");
        char const *name = data_.getTerm(term.function()).symbol();
        // for simplicity weights do not support +/2 and -/2
        if (strcmp(name, "*") == 0 && term.size() == 2) {
            requireWeight(*term.begin());
            requireWeight(*(term.begin() + 1));
        }
        else if (strcmp(name, "-") == 0 && strcmp(name, "+") == 0 && term.size() == 1) {
            requireWeight(*term.begin());
        }
        else {
            requireVariable(termId);
        }
    }
    return termId;
}

void FoundedOutput::collectVariablesWeightPrio(VariableSet &vars, Id_t termId) const {
    auto &&term = data_.getTerm(termId);
    if (term.type() == Theory_t::Compound) {
        require(term.isFunction(), "weight with optional priority expected");
        char const *name = data_.getTerm(term.function()).symbol();
        if (strcmp(name, "@") == 0 && term.size() == 2) {
            requireWeight(*term.begin());
            collectVariables(vars, *term.begin());
        }
        else {
            requireWeight(termId);
            collectVariables(vars, termId);
        }
    }
    else {
        collectVariables(vars, termId);
    }
}

void FoundedOutput::collectVariables(VariableSet &vars, Id_t termId) const {
    auto &&term = data_.getTerm(termId);
    switch (term.type()) {
        case Theory_t::Number: { break; }
        case Theory_t::Symbol: {
            vars.emplace(termId);
            break;
        }
        case Theory_t::Compound: {
            if (term.isFunction() && isOp(data_.getTerm(term.function()).symbol())) {
                for (auto &&t : term) {
                    collectVariables(vars, t);
                }
            }
            else {
                vars.emplace(requireVariable(termId));
            }
        }
    }
}

FoundedOutput::VariableSet FoundedOutput::collectVariables(Potassco::TheoryAtom const &atom) const {
    VariableSet vars;
    for (auto &&elemId : atom) {
        auto &&elem = requireEmptyCondition(elemId);
        require(elem.size() > 0, "not a valid constraint");
        collectVariables(vars, *elem.begin());
    }
    if (atom.guard()) {
        collectVariables(vars, *atom.rhs());
    }
    return vars;
}

FoundedOutput::Variable &FoundedOutput::mapVar(Id_t var) {
    // TODO: a variable should only receive an atom if necessary
    auto &&ret = varMap_.emplace(var, atoms_);
    if (ret.second) { ++atoms_; }
    return ret.first->second;
}

void FoundedOutput::rewriteConstraint(Gringo::Output::TheoryData &data, TheoryAtom const &atom) {
    // Outline:
    // - collect contained csp variables      { v1, ..., vn }
    // - associate an atom with each variable { a1, ..., an }
    // - let a be the atom associated with the theory atom A
    // - add rule: a :- A, a1, ..., an.
    //   where A is associated with a fresh atom a'
    require(atom.atom() > 0, "theory atoms must be associated with aspif atoms");
    VariableSet vars = collectVariables(atom);
    std::vector<WeightLit_t> body;
    body.reserve(vars.size() + 1);
    for (auto &&v : vars) {
        auto it = varMap_.find(v);
        if (it == varMap_.end()) {
            std::vector<Atom_t> head({atom.atom()});
            body.clear();
            body.push_back({lit(atom.atom()), 1});
            rule({Head_t::Disjunctive, {nullptr, 0}}, {Body_t::Normal, 1, toSpan(body)});
            return;
        }
        if (!it->second.defined) { body.push_back({lit(it->second.atom), 1}); }
    }

    body.push_back({lit(rewriteAtom(data, atom, true)), 1});
    std::vector<Atom_t> head({atom.atom()});
    rule({Head_t::Disjunctive, toSpan(head)}, {Body_t::Normal, static_cast<Weight_t>(body.size()), toSpan(body)});
}

Atom_t FoundedOutput::addSum(Gringo::Output::TheoryData &data, Id_t term, char const *rel, Id_t rhs) {
    Id_t t = rewriteTerm(data, term);
    Id_t elem = data.addElem({&t, 1}, {});
    auto &&newAtom = [&]() { return atoms_++; };
    auto &&ret = data.addAtom(
        newAtom,
        data.addTerm("sum"),
        {&elem, 1},
        data.addTerm(rel),
        rhs);
    return ret.first.atom();
}

FoundedOutput::LinearTerm FoundedOutput::combine(LinearTerm &&a, LinearTerm &&b, Op op) {
    if (op == Op::Mul) {
        if (!a.terms.empty() && !b.terms.empty()) {
            require(false, "not a linear term");
        }
        LinearTerm const &e = a.terms.empty() ? a : b;
        LinearTerm ret = a.terms.empty() ? std::move(b) : std::move(a);
        ret.fixed *= e.fixed;
        for (auto &&t : ret.terms) { t.second *= e.fixed; }
        return ret;
    }
    LinearTerm ret = std::move(a);
    for (auto &&t : b.terms) {
        ret.terms.emplace_back(t);
        if (op == Op::Sub) { ret.terms.back().second *= -1; }
    }
    ret.fixed += op == Op::Add ? b.fixed : -b.fixed;
    return ret;
}

FoundedOutput::LinearTerm FoundedOutput::parseLinearTerm(Id_t ti) {
    auto &&tp = data_.getTerm(ti);
    switch (tp.type()) {
        case Theory_t::Number: { return {tp.number()}; }
        case Theory_t::Symbol: { return {0, {{ti, 1}}}; }
        default: {
            if (tp.isFunction()) {
                char const *name = data_.getTerm(tp.function()).symbol();
                if (isOp(name)) {
                    Op op = Op::Add;
                    if (strcmp(name, "-") == 0) { op = Op::Sub; }
                    if (strcmp(name, "*") == 0) { op = Op::Mul; }
                    if (tp.size() == 2) {
                        return combine(parseLinearTerm(*tp.begin()), parseLinearTerm(*(tp.begin() + 1)), op);
                    }
                    else if (tp.size() == 1 && op != Op::Mul) {
                        return combine({0}, parseLinearTerm(*tp.begin()), op);
                    }
                }
                else {
                    return {0, {{requireVariable(ti), 1}}};
                }
            }
        }
    }
    require(false, "linear term expected");
    return {0, {}};
}

void FoundedOutput::rewriteDom(TheoryAtom const &atom) {
    constexpr char const *message = "invalid " ASSIGN " atom";
    Disjunction assign{atom.atom()};
    for (auto &&elemId : atom) {
        auto &&elem = requireEmptyCondition(elemId);
        require(elem.size() == 1, message);
        auto &&term = data_.getTerm(*elem.begin());
        require(term.type() == Theory_t::Compound && term.size() == 2 && term.isFunction(), message);
        auto &&name = data_.getTerm(term.function()).symbol();
        require(strcmp(name, ":=") == 0, message);
        auto &&var = requireVariable(*term.begin());
        auto &&rng = data_.getTerm(*(term.begin() + 1));
        Id_t termLeft, termRight;
        if (rng.type() == Theory_t::Compound && rng.size() == 2 && rng.isFunction() && strcmp(data_.getTerm(rng.function()).symbol(), "..") == 0) {
            termLeft = *rng.begin();
            termRight = *(rng.begin() + 1);
        }
        else {
            termLeft = termRight = *(term.begin() + 1);
        }
        assign.add(var, parseLinearTerm(termLeft), parseLinearTerm(termRight));
    }
    groupBy(assign.elems, std::equal_to<Disjunction::Elements::value_type>());
    assign_.emplace_back(std::move(assign));
}

void FoundedOutput::printAssign(Gringo::Output::TheoryData &data, Disjunction const &assign) {
    // Note: could be done with a simple vector as well but I am too lazy right now
    // TODO: this function should be made more general
    // - a variable dependency graph should be build here
    // - for the stratified part of the dependency graph, a better domain approximation should be calculated
    // - for example: &def { a=1..3 }. &def { b = 2 * a }.
    //   - a has domain {1, ..., 3}
    //   - b has domain {2, ..., 6} (or even better {2,4,6})
    // - variables that occur in non-trivial strongly connected components of the graph
    //   or depend on such components are passed to clingcon without an accompanying domain declaration
    //   as is done at the moment in the general case
    // TODO: furthermore, factual domain declarations could be detected.
    //   for such declarations no "founded" atoms have to be introduced because they are always guaranteed to be defined.
    //   this would allow for creating a more compact propositional representation.
    //   if all variables are defined by facts, this representation would be equivalent and as efficient as a standard constraint ASP program
    // TODO: the current implementation implements a polynomial translation
    //   it applies a tseitin-translation to the HT-formula to remove nested formulas in disjunctions and,
    //   afterward, uses strong-equivalence preserving rewritings to obtain a disjunctive logic program
    //   this translation introduces loops which very likely have a detrimental effect on the performance of constraint ASP solvers
    //   ideally, small disjunctions would simply be unfolded
    std::map<int, std::vector<std::unique_ptr<Define>>> domain;
    for (auto &&a : assign.elems) {
        auto &&dom = mapVar(a.var);
        if (a.left.constant() && a.right.constant()) {
            dom.extend(a.left.fixed, a.right.fixed);
            domain[a.var].emplace_back(new SimpleDefine(a.left.fixed, a.right.fixed));
        }
        else {
            dom.unbind();
            domain[a.var].emplace_back(new GeneralDefine(a.left, a.right));
        }
    }
    std::vector<Atom_t> head;
    WeightLit_t body = {lit(assign.atom), 1};
    for (auto &&ent : domain) {
        Variable &var = mapVar(ent.first);
        for (auto &&rng : ent.second) {
            Atom_t c = atoms_++;
            head.emplace_back(c);
            rng->encode(data, *this, ent.first, var, c);
        }
    }
    // c1, ..., cn :- a.
    rule({Head_t::Disjunctive, toSpan(head)}, {Body_t::Normal, 1, {&body, 1}});
}

void FoundedOutput::rewriteShow(TheoryAtom const &atom) {
    constexpr char const *message = "invalid show directive";
    for (auto &&elemId : atom) {
        auto &&elem = requireEmptyCondition(elemId);
        require(elem.size() == 1, message);
        auto &&term = data_.getTerm(*elem.begin());
        // NOTE: only slash notation is supported
        require(term.type() == Theory_t::Compound && term.size() == 2 && term.isFunction(), message);
        auto &&slash = data_.getTerm(term.function()).symbol();
        require(strcmp(slash, "/") == 0, message);
        auto &&name = data_.getTerm(*term.begin());
        require(name.type() == Theory_t::Symbol, message);
        auto &&arity = data_.getTerm(*(term.begin() + 1));
        require(arity.type() == Theory_t::Number, message);
        showTable_.emplace(name.symbol(), arity.number());
    }
}

void FoundedOutput::rewriteMinimize(Gringo::Output::TheoryData &data, TheoryAtom const &atom) {
    rewriteAtom(data, atom, false, [this](TheoryElement const &elem){
        require(elem.size() >= 1, "invalid minimize directive");
        VariableSet vars;
        collectVariablesWeightPrio(vars, *elem.begin());
        for (auto &&v : vars) {
            if (varMap_.find(v) == varMap_.end()) {
                return false;
            }
        }
        return true;
    });
}

Id_t FoundedOutput::rewriteTerm(Gringo::Output::TheoryData &data, Id_t termId) {
    auto &&term = data_.getTerm(termId);
    switch (term.type()) {
        case Theory_t::Number: { return data.addTerm(term.number()); }
        case Theory_t::Symbol: { return data.addTerm(term.symbol()); }
        case Theory_t::Compound: {
            std::vector<Id_t> terms;
            terms.reserve(term.size());
            for (auto &&t : term) {
                terms.emplace_back(rewriteTerm(data, t));
            }
            return term.isFunction()
                ? data.addTerm(rewriteTerm(data, term.function()), toSpan(terms))
                : data.addTerm(static_cast<Potassco::TupleType>(term.compound()), toSpan(terms));
        }
    }
    throw std::logic_error("must not happen");
}

Id_t FoundedOutput::addSum(Gringo::Output::TheoryData &data, LinearTerm const &term, char const *rel, Id_t rhs) {
    std::vector<Id_t> elems;
    Id_t te = data.addTerm(term.fixed);
    elems.emplace_back(data.addElem({&te, 1}, {}));
    for (auto &t : term.terms) {
        Id_t cv[2] = { data.addTerm(t.second), rewriteTerm(data, t.first) };
        te = data.addTerm(data.addTerm("*"), {cv, 2});
        elems.emplace_back(data.addElem({&te, 1}, {}));
    }
    auto &&newAtom = [&]() { return atoms_++; };
    auto &&ret = data.addAtom(
        newAtom,
        data.addTerm("sum"),
        toSpan(elems),
        data.addTerm(rel),
        rewriteTerm(data, rhs));
    return ret.first.atom();
}

Atom_t FoundedOutput::rewriteAtom(Gringo::Output::TheoryData &data, TheoryAtom const &atom, bool reMap) {
    return rewriteAtom(data, atom, reMap, [](TheoryElement const &) { return true; });
}

template <class ElemFilter>
Atom_t FoundedOutput::rewriteAtom(Gringo::Output::TheoryData &data, TheoryAtom const &atom, bool reMap, ElemFilter f) {
    std::vector<Id_t> elems;
    elems.reserve(atom.size());
    for (auto &elemId : atom) {
        auto &&elem = data_.getElement(elemId);
        if (f(elem)) {
            std::vector<Id_t> tuple;
            tuple.reserve(elem.size());
            for (auto &&term : elem) {
                tuple.emplace_back(rewriteTerm(data, term));
            }
            Gringo::Output::LitVec cond;
            if (elem.condition()) {
                cond.reserve(conditions_[elem.condition() - 1].size());
                for (auto &&lit : conditions_[elem.condition() - 1]) {
                    cond.emplace_back(Gringo::Output::LiteralId{
                        lit > 0 ? Gringo::NAF::POS : Gringo::NAF::NOT,
                        Gringo::Output::AtomType::Aux,
                        Potassco::atom(std::abs(lit)), 0});
                }
            }
            elems.emplace_back(data.addElem(toSpan(tuple), std::move(cond)));
        }
    }
    auto &&newAtom = [&]() {
        return atom.atom() && reMap
            ? atoms_++
            : atom.atom();
    };
    return (atom.guard()
        ? data.addAtom(newAtom, rewriteTerm(data, atom.term()), toSpan(elems), rewriteTerm(data, *atom.guard()), rewriteTerm(data, *atom.rhs()))
        : data.addAtom(newAtom, rewriteTerm(data, atom.term()), toSpan(elems))).first.atom();
}

void FoundedOutput::addDom(Gringo::Output::TheoryData &data, Id_t var, std::vector<std::pair<int, int>> const &def) {
    std::vector<Id_t> elems;
    for (auto &&d : def) {
        Id_t terms[2] = {data.addTerm(std::max(min_, d.first)), data.addTerm(std::min(max_, d.second))};
        Id_t tuple = data.addTerm(data.addTerm(".."), {terms, 2});
        elems.emplace_back(data.addElem({&tuple, 1}, {}));
    }
    auto &&newAtom = [&]() { return 0; };
    data.addAtom(
        newAtom,
        data.addTerm("dom"),
        toSpan(elems),
        data.addTerm("="),
        rewriteTerm(data, var));
}

void FoundedOutput::showVariable(Gringo::Output::TheoryData &data, Id_t varId, Variable &var, std::vector<Id_t> &elems) {
    bool show = showTable_.empty();
    if (!show) {
        auto &&var = data_.getTerm(varId);
        int arity = 0;
        char const *name;
        if (var.type() == Theory_t::Symbol) {
            name = var.symbol();
        }
        else {
            require(var.type() == Theory_t::Compound && var.isFunction(), "not a valid variable");
            name = data_.getTerm(var.function()).symbol();
            arity = var.size();
        }
        show = showTable_.find(std::make_pair(name, arity)) != showTable_.end();
    }
    if (show) {
        Id_t term = rewriteTerm(data, varId);
        Gringo::Output::LitVec cond;
        if (!var.defined) { cond.emplace_back(Gringo::Output::LiteralId{Gringo::NAF::POS, Gringo::Output::AtomType::Aux, var.atom, 0}); }
        elems.emplace_back(data.addElem({&term, 1}, std::move(cond)));
    }
}

void FoundedOutput::endStep() {
    TheoryData d;
    Gringo::Output::TheoryData data(d);
    for (auto &&atom : data_) {
        auto &&term = data_.getTerm(atom->term());
        if (term.type() == Theory_t::Symbol) {
            auto &&name = term.symbol();
            if      (strcmp(name, ASSIGN) == 0) { rewriteDom(*atom); }
            else if (strcmp(name, "show") == 0) { rewriteShow(*atom); }
        }
    }
    // TODO: detect defined variables
    for (auto &&assign : assign_) {
        Id_t var;
        if (assign.defines(*this, var)) {
            Variable &v = mapVar(var);
            if (!v.defined) { v.defined = true; }
        }
    }
    for (auto &&assign : assign_) {
        printAssign(data, assign);
    }
    std::vector<Id_t> elems;
    for (auto &&ent : varMap_) {
        showVariable(data, ent.first, ent.second, elems);
        if (!ent.second.defined) {
            if (ent.second.bounded()) {
                ent.second.extend(0, 0);
            }
            // :- not v, #sum {v} != 0.
            WeightLit_t body[2] = {{-lit(ent.second.atom), 1}, {lit(addSum(data, ent.first, "!=", data.addTerm(0))), 1}};
            rule({Head_t::Disjunctive, {nullptr, 0}}, {Body_t::Normal, 1, {body, 2}});
        }
    }
    data.addAtom(
        [&]() { return 0; },
        data.addTerm("show"),
        toSpan(elems));
    for (auto &&ent : varMap_) {
        // &dom { l1..r1; ...; ln..rn } = v.
        if (ent.second.bounded()) {
            addDom(data, ent.first, ent.second.domain);
        }
        else if (ent.second.bounded(min_, max_)) {
            addDom(data, ent.first, {{min_, max_}});
        }
    }
    for (auto &&atom : data_) {
        auto &&term = data_.getTerm(atom->term());
        if (term.type() == Theory_t::Symbol) {
            auto &&name = term.symbol();
            if      (strcmp(name, ASSIGN)     == 0) {                                 continue; }
            else if (strcmp(name, "show")     == 0) {                                 continue; }
            else if (strcmp(name, "sum")      == 0) { rewriteConstraint(data, *atom); continue; }
            else if (strcmp(name, "distinct") == 0) { rewriteConstraint(data, *atom); continue; }
            else if (strcmp(name, "minimize") == 0) { rewriteMinimize(data, *atom);   continue; }
        }
        rewriteAtom(data, *atom, false);
    }
    TheoryPrinter p(data, out_);
    for (auto &&atom : data.data()) {
        p.printTheoryAtom(*atom);
    }
    out_ << "0\n";
}

// }}}1
