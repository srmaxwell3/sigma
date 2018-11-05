#include <cstdio>
#include <map>
using std::map;
#include <sstream>
using std::ostringstream;
#include <string>
using std::string;
#include <utility>
using std::pair;
using std::make_pair;
#include <vector>
using std::vector;

enum NonTerminalOperator {
  opIllegal,

  opDot,

  opAdd,
  opSub,
  opMul,
  opDiv,

  eoNonTerminalOperator
};

// Opr < Nam < Lit

class Expr {
 public:
  Expr() { }
  virtual ~Expr() { }
  virtual bool operator<(Expr const &that) const { return false; }
  // virtual bool operator==(Expr const &that) const { return this == &that; }
  virtual bool IsCongruentTo(Expr const *that) const { return false; }
  virtual bool IsANonTerminal() const { return false; }
  virtual NonTerminalOperator Operator() const { return opIllegal; }
  virtual size_t NumberOfOperands() const { return 0; }
  virtual Expr const *Operand(size_t i) const { static Expr expr; return &expr; }
  virtual bool IsALiteral() const { return false; }
  virtual int Literal() const { return 0; }
  virtual bool IsAName() const { return false; }
  virtual string const &Name() const { static string empty; return empty; }
  virtual size_t Cost() const { return 0; }
  virtual string ToString() const { return "Expr()"; }
};

class Lit : public Expr {
 public:
  Lit(int _literal) : literal(_literal) {}

  // bool operator==(Expr const &that) const {
  //   if (Lit const *t = dynamic_cast<Lit const *>(&that)) {
  //     return literal == t->literal;
  //   }
  //   return false;
  // }
  bool operator<(Expr const &that) const {
    if (Lit const *t = dynamic_cast<Lit const *>(&that)) {
      return literal < t->literal;
    }
    return false;
  }

  bool IsCongruentTo(Expr const *that) const {
    if (Lit const *t = dynamic_cast<Lit const *>(that)) {
      return true;
    }
    return false;
  }

  bool IsALiteral() const { return true; }
  size_t Cost() const { return 1; }
  int Literal() const {
    return literal;
  }
  virtual string ToString() const {
    ostringstream result;
    result << literal;
    return result.str();
  }
 private:
  int literal;
};

class Nam : public Expr {
 public:
  Nam(string const &_name) : name(_name) {}

  // Nam < Lit
  bool operator<(Expr const &that) const {
    if (Lit const *t = dynamic_cast<Lit const *>(&that)) {
      return true;
    }
    if (Nam const *t = dynamic_cast<Nam const *>(&that)) {
      return name < t->name;
    }
    return false;
  }
  // bool operator==(Expr const &that) const {
  //   if (Nam const *t = dynamic_cast<Nam const *>(&that)) {
  //     return name == t->name;
  //   }
  //   return false;
  // }

  bool IsCongruentTo(Expr const *that) const {
    if (Nam const *t = dynamic_cast<Nam const *>(that)) {
      return name == t->name;
    }
    return false;
  }
  bool IsAName() const { return true; }
  size_t Cost() const { return 1; }
  string const &Name() const {
    return name;
  }
  virtual string ToString() const {
    ostringstream result;
    result << name;
    return result.str();
  }
 private:
  string name;
};

class Opr : public Expr {
 public:
  Opr(NonTerminalOperator _opr, vector<Expr *> _operands) : opr(_opr), operands(_operands) {}

  // Nam < Lit
  bool operator<(Expr const *that) const {
    if (Lit const *t = dynamic_cast<Lit const *>(that)) {
      return true;
    }
    if (Nam const *t = dynamic_cast<Nam const *>(that)) {
      return true;
    }
    if (Opr const *t = dynamic_cast<Opr const *>(that)) {
      if (opr < t->opr) {
        return true;
      } else if (opr == t->opr) {
        if (NumberOfOperands() < t->NumberOfOperands()) {
          return true;
        } else if (NumberOfOperands() == t->NumberOfOperands()) {
          for (size_t i = 0; i < NumberOfOperands(); i += 1) {
            if (operands[i] < t->operands[i]) {
              return true;
            }
            if (t->operands[i] < operands[i]) {
              return false;
            }
          }
        }
      }
    }
    return false;
  }
  // bool operator==(Expr const *that) const {
  //   if (Opr const *t = dynamic_cast<Opr const *>(that)) {
  //     if (opr == t->opr) {
  //       if (NumberOfOperands() == t->NumberOfOperands()) {
  //         for (size_t i = 0; i < NumberOfOperands(); i += 1) {
  //           if (operands[i] != t->operands[i]) {
  //             return false;
  //           }
  //         }
  //         return true;
  //       }
  //     }
  //   }
  //   return false;
  // }

  bool IsCongruentTo(Expr const *that) const {
    if (Opr const *t = dynamic_cast<Opr const *>(that)) {
      if (opr == t->opr) {
        if (NumberOfOperands() == t->NumberOfOperands()) {
          for (size_t i = 0; i < NumberOfOperands(); i += 1) {
            if (operands[i] != t->operands[i]) {
              return false;
            }
          }
          return true;
        }
      }
    }
    return false;
  }
  bool IsANonTerminal() const { return true; }
  size_t Cost() const {
    size_t cost = 1;
    for (auto const &operand : operands) {
      cost += operand->Cost();
    }
    return cost;
  }
  NonTerminalOperator Operator() const {
    return opr;
  }
  Expr *Operand(size_t i) const { return operands.at(i); }
  size_t NumberOfOperands() const { return operands.size(); }

 protected:
  NonTerminalOperator opr;
  vector<Expr *> operands;
};

class Dot : public Opr {
 public:
  Dot(Expr *_operand) : Opr(opDot, { _operand }) { }
  virtual string ToString() const {
    ostringstream result;
    result << "." << operands[0]->ToString();
    return result.str();
  }
};

class Add : public Opr {
 public:
  Add(vector<Expr *> _operands) : Opr(opAdd, _operands) { }
  virtual string ToString() const {
    ostringstream result;
    result << "+(" << operands[0]->ToString();
    for (size_t i = 1; i < NumberOfOperands(); i += 1) {
      result << "," << operands[i]->ToString();
    }
    result << ")";
    return result.str();
  }
};

class Mul : public Opr {
 public:
  Mul(vector<Expr *> _operands) : Opr(opMul, _operands) { }
  virtual string ToString() const {
    ostringstream result;
    result << "*(" << operands[0]->ToString();
    for (size_t i = 1; i < NumberOfOperands(); i += 1) {
      result << "," << operands[i]->ToString();
    }
    result << ")";
    return result.str();
  }
};

double F(Expr const *e1, Expr const *e2) {
  // fprintf(stdout, "F(e1 = %s, e2 = %s)\n", e1->ToString().c_str(), e2->ToString().c_str());

  // N is the set of non-terminals of E.
  // T is the set of terminals.
  // L is the set of literals, and
  // I the set of names.

  if (e1->IsANonTerminal() ^ e2->IsANonTerminal()) {
    return 1.0;
  }
  if (e1->IsANonTerminal()) {
    if (e1->Operator() != e2->Operator()) {
      return 1.0;
    }
    if (e1->NumberOfOperands() != e2->NumberOfOperands()) {
      return 1.0;
    }
    double s = 0.0;
    for (size_t i = 0; i < e1->NumberOfOperands(); i += 1) {
      s += F(e1->Operand(i), e2->Operand(i));
    }
    return s;
  }
  if (e1->IsALiteral() && e2->IsALiteral()) {
    if (e1->Literal() != e2->Literal()) {
      return 1.0;
    }
    return 0.0;
  }
  if (e1->IsAName() && e2->IsAName()) {
    if (e1->Name() != e2->Name()) {
      return 1.0;
    }
    return 0.0;
  }
  return 0.0;
}

double G(Expr const *e1, Expr const *e2) {
  // fprintf(stdout, "G(e1 = %s, e2 = %s)\n", e1->ToString().c_str(), e2->ToString().c_str());

  // N is the set of non-terminals of E.
  // T is the set of terminals.
  // L is the set of literals, and
  // I the set of names.

  if (e1->IsANonTerminal() ^ e2->IsANonTerminal()) {
    return 1.0;
  }
  if (e1->IsANonTerminal()) {
    if (e1->Operator() != e2->Operator()) {
      return 1.0;
    }
    if (e1->NumberOfOperands() != e2->NumberOfOperands()) {
      return 1.0;
    }
    double s = 0.0;
    for (size_t i = 0; i < e1->NumberOfOperands(); i += 1) {
      s += (1.0 / e1->NumberOfOperands()) * G(e1->Operand(i), e2->Operand(i));
    }
    return s;
  }
  if (e1->IsALiteral() && e2->IsALiteral()) {
    if (e1->Literal() != e2->Literal()) {
      return 1.0;
    }
    return 0.0;
  }
  if (e1->IsAName() && e2->IsAName()) {
    if (e1->Name() != e2->Name()) {
      return 1.0;
    }
    return 0.0;
  }
  return 0.0;
}

// Whenever the recursive subroutine S encounters a pair of dissimilar
// subexpressions of e and e' it calls the subroutine TRYPARMS. The
// subroutine TRYPARMS determines if a new parameter must be
// created. If so, it increments NPARMS by one and decreases COSTSAV
// by the estimated object code size of the parameter subexpression of
// e. When control returns to SIGMA from S, the variable NPARMS
// contains the number of parameters necessary to evaluate the pair e,
// e' by a common code sequence and COSTSAV contains an estimate of
// the size of the object code sequence sharable by the expressions.

size_t NParams;
size_t Cost;
size_t CostSav;
size_t M;

void TryParams
    (Expr const *e1,
     Expr const *e2,
     map<string, size_t> (&params)[2],
     map<pair<string, string>, size_t> &pairs
    )
{
  auto &e1Count = params[0][e1->ToString()];
  auto &e2Count = params[1][e2->ToString()];
  auto &pairCount = pairs[make_pair(e1->ToString(), e2->ToString())];

  // if (e1Count == 0 /* && e2Count == 0 */) {
  //   NParams += 1;
  // }
  e1Count += 1;
  CostSav -= e1->Cost();
  e2Count += 1;
  CostSav -= e2->Cost();
  pairCount += 1;
  NParams = pairs.size();
  M = params[0].size() + params[1].size();
  fprintf(stdout, "TryParams(%s [%lu, %lu], %s [%lu, %lu]): NParams = %lu, Cost = %lu, CostSav = %lu, M = %lu\n",
          e1->ToString().c_str(),
          e1Count,
          e1->Cost(),
          e2->ToString().c_str(),
          e2Count,
          e2->Cost(),
          NParams,
          Cost,
          CostSav,
          M
         );
}

void S
    (Expr const *e1,
     Expr const *e2,
     map<string, size_t> (&params)[2],
     map<pair<string, string>, size_t> &pairs
    )
{
  // fprintf(stdout, "S(e1 = %s, e2 = %s)\n", e1->ToString().c_str(), e2->ToString().c_str());

  // N is the set of non-terminals of E.
  // T is the set of terminals.
  // L is the set of literals, and
  // I the set of names.

  if (e1->IsANonTerminal() ^ e2->IsANonTerminal()) {
    TryParams(e1, e2, params, pairs);
    return;
  }
  if (e1->IsANonTerminal()) {
    if (e1->Operator() != e2->Operator()) {
      TryParams(e1, e2, params, pairs);
      return;
    }
    if (e1->NumberOfOperands() != e2->NumberOfOperands()) {
      TryParams(e1, e2, params, pairs);
      return;
    }
    if (!e1->IsCongruentTo(e2)) {
      for (size_t i = 0; i < e1->NumberOfOperands(); i += 1) {
        S(e1->Operand(i), e2->Operand(i), params, pairs);
      }
    }
    return;
  }
  if (e1->IsALiteral() && e2->IsALiteral()) {
    if (e1->Literal() != e2->Literal()) {
      TryParams(e1, e2, params, pairs);
    }
    return;
  }
  if (e1->IsAName() && e2->IsAName()) {
    if (e1->Name() != e2->Name()) {
      TryParams(e1, e2, params, pairs);
    }
    return;
  }
  TryParams(e1, e2, params, pairs);
}

double Sigma(Expr const *e1, Expr const *e2) {
  fprintf(stdout, "Sigma(e1 = %s, e2 = %s)\n", e1->ToString().c_str(), e2->ToString().c_str());

  // The subroutine S does a coordinated tree walk on the expressions
  // e and e' setting the variables NPARMS to the nurnber of
  // parameters and COSTSAV to the amount of code saved by the shared
  // code sequences. The subroutine TRYPARMS (not defined here)
  // increments NPARMS and decrements COSTSAV by e[cost] if a new
  // parameter must be created, e[cost] is the amount of code
  // necessary to evaluate the entire expression e.  e[count] is the
  // number of formally identical instances of this expression.

  if (e1->IsCongruentTo(e2)) {
    return 0.0;
  }

  CostSav = Cost = e1->Cost();
  map<string, size_t> params[2];
  map<pair<string, string>, size_t> pairs;
  S(e1, e2, params, pairs);
  for (auto const &p : params[0]) {
    fprintf(stdout, "    [0] %2lu: %s\n", p.second, p.first.c_str());
  }
  for (auto const &p : params[1]) {
    fprintf(stdout, "    [1] %2lu: %s\n", p.second, p.first.c_str());
  }
  for (auto const &p : pairs) {
    fprintf(stdout, "    [p] %2lu: %s; %s\n", p.second, p.first.first.c_str(), p.first.second.c_str());
  }

  // The final expression in the body of SIGMA requires some
  // explanation.  The numerator is the estimated cost in code size of
  // the overhead required to set up parameters (M * NPARMS), call (+
  // M), and return (+ 1) from a similarity-created subroutine.  The
  // denominator is the amount of code saved by replacing M - 1 of the
  // expressions with calls to a common sequence of code.  Hence if
  // SIGMA(e, e') < 1, then code size will be reduced by implementing
  // e and e' as calls on a common subroutine.  The application of the
  // similarity function SIGMA (more precisely its subroutine S)
  // partitions an expression into a body and a collection of
  // parameter expressions.  In subsequent discussions, body(e) refers
  // to the expression resulting from the removal of the parameter
  // nodes in e, and parms(e) refers to the set of sub-expressions
  // identified by S as parameters of e.

  return double(M * NParams + M + 1) / double((M - 1) * (Cost - CostSav));
}

int main(int argc, char *const argv[]) {
  // Example
  // Define:

  // e1: (.A + 1) * (.A + 2) * (.A + 3),
  Expr *e1 =
    new Mul({ new Add({ new Dot(new Nam("A")), new Lit(1) }),
              new Add({ new Dot(new Nam("A")), new Lit(2) }),
              new Add({ new Dot(new Nam("A")), new Lit(3) })
            }
           );
  // e2: (.B + 1) * (.B + 2) * (.B + 3),
  Expr *e2 =
    new Mul({ new Add({ new Dot(new Nam("B")), new Lit(1) }),
              new Add({ new Dot(new Nam("B")), new Lit(2) }),
              new Add({ new Dot(new Nam("B")), new Lit(3) })
            }
           );
  // e3: (.A + 1) * (.B + 2) * (.C + 3),
  Expr *e3 =
    new Mul({ new Add({ new Dot(new Nam("A")), new Lit(1) }),
              new Add({ new Dot(new Nam("B")), new Lit(2) }),
              new Add({ new Dot(new Nam("C")), new Lit(3) })
            }
           );

  // e4: .A + .B * (.C + .E),
  Expr *e4 =
    new Add({ new Dot(new Nam("A")),
              new Mul({ new Dot(new Nam("B")),
                        new Add({ new Dot(new Nam("C")),
                                  new Dot(new Nam("E"))
                                }
                               )
                      }
                     )
            }
           );
  // e5: .A + .B * (.D + .E),
  Expr *e5 =
      new Add({ new Dot(new Nam("A")),
                new Mul({ new Dot(new Nam("B")),
                          new Add({ new Dot(new Nam("D")),
                                    new Dot(new Nam("E"))
                                  }
                                 )
                        }
                       )
              }
             );

  // e6: .A +.B,
  Expr *e6 = new Add({ new Dot(new Nam("A")), new Dot(new Nam("B")) });
  // e7: .A + .C
  Expr *e7 = new Add({ new Dot(new Nam("A")), new Dot(new Nam("C")) });

  // The following table shows the values returned from F, G, and SIGMA.
  // f         F   G    SIGMA   (M * NParams + M + 1) / (M - 1) * CostSav
  // f(e1,e2) 3.0 0.5   0.625 = (2 * 1       + 2 + 1) / (2 - 1) * 8
  // f(e1,e3) 3.0 0.5   1.125 = (2 * 3       + 2 + 1) / (2 - 1) * 8
  // f(e4,e5) 1.0 0.125 1.5 =   (2 * 1       + 2 + 1) / (2 - 1) * 4
  // f(e6,e7) 1.0 0.5   2.5 =   (2 * 1       + 2 + 1) / (2 - 1) * 2

  // f         F   G    SIGMA   (M * NParams + M + 1) / (M - 1) * CostSav
  // f(e1,e2) 3.0 0.5   0.625 = (2 * 1       + 2 + 1) / (2 - 1) * 8
  // e1: (.A + 1) * (.A + 2) * (.A + 3),
  // e2: (.B + 1) * (.B + 2) * (.B + 3),
  fprintf(stdout,
          "f(e1, e2) %.3f %.3f %.3f",
          F(e1, e2),
          G(e1, e2),
          Sigma(e1, e2)
         );
  fprintf(stdout,
          " = (%lu * %lu       + %lu + 1) / (%lu - 1) * %lu\n\n",
          M,
          NParams,
          M,
          M,
          CostSav // (Cost - CostSav)
         );

  // f         F   G    SIGMA   (M * NParams + M + 1) / (M - 1) * CostSav
  // f(e1,e3) 3.0 0.5   1.325 = (2 * 3       + 2 + 1) / (2 - 1) * 8
  // e1: (.A + 1) * (.A + 2) * (.A + 3),
  // e3: (.A + 1) * (.B + 2) * (.C + 3),
  fprintf(stdout,
          "f(e1, e3) %.3f %.3f %.3f",
          F(e1, e3),
          G(e1, e3),
          Sigma(e1, e3)
         );
  fprintf(stdout,
          " = (%lu * %lu       + %lu + 1) / (%lu - 1) * %lu\n\n",
          M,
          NParams,
          M,
          M,
          CostSav // (Cost - CostSav)
         );

  // f         F   G    SIGMA   (M * NParams + M + 1) / (M - 1) * CostSav
  // f(e4,e5) 1.0 0.325 1.5 =   (2 * 1       + 2 + 1) / (2 - 1) * 4
  // e4: .A + .B * (.C + .E),
  // e5: .A + .B * (.D + .E),
  fprintf(stdout,
          "f(e4, e5) %.3f %.3f %.3f",
          F(e4, e5),
          G(e4, e5),
          Sigma(e4, e5)
         );
  fprintf(stdout,
          " = (%lu * %lu       + %lu + 1) / (%lu - 1) * %lu\n\n",
          M,
          NParams,
          M,
          M,
          CostSav // (Cost - CostSav)
         );

  // f         F   G    SIGMA   (M * NParams + M + 1) / (M - 1) * CostSav
  // f(e6,e7) 1.0 0.5   2.5 =   (2 * 1       + 2 + 1) / (2 - 1) * 2
  // e6: .A +.B,
  // e7: .A + .C
  fprintf(stdout,
          "f(e6, e7) %.3f %.3f %.3f",
          F(e6, e7),
          G(e6, e7),
          Sigma(e6, e7)
         );
  fprintf(stdout,
          " = (%lu * %lu       + %lu + 1) / (%lu - 1) * %lu\n\n",
          M,
          NParams,
          M,
          M,
          CostSav // (Cost - CostSav)
         );
  return 0;
}

