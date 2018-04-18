#include <Rcpp.h>
using namespace Rcpp;

struct SimpleAbstract{
  int a = 0;
};

struct SimpleDerived1 : SimpleAbstract {
  SimpleDerived1() : SimpleAbstract() { };
};

// ... imagine 16 more SimpleDerivedX structs

class BaseParent {	 // Declares the interface for the dispatcher
public:
  // Declare overloads for each kind of a node to dispatch
  virtual void DoSimpleThing() = 0;
};


template <typename SIMPLE_T>
struct BaseInterface : BaseParent {
  static_assert(std::is_base_of<SimpleAbstract, SIMPLE_T>::value, "Must have a valid simple type.");
  SIMPLE_T thing;

  BaseInterface() = default;
  BaseInterface(SIMPLE_T& specific_thing) : thing(specific_thing) {}
  virtual void DoSimpleThing() { std::cout << "Do simple thing" << std::endl;  }
  virtual void DoOtherThing() = 0; // pure virtual
};

template <typename SIMPLE_T>
struct Derived1 : public BaseInterface<SIMPLE_T> {
  Derived1() = default;
  Derived1(SIMPLE_T& specific_thing) : BaseInterface<SIMPLE_T>(specific_thing) {}
  void DoOtherThing() override { std::cout << "Do other thing" << std::endl; };
};

// template <typename SIMPLE_T> class Derived1<SIMPLE_T>::Derived1<SimpleDerived1>{};

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  SimpleDerived1 inst_thing;
  Derived1<SimpleDerived1> y = Derived1<SimpleDerived1>(inst_thing);
  y.DoSimpleThing();
  return(x);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
