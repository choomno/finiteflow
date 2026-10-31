#include <fflow/mp_functions.hh>
#include <fflow/primes.hh>
#include <cstdio>
#include <string>

#define CHECK(condition) do { if (!(condition)) { \
  std::fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
  return 1; } } while (0)

int main()
{
  using namespace fflow;
  const Mod mod(BIG_UINT_PRIMES[0]);
  const UInt point[] = {2};
  MPReconstructedPoly literal(1);
  literal.unserialize("1 1 1/2* 1 ");
  // Deserialization populates rational coefficients, not CRT residues.
  MPReconstructedPoly copied(literal);
  CHECK(literal.eval(point,mod)==1);
  CHECK(copied.eval(point,mod)==1);
  MPReconstructedPoly recopied(1);
  recopied.copy(copied);
  CHECK(recopied.eval(point,mod)==1);

  const std::string variables[] = {"x"};
  const std::string expression="(x+1)/(x-3)";
  MPReconstructedRatFun parsed(1);
  CHECK(parsed.parse(variables,1,expression.data(),expression.data()+expression.size())==SUCCESS);
  MPReconstructedRatFun parsed_copy(parsed);
  const UInt regular[] = {4};
  CHECK(parsed.eval(regular,mod)==5);
  CHECK(parsed_copy.eval(regular,mod)==5);
  std::puts("parsed polynomial and rational-function copies preserve exact coefficients");
  return 0;
}
