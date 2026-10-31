#include <fflow/capi.h>
#include <fflow/mp_functions.hh>
#include <fflow/primes.hh>
#include <fflow/mp_gcd.hh>
#include <cstdio>
#include <string>
#include <vector>

#define CHECK(condition) do { if (!(condition)) { \
  std::fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
  return false; } } while (0)

using namespace fflow;

namespace {

struct GraphGuard {
  FFGraph id;
  explicit GraphGuard(FFGraph value) : id(value) {}
  ~GraphGuard() { if (!ffIsError(id)) ffDeleteGraph(id); }
};

struct RatFunGuard {
  FFRatFunList * ptr;
  explicit RatFunGuard(FFRatFunList * value) : ptr(value) {}
  ~RatFunGuard() { if (ptr) ffFreeRatFun(ptr); }
};

FFRatFunList * constants(const std::vector<std::string> & values)
{
  std::vector<unsigned> sizes(values.size(), 1);
  std::vector<uint16_t> exponents(2*values.size(), 0);
  std::vector<const char *> coefficients;
  for (const auto & value : values) {
    coefficients.push_back(value.c_str());
    coefficients.push_back("1");
  }
  return ffNewRatFunList(1, values.size(), sizes.data(), sizes.data(),
                         coefficients.data(), exponents.data());
}

bool expect_value(FFGraph graph, const FFUInt * x, unsigned prime, UInt expected)
{
  FFUInt * result = ffEvaluateGraph(graph, x, prime);
  CHECK(result != nullptr);
  const UInt actual = result[0];
  ffFreeMemoryU64(result);
  CHECK(actual == expected);
  return true;
}

bool expect_failure(FFGraph graph, const FFUInt * x, unsigned prime)
{
  FFUInt * result = ffEvaluateGraph(graph, x, prime);
  if (result) ffFreeMemoryU64(result);
  CHECK(result == nullptr);
  return true;
}

bool polynomial_and_rational_evaluation()
{
  const UInt p = BIG_UINT_PRIMES[0];
  const Mod bad(p), good(BIG_UINT_PRIMES[1]);
  const std::string coefficient = "1/" + std::to_string(p);
  const MPRational scalar(1);
  const unsigned position = 0;
  CHECK(horner_map_mprational(&scalar, good, 1, &position, nullptr) == FAILED);
  CHECK(horner_map_mprational(nullptr, good, 0, nullptr, nullptr) == SUCCESS);
  MPReconstructedPoly polynomial(1);
  polynomial.unserialize("1 1 " + coefficient + "* 1 ");
  // Copy an ordinary serialized polynomial, whose CRT-residue vector is empty.
  MPReconstructedPoly copied(polynomial);
  const UInt x[1] = {2};
  CHECK(copied.eval(x, bad) == FAILED);
  CHECK(polynomial.eval(x, bad) == FAILED);
  CHECK(polynomial.eval(x, bad) == FAILED);
  const UInt value = polynomial.eval(x, good);
  CHECK(value != FAILED);
  CHECK(mul_mod(value, p % good.n(), good) == 2);
  CHECK(copied.eval(x, good) == value);
  CHECK(polynomial.eval(x, bad) == FAILED);
  CHECK(polynomial.monomialsMod(bad) == FAILED);
  CHECK(!polynomial.toHornerPtr(bad));
  CHECK(polynomial.toHornerPtr(good));

  MPReconstructedPoly one(1);
  one.to_one();
  MPReconstructedRatFun function(std::move(polynomial), std::move(one));
  CHECK(function.eval(x, bad) == FAILED);
  CHECK(function.eval(x, good) == value);
  HornerRatFunPtr converted;
  CHECK(function.toHornerRatFunPtr(converted, bad) == FAILED);
  CHECK(!converted.num() && !converted.den());
  CHECK(converted.required_workspace() == 0);
  UInt workspace[4], precomputed[1];
  CHECK(converted.eval(1, workspace, x, precomputed, bad) == FAILED);
  CHECK(function.toHornerRatFunPtr(converted, good) == SUCCESS);

  MPReconstructedRatFun zero(1);
  zero.to_zero();
  HornerRatFunPtr zero_horner;
  CHECK(zero.toHornerRatFunPtr(zero_horner, good) == SUCCESS);
  std::vector<UInt> zero_workspace(zero_horner.required_workspace());
  precomp_array_mul_shoup(x, 1, good, precomputed);
  CHECK(zero_horner.eval(1, zero_workspace.data(), x, precomputed, good) == 0);
  // A null Horner numerator is also the established representation of zero.
  zero_horner.num_ptr().reset();
  CHECK(zero_horner.eval(1, zero_workspace.data(), x, precomputed, good) == 0);

  MPReconstructedRatFun pole(1);
  const std::string vars[1] = {"x"};
  const std::string expression = "1/(x-3)";
  CHECK(pole.parse(vars, 1, expression.data(), expression.data()+expression.size()) == SUCCESS);
  MPReconstructedRatFun copied_pole(pole);
  const UInt singular[1] = {3}, regular[1] = {4};
  CHECK(pole.eval(singular, bad) == FAILED);
  CHECK(pole.eval(regular, bad) == 1);
  CHECK(copied_pole.eval(singular, bad) == FAILED);
  CHECK(copied_pole.eval(regular, bad) == 1);

  URationalFunction univariate(0, 1);
  univariate.num(0) = 1;
  univariate.den(0) = neg_mod(3, bad);
  univariate.den(1) = 1;
  CHECK(univariate.eval(3, bad) == FAILED);
  CHECK(univariate.eval(4, bad) == 1);
  SparseRationalFunction sparse(1);
  sparse.fromUnivariate(univariate);
  CHECK(sparse.eval(singular, bad) == FAILED);
  CHECK(sparse.eval(regular, bad) == 1);
  return true;
}

bool mapped_function_and_cache()
{
  const UInt p = BIG_UINT_PRIMES[0];
  const std::string coefficient = "1/" + std::to_string(p);
  RatFunGuard rf(constants({coefficient}));
  CHECK(rf.ptr != nullptr);
  FFNode input;
  GraphGuard graph(ffNewGraphWithInput(1, &input));
  const FFNode node = ffAlgRatFunEval(graph.id, input, rf.ptr);
  CHECK(!ffIsError(node));
  CHECK(ffSetOutputNode(graph.id, node) == FF_SUCCESS);
  const FFUInt x[1] = {2};
  CHECK(expect_failure(graph.id, x, 0));
  CHECK(expect_failure(graph.id, x, 0));
  const UInt expected = mul_inv(p % BIG_UINT_PRIMES[1], Mod(BIG_UINT_PRIMES[1]));
  CHECK(expect_value(graph.id, x, 1, expected));
  CHECK(expect_failure(graph.id, x, 0));
  CHECK(expect_value(graph.id, x, 1, expected));
  FFUInt * direct = ffEvaluateRatFunList(rf.ptr, x, 0);
  if (direct) ffFreeMemoryU64(direct);
  CHECK(direct == nullptr);
  return true;
}

bool constant_and_expression_nodes()
{
  const UInt p = BIG_UINT_PRIMES[0];
  const std::string coefficient = "1/" + std::to_string(p);
  const UInt expected = mul_inv(p % BIG_UINT_PRIMES[1], Mod(BIG_UINT_PRIMES[1]));
  {
    GraphGuard graph(ffNewGraph());
    const char * values[1] = {coefficient.c_str()};
    const FFNode node = ffAlgRatNumEval(graph.id, values, 1);
    CHECK(!ffIsError(node));
    CHECK(ffSetOutputNode(graph.id, node) == FF_SUCCESS);
    CHECK(expect_failure(graph.id, nullptr, 0));
    CHECK(expect_failure(graph.id, nullptr, 0));
    CHECK(expect_value(graph.id, nullptr, 1, expected));
    CHECK(expect_failure(graph.id, nullptr, 0));
  }
  {
    FFNode input;
    GraphGuard graph(ffNewGraphWithInput(1, &input));
    const std::string expression = coefficient + "+x";
    const char * vars[1] = {"x"}, * functions[1] = {expression.c_str()};
    const FFNode node = ffAlgRatExprEval(graph.id, input, vars, 1, nullptr, functions, 1);
    CHECK(!ffIsError(node));
    CHECK(ffSetOutputNode(graph.id, node) == FF_SUCCESS);
    const FFUInt x[1] = {2};
    CHECK(expect_failure(graph.id, x, 0));
    CHECK(expect_failure(graph.id, x, 0));
    CHECK(expect_value(graph.id, x, 1, add_mod(expected, 2, Mod(BIG_UINT_PRIMES[1]))));
    CHECK(expect_failure(graph.id, x, 0));
  }
  return true;
}

bool solver_coefficient_mapping(bool analytic, bool sparse)
{
  const UInt p = BIG_UINT_PRIMES[0];
  const std::string coefficient = "1/" + std::to_string(p);
  FFNode input = 0;
  GraphGuard graph(analytic ? ffNewGraphWithInput(1, &input) : ffNewGraph());
  RatFunGuard rf(analytic ? constants({coefficient, "1"}) : nullptr);
  FFNode node;
  const char * values[2] = {coefficient.c_str(), "1"};
  const unsigned nonzero_count[1] = {2}, columns[2] = {0, 1};
  const size_t indices[2] = {0, 1};
  if (analytic) {
    CHECK(rf.ptr != nullptr);
    node = sparse ? ffAlgAnalyticSparseLSolve(graph.id, input, 1, 1,
                    nonzero_count, columns, indices, rf.ptr, nullptr, 0)
                  : ffAlgAnalyticDenseLSolve(graph.id, input, 1, 1, rf.ptr, nullptr, 0);
  } else {
    node = sparse ? ffAlgNumericSparseLSolve(graph.id, 1, 1,
                    nonzero_count, columns, indices, values, 2, nullptr, 0)
                  : ffAlgNumericDenseLSolve(graph.id, 1, 1, values, nullptr, 0);
  }
  CHECK(!ffIsError(node));
  CHECK(ffSetOutputNode(graph.id, node) == FF_SUCCESS);
  FFLearnOptions learn = {};
  learn.prime_no = 1;
  CHECK(ffLearnEx(graph.id, learn) == FF_SUCCESS);
  const FFUInt x[1] = {2};
  const FFUInt * point = analytic ? x : nullptr;
  CHECK(expect_value(graph.id, point, 1, p % BIG_UINT_PRIMES[1]));
  CHECK(expect_failure(graph.id, point, 0));
  CHECK(expect_failure(graph.id, point, 0));
  CHECK(expect_value(graph.id, point, 1, p % BIG_UINT_PRIMES[1]));
  return true;
}

bool zero_parameter_matrices()
{
  FFNode input;
  GraphGuard g(ffNewGraphWithInput(0,&input));
  const unsigned sizes[] = {1,1,1,1,1,1};
  const char * coeffs[] = {"1","1","2","1","0","1","2","1","4","1","0","1"};
  FFRatFunList * rf=ffNewRatFunList(0,6,sizes,sizes,coeffs,nullptr);
  CHECK(rf);
  CHECK(ffIsError(ffAlgRatFunEval(g.id,input,rf)));
  CHECK(ffIsError(ffAlgAnalyticDenseLSolve(g.id,input,2,2,rf,nullptr,0)));
  const unsigned counts[] = {2,2}, cols[] = {0,1,0,1};
  const size_t indices[] = {0,1,3,4};
  CHECK(ffIsError(ffAlgAnalyticSparseLSolve(g.id,input,2,2,counts,cols,
                                          indices,rf,nullptr,0)));
  FFUInt * direct=ffEvaluateRatFunList(rf,nullptr,0);
  CHECK(direct);
  const FFUInt expected[] = {1,2,0,2,4,0};
  for (unsigned i=0; i<6; ++i) CHECK(direct[i]==expected[i]);
  ffFreeMemoryU64(direct);
  ffFreeRatFun(rf);
  const char * values[] = {"1","2","2","4"};
  const size_t index[] = {0,1,2,3};
  FFNode node=ffAlgNumericSparseLSolve(g.id,2,2,counts,cols,index,values,4,nullptr,0);
  CHECK(!ffIsError(node));
  CHECK(ffLSolveOnlyHomogeneous(g.id,node)==FF_SUCCESS);
  CHECK(ffSetOutputNode(g.id,node)==FF_SUCCESS);
  CHECK(ffLSolveValidateFullSystem(g.id,node,true)==FF_SUCCESS);
  CHECK(ffLearn(g.id)==FF_SUCCESS);
  FFUInt * output=ffEvaluateGraph(g.id,nullptr,0);
  CHECK(output && output[0]==ffPrimeNo(0)-2);
  ffFreeMemoryU64(output);
  FFRecOptions options = {}; options.max_primes=3; options.n_threads=1;
  char ** result=ffReconstructNumeric(g.id,options);
  CHECK(result && result[0]);
  mpq_t rational;
  mpq_init(rational);
  CHECK(mpq_set_str(rational,result[0],10)==0);
  mpq_canonicalize(rational);
  CHECK(mpq_cmp_si(rational,-2,1)==0);
  mpq_clear(rational);
  ffFreeCStrArray(result);
  return true;
}

} // namespace

int main()
{
  std::puts("Checking polynomial and rational evaluation...");
  if (!polynomial_and_rational_evaluation()) return 1;
  std::puts("Checking mapped functions...");
  if (!mapped_function_and_cache()) return 1;
  std::puts("Checking constants and expressions...");
  if (!constant_and_expression_nodes()) return 1;
  for (bool analytic : {false, true})
    for (bool sparse : {false, true})
      {
        std::printf("Checking %s %s solver...\n", analytic ? "analytic" : "numeric", sparse ? "sparse" : "dense");
        if (!solver_coefficient_mapping(analytic, sparse)) return 1;
      }
  if (!zero_parameter_matrices()) return 1;
  std::puts("Bad-modulus and rational-pole regression tests passed.");
  return 0;
}
