#include <fflow/capi.h>
#include <fflow/numeric_solver.hh>
#include <unistd.h>
#include <gmp.h>
#include <cstdio>
#include <string>
#include <vector>

#define CHECK(condition) do { if (!(condition)) { \
  std::fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
  return false; } } while (0)

namespace {
#include "solver_test_helpers.hh"

bool strict_bad_prime_and_relearning(bool sparse, bool strict = true)
{
  FFNode input;
  GraphGuard g(ffNewGraphWithInput(1,&input));
  const auto p0 = ffPrimeNo(0);
  FFNode node = solver(g.graph,input,sparse,"1+"+std::to_string(p0)+"*x");
  CHECK(!ffIsError(node));
  CHECK(ffLSolveOnlyHomogeneous(g.graph,node)==FF_SUCCESS);
  if (sparse) CHECK(ffLSolveSparseOutput(g.graph,node,false)==FF_SUCCESS);
  CHECK(ffSetOutputNode(g.graph,node)==FF_SUCCESS);
  if (strict) CHECK(ffLSolveValidateFullSystem(g.graph,node,true)==FF_SUCCESS);
  FFLearnOptions initial = {0};
  CHECK(ffLearnEx(g.graph,initial)==FF_SUCCESS);
  CHECK(ffLSolveNIndepVars(g.graph,node,0)==1);
  if (sparse && strict) {
    CHECK(ffLSolveMarkAndSweepEqs(g.graph,node)==FF_SUCCESS);
    CHECK(ffLSolveDeleteUnneededEqs(g.graph,node)==FF_ERROR);
  }
  const FFUInt point[] = {2};
  FFUInt * stale = ffEvaluateGraph(g.graph,point,1);
  if (stale) ffFreeMemoryU64(stale);
  if (strict) CHECK(stale==nullptr);
  else CHECK(stale!=nullptr);
  CHECK(ffResetLSolveLearning(g.graph)==FF_SUCCESS);
  FFUInt * unlearned = ffEvaluateGraph(g.graph,point,1);
  if (unlearned) ffFreeMemoryU64(unlearned);
  CHECK(unlearned==nullptr);
  FFGraphLearnOptions retry = {1,4,3};
  CHECK(ffLearnWithOptions(g.graph,retry)==FF_SUCCESS);
  CHECK(ffLSolveNIndepVars(g.graph,node,0)==0);
  CHECK(ffGraphNParsOut(g.graph)==0);
  FFUInt * result = ffEvaluateGraph(g.graph,point,1);
  CHECK(result!=nullptr);
  ffFreeMemoryU64(result);
  return true;
}

bool pruned_source_cannot_be_relearned()
{
  FFNode input;
  GraphGuard graph(ffNewGraphWithInput(1,&input));
  const FFNode node=solver(graph.graph,input,true,"2");
  CHECK(!ffIsError(node));
  CHECK(ffSetOutputNode(graph.graph,node)==FF_SUCCESS);
  CHECK(ffLearn(graph.graph)==FF_SUCCESS);
  CHECK(ffLSolveDeleteUnneededEqs(graph.graph,node)==FF_SUCCESS);
  CHECK(ffResetLSolveLearning(graph.graph)==FF_ERROR);
  return true;
}

bool ex_source_pruning_and_legacy_budgets()
{
  FFNode input;
  for (bool pruned : {false,true}) {
    GraphGuard ex(ffNewGraphWithInput(1,&input));
  const char * values[] = {"1",pruned ? "-1" : "1"};
  FFNode constants=ffAlgRatNumEval(ex.graph,values,2);
  CHECK(!ffIsError(constants));
  const char * vars[] = {"x"}, * functions[] = {"1"};
  FFRatFunList * rf=ffParseRatFun(vars,1,functions,1);
  CHECK(rf);
  FFNode inputs[] = {input,constants};
  const unsigned counts[] = {1}, cols[] = {0}, nweights[] = {2};
  const unsigned weights[] = {0,0,0,1};
  const size_t indices[] = {0,0};
  FFNode node=ffAlgAnalyticSparseLSolveEx(ex.graph,inputs,2,1,1,counts,cols,
                                         nweights,weights,indices,rf,nullptr,0);
  ffFreeRatFun(rf);
  CHECK(!ffIsError(node));
  CHECK(ffLSolveValidateFullSystem(ex.graph,node,true)==FF_ERROR);
  CHECK(ffLSolveValidateFullSystem(ex.graph,node,false)==FF_SUCCESS);
  CHECK(ffSetOutputNode(ex.graph,node)==FF_SUCCESS);
  CHECK(ffLearn(ex.graph)==FF_SUCCESS);
  // The weighted coefficient 1*1 + 1*(-1) was destructively deleted.
  CHECK(ffResetLSolveLearning(ex.graph)==(pruned ? FF_ERROR : FF_SUCCESS));
  if (!pruned) CHECK(ffLearn(ex.graph)==FF_SUCCESS);

  }

  GraphGuard legacy(ffNewGraphWithInput(1,&input));
  FFNode node=solver(legacy.graph,input,true,"2");
  CHECK(!ffIsError(node) && ffSetOutputNode(legacy.graph,node)==FF_SUCCESS);
  CHECK(ffLSolveValidateFullSystem(legacy.graph,node,false)==FF_SUCCESS);
  FFGraphLearnOptions insufficient = {0,0,1};
  CHECK(ffLearnWithOptions(legacy.graph,insufficient)==FF_ERROR);
  FFLearnOptions defaults = {0};
  CHECK(ffLearnEx(legacy.graph,defaults)==FF_SUCCESS);
  GraphGuard retained_prime(ffNewGraphWithInput(1,&input));
  node=solver(retained_prime.graph,input,true,"1/"+std::to_string(ffPrimeNo(0)));
  CHECK(!ffIsError(node) && ffSetOutputNode(retained_prime.graph,node)==FF_SUCCESS);
  FFGraphLearnOptions invalid_at_prime1 = {1,0,1};
  CHECK(ffLearnWithOptions(retained_prime.graph,invalid_at_prime1)==FF_ERROR);
  CHECK(ffLearn(retained_prime.graph)==FF_SUCCESS);
  return true;
}

bool exact_reciprocal(const FFRatFunList * rf)
{
  CHECK(ffRatFunListSize(rf)==1 && ffRatFunListNVars(rf)==2);
  CHECK(ffRatFunNumNTerms(rf,0)==1 && ffRatFunDenNTerms(rf,0)==2);
  char ** num = ffRatFunNumCoeffs(rf,0);
  char ** den = ffRatFunDenCoeffs(rf,0);
  uint16_t * ne = ffRatFunNumExponents(rf,0);
  uint16_t * de = ffRatFunDenExponents(rf,0);
  CHECK(num && den && ne && de && ne[0]==0 && ne[1]==0);
  mpq_t constant, value, expected;
  mpq_inits(constant,value,expected,nullptr);
  CHECK(mpq_set_str(constant,num[0],10)==0);
  mpq_canonicalize(constant);
  CHECK(mpq_sgn(constant)!=0);
  bool has_x=false, has_constant=false;
  for (unsigned i=0; i<2; ++i) {
    CHECK(mpq_set_str(value,den[i],10)==0);
    mpq_canonicalize(value);
    if (de[2*i]==1 && de[2*i+1]==0) {
      CHECK(mpq_equal(value,constant));
      has_x=true;
    } else {
      CHECK(de[2*i]==0 && de[2*i+1]==0);
      mpq_set_si(expected,-2,1);
      mpq_mul(expected,expected,constant);
      CHECK(mpq_equal(value,expected));
      has_constant=true;
    }
  }
  mpq_clears(constant,value,expected,nullptr);
  ffFreeCStrArray(num); ffFreeCStrArray(den);
  ffFreeMemoryU16(ne); ffFreeMemoryU16(de);
  CHECK(has_x && has_constant);
  return true;
}

bool shifts_and_reconstruction_reset()
{
  FFNode input;
  GraphGuard g(ffNewGraphWithInput(2,&input));
  const char * vars[] = {"x","y"}, * functions[] = {"1/(x-2)"};
  FFRatFunList * rf = ffParseRatFun(vars,2,functions,1);
  CHECK(rf!=nullptr);
  FFNode node = ffAlgRatFunEval(g.graph,input,rf);
  ffFreeRatFun(rf);
  CHECK(!ffIsError(node) && ffSetOutputNode(g.graph,node)==FF_SUCCESS);
  const FFUInt bad_shift[] = {2,0}, good_shift[] = {3,0};
  CHECK(ffSetReconstructionShift(g.graph,bad_shift,1)==FF_ERROR);
  const FFUInt invalid[] = {ffPrimeNo(ffNAvailablePrimes()-1),0};
  CHECK(ffSetReconstructionShift(g.graph,invalid,2)==FF_ERROR);
  CHECK(ffSetReconstructionShift(g.graph,bad_shift,2)==FF_SUCCESS);
  FFRecOptions options = {}; options.max_primes=2; options.max_deg=4; options.n_threads=1;
  FFRatFunList * result = nullptr;
  CHECK(ffReconstructFunction(g.graph,options,&result)!=FF_SUCCESS);
  CHECK(result==nullptr);
  CHECK(ffSetReconstructionShift(g.graph,good_shift,2)==FF_SUCCESS);
  CHECK(ffReconstructFunction(g.graph,options,&result)==FF_SUCCESS && result);
  CHECK(exact_reciprocal(result));
  ffFreeRatFun(result);
  CHECK(ffResetReconstruction(g.graph)==FF_SUCCESS);
  unsigned size=0;
  FFUInt * shift = ffGetReconstructionShift(g.graph,&size);
  CHECK(shift && size==2 && shift[0]==3 && shift[1]==0);
  ffFreeMemoryU64(shift);
  result=nullptr;
  CHECK(ffReconstructFunction(g.graph,options,&result)==FF_SUCCESS && result);
  CHECK(exact_reciprocal(result)); ffFreeRatFun(result);
  CHECK(ffSetReconstructionShift(g.graph,nullptr,0)==FF_SUCCESS);
  shift=ffGetReconstructionShift(g.graph,&size);
  CHECK(shift && size==2 && shift[0]!=3);
  ffFreeMemoryU64(shift);
  GraphGuard resized(ffNewGraphWithInput(2,&input));
  CHECK(ffSetReconstructionShift(resized.graph,good_shift,2)==FF_SUCCESS);
  shift=ffGetReconstructionShift(resized.graph,&size);
  CHECK(shift && size==2 && shift[0]==3);
  ffFreeMemoryU64(shift);
  CHECK(!ffIsError(ffSetGraphInput(resized.graph,3)));
  shift=ffGetReconstructionShift(resized.graph,&size);
  CHECK(shift && size==3 && shift[0]!=3);
  ffFreeMemoryU64(shift);
  CHECK(ffResetReconstruction(~unsigned(0))==FF_ERROR);
  FFLearnOptions invalid_prime = {ffNAvailablePrimes()};
  CHECK(ffLearnEx(g.graph,invalid_prime)==FF_ERROR);
  return true;
}
}

int main()
{
  for (bool sparse : {false,true}) for (bool strict : {false,true})
    if (!strict_bad_prime_and_relearning(sparse,strict)) return 1;
  if (!pruned_source_cannot_be_relearned()) return 1;
  if (!ex_source_pruning_and_legacy_budgets()) return 1;
  if (!shifts_and_reconstruction_reset()) return 1;
  std::puts("reconstruction-control regression groups passed");
  return 0;
}
