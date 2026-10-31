#include <fflow/capi.h>
#include <fflow/numeric_solver.hh>
#include <unistd.h>
#include <cstdio>
#include <string>
#include <vector>

#define CHECK(condition) do { if (!(condition)) { \
  std::fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
  return false; } } while (0)

namespace {
#include "solver_test_helpers.hh"

bool strict_bad_prime(bool sparse, bool strict = true)
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
  return true;
}

bool pruned_and_unsupported_requests()
{
  FFNode input;
  GraphGuard g(ffNewGraphWithInput(1,&input));
  FFNode node = solver(g.graph,input,true,"2");
  CHECK(!ffIsError(node));
  CHECK(ffSetOutputNode(g.graph,node)==FF_SUCCESS);
  CHECK(ffLearn(g.graph)==FF_SUCCESS);
  CHECK(ffLSolveDeleteUnneededEqs(g.graph,node)==FF_SUCCESS);

  GraphGuard partial(ffNewGraphWithInput(1,&input));
  const unsigned needed[] = {0};
  node = solver(partial.graph,input,true,"2",needed,1);
  CHECK(ffLSolveValidateFullSystem(partial.graph,node,true)==FF_ERROR);

  GraphGuard particular(ffNewGraphWithInput(1,&input));
  node = solver(particular.graph,input,true,"2");
  CHECK(ffLSolveOnlyNonHomogeneous(particular.graph,node)==FF_SUCCESS);
  CHECK(ffLSolveValidateFullSystem(particular.graph,node,true)==FF_ERROR);

  GraphGuard strict(ffNewGraphWithInput(1,&input));
  node = solver(strict.graph,input,true,"2");
  CHECK(ffLSolveValidateFullSystem(strict.graph,node,true)==FF_SUCCESS);
  CHECK(ffLSolveOnlyNonHomogeneous(strict.graph,node)==FF_ERROR);
  CHECK(ffLSolveSparseOutputWithMaxCol(strict.graph,node,1,false,false)==FF_ERROR);
  return true;
}

bool legacy_layout_cannot_claim_original_source()
{
  fflow::NumericSparseSolver solver;
  fflow::NumericSparseSolverData data;
  const unsigned needed[] = {0};
  solver.nparsin.resize(1); solver.nparsin[0]=0;
  solver.init(0,1,needed,1,data);
  char filename[]="/tmp/finiteflow-layout-XXXXXX";
  const int fd=mkstemp(filename);
  CHECK(fd>=0); close(fd);
  CHECK(solver.dump_info(&data,filename)==fflow::SUCCESS);
  CHECK(solver.validate_full_system()==fflow::SUCCESS);
  CHECK(solver.load_info(&data,filename)==fflow::FAILED);
  CHECK(solver.full_system_validation());
  CHECK(solver.validate_full_system(false)==fflow::SUCCESS);
  CHECK(solver.load_info(&data,filename)==fflow::SUCCESS);
  CHECK(!solver.source_available());
  CHECK(solver.restart_learning(&data)==fflow::FAILED);
  CHECK(solver.validate_full_system()==fflow::FAILED);
  CHECK(std::remove(filename)==0);
  return true;
}

bool dropped_affine_row_is_detected()
{
  for (bool sparse : {false,true}) for (bool homogeneous : {false,true})
    for (bool strict : {false,true}) for (bool sparse_output : {false,true}) {
      if (!sparse && sparse_output) continue;
      FFNode input;
      GraphGuard g(ffNewGraphWithInput(1,&input));
      const std::string last="1+"+std::to_string(ffPrimeNo(0))+"*x";
      const char * vars[] = {"x"};
      FFNode node;
      if (sparse) {
        const char * functions[] = {"1",last.c_str()};
        FFRatFunList * rf=ffParseRatFun(vars,1,functions,2);
        CHECK(rf);
        const size_t indices[] = {0,0,1,0};
        const unsigned counts[] = {2,2}, cols[] = {0,1,0,1};
        node=ffAlgAnalyticSparseLSolve(g.graph,input,2,1,counts,cols,indices,rf,nullptr,0);
        ffFreeRatFun(rf);
        CHECK(ffLSolveSparseOutput(g.graph,node,sparse_output)==FF_SUCCESS);
      } else {
        const char * functions[] = {"1","1",last.c_str(),"1"};
        FFRatFunList * rf=ffParseRatFun(vars,1,functions,4);
        CHECK(rf);
        node=ffAlgAnalyticDenseLSolve(g.graph,input,2,1,rf,nullptr,0);
        ffFreeRatFun(rf);
      }
      CHECK(!ffIsError(node));
      if (homogeneous) CHECK(ffLSolveOnlyHomogeneous(g.graph,node)==FF_SUCCESS);
      CHECK(ffSetOutputNode(g.graph,node)==FF_SUCCESS);
      if (strict) CHECK(ffLSolveValidateFullSystem(g.graph,node,true)==FF_SUCCESS);
      CHECK(ffLearn(g.graph)==FF_SUCCESS);
      const FFUInt point[] = {2};
      FFUInt * output=ffEvaluateGraph(g.graph,point,1);
      CHECK(strict ? output==nullptr : output!=nullptr);
      if (output) ffFreeMemoryU64(output);
    }
  return true;
}

bool weighted_strict_mode_is_unsupported()
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

  }

  return true;
}

bool strict_valid_output_layouts()
{
  for (bool sparse : {false,true}) {
    for (bool homogeneous : {false,true}) {
      FFNode input;
      GraphGuard g(ffNewGraphWithInput(1,&input));
      const char * vars[] = {"x"};
      const char * dense[] = {"1","x","0","1"};
      FFRatFunList * rf=ffParseRatFun(vars,1,dense,4);
      CHECK(rf);
      FFNode node;
      if (sparse) {
        // x0 + x*x1 = 1; x2 is a free zero column.
        const size_t indices[] = {0,1,3};
        FFIdxRatFunList * idx=ffMoveRatFunToIdx(rf,indices,3);
        CHECK(idx);
        const unsigned counts[] = {3}, cols[] = {0,1,3};
        node=ffAlgAnalyticSparseLSolveIdx(g.graph,input,1,3,counts,cols,idx,nullptr,0);
        ffFreeIdxRatFun(idx);
        CHECK(ffLSolveSparseOutput(g.graph,node,true)==FF_SUCCESS);
      } else {
        node=ffAlgAnalyticDenseLSolve(g.graph,input,1,3,rf,nullptr,0);
      }
      ffFreeRatFun(rf);
      CHECK(!ffIsError(node));
      if (homogeneous) CHECK(ffLSolveOnlyHomogeneous(g.graph,node)==FF_SUCCESS);
      CHECK(ffSetOutputNode(g.graph,node)==FF_SUCCESS);
      CHECK(ffLSolveValidateFullSystem(g.graph,node,true)==FF_SUCCESS);
      CHECK(ffLearn(g.graph)==FF_SUCCESS);
      const FFUInt point[] = {2};
      FFUInt * result=ffEvaluateGraph(g.graph,point,0);
      CHECK(result);
      CHECK(result[0]==ffPrimeNo(0)-2);
      if (!homogeneous) {
        const unsigned offset=sparse?1:2;
        CHECK(result[offset]==1);
      }
      ffFreeMemoryU64(result);
    }
  }
  return true;
}
}

int main()
{
  for (bool sparse : {false,true}) for (bool strict : {false,true})
    if (!strict_bad_prime(sparse,strict)) return 1;
  if (!pruned_and_unsupported_requests()) return 1;
  if (!legacy_layout_cannot_claim_original_source()) return 1;
  if (!dropped_affine_row_is_detected()) return 1;
  if (!weighted_strict_mode_is_unsupported()) return 1;
  if (!strict_valid_output_layouts()) return 1;
  std::puts("strict original-system regression groups passed");
  return 0;
}
