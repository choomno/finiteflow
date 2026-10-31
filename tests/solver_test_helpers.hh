#ifndef FFLOW_SOLVER_TEST_HELPERS_HH
#define FFLOW_SOLVER_TEST_HELPERS_HH

struct GraphGuard {
  FFGraph graph;
  explicit GraphGuard(FFGraph g) : graph(g) {}
  ~GraphGuard() { ffDeleteGraph(graph); }
};

FFNode solver(FFGraph graph, FFNode input, bool sparse,
              const std::string & last, const unsigned * needed = nullptr,
              unsigned nneeded = 0)
{
  const char * vars[] = {"x"};
  FFNode node;
  if (sparse) {
    const char * functions[] = {"1",last.c_str()};
    FFRatFunList * rf = ffParseRatFun(vars,1,functions,2);
    if (!rf) return ~unsigned(0);
    const size_t indexes[] = {0,0,0,1};
    FFIdxRatFunList * idx = ffMoveRatFunToIdx(rf,indexes,4);
    ffFreeRatFun(rf);
    if (!idx) return ~unsigned(0);
    const unsigned counts[] = {2,2}, cols[] = {0,1,0,1};
    node = ffAlgAnalyticSparseLSolveIdx(graph,input,2,2,counts,cols,idx,needed,nneeded);
    ffFreeIdxRatFun(idx);
  } else {
    const char * functions[] = {"1","1","0","1",last.c_str(),"0"};
    FFRatFunList * rf = ffParseRatFun(vars,1,functions,6);
    if (!rf) return ~unsigned(0);
    node = ffAlgAnalyticDenseLSolve(graph,input,2,2,rf,needed,nneeded);
    ffFreeRatFun(rf);
  }
  return node;
}


#endif
