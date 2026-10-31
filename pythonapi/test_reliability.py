"""Small installed-binding regressions; no research fixtures or external solvers."""
from fractions import Fraction

import fflow as ff
import pytest


def test_numeric_strict_kernel_and_relearning():
    graph, _ = ff.NewGraphWithInput(0)
    try:
        node = ff.AlgNumericSparseLSolve(graph, 2, [[0, 1], [0, 1]], ['1', '2', '2', '4'])
        ff.LSolveOnlyHomogeneous(graph, node)
        ff.LSolveSparseOutput(graph, node, False)
        ff.LSolveValidateFullSystem(graph, node, True)
        ff.SetOutputNode(graph, node)
        ff.LearnEx(graph, prime_no=0, n_singular=4, n_learning=3)
        assert list(map(Fraction, ff.ReconstructNumeric(graph, max_primes=3))) == [Fraction(-2)]
        ff.ResetLSolveLearning(graph)
        ff.LearnEx(graph, prime_no=1, n_singular=4, n_learning=3)
        assert list(map(Fraction, ff.ReconstructNumeric(graph, start_mod=1, max_primes=3))) == [Fraction(-2)]
    finally:
        ff.DeleteGraph(graph)


def test_shift_reset_preserves_explicit_shift_and_can_restore_default():
    graph, _ = ff.NewGraphWithInput(2)
    try:
        default = ff.GetReconstructionShift(graph)
        assert len(default) == 2
        ff.SetReconstructionShift(graph, [1, 2])
        assert ff.GetReconstructionShift(graph) == [1, 2]
        ff.ResetReconstruction(graph)
        assert ff.GetReconstructionShift(graph) == [1, 2]
        ff.SetReconstructionShift(graph, None)
        assert ff.GetReconstructionShift(graph) == default
    finally:
        ff.DeleteGraph(graph)


@pytest.mark.parametrize('shift', [[True], [-1], [1, 2], [2**64 - 1]])
def test_invalid_shift_is_rejected_without_changing_previous_state(shift):
    graph, _ = ff.NewGraphWithInput(1)
    try:
        ff.SetReconstructionShift(graph, [3])
        with pytest.raises((ValueError, ff.FFlowError)):
            ff.SetReconstructionShift(graph, shift)
        assert ff.GetReconstructionShift(graph) == [3]
    finally:
        ff.DeleteGraph(graph)
