"""
Docking trajectory entry point.

As of Level 2 the trajectory is produced by a *real* rigid-body docking engine
(see ``docking/engine.py``): a Vina-style scoring function searched with a
multi-start simulated-annealing optimizer. The ligand genuinely explores, is
rejected for clashes, reorients, forms interactions, and converges on the
lowest-energy pose the engine finds — and that real search path *is* the movie.

This module keeps the historical ``generate_trajectory`` signature so the rest
of the app is unaffected; it simply delegates to the engine.
"""
from __future__ import annotations

from docking import engine

# re-exported so callers that introspect phases keep working
PHASES = engine.PHASES


def generate_trajectory(record, n_frames=1000, seed=7, beats=None, exhaustiveness=8):
    """Dock ``record`` and return a viewer-ready trajectory dict (see engine.dock)."""
    return engine.dock(record, n_frames=n_frames, seed=seed, beats=beats,
                       exhaustiveness=exhaustiveness)
