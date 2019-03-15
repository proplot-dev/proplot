"""
A collection of useful variables.
"""

__all__ = ['PARROT_STATE', 'FUNNY_WALK_STEPS']


PARROT_STATE = 'dead'
"""The global state of the parrot."""

FUNNY_WALK_STEPS = [['left', 'right'],
                    ['left', 'jump', 'right', 'jump'],
                    ['swim']]
"""List of different possible walk.

Each item contains a list of steps.
"""

# A variable not in __all__ should not be propagated.
NOTHING_HAPPENS = 0

# Even if it has a docstring
REALLY_NOTHING = 1
"""Really nothing."""
