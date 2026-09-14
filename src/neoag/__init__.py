"""neoag — neoantigen prediction from SNVs and Indels.

Modern Python 3 re-implementation of the pipeline originally written by
Giuseppe Rospo (method: Germano et al., Nature 2017). The scientific logic
(peptide window generation, frameshift translation, netMHC consensus) is
preserved; only the implementation has been modernized. See CITATIONS.md.
"""

__version__ = "0.1.0"
