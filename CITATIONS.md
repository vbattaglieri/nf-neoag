# Citations & attribution

## Original method

The neoantigen-prediction **method** and its original reference implementation
(a set of Python 2 / bash scripts, not distributed with this repository) were
developed by **Giuseppe Rospo**, a former member of our group.
**All credit for the original pipeline and methodology goes to him.**

This repository is an independent **Python 3 reimplementation** of that method,
wrapped as a **Nextflow pipeline**. The goal is portability and reproducibility;
the scientific logic is preserved, with two small documented deviations (see the
README) that correct legacy behaviour on current reference files.

### Please cite

- Germano G, *et al.* **Inactivation of DNA repair triggers neoantigen
  generation and impairs tumour growth.** *Nature* (2017).

<!-- Add further publications describing / using this pipeline below: -->
<!-- - Author X, et al. Title. Journal (Year). -->

## Tools

- **netMHC 4.0** — Andreatta M & Nielsen M, *Bioinformatics* (2016);
  Nielsen M et al., *Protein Sci* (2003). (License-restricted; not distributed here.)
- **Nextflow** — Di Tommaso P et al., *Nat Biotechnol* (2017).
