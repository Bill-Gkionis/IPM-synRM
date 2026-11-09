# FEMM 5.0 against FEMM 4.2: the report is in the project notes

The documentation of the motor project lives in the project notes only (decided 2026-11-05). The report that used to be this page is two notes there:

- **4.2 vs 5.0**,
  whether FEMM 5.0 gives the same answers as FEMM 4.2. Every study, every flag,
  the node counts, how Python drives each FEMM, the open questions, Self-check.
- **M19**,
  the steel, FEMM 5.0's "M-19 14 mil" against FEMM 4.2's "M-19 Steel".

What stays in this folder, because it is code or data: the table `femm5_vs_femm42.csv`; the scripts `compare.py`, `mesh_checks.py`, `m19_femm5.py` and
`provenance.py` (each docstring says how to run it and what numbers to expect); and the figures' originals, `out/validation/mesh/mesh_convergence.png` and
`out/validation/m19_curves.png`, with a copy of each beside the scripts that make them (`mesh_convergence.png`, `m19_curves.png`). `femm42-out`, next to this
folder, is FEMM 4.2's run, read only (its SHA-256 list is `femm42-out.sha256`).
