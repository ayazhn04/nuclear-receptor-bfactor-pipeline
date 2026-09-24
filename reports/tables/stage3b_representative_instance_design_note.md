# Design Note: Representative-Instance Selection Granularity (Future Stage)

**Status: descriptive design note only. No representative instance has been
selected. No final completeness threshold has been chosen. This note exists
to record a structural constraint discovered during Stage 3B.2, to be
honored whenever a future stage performs representative-instance selection.**

## The constraint

Stage 3B's coordinate-completeness audit is keyed at the level of

    (uniprot_id, pdb_id, polymer_entity_id, instance_id, model_num)

A single PDB entry frequently contains multiple polymer entities (e.g. a
receptor construct plus a separate co-crystallized peptide or a different
receptor's LBD in a heterodimer), and a single polymer entity is frequently
present as multiple structural instances (chains) within one asymmetric
unit. Stage 3B.2's instance-heterogeneity audit
(`data/manifests/stage3b_instance_coverage_heterogeneity.csv`) found 1,012
receptor × PDB × polymer-entity candidates with more than one instance, of
which 111 show a coverage range of 0.05 or more between their best- and
worst-observed instance (max observed range: 0.197).

This means **representative-instance selection cannot be performed once per
PDB entry.** Two different chains of the same PDB entry can belong to two
different receptors (in a heterodimer structure) or can be the same receptor
construct with materially different observed coverage (as the heterogeneity
audit shows). Selecting "one chain per PDB" would silently conflate
different receptors or discard a more-complete instance in favor of an
arbitrary one.

## Required granularity for future representative-instance selection

Any future stage that selects a representative instance MUST make that
selection independently within each

    (project receptor, PDB entry, receptor's polymer entity)

triple — i.e., once per row of `stage3b_instance_coverage_heterogeneity.csv`
/ `stage3b_any_all_threshold_comparison.csv`, never once per PDB ID alone.

## Required ordering of future decisions

The following decisions are logically and causally prior to any
representative-instance selection, and must each be made and documented
explicitly before it, not folded into it implicitly:

1. Final coordinate-completeness threshold (not yet selected; Stage 3B.2
   provides only a descriptive sensitivity comparison across
   100/99/95/90/80%, including the "any instance passes" vs. "all instances
   pass" distinction in
   `data/manifests/stage3b_any_all_threshold_comparison.csv`).
2. Representative-instance selection (per the granularity above).
3. Normalized B-factor computation.
4. Clustering / structural comparison across receptors.
5. Ligand-bound vs. apo (ligand-state) effect analysis.
6. Any statistical significance testing.

Each of these depends on the ones before it; performing them out of order
(e.g. normalizing B-factors before a representative instance is fixed, or
selecting a representative instance before a completeness threshold is
fixed) would make the analysis irreproducible and its inputs ambiguous.
Stage 3B and Stage 3B.2 deliberately stop before step 1.
