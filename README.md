# Protein CAD

A CAD-style editor for drawing protein backbones to use as templates for inverse folding
(ProteinMPNN and similar). Extrude α-helices from a sketch plane, lay β-sheets on it, join chain ends
with loops or draw a loop freehand, select and move or rotate regions with a gizmo, mutate residues,
edit numbering and torsions, then export PDB or mmCIF.

Runs entirely in the browser (WebGL2). Nothing is uploaded and there is nothing to install.

## Start

Open `index.html`, or the single-file build `dist/protein_cad.html`, in Chrome, Edge or Firefox.
To rebuild the single file after editing: `python tools/build_single_file.py`.

The model is saved in the browser as you work and comes back when you reopen the page; **New**
starts over. Export to keep a copy outside the browser.

## Workflow

1. **Sketch plane.** The grid is the sketch plane, XY with Z up to start with. The Sketch plane
   panel (shown with the Helix, Sheet and Draw tools) sets it to XY / YZ / ZX, to face the camera,
   to a best fit through the selected residues (for example a sheet), or perpendicular to the
   selection's long axis at its C-terminal end (to stack or continue from a helix). Use **Offset**
   to slide it along its normal and **Snap** to set the grid step.
2. **Helix (H).** Click the plane for a helix of the default length, or press and drag along the
   dashed normal: the drag sets the length (1.54 Å per residue) and whether the helix goes above or
   below the plane. *On plane* chooses whether the N- or the C-terminus sits on the plane, and
   *Phase* rotates the helix about its axis. The 10 Å footprint circle is a typical
   helix–helix packing distance, so circles that touch give a sensible bundle.
3. **Sheet (B).** Drag across the plane to draw the first strand (N→C along the stroke). The other
   strands are laid beside it 4.8 Å apart, antiparallel or parallel, with optional twist and
   register shift. Each strand is its own chain.
4. **Connect (C).** Either select the last residue of one chain and the first residue of another,
   or select whole chains in the order they should be joined (Ctrl-click them in the chain list),
   then press Connect. Selecting two residues of the *same* chain rebuilds the stretch between them.
   Loop length is automatic (the shortest loop that closes without clashes) unless you set it.
5. **Draw loop (D).** Draw a path. Start on a C-terminus (▲) to grow from it and/or end on an
   N-terminus (●) to join it; without either, the loop becomes a new chain. *Spacing* sets Å of
   drawing per residue (3.8 is fully extended).
6. **Refine (R).** Select a region and press Refine to regularise it: bond lengths, angles and
   peptide units go back to ideal values, overlapping atoms are pushed apart, and small chain breaks
   are closed. Use it after dragging residues around, after closing a gap by hand, or on anything
   imported that you have edited. Residues outside the selection stay fixed and act as context, so
   **select a margin of a few residues either side of the damage**. *Check geometry* reports what is
   wrong without changing anything.
7. **Edit.** Select residues by clicking (Ctrl toggles, Shift extends along the chain, double-click
   takes the whole helix / strand / loop, Shift-drag draws a box), or in the sequence strip. Move
   and rotate with the gizmo (Ctrl snaps; **L** switches between world and local axes, where local
   X′ is the long axis of a helix) or with numbers in the inspector. A single residue can be mutated,
   renumbered, moved, and have φ / ψ / ω set (choosing which side of the chain moves).
8. **Templates.** *Backbone*, *Backbone → GLY*, *Poly-ALA* and *Poly-GLY* in Chain tools (whole
   chains) or in the inspector (selection). Export → PDB or mmCIF, optionally backbone atoms only.

Mouse: right-drag orbits in any tool, left-drag on empty space orbits in the Select tool,
middle-drag or Shift + right-drag pans, wheel zooms. Press **?** for the full key list.

## How the geometry is built

**Ideal backbone.** Residues are placed by NeRF (Parsons et al. 2005) from Engh & Huber bond
lengths and angles: N–CA 1.458, CA–C 1.525, C–N 1.329, C=O 1.231 Å; N–CA–C 111.2°, CA–C–N 116.2°,
C–N–CA 121.7°; ω = 180°. Helices use φ/ψ = −57.8/−47.0 (α), −49/−26 (3₁₀) or −57/−70 (π); strands
−139/135 (antiparallel) or −119/113 (parallel).

**Placing elements.** The axis of a helix or strand is the screw axis of the rigid motion that maps
one residue's frame onto the next. That is exact for a regular segment, unlike a PCA of the CA atoms,
which tilts a few degrees for a non-integer number of turns. Elements are rotated so that this axis
lies along the plane normal (helix) or along the drawn stroke in the plane (strand). Strands are
turned about their own axis so that the residues facing each other across the sheet point their
side chains the same way.

**Side chains.** CB is placed so that it matches the ProteinMPNN virtual-CB formula to 0.003 Å.
Side chains are built from internal coordinates with a common rotamer for each residue type
(Lovell et al. 2000 modal χ values). Proline uses a narrower N–CA–CB angle and an exo pucker so that
its ring closes at 1.47 Å. Mutating a residue to the same type keeps its current χ angles.

**Loops.** Loops are solved with cyclic coordinate descent (Canutescu & Dunbrack 2003) generalised
to many weighted targets. Each rotatable bond is turned to the angle that best brings every target
downstream of it onto its goal, which has a closed form. Non-glycine φ/ψ are held inside the αR, β,
PPII and αL regions, and glycine is held anywhere except near φ = 0. αL costs a small penalty: it is
what tight turns need, and an inverse-folding model will put Gly, Asn or Asp there. Each loop is
solved from 10–30 random starts. Solutions are scored on closure, fit to the drawing, backbone
clashes and αL use; clash cut-offs sit below backbone H-bond distances, so tight turns survive.
For **Connect** the targets are the N, CA and C of the residue being joined. For a **drawn loop**
the path is smoothed, resampled at the chosen spacing and used as CA targets.

What to expect from loops:

- Connecting adjacent antiparallel strands gives 2-residue β-hairpins (closure ≈ 0.07 Å, no
  clashes), with either ALA or GLY loop residues.
- A helix-to-helix connection usually closes to a few hundredths of an Å. The status bar reports
  the closure and any close contacts after every connection. Above 0.5 Å it says the loop is not
  closed; set a longer loop.
- **Drawn loops follow the drawing loosely near an anchored end.** The first CA out of a chain end
  always sits within about 20° of that residue's CA→C bond, so a path that leaves a helix at a sharp
  angle cannot be followed for the first two residues. To bend toward the drawing, the start
  residue's own φ may change. A loop drawn from a free end or grown from one C-terminus sits within
  about 1–2 Å of the path away from the anchor. A loop that must also *join* an N-terminus is
  solved by fitting each half from its own end and closing a short 3–5 residue bridge in the middle.
  It closes reliably but can sit 3–5 Å off an awkward drawing. If the fit is poor, draw a path
  that leaves and enters along the chain direction, or use Connect.

**Refinement.** Refine minimises a restrained target function over the Cartesian coordinates of the
selected residues, with everything else held fixed. Every term is a harmonic restraint on an
interatomic distance, apart from one signed-volume term per CA. A bond *angle* is restrained as the
1-3 distance its ideal value implies, and the peptide torsion as the 1-4 CA(i)···CA(i+1) distance
(3.80 Å trans, 2.77 cis, whichever isomer the bond already has); written that way the gradients are
simple and have no singularities, unlike angle and torsion derivatives. The terms are:

| Term | What it does |
|---|---|
| Bonds and angles | N–CA, CA–C, C–N, C=O, CA–CB and the 1-3 distances, to Engh & Huber values |
| Peptide unit | CA···CA and O···CA across the bond: keeps it planar and trans (or cis) |
| Chirality | the signed volume (N−CA)·((C−CA)×(CB−CA)), so CA centres stay L |
| Side chains | all internal distances, taken from a reference residue at the same χ, so the side chain travels with the backbone and keeps its rotamer (or is repacked — below) |
| Links | disulfides pulled to 2.03 Å; a ligand or metal within 2.4 Å of a residue held at its current distance |
| Soft spheres | atoms closer than 0.88 × the sum of their van der Waals radii are pushed apart (N···O capped at 2.80 Å so hydrogen bonds survive) |
| Conformation | two weak 1-4 distances per peptide that pin φ and ψ near their current values |
| φ | a one-sided push out of the band around φ = 0 that no non-glycine residue can occupy |
| Positions | a weak restraint of every atom to where it started |

Minimised with L-BFGS and an Armijo line search, then the restraints are rebuilt from the result and
it runs again, up to three times — a φ forced into the forbidden band while a gap was being closed
works its way out on the second pass. Typical results: 0.6 Å of random noise on a 10-residue helix
comes back to 0.001 Å bond RMS and no clashes in under 60 ms; crambin refines in ~180 ms, moving
0.15 Å RMS from the crystal coordinates with its three disulfides intact.

What it will and will not do:

- It **regularises, it does not refold.** φ and ψ are held close to what they are. If the
  conformation itself is wrong, fix it with the torsion fields, Connect, or Draw loop. Repacking
  changes χ only, never the backbone.
- It closes chain breaks **up to the "Close gaps ≤" distance** (4 Å by default) by pulling the two
  sides together. Wider gaps are left alone and reported — build a loop there with Connect instead.
- **Side chains** are set by *Side chains* in the panel: *Rebuild at same χ* (the default) keeps the
  rotamer and just regularises the geometry, *Repack* searches χ (below), *Leave alone* holds the
  side chains exactly as they are.
- *Hold shape* trades faithfulness to your positions against ideal geometry; *Clash target* sets how
  hard atoms push apart. Refining an already-clean region is a no-op (it moves ~0.01 Å).

**Repacking side chains.** With *Side chains › Repack*, each minimisation pass is preceded by a
search over the selected residues' χ angles for the rotamers that fit. Rotamers are the staggered
positions about each sp3 bond plus the ring and amide orientations each residue type actually uses
(3 for Ser through 81 for Arg), expanded combinatorially and pruned of combinations that hit their
own backbone; two-fold symmetric groups — Phe and Tyr rings, Asp and Glu carboxylates — only span
180°. Scoring is the same soft-sphere function as refinement, plus a small χ1 preference as a
tie-break. The assignment is found by going round the sites repeatedly, picking each residue's best
rotamer given its neighbours, from four starting points.

This is **not** a backbone-dependent library with real frequencies, and there is no attractive term
— no packing energy, no hydrogen bonds. It resolves sterics, which is what a template needs, since
an inverse-folding model is about to choose the sequence anyway.

Three things it is careful about:

- **Only the selection moves.** Side chains outside it are fixed environment, so a clash with a
  fixed neighbour is resolved by moving the selected residue.
- **Committed side chains are left alone**: one end of a disulfide, or anything within bonding
  distance of a ligand or metal. The packer has no notion of a bond, so an S–S at 2.03 Å would read
  as a terrible clash and it would cheerfully rotate out of it. Both ends of all three crambin
  disulfides are detected and held.
- **A rotamer is only changed if it is a real improvement** (by more than a fixed margin), so a
  clash-free side chain you placed deliberately is not swapped for a marginally better score.
  Repacking an already-clean structure changes nothing.

What this buys you is a backbone that does not have to distort to make room. Two poly-Trp helices
9 Å apart, refined to the same clash-free result: the backbone moves **1.10 Å RMS** with rigid side
chains, **0.40 Å** when they can repack. For poly-Arg it is 0.42 Å against 0.10 Å.

**Warnings in the view.** Chain breaks (C–N longer than 2.0 Å) are drawn as dashed red lines with
their length. Contacts closer than 2.5 Å between backbone and CB atoms of non-neighbouring residues
are marked with orange crosses. Both are counted in the status bar.

## Files

- **Import.** `.pdb`, `.ent` and `.cif` / `.mmcif`. Open replaces the model; Import adds to it,
  renaming chain IDs that are already used. The first model and first alternate location are
  kept. Hydrogens and waters are dropped, and MSE, HID/HIE/HIP, CYX and similar names are folded
  into the standard 20. Residues with N, CA and C go into chains; everything else becomes a ligand,
  which can be shown, selected, moved and exported. Secondary structure is read from HELIX / SHEET
  or `_struct_conf` / `_struct_sheet_range` when present, and otherwise assigned from φ/ψ.
- **Export.** PDB with HELIX / SHEET records, one TER per chain, ligands as HETATM; or mmCIF
  (`_atom_site` with both label and auth fields). PDB keeps only the first character of a chain ID,
  so use mmCIF for longer IDs.

## Tests

`python tools/serve.py`, then open http://127.0.0.1:8770/tests/core_test.html. It runs 46 checks:
NeRF round trips, helix handedness, rise and H-bonds, CB chirality, all 20 side chains, Kabsch,
element placement, sheet pairing, torsion edits, chain extension, loop closure and Ramachandran
regions, β-hairpins, drawn-loop fits, PDB / mmCIF round trips including 1CRN, and refinement
(analytic gradients against finite differences, noise removal for every residue type, rotamers
kept, context held fixed, gap closing, backbone-only models, φ correction, and 1CRN refined without
losing its fold or its disulfides) and repacking (the rotamer library, clashes cleared with less
backbone movement than rigid side chains need, disulfides left alone, and a clean structure left
untouched).

## Not done

- Refinement is geometric: ideal bond geometry, sterics and the restraints listed above. There is no
  force field, no electrostatics, no solvation and no attractive packing term, and the rotamer
  library is backbone-independent. It will not tell you whether a design is foldable.
- No hydrogen-bond-aware sheet registry: strands are placed geometrically, and pairing quality comes
  from the ideal geometry only.
- No export of design masks (fixed / designable positions for ProteinMPNN); the template is the
  structure alone.

## Licence

MIT — see [LICENSE](LICENSE). Use it, change it, fork it.

The test structure `tests/1crn.pdb` is [PDB entry 1CRN](https://www.rcsb.org/structure/1CRN),
redistributed from the public RCSB archive.
