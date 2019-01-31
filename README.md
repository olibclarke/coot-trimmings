# coot-trimmings
Python customizations for the macromolecular model building software Coot.

Copy python file (`oli_custom.py`) to the `~/.coot-preferences/` directory (hidden dir, copy on command line, e.g `cp oli_custom.py ~/.coot-preferences/`) and restart coot.

You should see a new menu ("Custom") and a bunch of new key bindings, as well as a couple of new toolbar buttons (e.g. "sequence context").

(*There are some non-default settings which I like but you may not! Check the bottom of this readme for details and tweak as you prefer*)

Email me if you run into any trouble (olibclarke at gmail dot com).

I haven't tested most of these on WinCoot, or on older versions of Coot. Most testing done with Coot 0.8.8 or later.

Also, many of the functions here haven't been tested on molecules with insertion codes. If this is important to you, let me know and I'll fix.

# Custom keybindings
Ordered approximately by interestingness/usefulness. YMMV.

NEW "G": Toggle display of the active map between a local mesh and a global surface.

NEW "L": Set contour level of scrollable map in sigma by entering a new value.

NEW "Y"/"T": Cycle phi/psi angles of terminal residue. Useful when manually building in conjunction with "y"

"M": Mutate active residue by entered single letter code (case-insensitive).

"?": Display only active model; If only active model displayed, cycle display of models.

"~": Display only active map; If only active map displayed, cycle display of maps.

"\`": Toggle display of all maps.

"/": Toggle display of all models.

"]"/"[": Cycle representation mode forward/back for active model.

"}"/"{": Cycle symmetry representation mode forward/back for active model.

"R": Cycle through rotamers for active residue.

"H": Toggle hide/display of modelling toolbar.

"h": Place helix here, with H-bond restraints added.

"Tab": Clear pending picks.

"!"..."(": Set active map contour to 1,2,3...9 sigma.

">"/"<": Next/previous residue.

"Q": Save and overwrite active model (makes backup in case of accidents).

"A": Real space refine zone (click start and end of zone)

"Z": Clear distances and labels.

"P": Place atom at pointer.

"q": Pepflip active residue.

"a": Auto-refine zone of 10 residues centered on active residue.

"J": Jiggle fit active residue.

"D": Toggle display of environment distances.

"X": Delete active residue.

"K": Kill sidechain of active residue.

"k": Fill sidechain of active residue.

"w": Place water without refinement.

"W": Place water with refinement.

"y": Add terminal residue.

"r": Refine three residues centered on active residue.

"V": Undo symmetry view.

"z": Undo for active model.

"x": Redo for active model.

"O": Go to equivalent residue on NCS master chain.

"|":/"\_": Increase/decrease active map level by 0.5 sigma.

# Custom menu items
(Very incomplete list - these are some highlights)

## _Display_
* Colour by rotamer probability/missing atomsl; hydrophobics/polars; +ve/-ve charge; entered subset of residues
* Highlight chainbreaks with dotted lines; red >50 residues, orange 15-50 residues missing, gray <15 residues missing
* Color active segment (covalently connected polymer segment).
* Open current view in UCSF chimera (requires chimera in PATH): Hopefully does what it says on the box, including changing view and display of maps. Useful as a starting point for making density figures.

## _Renumber_
* Renumber active chain by current res: Adjusts sequence numbering of chain so that active residue matches entered number.
* Renumber active segment (contiguous stretch of sequence, bounded by chain breaks) by current res: Adjusts sequence numbering of segment so that active residue matches entered number. Checks for overlapping numbering.

## _Settings_
* Auto-scale B-factor coloring
* Change default B for new residues to mean B of active model.

## _Build_
* Forced addition of terminal residue: adds residue, ignoring map, for when you disagree with Coot's assessment of residue placeability.
* Rebuild backbone: Uses db_mainchain to rebuild the selected zone based on similar sequence fragments in a database of high resolution structures.
* Make alkyl chain of length n: Makes alkyl chain and jiggle-fits to map. Useful for preliminary modelling of lipids/detergents.
* Grow helix, grow strand etc. Grows selected helix/strand by entered number of residues assuming ideal geometry.

## _Mutate_
* Mutate range to ALA
* Mutate range to UNK
* Mutate Mets to MSE and vice-versa
* Mutate active chain to template sequence: Mutates active chain to sequence pasted into textbox, assuming numbering the same - i.e. first residue in pasted sequence is residue 1 in numbering scheme of model. Filters non-sequence (e.g. formatting) characters from pasted sequence. Once sequence is pasted once, it is remembered after, to speed up multiple register adjustments.

## _Copy_
* Cut/copy chain, segment, selected fragment, etc.

## _Delete_
* Delete active chain
* Delete active segment
* Delete sidechains in range
* Delete hydrogens in active molecule

## _Merge_
* Merge two molecules or chains by clicking.

## _Maps_
* Sharpen by entered B-factor. Absolute not relative; 0 always returns to original map. (MTZ only!)
* Change hi-res-limit for map; creates low-pass filtered version of active map (MTZ only!)
* Go to center of active map.
* Set refinement map to active (scrollable) map.



# Toolbar buttons
* Measure distance
* Toggle display of symmetry copies
* Sequence context (displays local sequence before/after active residue)
* Accept RSR (forces acceptance of current real space refine results)

# Non-default setttings

`set_symmetry_colour(255,35,0)`

Makes sym copies brighter (yellow), rather than grey.


`set_refine_max_residues(100)`

Increases max number of residues included in a refinement.

`set_rotamer_search_mode(ROTAMERSEARCHLOWRES)`

Use "backrub" rotamers (better at low res).


`set_refine_ramachandran_angles(1)`

Turn ramachandran restraints on by default during real space refinement


`set_map_sampling_rate(2.0)`

Increase map  sampling a bit for nicer looking maps.


`allow_duplicate_sequence_numbers()`

So Coot doesn't choke on the occasional entry in the PDB which has duplicate residue numbers.


`set_add_terminal_residue_n_phi_psi_trials(1000)`

Increase number of trials when adding a terminal residue


`set_refinement_drag_elasticity(0.5)`

Change RSR back to old "rubber-banding" behaviour.


`set_matrix(20.0)`

Change refinement weighting to be more appropriarte for low res (more weighting on geometry).


`set_smooth_scroll_flag(0)`

Make scrolling faster.


`set_active_map_drag_flag(0)`

Update mag after dragging - better for large map radius.


`set_show_environment_distances(0)
set_show_environment_distances_bumps(0)
set_show_environment_distances_h_bonds(1)`

Don't show environment distances by default, and only show H-bonds when shown.


`set_environment_distances_distance_limits(2.1,3.2)`

Set distance limits for environment distances.


`set_map_radius(20)`

Set default map radius.


`set_default_temperature_factor_for_new_atoms(50.0)`

Increase default B for new atoms.


`set_add_terminal_residue_do_post_refine(1)
set_terminal_residue_do_rigid_body_refine(0)`

Post-refine when adding term res and don't rigid body fit  (better for low res).


`set_mutate_auto_fit_do_post_refine(1)`

Real-space refine after mutating residue.


`set_scroll_by_wheel_mouse(0)`

Turn scroll-wheel adjustment of contouring off - better for trackpads, large maps.


`set_symmetry_size(30)`

Set default symmetry radius to 30 Å.


`set_nomenclature_errors_on_read("ignore")`

Turn off nomenclature errors, because they are annoying and omnipresent.

