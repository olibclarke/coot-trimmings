# coot-trimmings
Python customizations for the macromolecular model building software Coot.

Copy to ~/.coot-preferences (hidden dir, copy on command line) and restart coot.

You should see a new menu ("Custom") and a bunch of new key bindings, as well as a couple of new toolbar buttons (e.g. "sequence context").

Email me if you run into any trouble (olibclarke at gmail dot com).

I haven't tested most of these on WinCoot, or on older versions of Coot. Most testing done with Coot 0.8.8 or later.

Also, many of the functions here haven't been tested on molecules with insertion codes. If this is important to you, let me know and I'll fix.

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


# Custom keybindings
To be added

# Custom menu items
To be added

# Toolbar buttons
* Measure distance
* Toggle display of symmetry copies
* Sequence context (displays local sequence before/after active residue)
* Accept RSR (forces acceptance of current real space refine results)

