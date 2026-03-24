# coot-trimmings
Python customizations for the macromolecular model building software Coot.

Copy python file (`coot_trimmings.py`) to the `~/.coot-preferences/` directory (hidden dir, copy on command line, e.g `cp coot_trimmings.py ~/.coot-preferences/`) and restart coot.

You should see a new menu ("Custom") and a bunch of new key bindings, as well as a couple of new toolbar buttons (e.g. "sequence context").

(If you see an errror like `<!DOCTYPE html> ^ SyntaxError: invalid syntax`, you have downloaded the webpage rather than the actual Python script. Make sure you clicked "download raw file" when downloading.)

(*NOTE: There are some non-default settings (e.g. I don't like using the scroll-wheel for changing map contours) which I like but you may not! Check the bottom of this readme for details and tweak as you prefer.*)

Email me if you run into any trouble (olibclarke at gmail dot com).

I haven't tested most of these on WinCoot, or on older versions of Coot. Most testing done with Coot 0.9 or later (_Doesn't work with Coot 1.1! But I have new scripts for that here: https://github.com/olibclarke/coot1-trimmings/tree/main_).

Also, many of the functions here haven't been tested on molecules with insertion codes. If this is important to you, let me know and I'll fix.

**_Note: I turn off scrolling contour with mouse wheel by default (`set_scroll_by_wheel_mouse(0)`), as it doesn't play well with Mac trackpads - if you prefer it, you can change it back to `set_scroll_by_wheel_mouse(1)` in the script._**

## Example use case/Tutorial
1. Load up an EM-map and model (e.g. `8HEZ` & `EMD-34705`)
2. Adjust the contour - `Shift-1` - `Shift-9` set the level from `1` to `9` sigma respectively, and `=` and `-` nudge the contour level up and down.
3. Zoom in - you might notice the mesh is a bit coarse, which can make fine details difficult to interpret. This is a 2.8 Å structure, but was solved at 1.1 Å/pixel. Let's try resampling and brightening things up a bit to make it easier to see. To do this use "Custom -> Maps... -> Resample active EM map to 0.5 A/pixel". (_Note: there are some smarts here - if the map is an X-ray map, it will not be resampled, but color changes etc will still happen; if it is an EM map but is already finer than 0.5 Å/pix, the sampling will not be changed; if it is coarser than 0.5 Å/pix, it will be resampled to 0.5 Å/pix and the original map will be closed._)
4. Use `;` and `'` to adjust the radius of the map in local mode.
5. Sometimes, you want to quickly toggle the model off to better see the map, or vice-versa. To toggle the model on/off, press `/`; to switch between multiple models, press `Shift-/`. To toggle the map, press `` ` ``; to switch between multiple maps, press `Shift-``.
6. Try cycling the representation of the model - press the `[` and `]` keys to move back and forward between the different modes (e.g. CA-only and all atom).
7. Try cylinder refine - center on an atom, and press `a`. It will refine +/- 5 residues around the center, and all segments within 4 Å of the primary range.
8. Try quick refine zone - press `A`, then click two atoms to define a range and refine.
9. Try copy/pasting a ligand - center on a ligand or solvent molecule, and press `C`. A status message in the lower bar should indicate that the molecule has been copied. Move somewhere else, and press `V`. A copy should be pasted and merged into the active molecule. You can repeat this operation multiple times.
10. Try some of the custom coloring modes - try coloring by density fit, ramachandran outliers or clash score. If you cycle through the representation modes (`[` and `]`) you will see that there are custom coloring modes for both all atom and CA-only representation.

## Custom menu

The script creates a `Custom` menu with these entries:

- `Custom keybindings...`
- `Display...`
- `Colour...`
- `Fit...`
- `Renumber...`
- `Settings...`
- `Build...`
- `Mutate...`
- `Copy...`
- `Delete...`
- `Merge...`
- `Maps...`

Some notable menu functions:

- `Colour...`
  - color by rotamer outliers and missing atoms
  - color by hydrophobic/polars
  - color by charge
  - color by Ramachandran outliers
  - color by density fit
  - color by NCS difference
  - highlight chain breaks
- `Fit...`
  - fit chains and segments to map
  - local cylinder refine
  - smart self restraints with user-entered cutoff
- `Settings...`
  - auto-scale B-factor colouring
  - set B for new atoms from active molecule mean B
  - **set proportional editing radius** 
- `Copy...`
  - copy/cut chains and segments
  - **copy active ligand/ion/solvent**
  - **paste copied ligand/ion/solvent** 
- `Maps...`
  - map sharpening / blurring helpers
  - low-pass filtered map generation
  - resample active EM map to `0.5 A/pixel`
- `Display...`
  - Load molecular symmetry copies from file metadata

## Proportional editing radius explanation
Changing the proportional editing radius changes the "sphere of influence" when dragging during interactive dragged refinement. A larger proportional editing radius means one can drag a larger region without disturibng local geometry - e.g. think moving a domain, or adjusting the register of a helix. There is a mouse interface for this currently (Ctrl-scroll during dragged refinement) but it can be fragile on some systems. 

I have added a dialog (_Custom...Settings...Set Proportional Editing Radius_) to directly set this radius; However, because there is no way in the API to get the current radius, it assumes that you are not using the Ctrl-scroll method for adjusting the radius. If you are, unexpected results may ensue.

## Custom keybindings

These are the main custom shortcuts defined by the script:

### Display and navigation

| Key | Action |
| --- | --- |
| `G` | Toggle the active map between local mesh view and global solid view |
| `?` | Show only the active model, or cycle displayed models |
| `~` | Show only the active map, or cycle displayed maps |
| `` ` `` | Toggle display of all maps |
| `/` | Toggle display of all model molecules |
| `[` | Cycle model representation mode |
| `]` | Cycle model representation mode backward |
| `{` | Cycle symmetry representation mode |
| `}` | Cycle symmetry representation mode backward |
| `>` | Go to next residue in the current chain |
| `<` | Go to previous residue in the current chain |
| `o` | Go to the next NCS-related chain |
| `O` | Go to the NCS master chain |
| `D` | Toggle environment distances |
| `Z` | Clear labels and distances |
| `v` | Undo symmetry view, only if symmetry is currently shown |
| `Tab` | Clear pending picks |
| `b` | Go to the nearest density peak near the current rotation centre |

### Map control

| Key | Action |
| --- | --- |
| `!` to `(` | Set the current map to `1` to `9` sigma |
| `=` | Increase current map contour |
| `-` | Decrease current map contour |
| `;` | Decrease map radius |
| `'` | Increase map radius |

### Building and editing

| Key | Action |
| --- | --- |
| `h` | Place helix here |
| `m` | Measure distance |
| `w` | Place water |
| `W` | Place water and refine |
| `y` | Add terminal residue |
| `Y` | Cycle terminal-residue phi |
| `T` | Cycle terminal-residue psi |
| `M` | Mutate active residue by one-letter code |
| `q` | Flip peptide |
| `X` | Delete active residue |
| `K` | Delete active sidechain |
| `k` | Fill partial sidechain |
| `C` | Copy active ligand/ion/solvent |
| `V` | Paste copied ligand/ion/solvent at the pointer |

### Refinement and fitting

| Key | Action |
| --- | --- |
| `A` | Refine the clicked residue range |
| `a` | Local cylinder refine around the active residue or ligand |
| `r` | Refine triple around the active residue |
| `J` | Jiggle-fit the active non-polymer residue |
| `R` | Cycle rotamers for the active residue |
| `g` | Generate smart local extra restraints for the active model |

### Saving and history

| Key | Action |
| --- | --- |
| `Q` | Save and overwrite the active model |
| `z` | Undo |
| `x` | Redo |

## Smart local extra restraints

The `g` shortcut and the `Custom -> Fit -> Smart self restrain active mol...` menu item both generate extra Geman-McClure restraints.

Current behavior:

- same-chain local restraints are generated for inter-residue contacts within the chosen distance cutoff
- residue pairs more than `10` positions apart in true chain order are excluded from the local pass
- long-range backbone `N···O/OXT` contacts are retained
- existing extra restraints on the active molecule are cleared first

The hotkey uses a fixed cutoff of `3.7 A`.

The menu item lets the user enter the cutoff explicitly.

## Smart ligand copy / paste

The script includes a simple ligand/ion/solvent copy buffer.

Copy:

- works on the active non-polymer residue
- stores a hidden copied template molecule

Paste:

- duplicates the stored template
- translates it to the current pointer / rotation centre
- merges it into the active model
- supports repeated paste operations

This is available both from hotkeys (`C` / `V`) and from `Custom -> Copy...`.

## Colouring modes

Several custom colouring modes use Coot's user-defined residue colouring.

Current built-in custom colouring options include:

- rotamer outliers and missing atoms
- hydrophobic / polar residue classes
- charge colouring
- protein vs nucleic-acid colouring
- water colouring
- Ramachandran outliers
- density fit
- NCS difference
- Clash score

Details:

- Ramachandran colouring
  - red = outlier
  - orange = allowed/disfavored
- Density-fit colouring
  - spectral colouring
  - blue = model/map correlation `CC 1.0`
  - red = model/map correlation `CC 0.0`
- NCS-difference colouring
  - spectral colouring on a fixed `0 to 2 A` scale
  - blue = low NCS difference
  - red = high NCS difference
  - values above `2 A` saturate at the red end
- Clash score coloring
  - spectral colouring by per-residue maximum clash overlap
  - blue = low clash
  - red = high clash  

## Molecular symmetry
Load molecular symmetry copies from file metadata reads deposited biological-assembly operators from the active model file and displays the corresponding molecular symmetry copies in Coot. For PDB files it uses REMARK 350 BIOMT; for mmCIF files it uses _pdbx_struct_oper_list together with _pdbx_struct_assembly_gen. This is useful for regenerating full symmetric assemblies from a single deposited chain or asymmetric unit.

if _"Load molecular symmetry copies from file metadata"_ does not find usable BIOMT / mmCIF assembly operators, it opens a manual entry prompt so you can enter a point group such as C13.

It assumes:
- symmetry center = unit-cell center
- principal Cn axis = z


## EM map resampling / restyling

The menu item:

- `Custom -> Maps... -> Resample active EM map to 0.5 A/pixel`

does the following:

- leaves difference maps unchanged
- reads the CCP4/MRC header to determine grid spacing
- if the spacing is already `<= 0.5 A/pixel`, it simply restyles the current map (changes color to a default blue)
- otherwise resamples to `0.5 A/pixel`, preserves the current contour sigma, makes the new map active, and closes the old map

## Toolbar buttons
* Measure distance
* Toggle display of symmetry copies
* Sequence context (displays local sequence before/after active residue)
* Accept RSR (forces acceptance of current real space refine results)

## Non-default setttings

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

