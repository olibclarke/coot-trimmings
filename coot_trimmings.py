#Oli's Coot customizations:
#Usage: Copy to ~/.coot-preferences
#When Coot is restarted, a new "Custom" menu will appear, with some new shortcuts for various model-building tasks.
#Also there will be a bunch of new keyboard shortcuts - check "Extensions...Settings...Key bindings" for details.

import math

try:
  import gap as _plain_gap_module
except:
  _plain_gap_module=None

MOLECULAR_SYMMETRY_SESSION_STATE={}

#****Settings****
#Make symmetry copies a brighter color
set_symmetry_colour(255,35,0)

#turn off scrolling to change contour (for trackpads)
set_scroll_by_wheel_mouse(0)

#Consider clashes when autofitting rotamers.
set_auto_fit_best_rotamer_clash_flag(1)

#Increase the default limit for the max number of residues to refine
set_refine_max_residues(100)

#Sets "Backrub" rotamers as default (best at low res)
set_rotamer_search_mode(ROTAMERSEARCHLOWRES)

#Keep ramachandran restraints off by default
set_refine_ramachandran_angles(0)

#Use finer map sampling
set_map_sampling_rate(3.0)

#Allow duplicate sequence numbers (otherwise some PDBs won't load)
allow_duplicate_sequence_numbers()
#Increase number of trials for add terminal residue
set_add_terminal_residue_n_phi_psi_trials(1000)

#set default refinement weighting for low resolution
set_matrix(20.0)

#Make scrolling faster
set_smooth_scroll_flag(0)

#Update map after dragging - improves speed with large map radius
set_active_map_drag_flag(0)

#Default to not showing environment distances, and only showing h-bonds if shown
set_show_environment_distances(0)
set_show_environment_distances_bumps(0)
set_show_environment_distances_h_bonds(1)

#Set distance limits for environment distances
set_environment_distances_distance_limits(2.1,3.2)

#Set default map radius
set_map_radius(20)
set_map_radius_em(20)

#Set default bond thickness
set_default_bond_thickness(3)

#Use variable width bonds
try:
  set_use_variable_bond_thickness(1)
  set_default_bond_thickness(7)
except:
  print("Your coot is a bit old... consider upgrading...")

#Increase default B for new atoms
set_default_temperature_factor_for_new_atoms(50.0)

#Post-refine when adding term res and don't rigid body fit  (better for low res)
set_add_terminal_residue_do_post_refine(1)
set_terminal_residue_do_rigid_body_refine(0)

#real space refine after mutating residue
set_mutate_auto_fit_do_post_refine(1)

#Set symmetry radius to 30 A
set_symmetry_size(30)

#Ignore nomenclature errors
set_nomenclature_errors_on_read("ignore")

def _safe_low_density_average(imol_map, imol, chain_id, start_resno, stop_resno):
  map_coords=[]
  map_density=[]
  for resno in range(start_resno, stop_resno + 1):
    atom_ls=residue_info(imol, chain_id, resno, "")
    if not atom_ls:
      continue
    for atom in atom_ls:
      if atom[0][0] in [' N  ', ' CA ', ' CB ', ' C  ', ' O  ']:
        map_coords.append(atom[2])
  for [x, y, z] in map_coords:
    map_density.append(density_at_point(imol_map, x, y, z))
  if not map_density:
    return 0.0
  map_density.sort()
  cut_off=len(map_density) // 5
  if (cut_off <= 0):
    cut_off=1
  map_average=sum(map_density[0:cut_off]) / float(cut_off)
  return map_average

def _patch_gap_low_density_average():
  if _plain_gap_module is not None:
    _plain_gap_module.low_density_average=_safe_low_density_average
  try:
    fit_gap.func_globals["low_density_average"]=_safe_low_density_average
  except:
    pass
  return True

_patch_gap_low_density_average()

#****Modules****
add_module_cryo_em_gui()
add_module_refine()
#add_module_restraints()
add_module_carbohydrate_gui()

#****Keybindings and toolbar buttons****
# Measure distance shortcut
coot_toolbar_button("Measure distance", 
"do_distance_define()", icon_name="geom.svg")

#Toggle display of symmetry copies
coot_toolbar_button("Sym?", 
"set_show_symmetry_master(not get_show_symmetry())", 
icon_name="cell+symm.svg")



#Place helix here
add_key_binding("Place helix here","h",
lambda: place_helix_here())

#Toggle display of modelling toolbar (assumes initial state is shown)
add_key_binding("Toggle toolbar display","H",
lambda: toggle_toolbar_display())

#Toggle global display of map
add_key_binding("Toggle global view of map","G",
lambda: toggle_global_map_view())


#Quicksave active mol (overwrite orig)
add_key_binding("Save and overwrite active model","Q",
lambda: quicksave_active())

#Refine zond (click two atoms)
add_key_binding("Refine zone","A",
lambda: refine_residues_click())

#Bring up place atom at pointer dialog
add_key_binding("Place atom at pointer","P",
lambda: place_atom_at_pointer())


#Flip peptide.
add_key_binding("Flip peptide","q",
lambda: pepflip_active_residue())

#Local cylinder refinement around the active residue or ligand
add_key_binding("Auto refine zone","a",
lambda: auto_refine())

#Jiggle fit active non-polymer residue
add_key_binding("Jiggle Fit","J",
lambda: jiggle_fit_active_non_polymer_residue())

#Clear pending picks
add_key_binding("Clear Pending Picks","Tab",
lambda: clear_pending_picks())

#Toggle environment distances
add_key_binding("Toggle environment distances","D",
lambda: toggle_env_dist()) 

#Delete active residue
add_key_binding("Delete this residue","X",
lambda: using_active_atom(delete_residue,"aa_imol","aa_chain_id","aa_res_no","aa_ins_code"))

#Delete sidechain
add_key_binding("Kill Sidechain","K",
lambda: using_active_atom(delete_residue_sidechain,"aa_imol","aa_chain_id","aa_res_no","aa_ins_code",0))

#Fill sidechain
add_key_binding("Fill Sidechain","k",
lambda: using_active_atom(fill_partial_residue,"aa_imol","aa_chain_id","aa_res_no","aa_ins_code"))

#place water without refinement
add_key_binding("Add Water","w",
lambda: place_water_in_active_molecule())

#place water with refinement
add_key_binding("Add Water +","W",
lambda: add_water_and_refine())

#add terminal residue
add_key_binding("Add terminal residue","y",
lambda: add_term_shortcut())
 
#add terminal residue
add_key_binding("Cycle terminus phi","Y",
lambda: cycle_residue_phi())

#add terminal residue
add_key_binding("cycle terminus psi","T",
lambda: cycle_residue_psi())

 
#Refine active residue
add_key_binding("Refine Triple","r",
lambda: key_binding_refine_triple())

#Undo symmetry view only when symmetry is currently shown
def undo_symmetry_view_if_shown():
  if get_show_symmetry():
    return undo_symmetry_view()
  add_status_bar_text("Symmetry atoms are not currently shown")

#Undo Symm view
add_key_binding("Undo Symmetry View", "v",
lambda: undo_symmetry_view_if_shown())

#Cycle through rotamers for active reidue with 'R"
add_key_binding("Cycle rotamers","R",
lambda: cycle_rotamers())

#Undo function for keybinding. Undoes last change to active mol.
add_key_binding("Undo","z",
lambda: undo_visible())

#Redo function for keybinding. redoes last change to active mol.
add_key_binding("Redo","x",
lambda: redo_visible())


add_key_binding("Cycle representation mode forward","[",
lambda: cycle_rep_up_active())

add_key_binding("Cycle representation mode back","]",
lambda: cycle_rep_down_active())

add_key_binding("Cycle  symm representation mode forward","{",
lambda: cycle_symm_up_active())

add_key_binding("Cycle  symm representation mode back","}",
lambda: cycle_symm_down_active())

add_key_binding("Toggle map display","`",
lambda: toggle_map_display())

add_key_binding("Toggle mol display","/",
lambda: toggle_mol_display())

#Clear distances/labels
add_key_binding("Clear distances and labels","Z",
lambda: clear_distances_and_labels())

#Shortcut to set map level in sigma (useful for EM maps)
add_key_binding("Set map contour in sigma","L",
lambda: set_map_level_quickly())

#Show local sequence context for active residue
coot_toolbar_button("Sequence context",
"sequence_context()", icon_name="")


#Force accept real space refinemnt
coot_toolbar_button("Accept RSR",
"accept_regularizement()", icon_name="")

#mutate active residue to entered residue code (upper or lower case single-letter)
add_key_binding("Mutate by single letter code","M",
lambda: mutate_by_entered_code())


#Bind next_res() and prev_res() to ">" and "<"
add_key_binding("Next residue in chain",">",
lambda: next_res())
add_key_binding("Prev residue in chain","<",
lambda: prev_res())

#The nine key bindings elow allow easy setting of map
#level by rmsd - shift + any single digit integer
#sets the currently scrollable map to that level
# in rmsd. Useful when on a laptop with touchpad,
#when changing the contour using the scrollwheel is
#not practical and using +/- is too slow.
add_key_binding("Map to 1 sigma","!",
lambda: set_contour_level_in_sigma(scroll_wheel_map(),1))

add_key_binding("Map to 2 sigma","@",
lambda: set_contour_level_in_sigma(scroll_wheel_map(),2))

add_key_binding("Map to 3 sigma","#",
lambda: set_contour_level_in_sigma(scroll_wheel_map(),3))

add_key_binding("Map to 4 sigma","$",
lambda: set_contour_level_in_sigma(scroll_wheel_map(),4))

add_key_binding("Map to 5 sigma","%",
lambda: set_contour_level_in_sigma(scroll_wheel_map(),5))

add_key_binding("Map to 6 sigma","^",
lambda: set_contour_level_in_sigma(scroll_wheel_map(),6))

add_key_binding("Map to 7 sigma","&",
lambda: set_contour_level_in_sigma(scroll_wheel_map(),7))

add_key_binding("Map to 8 sigma","*",
lambda: set_contour_level_in_sigma(scroll_wheel_map(),8))

add_key_binding("Map to 9 sigma","(",
lambda: set_contour_level_in_sigma(scroll_wheel_map(),9))

add_key_binding("Map plus 0.5 sigma","|",
lambda: step_map_coarse_up(scroll_wheel_map()))

add_key_binding("Map minus 0.5 sigma","_",
lambda: step_map_coarse_down(scroll_wheel_map()))

#Undisplay all models except the active one.
#If only one model is displayed, cycle through
#all available models.

add_key_binding("Display only the active model","?",
lambda: display_only_active())  



#Undisplay all maps except the active one.
#If only one map is displayed, cycle through
#all available models.
add_key_binding("Display only the active map","~",
lambda: display_only_active_map())

#Go to equivalent residue on next NCS chain / master chain

def goto_next_ncs_chain():
  try:
    return skip_to_next_ncs_chain("forward")
  except:
    _status_message("Unable to skip to next NCS chain")
    return None

def goto_ncs_master():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  master_chain_id=_ncs_master_chain_id(mol_id)
  if not master_chain_id:
    _status_message("No NCS master chain found")
    return None
  start_chain_id=go_to_atom_chain_id()
  if start_chain_id == master_chain_id:
    return None
  visited_chain_ids={}
  while True:
    current_chain_id=go_to_atom_chain_id()
    if current_chain_id == master_chain_id:
      return None
    if visited_chain_ids.get(current_chain_id, 0):
      _status_message("Unable to reach NCS master chain")
      return None
    visited_chain_ids[current_chain_id]=1
    goto_next_ncs_chain()
    if go_to_atom_chain_id() == current_chain_id:
      _status_message("Unable to reach NCS master chain")
      return None
    if go_to_atom_chain_id() == start_chain_id:
      _status_message("Unable to reach NCS master chain")
      return None

add_key_binding("Go to next NCS chain","o",
lambda: goto_next_ncs_chain())

add_key_binding("Go to NCS master chain","O",
lambda: goto_ncs_master())

add_key_binding("Go to nearest density peak","b",
lambda: go_to_nearest_density_peak())

add_key_binding("Smart copy active non-polymer residue","C",
lambda: smart_copy_active_non_polymer_residue())

add_key_binding("Smart paste copied non-polymer residue","V",
lambda: smart_paste_copied_non_polymer_residue())

add_key_binding("Generate smart local extra restraints","g",
lambda: confirm_generate_smart_local_extra_restraints())

add_key_binding("Decrease map radius",";",
lambda: decrease_map_radius_with_status())

add_key_binding("Increase map radius","'",
lambda: increase_map_radius_with_status())

  
#****Misc. functions (for keybindings and scripting****
SMART_COPY_TEMPLATE_IMOL=None
SMART_COPY_SOURCE_CENTRE=None
SMART_COPY_RESIDUE_NAME=None
MAP_GLOBAL_VIEW_SETTINGS={}
EM_GLOBAL_MAP_COLOUR=(0.23137254901960785, 0.4549019607843137, 1.0)
MAX_EM_RESAMPLE_TARGET_VOXELS=250000000
PROPORTIONAL_EDITING_RADIUS=1.0

POLYMER_RESIDUE_NAMES=set([
  "A","C","G","U","T","DA","DC","DG","DT",
  "PSU","H2U","5MU","OMU","4SU",
  "1MA","2MA","6MA","MIA",
  "5MC","OMC",
  "1MG","2MG","7MG","M2G","OMG","YG",
  "ALA","UNK","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
  "LEU","LYS","MET","MSE","PHE","PRO","SER","THR","TRP","TYR","VAL",
  "P1L","SEP","TPO","PTR","CSO","CME","MLY","MLZ","HYP","KCX","ALY",
  "FME","PYL","SEC"
])

def _status_message(message):
  try:
    add_status_bar_text(message)
  except:
    print message

def _active_residue_or_status():
  residue=active_residue()
  if residue:
    return residue
  _status_message("No active residue")
  return None

def _active_molecule_or_status():
  residue=_active_residue_or_status()
  if not residue:
    return None
  return residue[0]

def _positive_int_from_entry(value, dialog_label, minimum=1):
  try:
    parsed=int(value)
  except:
    info_dialog("%s must be a whole number." % dialog_label)
    return None
  if parsed < minimum:
    if minimum == 1:
      info_dialog("%s must be at least 1." % dialog_label)
    else:
      info_dialog("%s must be at least %d." % (dialog_label, minimum))
    return None
  return parsed

_NCS_MASTER_CHAIN_ID_BUILTIN=globals().get("ncs_master_chain_id")

def _ncs_master_chain_id(imol):
  if callable(_NCS_MASTER_CHAIN_ID_BUILTIN):
    try:
      return _NCS_MASTER_CHAIN_ID_BUILTIN(imol)
    except:
      pass
  for helper_name in ["ncs_master_chain_id_py", "ncs_chain_ids"]:
    helper=globals().get(helper_name)
    if not callable(helper):
      continue
    try:
      result=helper(imol)
    except:
      continue
    if isinstance(result, (list, tuple)) and result:
      first=result[0]
      if isinstance(first, (list, tuple)) and first:
        return first[0]
      return first
    if result:
      return result
  return None

def _ncs_master_residue_spec(imol, chain_id, resno, ins_code):
  master_chain_id=_ncs_master_chain_id(imol)
  if not master_chain_id:
    return None
  if chain_id == master_chain_id:
    return [master_chain_id, resno, ins_code]
  try:
    diffs=ncs_chain_differences(imol, master_chain_id)
  except:
    diffs=False
  if not diffs:
    return None
  for i in range(0, len(diffs), 3):
    try:
      peer_chain_id=diffs[i]
      target_chain_id=diffs[i+1]
      residue_diffs=diffs[i+2]
    except:
      continue
    if peer_chain_id != chain_id or target_chain_id != master_chain_id:
      continue
    if not isinstance(residue_diffs, list):
      continue
    for residue_diff in residue_diffs:
      if not isinstance(residue_diff, list) or len(residue_diff) < 2:
        continue
      peer_residue=residue_diff[0]
      target_residue=residue_diff[1]
      if (isinstance(peer_residue, list) and len(peer_residue) >= 2 and
          peer_residue[0] == resno and peer_residue[1] == ins_code and
          isinstance(target_residue, list) and len(target_residue) >= 2):
        return [master_chain_id, target_residue[0], target_residue[1]]
  return [master_chain_id, resno, ins_code]

def _set_go_to_residue_atom(chain_id, resno, ins_code, atom_name, alt_conf):
  full_helper=globals().get("set_go_to_atom_chain_residue_atom_name_full")
  if callable(full_helper):
    try:
      return full_helper(chain_id, resno, ins_code, atom_name, alt_conf)
    except:
      pass
  return set_go_to_atom_chain_residue_atom_name(chain_id, resno, atom_name)

def _pick_existing_target_atom(imol, chain_id, resno, ins_code, preferred_atom_name, preferred_alt_conf):
  residue_atoms=residue_info_py(imol, chain_id, resno, ins_code) or []
  parsed_atoms=[]
  for atom in residue_atoms:
    parsed_atom=_parsed_atom_record(atom)
    if parsed_atom:
      parsed_atoms.append(parsed_atom)
  if not parsed_atoms:
    return (preferred_atom_name, preferred_alt_conf)

  stripped_preferred=preferred_atom_name.strip()
  for atom in parsed_atoms:
    if atom["name"].strip() == stripped_preferred and atom["alt_conf"] == preferred_alt_conf:
      return (atom["name"], atom["alt_conf"])
  for atom in parsed_atoms:
    if atom["name"].strip() == stripped_preferred:
      return (atom["name"], atom["alt_conf"])
  for fallback_name in ["CA", "P"]:
    for atom in parsed_atoms:
      if atom["name"].strip() == fallback_name:
        return (atom["name"], atom["alt_conf"])
  return (parsed_atoms[0]["name"], parsed_atoms[0]["alt_conf"])

def _scrollable_map_or_status():
  map_id=scroll_wheel_map()
  if map_id!=-1 and map_id in map_molecule_list():
    return map_id
  _status_message("No active map")
  return None

def _residue_is_polymer(mol_id, chain_id, resno, ins_code):
  if does_residue_exist_p(mol_id, chain_id, resno, ins_code)==0:
    return False
  return residue_name(mol_id, chain_id, resno, ins_code) in POLYMER_RESIDUE_NAMES

def _rotation_centre_xyz():
  return [rotation_centre_position(0), rotation_centre_position(1), rotation_centre_position(2)]

def _distance_sq(point_1, point_2):
  dx=point_1[0]-point_2[0]
  dy=point_1[1]-point_2[1]
  dz=point_1[2]-point_2[2]
  return dx*dx + dy*dy + dz*dz

def _normalize_vector(vector):
  length=math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2])
  if length < 0.0000001:
    return None
  return [vector[0]/length, vector[1]/length, vector[2]/length]

def _rotation_matrix_about_axis(axis, theta):
  normalized_axis=_normalize_vector(axis)
  if not normalized_axis:
    return None
  ux=normalized_axis[0]
  uy=normalized_axis[1]
  uz=normalized_axis[2]
  cos_theta=math.cos(theta)
  sin_theta=math.sin(theta)
  one_minus_cos=1.0-cos_theta
  return [
    cos_theta + ux*ux*one_minus_cos,
    ux*uy*one_minus_cos - uz*sin_theta,
    ux*uz*one_minus_cos + uy*sin_theta,
    uy*ux*one_minus_cos + uz*sin_theta,
    cos_theta + uy*uy*one_minus_cos,
    uy*uz*one_minus_cos - ux*sin_theta,
    uz*ux*one_minus_cos - uy*sin_theta,
    uz*uy*one_minus_cos + ux*sin_theta,
    cos_theta + uz*uz*one_minus_cos
  ]

def _zeroify_rotation_matrix(matrix):
  if not matrix:
    return None
  zeroified=[]
  for value in matrix:
    if abs(value) < 0.0000001:
      zeroified.append(0.0)
    elif abs(value-1.0) < 0.0000001:
      zeroified.append(1.0)
    elif abs(value+1.0) < 0.0000001:
      zeroified.append(-1.0)
    elif abs(value-0.5) < 0.0000001:
      zeroified.append(0.5)
    elif abs(value+0.5) < 0.0000001:
      zeroified.append(-0.5)
    else:
      zeroified.append(value)
  return zeroified

def _rotation_c2_xy(theta):
  return [
    math.cos(2.0*theta), math.sin(2.0*theta), 0.0,
    math.sin(2.0*theta), -math.cos(2.0*theta), 0.0,
    0.0, 0.0, -1.0
  ]

def _icosahedral_phi():
  return 0.5*(1.0+math.sqrt(5.0))

def _icosahedral_axes_c2():
  phi=_icosahedral_phi()
  axes=[]
  for sign_1 in [1.0, -1.0]:
    for sign_2 in [1.0, -1.0]:
      axes.append([0.0, sign_1, sign_2*phi])
      axes.append([sign_1, sign_2*phi, 0.0])
      axes.append([sign_2*phi, 0.0, sign_1])
  return axes

def _icosahedral_axes_c3():
  axes=[]
  for x in [1.0, -1.0]:
    for y in [1.0, -1.0]:
      for z in [1.0, -1.0]:
        axes.append([x, y, z])
  return axes

def _icosahedral_axes_c5():
  phi=_icosahedral_phi()
  axes=[]
  for sign_1 in [1.0, -1.0]:
    for sign_2 in [1.0, -1.0]:
      axes.append([0.0, sign_1/phi, sign_2*phi])
      axes.append([sign_1/phi, sign_2*phi, 0.0])
      axes.append([sign_2*phi, 0.0, sign_1/phi])
  return axes

def _unique_rotation_matrices(matrices):
  unique=[]
  seen={}
  for matrix in matrices:
    zeroified=_zeroify_rotation_matrix(matrix)
    key=tuple([round(value, 6) for value in zeroified])
    if not seen.has_key(key):
      seen[key]=True
      unique.append(zeroified)
  return unique

def _generate_point_group_rotation_matrices(point_group_symbol):
  symbol=point_group_symbol.strip().upper()
  if not symbol:
    return None

  if symbol.startswith("C") and len(symbol) > 1 and symbol[1:].isdigit():
    order=int(symbol[1:])
    if order < 2:
      return None
    matrices=[]
    for index in range(1, order):
      matrices.append(_rotation_matrix_about_axis([0.0, 0.0, 1.0], (2.0*math.pi*index)/float(order)))
    return _unique_rotation_matrices(matrices)

  if symbol.startswith("D") and len(symbol) > 1 and symbol[1:].isdigit():
    order=int(symbol[1:])
    if order < 2:
      return None
    matrices=[]
    for index in range(1, order):
      matrices.append(_rotation_matrix_about_axis([0.0, 0.0, 1.0], (2.0*math.pi*index)/float(order)))
    for index in range(order):
      matrices.append(_rotation_c2_xy((math.pi*index)/float(order)))
    return _unique_rotation_matrices(matrices)

  if symbol == "T":
    axes=[[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0],
          [1.0,1.0,1.0], [-1.0,-1.0,1.0], [1.0,-1.0,-1.0], [-1.0,1.0,-1.0]]
    matrices=[]
    for axis in axes[:3]:
      matrices.append(_rotation_matrix_about_axis(axis, math.pi))
    for axis in axes[3:]:
      matrices.append(_rotation_matrix_about_axis(axis, 2.0*math.pi/3.0))
      matrices.append(_rotation_matrix_about_axis(axis, 4.0*math.pi/3.0))
    return _unique_rotation_matrices(matrices)

  if symbol == "O":
    axes_c4=[[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
    axes_c3=[[1.0,1.0,1.0], [-1.0,-1.0,1.0], [1.0,-1.0,-1.0], [-1.0,1.0,-1.0]]
    axes_c2=[[1.0,1.0,0.0], [1.0,-1.0,0.0], [1.0,0.0,1.0],
             [1.0,0.0,-1.0], [0.0,1.0,1.0], [0.0,1.0,-1.0]]
    matrices=[]
    for axis in axes_c4:
      matrices.append(_rotation_matrix_about_axis(axis, math.pi/2.0))
      matrices.append(_rotation_matrix_about_axis(axis, math.pi))
      matrices.append(_rotation_matrix_about_axis(axis, 3.0*math.pi/2.0))
    for axis in axes_c3:
      matrices.append(_rotation_matrix_about_axis(axis, 2.0*math.pi/3.0))
      matrices.append(_rotation_matrix_about_axis(axis, 4.0*math.pi/3.0))
    for axis in axes_c2:
      matrices.append(_rotation_matrix_about_axis(axis, math.pi))
    return _unique_rotation_matrices(matrices)

  if symbol == "I":
    matrices=[]
    for axis in _icosahedral_axes_c2():
      matrices.append(_rotation_matrix_about_axis(axis, math.pi))
    for axis in _icosahedral_axes_c3():
      matrices.append(_rotation_matrix_about_axis(axis, 2.0*math.pi/3.0))
      matrices.append(_rotation_matrix_about_axis(axis, 4.0*math.pi/3.0))
    for axis in _icosahedral_axes_c5():
      matrices.append(_rotation_matrix_about_axis(axis, 2.0*math.pi/5.0))
      matrices.append(_rotation_matrix_about_axis(axis, 4.0*math.pi/5.0))
      matrices.append(_rotation_matrix_about_axis(axis, 6.0*math.pi/5.0))
      matrices.append(_rotation_matrix_about_axis(axis, 8.0*math.pi/5.0))
    return _unique_rotation_matrices(matrices)

  return None

def _strict_ncs_translation_about_origin(rotation_matrix, origin_xyz):
  ox=origin_xyz[0]
  oy=origin_xyz[1]
  oz=origin_xyz[2]
  rx=rotation_matrix[0]*ox + rotation_matrix[1]*oy + rotation_matrix[2]*oz
  ry=rotation_matrix[3]*ox + rotation_matrix[4]*oy + rotation_matrix[5]*oz
  rz=rotation_matrix[6]*ox + rotation_matrix[7]*oy + rotation_matrix[8]*oz
  return [ox-rx, oy-ry, oz-rz]

def _cell_centre_xyz(mol_id):
  cell_info=cell(mol_id)
  return [0.5*float(cell_info[0]), 0.5*float(cell_info[1]), 0.5*float(cell_info[2])]

def _apply_point_group_molecular_symmetry(mol_id, point_group_symbol, origin_xyz):
  global MOLECULAR_SYMMETRY_SESSION_STATE
  normalized_symbol=point_group_symbol.strip().upper()
  matrices=_generate_point_group_rotation_matrices(normalized_symbol)
  if not matrices:
    info_dialog("Unsupported symmetry group.\n\nSupported forms are Cn, Dn, T, O and I.\nExamples: C4, C13, D7, O, I")
    return None

  prior_state=MOLECULAR_SYMMETRY_SESSION_STATE.get(("metadata-display", mol_id))
  if prior_state:
    set_show_symmetry_master(1)
    _status_message("Molecular symmetry copies already loaded")
    return None

  added_count=0
  for rotation_matrix in matrices:
    if _matrix_is_close_to_identity(rotation_matrix, [0.0, 0.0, 0.0]):
      continue
    if _add_molecular_symmetry_operator(mol_id, rotation_matrix, origin_xyz):
      added_count=added_count+1
  if added_count == 0:
    info_dialog("Only the identity operator was generated for %s" % normalized_symbol)
    return None
  try:
    set_show_symmetry_molecule(mol_id, 1)
  except:
    pass
  set_show_symmetry_master(1)
  session_state={
    "file_name": molecule_name(mol_id),
    "parser": "Manual point group about cell centre",
    "operators_added": added_count,
    "symbol": normalized_symbol,
    "origin": tuple(origin_xyz)
  }
  MOLECULAR_SYMMETRY_SESSION_STATE[("metadata-display", mol_id)]=session_state
  _status_message("Loaded %s molecular symmetry copies for %s" % (added_count, normalized_symbol))
  return session_state

def _manual_point_group_entry_callback_for_cell_centre(mol_id):
  def _callback(text):
    if not text:
      info_dialog("Enter a symmetry group such as C4 or C13")
      return None
    origin_xyz=_cell_centre_xyz(mol_id)
    return _apply_point_group_molecular_symmetry(mol_id, text, origin_xyz)
  return _callback

def _prompt_manual_point_group_about_cell_centre(mol_id):
  generic_single_entry("No assembly operators were found.\n\nEnter a point group (e.g. C13).\nAssumption: symmetry centre at unit-cell centre,\nwith the principal Cn axis along z.",
                       "C13",
                       "Manual molecular symmetry",
                       _manual_point_group_entry_callback_for_cell_centre(mol_id))

def _apply_point_group_strict_ncs(mol_id, point_group_symbol):
  global MOLECULAR_SYMMETRY_SESSION_STATE
  normalized_symbol=point_group_symbol.strip().upper()
  matrices=_generate_point_group_rotation_matrices(normalized_symbol)
  if not matrices:
    info_dialog("Unsupported symmetry group.\n\nSupported forms are Cn, Dn, T, O and I.\nExamples: C4, C16, D7, O, I")
    return None

  origin_xyz=_rotation_centre_xyz()
  state_key=mol_id
  prior_state=MOLECULAR_SYMMETRY_SESSION_STATE.get(state_key)
  origin_key=(round(origin_xyz[0], 3), round(origin_xyz[1], 3), round(origin_xyz[2], 3))

  if prior_state:
    if prior_state["symbol"] == normalized_symbol and prior_state["origin"] == origin_key:
      set_show_strict_ncs(mol_id, 1)
      _status_message("Molecular symmetry already active for %s" % normalized_symbol)
      return matrices
    info_dialog("Molecular symmetry has already been activated for this molecule in this session.\n\nOld Coot does not expose a way to clear strict NCS matrices, so to change the symmetry group or origin you should reload the molecule.")
    return None

  for rotation_matrix in matrices:
    translation=_strict_ncs_translation_about_origin(rotation_matrix, origin_xyz)
    add_strict_ncs_matrix(mol_id, "A", "A",
                          rotation_matrix[0], rotation_matrix[1], rotation_matrix[2],
                          rotation_matrix[3], rotation_matrix[4], rotation_matrix[5],
                          rotation_matrix[6], rotation_matrix[7], rotation_matrix[8],
                          translation[0], translation[1], translation[2])

  set_show_strict_ncs(mol_id, 1)
  set_show_symmetry_master(1)
  MOLECULAR_SYMMETRY_SESSION_STATE[state_key]={
    "symbol": normalized_symbol,
    "origin": origin_key,
    "count": len(matrices)
  }
  _status_message("Activated %s molecular symmetry (%s operators)" % (normalized_symbol, len(matrices)))
  return matrices

def activate_molecular_symmetry_from_entry(text):
  residue=_active_residue_or_status()
  if not residue:
    return None
  if not text:
    info_dialog("Enter a symmetry group such as C4, C16, D7, O or I")
    return None
  return _apply_point_group_strict_ncs(residue[0], text)

def prompt_activate_molecular_symmetry():
  generic_single_entry("Point group symmetry about current rotation centre (e.g. C4, C16, D7, O, I)",
                       "C4",
                       "Activate molecular symmetry",
                       activate_molecular_symmetry_from_entry)

def _matrix_is_close_to_identity(rotation_matrix, translation_vector):
  if not rotation_matrix or not translation_vector:
    return False
  identity=[1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]
  for index in range(9):
    if abs(rotation_matrix[index]-identity[index]) > 0.00001:
      return False
  for value in translation_vector:
    if abs(value) > 0.00001:
      return False
  return True

def _active_molecule_source_file_or_status():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  file_name=molecule_name(mol_id)
  if not file_name:
    info_dialog("Active molecule does not have an associated file name")
    return None
  return (mol_id, file_name)

def _parse_pdb_biomt_operators(file_name):
  import os
  import re
  if not os.path.isfile(file_name):
    return None
  file_handle=open(file_name, "r")
  try:
    lines=file_handle.readlines()
  finally:
    file_handle.close()

  current_chain_ids=None
  operators={}
  biomt_row_regex=re.compile(r"^REMARK 350\s+BIOMT([123])\s+(\d+)\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)")

  for raw_line in lines:
    line=raw_line.rstrip("\n")
    if line.startswith("REMARK 350") and ("APPLY THE FOLLOWING TO CHAINS:" in line or "AND CHAINS:" in line):
      chain_text=line.split("CHAINS:",1)[1].strip()
      chain_text=chain_text.replace(";", "")
      chain_ids=[piece.strip() for piece in chain_text.split(",") if piece.strip()]
      if not current_chain_ids:
        current_chain_ids=[]
      for chain_id in chain_ids:
        if chain_id not in current_chain_ids:
          current_chain_ids.append(chain_id)
      continue
    match=biomt_row_regex.match(line)
    if not match:
      continue
    row_index=int(match.group(1))
    operator_id=int(match.group(2))
    values=[float(match.group(3)), float(match.group(4)), float(match.group(5)), float(match.group(6))]
    if not operators.has_key(operator_id):
      operators[operator_id]={"rows": {}, "chain_ids": list(current_chain_ids or [])}
    if current_chain_ids:
      operators[operator_id]["chain_ids"]=list(current_chain_ids)
    operators[operator_id]["rows"][row_index]=values

  parsed=[]
  for operator_id in sorted(operators.keys()):
    operator=operators[operator_id]
    rows=operator["rows"]
    if not (rows.has_key(1) and rows.has_key(2) and rows.has_key(3)):
      continue
    rotation_matrix=[
      rows[1][0], rows[1][1], rows[1][2],
      rows[2][0], rows[2][1], rows[2][2],
      rows[3][0], rows[3][1], rows[3][2]
    ]
    translation_vector=[rows[1][3], rows[2][3], rows[3][3]]
    parsed.append({
      "source": "BIOMT",
      "id": operator_id,
      "chain_ids": operator["chain_ids"],
      "rotation": _zeroify_rotation_matrix(rotation_matrix),
      "translation": translation_vector
    })
  return parsed

def _tokenize_cif_row_text(text):
  import shlex
  lexer=shlex.shlex(text, posix=True)
  lexer.whitespace_split=True
  lexer.commenters=''
  return list(lexer)

def _parse_cif_loop_rows(lines, start_index):
  headers=[]
  row_lines=[]
  index=start_index
  while index < len(lines):
    stripped=lines[index].strip()
    if stripped.startswith("_"):
      headers.append(stripped)
      index=index+1
      continue
    break
  while index < len(lines):
    stripped=lines[index].strip()
    if not stripped:
      index=index+1
      continue
    if stripped == "#":
      index=index+1
      break
    if stripped.startswith("loop_") or stripped.startswith("_"):
      break
    row_lines.append(stripped)
    index=index+1
  tokens=_tokenize_cif_row_text(" ".join(row_lines))
  if not headers:
    return ([], [], index)
  n_columns=len(headers)
  rows=[]
  current=[]
  for token in tokens:
    current.append(token)
    if len(current) == n_columns:
      rows.append(current)
      current=[]
  return (headers, rows, index)

def _parse_cif_struct_ncs_oper_operators(file_name):
  import os
  if not os.path.isfile(file_name):
    return None
  file_handle=open(file_name, "r")
  try:
    lines=file_handle.readlines()
  finally:
    file_handle.close()

  index=0
  while index < len(lines):
    if lines[index].strip() != "loop_":
      index=index+1
      continue
    headers, rows, next_index=_parse_cif_loop_rows(lines, index+1)
    index=next_index
    if "_struct_ncs_oper.id" not in headers:
      continue
    header_index={}
    for header_position in range(len(headers)):
      header_index[headers[header_position]]=header_position
    required_headers=[
      "_struct_ncs_oper.id",
      "_struct_ncs_oper.matrix[1][1]",
      "_struct_ncs_oper.matrix[1][2]",
      "_struct_ncs_oper.matrix[1][3]",
      "_struct_ncs_oper.matrix[2][1]",
      "_struct_ncs_oper.matrix[2][2]",
      "_struct_ncs_oper.matrix[2][3]",
      "_struct_ncs_oper.matrix[3][1]",
      "_struct_ncs_oper.matrix[3][2]",
      "_struct_ncs_oper.matrix[3][3]",
      "_struct_ncs_oper.vector[1]",
      "_struct_ncs_oper.vector[2]",
      "_struct_ncs_oper.vector[3]"
    ]
    missing_required=[header for header in required_headers if not header_index.has_key(header)]
    if missing_required:
      return None
    parsed=[]
    for row in rows:
      rotation_matrix=[
        float(row[header_index["_struct_ncs_oper.matrix[1][1]"]]),
        float(row[header_index["_struct_ncs_oper.matrix[1][2]"]]),
        float(row[header_index["_struct_ncs_oper.matrix[1][3]"]]),
        float(row[header_index["_struct_ncs_oper.matrix[2][1]"]]),
        float(row[header_index["_struct_ncs_oper.matrix[2][2]"]]),
        float(row[header_index["_struct_ncs_oper.matrix[2][3]"]]),
        float(row[header_index["_struct_ncs_oper.matrix[3][1]"]]),
        float(row[header_index["_struct_ncs_oper.matrix[3][2]"]]),
        float(row[header_index["_struct_ncs_oper.matrix[3][3]"]])
      ]
      translation_vector=[
        float(row[header_index["_struct_ncs_oper.vector[1]"]]),
        float(row[header_index["_struct_ncs_oper.vector[2]"]]),
        float(row[header_index["_struct_ncs_oper.vector[3]"]])
      ]
      parsed.append({
        "source": "_struct_ncs_oper",
        "id": row[header_index["_struct_ncs_oper.id"]],
        "rotation": _zeroify_rotation_matrix(rotation_matrix),
        "translation": translation_vector
      })
    return parsed
  return None

def _parse_cif_pdbx_struct_oper_list(file_name):
  import os
  if not os.path.isfile(file_name):
    return None
  file_handle=open(file_name, "r")
  try:
    lines=file_handle.readlines()
  finally:
    file_handle.close()

  index=0
  while index < len(lines):
    if lines[index].strip() != "loop_":
      index=index+1
      continue
    headers, rows, next_index=_parse_cif_loop_rows(lines, index+1)
    index=next_index
    if "_pdbx_struct_oper_list.id" not in headers:
      continue
    header_index={}
    for header_position in range(len(headers)):
      header_index[headers[header_position]]=header_position
    required_headers=[
      "_pdbx_struct_oper_list.id",
      "_pdbx_struct_oper_list.matrix[1][1]",
      "_pdbx_struct_oper_list.matrix[1][2]",
      "_pdbx_struct_oper_list.matrix[1][3]",
      "_pdbx_struct_oper_list.matrix[2][1]",
      "_pdbx_struct_oper_list.matrix[2][2]",
      "_pdbx_struct_oper_list.matrix[2][3]",
      "_pdbx_struct_oper_list.matrix[3][1]",
      "_pdbx_struct_oper_list.matrix[3][2]",
      "_pdbx_struct_oper_list.matrix[3][3]",
      "_pdbx_struct_oper_list.vector[1]",
      "_pdbx_struct_oper_list.vector[2]",
      "_pdbx_struct_oper_list.vector[3]"
    ]
    missing_required=[header for header in required_headers if not header_index.has_key(header)]
    if missing_required:
      return None
    parsed={}
    for row in rows:
      parsed[row[header_index["_pdbx_struct_oper_list.id"]]]={
        "source": "_pdbx_struct_oper_list",
        "id": row[header_index["_pdbx_struct_oper_list.id"]],
        "rotation": _zeroify_rotation_matrix([
          float(row[header_index["_pdbx_struct_oper_list.matrix[1][1]"]]),
          float(row[header_index["_pdbx_struct_oper_list.matrix[1][2]"]]),
          float(row[header_index["_pdbx_struct_oper_list.matrix[1][3]"]]),
          float(row[header_index["_pdbx_struct_oper_list.matrix[2][1]"]]),
          float(row[header_index["_pdbx_struct_oper_list.matrix[2][2]"]]),
          float(row[header_index["_pdbx_struct_oper_list.matrix[2][3]"]]),
          float(row[header_index["_pdbx_struct_oper_list.matrix[3][1]"]]),
          float(row[header_index["_pdbx_struct_oper_list.matrix[3][2]"]]),
          float(row[header_index["_pdbx_struct_oper_list.matrix[3][3]"]])
        ]),
        "translation": [
          float(row[header_index["_pdbx_struct_oper_list.vector[1]"]]),
          float(row[header_index["_pdbx_struct_oper_list.vector[2]"]]),
          float(row[header_index["_pdbx_struct_oper_list.vector[3]"]])
        ]
      }
    return parsed
  return None

def _parse_cif_pdbx_struct_assembly_gen_rows(file_name):
  import os
  if not os.path.isfile(file_name):
    return None
  file_handle=open(file_name, "r")
  try:
    lines=file_handle.readlines()
  finally:
    file_handle.close()

  index=0
  while index < len(lines):
    if lines[index].strip() != "loop_":
      index=index+1
      continue
    headers, rows, next_index=_parse_cif_loop_rows(lines, index+1)
    index=next_index
    if "_pdbx_struct_assembly_gen.assembly_id" not in headers:
      continue
    header_index={}
    for header_position in range(len(headers)):
      header_index[headers[header_position]]=header_position
    required_headers=[
      "_pdbx_struct_assembly_gen.assembly_id",
      "_pdbx_struct_assembly_gen.oper_expression",
      "_pdbx_struct_assembly_gen.asym_id_list"
    ]
    missing_required=[header for header in required_headers if not header_index.has_key(header)]
    if missing_required:
      return None
    parsed=[]
    for row in rows:
      parsed.append({
        "assembly_id": row[header_index["_pdbx_struct_assembly_gen.assembly_id"]],
        "oper_expression": row[header_index["_pdbx_struct_assembly_gen.oper_expression"]],
        "asym_id_list": row[header_index["_pdbx_struct_assembly_gen.asym_id_list"]]
      })
    return parsed

  simple_values={}
  for raw_line in lines:
    stripped=raw_line.strip()
    if not stripped.startswith("_pdbx_struct_assembly_gen."):
      continue
    pieces=stripped.split(None, 1)
    if len(pieces) < 2:
      continue
    simple_values[pieces[0]]=pieces[1]
  if simple_values.has_key("_pdbx_struct_assembly_gen.assembly_id") and \
     simple_values.has_key("_pdbx_struct_assembly_gen.oper_expression") and \
     simple_values.has_key("_pdbx_struct_assembly_gen.asym_id_list"):
    return [{
      "assembly_id": simple_values["_pdbx_struct_assembly_gen.assembly_id"],
      "oper_expression": simple_values["_pdbx_struct_assembly_gen.oper_expression"],
      "asym_id_list": simple_values["_pdbx_struct_assembly_gen.asym_id_list"]
    }]
  return None

def _expand_operator_token_list(expression_text):
  tokens=[]
  for piece in expression_text.split(","):
    token=piece.strip()
    if not token or token == "?":
      continue
    if "-" in token:
      range_parts=token.split("-", 1)
      try:
        start_value=int(range_parts[0])
        end_value=int(range_parts[1])
      except:
        tokens.append(token)
        continue
      step=1
      if end_value < start_value:
        step=-1
      for value in range(start_value, end_value+step, step):
        tokens.append(str(value))
    else:
      tokens.append(token)
  return tokens

def _expand_cif_oper_expression(oper_expression):
  import re
  if not oper_expression:
    return []
  expression=oper_expression.replace(" ", "")
  grouped_parts=re.findall(r"\(([^()]*)\)", expression)
  if not grouped_parts:
    grouped_parts=[expression]
  expanded_groups=[_expand_operator_token_list(group_text) for group_text in grouped_parts]
  if not expanded_groups:
    return []
  sequences=[[]]
  for group_tokens in expanded_groups:
    new_sequences=[]
    for sequence in sequences:
      for token in group_tokens:
        new_sequences.append(sequence+[token])
    sequences=new_sequences
  return sequences

def _matrix_multiply_3x3(left_matrix, right_matrix):
  return [
    left_matrix[0]*right_matrix[0] + left_matrix[1]*right_matrix[3] + left_matrix[2]*right_matrix[6],
    left_matrix[0]*right_matrix[1] + left_matrix[1]*right_matrix[4] + left_matrix[2]*right_matrix[7],
    left_matrix[0]*right_matrix[2] + left_matrix[1]*right_matrix[5] + left_matrix[2]*right_matrix[8],
    left_matrix[3]*right_matrix[0] + left_matrix[4]*right_matrix[3] + left_matrix[5]*right_matrix[6],
    left_matrix[3]*right_matrix[1] + left_matrix[4]*right_matrix[4] + left_matrix[5]*right_matrix[7],
    left_matrix[3]*right_matrix[2] + left_matrix[4]*right_matrix[5] + left_matrix[5]*right_matrix[8],
    left_matrix[6]*right_matrix[0] + left_matrix[7]*right_matrix[3] + left_matrix[8]*right_matrix[6],
    left_matrix[6]*right_matrix[1] + left_matrix[7]*right_matrix[4] + left_matrix[8]*right_matrix[7],
    left_matrix[6]*right_matrix[2] + left_matrix[7]*right_matrix[5] + left_matrix[8]*right_matrix[8]
  ]

def _compose_affine_operator_sequence(operator_ids, operator_lookup):
  rotation_matrix=[1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]
  translation_vector=[0.0,0.0,0.0]
  for operator_id in operator_ids:
    if not operator_lookup.has_key(operator_id):
      return None
    operator=operator_lookup[operator_id]
    rotation_matrix=_matrix_multiply_3x3(operator["rotation"], rotation_matrix)
    rotated_translation=_matrix_vector_product_3x3(operator["rotation"], translation_vector)
    translation_vector=[
      rotated_translation[0]+operator["translation"][0],
      rotated_translation[1]+operator["translation"][1],
      rotated_translation[2]+operator["translation"][2]
    ]
  return {
    "source": "_pdbx_struct_assembly_gen",
    "id": ",".join(operator_ids),
    "rotation": _zeroify_rotation_matrix(rotation_matrix),
    "translation": translation_vector
  }

def _parse_cif_assembly_operators(file_name, target_chain_ids):
  assembly_rows=_parse_cif_pdbx_struct_assembly_gen_rows(file_name)
  operator_lookup=_parse_cif_pdbx_struct_oper_list(file_name)
  if not assembly_rows or not operator_lookup:
    return None

  matching_rows=[]
  matching_assembly_id=None
  for row in assembly_rows:
    row_chain_ids=[piece.strip() for piece in row["asym_id_list"].split(",") if piece.strip()]
    if target_chain_ids and not [chain_id for chain_id in target_chain_ids if chain_id in row_chain_ids]:
      continue
    if matching_assembly_id is None:
      matching_assembly_id=row["assembly_id"]
    if row["assembly_id"] == matching_assembly_id:
      matching_rows.append(row)
  if not matching_rows:
    matching_assembly_id=assembly_rows[0]["assembly_id"]
    matching_rows=[row for row in assembly_rows if row["assembly_id"] == matching_assembly_id]

  parsed=[]
  seen_ids={}
  for row in matching_rows:
    operator_sequences=_expand_cif_oper_expression(row["oper_expression"])
    for operator_ids in operator_sequences:
      composed_operator=_compose_affine_operator_sequence(operator_ids, operator_lookup)
      if not composed_operator:
        continue
      if seen_ids.has_key(composed_operator["id"]):
        continue
      seen_ids[composed_operator["id"]]=True
      parsed.append(composed_operator)
  return parsed

def _molecule_chain_ids_for_symmetry(mol_id):
  try:
    chains=chain_ids(mol_id)
  except:
    chains=[]
  if not chains:
    residue=_active_residue_or_status()
    if residue and residue[0] == mol_id:
      return [residue[1]]
  return chains

def _active_chain_id_for_symmetry_source(mol_id):
  residue=_active_residue_or_status()
  if residue and residue[0] == mol_id:
    return residue[1]
  chains=_molecule_chain_ids_for_symmetry(mol_id)
  if len(chains) == 1:
    return chains[0]
  return None

def _matrix_vector_product_3x3(rotation_matrix, vector):
  return [
    rotation_matrix[0]*vector[0] + rotation_matrix[1]*vector[1] + rotation_matrix[2]*vector[2],
    rotation_matrix[3]*vector[0] + rotation_matrix[4]*vector[1] + rotation_matrix[5]*vector[2],
    rotation_matrix[6]*vector[0] + rotation_matrix[7]*vector[1] + rotation_matrix[8]*vector[2]
  ]

def _rotation_axis_from_matrix(rotation_matrix):
  axis=[
    rotation_matrix[7]-rotation_matrix[5],
    rotation_matrix[2]-rotation_matrix[6],
    rotation_matrix[3]-rotation_matrix[1]
  ]
  axis_length=math.sqrt(_distance_sq(axis, [0.0, 0.0, 0.0]))
  if axis_length > 0.000001:
    return [component/axis_length for component in axis]
  rows=[
    [rotation_matrix[0]+1.0, rotation_matrix[1], rotation_matrix[2]],
    [rotation_matrix[3], rotation_matrix[4]+1.0, rotation_matrix[5]],
    [rotation_matrix[6], rotation_matrix[7], rotation_matrix[8]+1.0]
  ]
  best_row=None
  best_norm=0.0
  for row in rows:
    row_norm=math.sqrt(_distance_sq(row, [0.0, 0.0, 0.0]))
    if row_norm > best_norm:
      best_row=row
      best_norm=row_norm
  if best_row and best_norm > 0.000001:
    return [component/best_norm for component in best_row]
  return [0.0, 0.0, 1.0]

def _solve_3x3_linear_system(matrix_3x3, rhs_vector):
  augmented=[
    [float(matrix_3x3[0][0]), float(matrix_3x3[0][1]), float(matrix_3x3[0][2]), float(rhs_vector[0])],
    [float(matrix_3x3[1][0]), float(matrix_3x3[1][1]), float(matrix_3x3[1][2]), float(rhs_vector[1])],
    [float(matrix_3x3[2][0]), float(matrix_3x3[2][1]), float(matrix_3x3[2][2]), float(rhs_vector[2])]
  ]
  for pivot_index in range(3):
    best_row=pivot_index
    best_value=abs(augmented[pivot_index][pivot_index])
    for row_index in range(pivot_index+1, 3):
      value=abs(augmented[row_index][pivot_index])
      if value > best_value:
        best_row=row_index
        best_value=value
    if best_value < 0.00000001:
      return None
    if best_row != pivot_index:
      augmented[pivot_index], augmented[best_row]=augmented[best_row], augmented[pivot_index]
    pivot_value=augmented[pivot_index][pivot_index]
    for column_index in range(pivot_index, 4):
      augmented[pivot_index][column_index]=augmented[pivot_index][column_index]/pivot_value
    for row_index in range(3):
      if row_index == pivot_index:
        continue
      scale=augmented[row_index][pivot_index]
      if abs(scale) < 0.00000001:
        continue
      for column_index in range(pivot_index, 4):
        augmented[row_index][column_index]=augmented[row_index][column_index]-scale*augmented[pivot_index][column_index]
  return [augmented[0][3], augmented[1][3], augmented[2][3]]

def _origin_from_rotation_and_translation(rotation_matrix, translation_vector):
  axis=_rotation_axis_from_matrix(rotation_matrix)
  affine_rows=[
    [1.0-rotation_matrix[0],   -rotation_matrix[1],   -rotation_matrix[2]],
    [  -rotation_matrix[3], 1.0-rotation_matrix[4],   -rotation_matrix[5]],
    [  -rotation_matrix[6],   -rotation_matrix[7], 1.0-rotation_matrix[8]]
  ]
  affine_rhs=[translation_vector[0], translation_vector[1], translation_vector[2]]
  origin_vector=None
  for row_pair in ((0,1), (0,2), (1,2)):
    system_matrix=[
      affine_rows[row_pair[0]],
      affine_rows[row_pair[1]],
      [axis[0], axis[1], axis[2]]
    ]
    system_rhs=[affine_rhs[row_pair[0]], affine_rhs[row_pair[1]], 0.0]
    trial_origin=_solve_3x3_linear_system(system_matrix, system_rhs)
    if trial_origin is None:
      continue
    origin_vector=trial_origin
    break
  if origin_vector is None:
    return None
  residual_matrix=[
    1.0-rotation_matrix[0],   -rotation_matrix[1],   -rotation_matrix[2],
      -rotation_matrix[3], 1.0-rotation_matrix[4],   -rotation_matrix[5],
      -rotation_matrix[6],   -rotation_matrix[7], 1.0-rotation_matrix[8]
  ]
  predicted_translation=_matrix_vector_product_3x3(residual_matrix, origin_vector)
  residual=[
    predicted_translation[0]-translation_vector[0],
    predicted_translation[1]-translation_vector[1],
    predicted_translation[2]-translation_vector[2]
  ]
  if math.sqrt(_distance_sq(residual, [0.0, 0.0, 0.0])) > 0.01:
    return None
  return origin_vector

def _add_molecular_symmetry_operator(mol_id, rotation_matrix, origin_vector):
  try:
    add_molecular_symmetry(mol_id,
                           rotation_matrix[0], rotation_matrix[1], rotation_matrix[2],
                           rotation_matrix[3], rotation_matrix[4], rotation_matrix[5],
                           rotation_matrix[6], rotation_matrix[7], rotation_matrix[8],
                           origin_vector[0], origin_vector[1], origin_vector[2])
    return True
  except NameError:
    pass
  try:
    import _coot
    _coot.add_molecular_symmetry(mol_id,
                                 rotation_matrix[0], rotation_matrix[1], rotation_matrix[2],
                                 rotation_matrix[3], rotation_matrix[4], rotation_matrix[5],
                                 rotation_matrix[6], rotation_matrix[7], rotation_matrix[8],
                                 origin_vector[0], origin_vector[1], origin_vector[2])
    return True
  except:
    return False

def _load_display_molecular_symmetry_from_metadata():
  global MOLECULAR_SYMMETRY_SESSION_STATE
  active_info=_active_molecule_source_file_or_status()
  if not active_info:
    return None
  mol_id=active_info[0]
  file_name=active_info[1]
  lower_name=file_name.lower()

  prior_state=MOLECULAR_SYMMETRY_SESSION_STATE.get(("metadata-display", mol_id))
  if prior_state:
    set_show_symmetry_master(1)
    _status_message("Molecular symmetry copies already loaded")
    return prior_state

  if lower_name.endswith(".pdb") or lower_name.endswith(".ent"):
    operators=_parse_pdb_biomt_operators(file_name)
    parser_name="PDB BIOMT"
  elif lower_name.endswith(".cif") or lower_name.endswith(".mmcif"):
    target_chain_ids=_molecule_chain_ids_for_symmetry(mol_id)
    operators=_parse_cif_assembly_operators(file_name, target_chain_ids)
    parser_name="mmCIF assembly metadata"
    if not operators:
      operators=_parse_cif_struct_ncs_oper_operators(file_name)
      parser_name="mmCIF _struct_ncs_oper"
  else:
    info_dialog("Only PDB and mmCIF files are currently supported for molecular symmetry loading")
    return None

  if not operators:
    _prompt_manual_point_group_about_cell_centre(mol_id)
    return None

  target_chain_ids=_molecule_chain_ids_for_symmetry(mol_id)
  if not target_chain_ids:
    info_dialog("Could not determine chain IDs for the active molecule")
    return None

  added_count=0
  for operator in operators:
    rotation_matrix=operator["rotation"]
    translation_vector=operator["translation"]
    if _matrix_is_close_to_identity(rotation_matrix, translation_vector):
      continue
    if operator.get("source") == "BIOMT" or operator.get("source") == "_pdbx_struct_assembly_gen":
      origin_vector=_origin_from_rotation_and_translation(rotation_matrix, translation_vector)
      if origin_vector is not None:
        if _add_molecular_symmetry_operator(mol_id, rotation_matrix, origin_vector):
          added_count=added_count+1
          continue

  if added_count == 0:
    info_dialog("Only identity symmetry operators were found for the active molecule")
    return None

  try:
    set_show_symmetry_molecule(mol_id, 1)
  except:
    pass
  set_show_symmetry_master(1)
  session_state={
    "file_name": file_name,
    "parser": parser_name,
    "operators_added": added_count
  }
  MOLECULAR_SYMMETRY_SESSION_STATE[("metadata-display", mol_id)]=session_state
  _status_message("Loaded %s molecular symmetry operators from %s" % (added_count, parser_name))
  return session_state

def decrease_map_radius_with_status():
  current_radius=get_map_radius()
  new_radius=current_radius-2.0
  if new_radius < 2.0:
    new_radius=2.0
  set_map_radius(new_radius)
  try:
    set_map_radius_em(new_radius)
  except:
    pass
  _status_message("Map radius %.1f A" % new_radius)

def increase_map_radius_with_status():
  current_radius=get_map_radius()
  new_radius=current_radius+2.0
  if new_radius > 1000.0:
    new_radius=1000.0
  set_map_radius(new_radius)
  try:
    set_map_radius_em(new_radius)
  except:
    pass
  _status_message("Map radius %.1f A" % new_radius)

def stepped_sphere_refine_active_chain():
  residue=_active_residue_or_status()
  if not residue:
    return None
  return stepped_sphere_refine(residue[0], residue[1])

def color_rotamer_outliers_and_missing_atoms_for_active_molecule():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  return color_rotamer_outliers_and_missing_atoms(mol_id)

def color_polars_and_hphobs_for_active_molecule():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  return color_polars_and_hphobs(mol_id)

def color_by_charge_for_active_molecule():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  return color_by_charge(mol_id)

def color_protein_na_for_active_molecule():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  return color_protein_na(mol_id)

def color_waters_for_active_molecule():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  return color_waters(mol_id)

def _read_ccp4_mrc_grid_info(file_name):
  import os
  import struct
  if not file_name or not os.path.isfile(file_name):
    return None
  header_file=open(file_name, "rb")
  try:
    header=header_file.read(224)
  finally:
    header_file.close()
  if len(header) < 52:
    return None
  parsed=None
  for endian in ["<", ">"]:
    try:
      ints=struct.unpack(endian+"10i", header[:40])
      floats=struct.unpack(endian+"3f", header[40:52])
    except struct.error:
      continue
    mx=ints[7]
    my=ints[8]
    mz=ints[9]
    xlen=floats[0]
    ylen=floats[1]
    zlen=floats[2]
    if mx > 0 and my > 0 and mz > 0 and xlen > 0.0 and ylen > 0.0 and zlen > 0.0:
      parsed={
        "spacing_xyz": (xlen/float(mx), ylen/float(my), zlen/float(mz)),
        "grid_xyz": (mx, my, mz),
        "cell_xyz": (xlen, ylen, zlen)
      }
      break
  if not parsed:
    return None
  return parsed

def _read_ccp4_mrc_grid_spacing(file_name):
  grid_info=_read_ccp4_mrc_grid_info(file_name)
  if not grid_info:
    return None
  return min(grid_info["spacing_xyz"])

def _restyle_em_map_mesh(map_id, contour_sigma):
  if map_is_difference_map(map_id)!=0:
    return None
  set_draw_solid_density_surface(map_id,0)
  set_draw_map_standard_lines(map_id,1)
  set_map_colour(map_id,EM_GLOBAL_MAP_COLOUR[0],EM_GLOBAL_MAP_COLOUR[1],EM_GLOBAL_MAP_COLOUR[2])
  set_contour_level_in_sigma(map_id, contour_sigma)

def place_water_in_active_molecule():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  try:
    set_pointer_atom_molecule(mol_id)
  except:
    pass
  try:
    set_go_to_atom_molecule(mol_id)
  except:
    pass
  place_typed_atom_at_pointer("Water")
  return mol_id

def add_water_and_refine():
  mol_id=place_water_in_active_molecule()
  if mol_id is None:
    return None
  if imol_refinement_map()==-1:
    _status_message("Water added; no refinement map set")
    return None
  sphere_refine()
  accept_moving_atoms()
  return mol_id

def resample_active_map_for_em_half_angstrom():
  map_id=_scrollable_map_or_status()
  if map_id is None:
    return None
  if map_is_difference_map(map_id)!=0:
    _status_message("Difference maps are left unchanged")
    return None
  contour_sigma=get_contour_level_in_sigma(map_id)
  map_file_name=molecule_name(map_id)
  grid_info=_read_ccp4_mrc_grid_info(map_file_name)
  if grid_info:
    grid_spacing=min(grid_info["spacing_xyz"])
  else:
    grid_spacing=None
  if grid_spacing is None:
    _restyle_em_map_mesh(map_id, contour_sigma)
    _status_message("Grid spacing could not be determined; restyled active map")
    return None
  if grid_spacing <= 0.5:
    _restyle_em_map_mesh(map_id, contour_sigma)
    _status_message("Map already finer than 0.5 A/pixel; restyled active map")
    return map_id
  resample_factor=grid_spacing/0.5
  if grid_info:
    current_grid=grid_info["grid_xyz"]
    target_grid=[
      int(math.ceil(current_grid[0]*resample_factor)),
      int(math.ceil(current_grid[1]*resample_factor)),
      int(math.ceil(current_grid[2]*resample_factor))
    ]
    target_voxels=target_grid[0]*target_grid[1]*target_grid[2]
    if target_voxels > MAX_EM_RESAMPLE_TARGET_VOXELS:
      _restyle_em_map_mesh(map_id, contour_sigma)
      info_dialog("Refusing to resample this map to 0.5 A/pixel\nbecause the target grid would be too large.\n\nTarget grid: %d x %d x %d\nTarget size: %.1f million voxels\n\nThe active map has been restyled instead.\nIncrease MAX_EM_RESAMPLE_TARGET_VOXELS in the script\nif you really want to allow this." % (target_grid[0], target_grid[1], target_grid[2], target_voxels/1000000.0))
      _status_message("Target resampled map too large; restyled active map")
      return None
  old_refinement_map=(imol_refinement_map()==map_id)
  try:
    new_map_id=sharpen_blur_map_with_resampling(map_id, 0.0, resample_factor)
  except:
    new_map_id=-1
  if new_map_id not in molecule_number_list():
    _restyle_em_map_mesh(map_id, contour_sigma)
    _status_message("Resampling failed; restyled active map")
    return None
  _restyle_em_map_mesh(new_map_id, contour_sigma)
  set_scroll_wheel_map(new_map_id)
  set_scrollable_map(new_map_id)
  if old_refinement_map:
    set_imol_refinement_map(new_map_id)
  close_molecule(map_id)
  _status_message("Created EM-style 0.5 A/pixel map")
  return new_map_id

def _capture_view_state():
  return {
    "rotation_centre": _rotation_centre_xyz(),
    "zoom": zoom_factor(),
    "quaternion": [get_view_quaternion_internal(0), get_view_quaternion_internal(1),
                   get_view_quaternion_internal(2), get_view_quaternion_internal(3)]
  }

def _restore_view_state(view_state):
  set_rotation_centre(*view_state["rotation_centre"])
  set_view_quaternion(*view_state["quaternion"])
  set_zoom(view_state["zoom"])

def _close_smart_copy_template():
  global SMART_COPY_TEMPLATE_IMOL
  if SMART_COPY_TEMPLATE_IMOL in molecule_number_list():
    close_molecule(SMART_COPY_TEMPLATE_IMOL)
  SMART_COPY_TEMPLATE_IMOL=None

def _set_smart_copy_template(imol):
  global SMART_COPY_TEMPLATE_IMOL
  SMART_COPY_TEMPLATE_IMOL=imol
  try:
    set_mol_displayed(imol,0)
    set_mol_active(imol,0)
  except:
    pass

def _ensure_non_polymer_restraints_loaded_from_molecule(source_imol):
  if source_imol not in model_molecule_list():
    return False
  transferred_any=False
  seen_comp_ids={}
  for residue_spec in all_residues_sans_water(source_imol) or []:
    chain_id=residue_spec_to_chain_id(residue_spec)
    resno=residue_spec_to_res_no(residue_spec)
    ins_code=residue_spec_to_ins_code(residue_spec)
    if chain_id is False or resno is False or ins_code is False:
      continue
    if _residue_is_polymer(source_imol, chain_id, resno, ins_code):
      continue
    comp_id=residue_name(source_imol, chain_id, resno, ins_code)
    if not comp_id or comp_id in seen_comp_ids:
      continue
    seen_comp_ids[comp_id]=True
    try:
      restraints=monomer_restraints_for_molecule_py(comp_id, source_imol)
    except:
      restraints=False
    if isinstance(restraints, dict):
      try:
        set_monomer_restraints_py(comp_id, restraints)
        transferred_any=True
        continue
      except:
        pass
    try:
      add_dictionary_from_residue(source_imol, chain_id, resno, ins_code)
      transferred_any=True
    except:
      pass
  return transferred_any

def _smart_copy_atom_selection(chain_id, resno, ins_code):
  if ins_code:
    return None
  return "//%s/%s" %(chain_id,resno)

def _sorted_residue_specs(residue_specs):
  return sorted([list(spec) for spec in residue_specs],
                key=lambda spec: (spec[0], spec[1], spec[2]))

def _fill_short_polymer_gaps(molecule_id, polymer_specs, max_gap_to_fill):
  expanded_specs=set(polymer_specs)
  polymer_by_chain={}
  for chain_id, residue_no, ins_code in polymer_specs:
    if ins_code=="":
      polymer_by_chain.setdefault(chain_id,set()).add(residue_no)
  for chain_id, residue_numbers in polymer_by_chain.items():
    for _mol_id, segment_chain_id, segment_start, segment_end in segment_list_chain(molecule_id, chain_id):
      residues_in_segment=sorted([residue_no for residue_no in residue_numbers
                                  if segment_start <= residue_no <= segment_end])
      for index in range(len(residues_in_segment)-1):
        left_residue=residues_in_segment[index]
        right_residue=residues_in_segment[index+1]
        gap_size=right_residue-left_residue-1
        if gap_size <= 0 or gap_size > max_gap_to_fill:
          continue
        for fill_residue in range(left_residue+1, right_residue):
          if residue_exists_qm(molecule_id, segment_chain_id, fill_residue, ""):
            expanded_specs.add((segment_chain_id, fill_residue, ""))
  return expanded_specs

def _expand_polymer_contact_windows(molecule_id, polymer_specs, window_size):
  expanded_specs=set()
  polymer_by_chain={}
  for chain_id, residue_no, ins_code in polymer_specs:
    if ins_code=="":
      polymer_by_chain.setdefault(chain_id,set()).add(residue_no)
    else:
      expanded_specs.add((chain_id, residue_no, ins_code))
  for chain_id, residue_numbers in polymer_by_chain.items():
    fpr=first_polymer_residue(molecule_id, chain_id)
    lpr=last_polymer_residue(molecule_id, chain_id)
    if fpr==-10000 or lpr==-10000:
      continue
    for residue_no in residue_numbers:
      res_start=max(fpr, residue_no-window_size)
      res_end=min(lpr, residue_no+window_size)
      for resn in range(res_start, res_end+1):
        if residue_exists_qm(molecule_id, chain_id, resn, ""):
          expanded_specs.add((chain_id, resn, ""))
  return expanded_specs

def _prune_small_polymer_fragments(molecule_id, polymer_specs, minimum_fragment_size):
  kept_specs=set()
  polymer_by_chain={}
  for chain_id, residue_no, ins_code in polymer_specs:
    if ins_code=="":
      polymer_by_chain.setdefault(chain_id,set()).add(residue_no)
  for chain_id, residue_numbers in polymer_by_chain.items():
    for _mol_id, segment_chain_id, segment_start, segment_end in segment_list_chain(molecule_id, chain_id):
      residues_in_segment=sorted([residue_no for residue_no in residue_numbers
                                  if segment_start <= residue_no <= segment_end])
      if not residues_in_segment:
        continue
      cluster=[residues_in_segment[0]]
      for residue_no in residues_in_segment[1:]:
        if residue_no == cluster[-1] + 1:
          cluster.append(residue_no)
        else:
          if len(cluster) >= minimum_fragment_size:
            for cluster_residue in cluster:
              kept_specs.add((segment_chain_id, cluster_residue, ""))
          cluster=[residue_no]
      if len(cluster) >= minimum_fragment_size:
        for cluster_residue in cluster:
          kept_specs.add((segment_chain_id, cluster_residue, ""))
  return kept_specs

def smart_copy_active_non_polymer_residue():
  global SMART_COPY_SOURCE_CENTRE
  global SMART_COPY_RESIDUE_NAME
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  resno=residue[2]
  ins_code=residue[3]
  if _residue_is_polymer(mol_id, ch_id, resno, ins_code):
    _status_message("Active residue is polymer; nothing copied")
    return None
  atom_selection=_smart_copy_atom_selection(ch_id, resno, ins_code)
  if atom_selection is None:
    _status_message("Smart copy does not yet support insertion-code residues")
    return None
  view_state=_capture_view_state()
  recentre_state=recentre_on_read_pdb()
  copied_imol=-1
  try:
    set_recentre_on_read_pdb(0)
    copied_imol=new_molecule_by_atom_selection(mol_id, atom_selection)
  finally:
    set_recentre_on_read_pdb(recentre_state)
    _restore_view_state(view_state)
  if copied_imol==-1:
    _status_message("Unable to copy active residue")
    return None
  _close_smart_copy_template()
  _set_smart_copy_template(copied_imol)
  SMART_COPY_SOURCE_CENTRE=residue_centre_py(mol_id, ch_id, resno, ins_code)
  SMART_COPY_RESIDUE_NAME=residue_name(mol_id, ch_id, resno, ins_code)
  _status_message("Copied %s for smart paste" % SMART_COPY_RESIDUE_NAME)

def smart_paste_copied_non_polymer_residue():
  if SMART_COPY_TEMPLATE_IMOL not in molecule_number_list():
    _status_message("No copied ligand, ion, or water available")
    return None
  residue=active_residue()
  if residue:
    target_mol_id=residue[0]
  else:
    target_mol_id=go_to_atom_molecule_number()
    if target_mol_id not in model_molecule_list():
      _status_message("No active model for smart paste")
      return None
  if SMART_COPY_SOURCE_CENTRE is None:
    _status_message("Copied residue centre is unavailable")
    return None
  paste_imol=copy_molecule(SMART_COPY_TEMPLATE_IMOL)
  if paste_imol==-1:
    _status_message("Unable to prepare pasted residue")
    return None
  pointer_position=_rotation_centre_xyz()
  dx=pointer_position[0]-SMART_COPY_SOURCE_CENTRE[0]
  dy=pointer_position[1]-SMART_COPY_SOURCE_CENTRE[1]
  dz=pointer_position[2]-SMART_COPY_SOURCE_CENTRE[2]
  try:
    _ensure_non_polymer_restraints_loaded_from_molecule(SMART_COPY_TEMPLATE_IMOL)
    translate_molecule_by(paste_imol, dx, dy, dz)
    merge_molecules([paste_imol], target_mol_id)
  finally:
    if paste_imol in molecule_number_list():
      close_molecule(paste_imol)
  _status_message("Pasted %s at pointer" % (SMART_COPY_RESIDUE_NAME or "copied residue"))

def go_to_nearest_density_peak():
  map_id=_scrollable_map_or_status()
  if map_id is None:
    return None
  centre=_rotation_centre_xyz()
  def density_here(point):
    try:
      return density_at_point(map_id, point[0], point[1], point[2])
    except:
      return None

  peak_point=[centre[0], centre[1], centre[2]]
  peak_density=density_here(peak_point)
  if peak_density is None:
    _status_message("Unable to evaluate map density at the current position")
    return None

  for step_size in [0.5, 0.25, 0.1, 0.05]:
    for step_index in range(40):
      best_point=peak_point
      best_density=peak_density
      for dx in [-step_size, 0.0, step_size]:
        for dy in [-step_size, 0.0, step_size]:
          for dz in [-step_size, 0.0, step_size]:
            if dx==0.0 and dy==0.0 and dz==0.0:
              continue
            trial_point=[peak_point[0]+dx, peak_point[1]+dy, peak_point[2]+dz]
            trial_density=density_here(trial_point)
            if trial_density is None:
              continue
            if trial_density > best_density:
              best_point=trial_point
              best_density=trial_density
      if best_point == peak_point:
        break
      peak_point=best_point
      peak_density=best_density

  if _distance_sq(peak_point, centre) > 4.0*4.0:
    _status_message("No nearby density peak within 4 A")
    return None

  set_rotation_centre(peak_point[0], peak_point[1], peak_point[2])
  return peak_point

def _parsed_atom_record(atom):
  position=residue_atom_to_position(atom)
  if not position:
    return None
  return {
    "name": residue_atom_to_atom_name(atom),
    "alt_conf": residue_atom2alt_conf(atom),
    "position": position
  }

def _residue_record_for_restraints(mol_id, residue_spec):
  chain_id, resno, ins_code=residue_spec
  residue_atoms=residue_info_py(mol_id, chain_id, resno, ins_code) or []
  parsed_atoms=[]
  backbone_atoms=[]
  for atom in residue_atoms:
    parsed_atom=_parsed_atom_record(atom)
    if not parsed_atom:
      continue
    parsed_atoms.append(parsed_atom)
    atom_name=parsed_atom["name"].strip()
    if atom_name in ["N", "O", "OXT"]:
      backbone_atoms.append(parsed_atom)
  if not parsed_atoms:
    return None
  xs=[atom["position"][0] for atom in parsed_atoms]
  ys=[atom["position"][1] for atom in parsed_atoms]
  zs=[atom["position"][2] for atom in parsed_atoms]
  return {
    "spec": residue_spec,
    "atoms": parsed_atoms,
    "backbone_atoms": backbone_atoms,
    "bbox": (min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)),
    "polymer": _residue_is_polymer(mol_id, chain_id, resno, ins_code)
  }

def _generate_smart_local_extra_restraints_for_mol(mol_id, distance_cutoff=3.7):
  if mol_id not in model_molecule_list():
    _status_message("No active model")
    return None
  distance_cutoff_sq=distance_cutoff*distance_cutoff
  max_sequence_separation=10
  restraint_esd=0.05
  sequence_index_by_residue={}
  ordered_residue_specs_by_chain={}
  next_sequence_index_by_chain={}
  residue_serial_entries=all_residues_with_serial_numbers(mol_id) or []
  for residue_entry in residue_serial_entries:
    if not residue_entry or len(residue_entry) < 4:
      continue
    residue_spec=residue_entry[1:]
    chain_id=residue_spec_to_chain_id(residue_spec)
    resno=residue_spec_to_res_no(residue_spec)
    ins_code=residue_spec_to_ins_code(residue_spec)
    if chain_id is False or resno is False or ins_code is False:
      continue
    sequence_index_by_residue[(chain_id,resno,ins_code)]=next_sequence_index_by_chain.get(chain_id,0)
    next_sequence_index_by_chain[chain_id]=next_sequence_index_by_chain.get(chain_id,0)+1
    ordered_residue_specs_by_chain.setdefault(chain_id,[]).append((chain_id,resno,ins_code))
  residue_records={}
  for chain_id, residue_specs in ordered_residue_specs_by_chain.items():
    for residue_spec in residue_specs:
      residue_records[residue_spec]=_residue_record_for_restraints(mol_id, residue_spec)

  def atom_distance_sq(atom_1, atom_2):
    dx=atom_1["position"][0]-atom_2["position"][0]
    dy=atom_1["position"][1]-atom_2["position"][1]
    dz=atom_1["position"][2]-atom_2["position"][2]
    return dx*dx + dy*dy + dz*dz

  def residue_records_can_contact(record_1, record_2, cutoff):
    bbox_1=record_1["bbox"]
    bbox_2=record_2["bbox"]
    if bbox_1[1] + cutoff < bbox_2[0] or bbox_2[1] + cutoff < bbox_1[0]:
      return False
    if bbox_1[3] + cutoff < bbox_2[2] or bbox_2[3] + cutoff < bbox_1[2]:
      return False
    if bbox_1[5] + cutoff < bbox_2[4] or bbox_2[5] + cutoff < bbox_1[4]:
      return False
    return True

  def sequence_separation(residue_spec_1, residue_spec_2):
    if residue_spec_1[0] != residue_spec_2[0]:
      return None
    index_1=sequence_index_by_residue.get(residue_spec_1)
    index_2=sequence_index_by_residue.get(residue_spec_2)
    if index_1 is None or index_2 is None:
      return None
    return abs(index_1-index_2)

  def ordered_pair_key(spec_1, spec_2):
    if spec_2 < spec_1:
      return (spec_2, spec_1)
    return (spec_1, spec_2)

  def is_backbone_hbond_candidate(atom_name_1, atom_name_2):
    atom_name_1=atom_name_1.strip()
    atom_name_2=atom_name_2.strip()
    return ((atom_name_1=="N" and atom_name_2 in ["O","OXT"])
            or (atom_name_2=="N" and atom_name_1 in ["O","OXT"]))

  added_restraints=[0]
  seen_restraints=set()
  processed_backbone_pairs=set()
  delete_all_extra_restraints(mol_id)

  def add_restraints_between_records(record_1, record_2, backbone_only):
    if not record_1 or not record_2:
      return
    if backbone_only and not (record_1["polymer"] and record_2["polymer"]):
      return
    if not residue_records_can_contact(record_1, record_2, distance_cutoff):
      return
    residue_pair_key=ordered_pair_key(record_1["spec"], record_2["spec"])
    if backbone_only:
      atoms_1=record_1["backbone_atoms"]
      atoms_2=record_2["backbone_atoms"]
    else:
      atoms_1=record_1["atoms"]
      atoms_2=record_2["atoms"]
    for atom_1 in atoms_1:
      atom_name_1=atom_1["name"]
      alt_conf_1=atom_1["alt_conf"]
      for atom_2 in atoms_2:
        atom_name_2=atom_2["name"]
        alt_conf_2=atom_2["alt_conf"]
        if backbone_only and not is_backbone_hbond_candidate(atom_name_1, atom_name_2):
          continue
        restraint_key=(residue_pair_key, (atom_name_1, alt_conf_1), (atom_name_2, alt_conf_2))
        reverse_restraint_key=(residue_pair_key, (atom_name_2, alt_conf_2), (atom_name_1, alt_conf_1))
        if restraint_key in seen_restraints or reverse_restraint_key in seen_restraints:
          continue
        distance_sq=atom_distance_sq(atom_1, atom_2)
        if distance_sq >= distance_cutoff_sq:
          continue
        distance=math.sqrt(distance_sq)
        add_extra_geman_mcclure_restraint(
          mol_id,
          record_1["spec"][0], record_1["spec"][1], record_1["spec"][2], atom_name_1, alt_conf_1,
          record_2["spec"][0], record_2["spec"][1], record_2["spec"][2], atom_name_2, alt_conf_2,
          distance, restraint_esd)
        seen_restraints.add(restraint_key)
        added_restraints[0]=added_restraints[0]+1

  processed_local_pairs=set()
  for residue_spec_1 in residue_records.keys():
    record_1=residue_records.get(residue_spec_1)
    if not record_1:
      continue
    nearby_residue_specs=residues_near_residue(mol_id, list(residue_spec_1), distance_cutoff) or []
    for nearby_spec in nearby_residue_specs:
      if not nearby_spec or len(nearby_spec) < 3:
        continue
      residue_spec_2=(nearby_spec[0], nearby_spec[1], nearby_spec[2])
      if residue_spec_1 == residue_spec_2:
        continue
      pair_key=ordered_pair_key(residue_spec_1, residue_spec_2)
      if pair_key in processed_local_pairs:
        continue
      processed_local_pairs.add(pair_key)
      separation=sequence_separation(residue_spec_1, residue_spec_2)
      if separation is None or separation > max_sequence_separation:
        continue
      add_restraints_between_records(record_1, residue_records.get(residue_spec_2), False)

  for residue_spec_1 in residue_records.keys():
    record_1=residue_records.get(residue_spec_1)
    if not record_1 or not record_1["polymer"]:
      continue
    nearby_residue_specs=residues_near_residue(mol_id, list(residue_spec_1), distance_cutoff) or []
    for nearby_spec in nearby_residue_specs:
      if not nearby_spec or len(nearby_spec) < 3:
        continue
      residue_spec_2=(nearby_spec[0], nearby_spec[1], nearby_spec[2])
      if residue_spec_1 == residue_spec_2:
        continue
      pair_key=ordered_pair_key(residue_spec_1, residue_spec_2)
      if pair_key in processed_backbone_pairs:
        continue
      processed_backbone_pairs.add(pair_key)
      separation=sequence_separation(residue_spec_1, residue_spec_2)
      if separation is not None and separation <= max_sequence_separation:
        continue
      add_restraints_between_records(record_1, residue_records.get(residue_spec_2), True)

  set_show_extra_restraints(mol_id,0)
  set_show_extra_restraints(mol_id,1)
  _status_message("Smart local extra restraints generated: %s added" % added_restraints[0])
  print "Smart local extra restraints generated: %s added" % added_restraints[0]
  return added_restraints[0]

def generate_smart_local_extra_restraints():
  residue=_active_residue_or_status()
  if not residue:
    return None
  return _generate_smart_local_extra_restraints_for_mol(residue[0])

def generate_smart_local_extra_restraints_with_cutoff(distance_cutoff):
  residue=_active_residue_or_status()
  if not residue:
    return None
  try:
    parsed_cutoff=float(distance_cutoff)
  except:
    info_dialog("Minimum interatomic distance must be a number")
    return None
  if parsed_cutoff <= 0.0:
    info_dialog("Minimum interatomic distance must be greater than 0")
    return None
  return _generate_smart_local_extra_restraints_for_mol(residue[0], parsed_cutoff)

def prompt_generate_smart_local_extra_restraints():
  generic_single_entry("Minimum interatomic distance for smart restraints (A)",
                       "3.7",
                       "Generate smart self restraints",
                       generate_smart_local_extra_restraints_with_cutoff)

def confirm_generate_smart_local_extra_restraints():
  proceed=yes_no_dialog(
    "Generate smart local extra restraints for the active model?\n\nThis may be slow for large structures.",
    "Generate restraints")
  if not proceed:
    return None
  return generate_smart_local_extra_restraints()

def jiggle_fit_active_non_polymer_residue():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  resno=residue[2]
  ins_code=residue[3]
  if _residue_is_polymer(mol_id, ch_id, resno, ins_code):
    _status_message("Jiggle fit only applies to non-polymer residues")
    return None
  _status_message("Jiggle fitting active non-polymer residue")
  return fit_to_map_by_random_jiggle(mol_id, ch_id, resno, ins_code, 100, 1.0)

def decrease_proportional_editing_radius_with_status():
  global PROPORTIONAL_EDITING_RADIUS
  if PROPORTIONAL_EDITING_RADIUS <= 1.0:
    PROPORTIONAL_EDITING_RADIUS=1.0
    _status_message("Proportional editing radius: %.1f A" % PROPORTIONAL_EDITING_RADIUS)
    return None
  decrease_proportional_editing_radius()
  PROPORTIONAL_EDITING_RADIUS=max(1.0, PROPORTIONAL_EDITING_RADIUS-1.0)
  _status_message("Proportional editing radius: %.1f A" % PROPORTIONAL_EDITING_RADIUS)

def increase_proportional_editing_radius_with_status():
  global PROPORTIONAL_EDITING_RADIUS
  if PROPORTIONAL_EDITING_RADIUS >= 1000.0:
    PROPORTIONAL_EDITING_RADIUS=1000.0
    _status_message("Proportional editing radius: %.1f A" % PROPORTIONAL_EDITING_RADIUS)
    return None
  increase_proportional_editing_radius()
  PROPORTIONAL_EDITING_RADIUS=min(1000.0, PROPORTIONAL_EDITING_RADIUS+1.0)
  _status_message("Proportional editing radius: %.1f A" % PROPORTIONAL_EDITING_RADIUS)

def set_proportional_editing_radius():
  current_radius_string="%.1f" % PROPORTIONAL_EDITING_RADIUS
  def set_proportional_editing_radius_from_text(text):
    global PROPORTIONAL_EDITING_RADIUS
    try:
      requested_radius=float(text)
    except ValueError:
      info_dialog("Radius must be a number")
      return None
    if requested_radius < 1.0 or requested_radius > 1000.0:
      info_dialog("Proportional editing radius must be between 1 and 1000 A")
      return None
    target_radius=float(int(round(requested_radius)))
    while PROPORTIONAL_EDITING_RADIUS < target_radius:
      increase_proportional_editing_radius()
      PROPORTIONAL_EDITING_RADIUS=PROPORTIONAL_EDITING_RADIUS+1.0
    while PROPORTIONAL_EDITING_RADIUS > target_radius:
      decrease_proportional_editing_radius()
      PROPORTIONAL_EDITING_RADIUS=max(1.0, PROPORTIONAL_EDITING_RADIUS-1.0)
    if abs(target_radius-requested_radius) > 0.001:
      _status_message("Proportional editing radius set to %.1f A (rounded to 1 A steps)" % PROPORTIONAL_EDITING_RADIUS)
    else:
      _status_message("Proportional editing radius set to %.1f A" % PROPORTIONAL_EDITING_RADIUS)
  generic_single_entry("New proportional editing radius (1 A steps)", current_radius_string,
                       "Set proportional editing radius", set_proportional_editing_radius_from_text)

def display_only_active_map():
  active_map=scroll_wheel_map()
  displayed_maps=[map_id for map_id in map_molecule_list() if map_is_displayed(map_id)]
  if len(displayed_maps)==1 and active_map!=displayed_maps[0]:
    set_scroll_wheel_map(displayed_maps[0])
    set_scrollable_map(displayed_maps[0])
    return None
  if not scroll_wheel_map() in map_molecule_list():
    for map_id in map_molecule_list():
      if (map_is_displayed(map_id)==1) and (map_id!=active_map):
        set_scroll_wheel_map(map_id)
        set_scrollable_map(map_id)
      else:
        set_scroll_wheel_map(map_molecule_list()[0])
        set_scrollable_map(map_molecule_list()[0])
    active_map=scroll_wheel_map()
  if not map_is_displayed(active_map):
    set_map_displayed(active_map,1)
  displayed_maps_count=0
  for map_id in map_molecule_list():
    displayed_maps_count=displayed_maps_count+map_is_displayed(map_id)
    if (map_is_displayed(map_id)==1) and (map_id!=active_map):
      set_map_displayed(map_id,0)
    if map_is_displayed(map_id):
      displayed_map=map_id
  if displayed_maps_count==1:
    index_displayed=map_molecule_list().index(active_map)
    try:
      next_map=map_molecule_list()[index_displayed+1]
    except IndexError:
      next_map=map_molecule_list()[0]
    set_map_displayed(active_map,0)
    set_map_displayed(next_map,1)
    set_scroll_wheel_map(next_map)
    set_scrollable_map(next_map)
  for map_id in map_molecule_list():
    if map_is_displayed(map_id):
      set_scrollable_map(map_id)
      set_scroll_wheel_map(map_id) #New

def _ensure_single_displayed_map_is_scrollable():
  displayed_maps=[map_id for map_id in map_molecule_list() if map_is_displayed(map_id)]
  if len(displayed_maps) == 1:
    set_scroll_wheel_map(displayed_maps[0])
    set_scrollable_map(displayed_maps[0])
    return displayed_maps[0]
  return None

def hide_active_mol():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  set_mol_displayed(mol_id,0)

def display_only_active():
  models=model_molecule_list()
  if not models:
    return None

  residue=_active_residue_or_status()
  if residue:
    mol_id_active=residue[0]
  else:
    try:
      mol_id_active=go_to_atom_molecule_number()
    except:
      mol_id_active=-1
    if mol_id_active not in models:
      displayed_models=[mol_id for mol_id in models if mol_is_displayed(mol_id)]
      if displayed_models:
        mol_id_active=displayed_models[0]
      else:
        mol_id_active=models[0]

  displayed_models=[mol_id for mol_id in models if mol_is_displayed(mol_id)]
  displayed_mols_count=len(displayed_models)

  if displayed_mols_count==1:
    displayed_mol=displayed_models[0]
    index_displayed=models.index(displayed_mol)
    try:
      next_mol=models[index_displayed+1]
    except IndexError:
      next_mol=models[0]
    set_mol_displayed(displayed_mol,0)
    set_mol_displayed(next_mol,1)
    return None

  for mol_id in models:
    if mol_id == mol_id_active:
      if mol_is_displayed(mol_id)==0:
        set_mol_displayed(mol_id,1)
    else:
      if mol_is_displayed(mol_id)==1:
        set_mol_displayed(mol_id,0)

#def user_defined_add_3_10_helix_restraints():
#  def make_restr(*args):
#    spec_1 = args[0]
#    spec_2 = args[1]
#    chain_id_1 = spec_1[2]
#    chain_id_2 = spec_2[2]
#    res_no_1   = spec_1[3]
#    res_no_2   = spec_2[3]
#    imol       = spec_1[1]
#    if (chain_id_1 == chain_id_2):
#      # if backwards, swap 'em
#      if res_no_2 < res_no_1:
#        tmp = res_no_1
#        res_no_1 = res_no_2
#        res_no_2 = tmp
#      for rn in range(res_no_1, res_no_2 - 2):
#        if (rn + 3 <= res_no_2):
#          add_extra_bond_restraint(imol, chain_id_1, rn    , "", " O  ", "", chain_id_1, rn + 3, "", " N  ", "", 2.91, 0.035)
#  user_defined_click(2, make_restr)
    
def step_map_coarse_up(mol_id):
  current_level=get_contour_level_in_sigma(mol_id)
  if (current_level >= 0.5) and (current_level <= 10.0):
    new_level=current_level+0.5
  elif (current_level<0.5):
    new_level=0.5
  elif (current_level>10.0):
    new_level=10.0
  set_contour_level_in_sigma(mol_id,new_level)

def step_map_coarse_down(mol_id):
  current_level=get_contour_level_in_sigma(mol_id)
  if (current_level >= 0.5) and (current_level <= 10.0):
    new_level=current_level-0.5
  elif (current_level<0.5):
    new_level=0.5
  elif (current_level>10.0):
    new_level=10.0
  set_contour_level_in_sigma(mol_id,new_level)

def toggle_global_map_view():
  map_id=_scrollable_map_or_status()
  if map_id is None:
    return None
  if map_is_difference_map(map_id)!=0:
    _status_message("Global map view does not apply to difference maps")
    return None
  state=MAP_GLOBAL_VIEW_SETTINGS.get(map_id)
  if not state or not state.get("enabled"):
    try:
      current_colour=map_colour_components(map_id)
    except:
      current_colour=get_map_colour(map_id)
    MAP_GLOBAL_VIEW_SETTINGS[map_id]={
      "enabled": 1,
      "colour": current_colour,
      "radius": get_map_radius()
    }
    max_radius=map_cell(map_id)[0]/2.0
    set_draw_solid_density_surface(map_id,1)
    set_draw_map_standard_lines(map_id,0)
    try:
      set_solid_density_surface_opacity(map_id,0.8)
    except:
      pass
    set_flat_shading_for_solid_density_surface(0)
    set_map_colour(map_id,EM_GLOBAL_MAP_COLOUR[0],EM_GLOBAL_MAP_COLOUR[1],EM_GLOBAL_MAP_COLOUR[2])
    set_map_radius(max_radius)
    try:
      set_map_radius_em(max_radius)
    except:
      pass
  else:
    set_draw_solid_density_surface(map_id,0)
    set_draw_map_standard_lines(map_id,1)
    if state.has_key("colour"):
      saved_colour=state.get("colour")
    else:
      saved_colour=get_map_colour(map_id)
    set_map_colour(map_id,saved_colour[0],saved_colour[1],saved_colour[2])
    restored_radius=state.get("radius",20.0)
    set_map_radius(restored_radius)
    try:
      set_map_radius_em(restored_radius)
    except:
      pass
    MAP_GLOBAL_VIEW_SETTINGS[map_id]["enabled"]=0
    
  
#Go to next residue in current polymer chain.
def next_res():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  resn=residue[2]
  atom_name=residue[4]
  sn=get_sn_from_resno(mol_id,ch_id,resn)
  next_sn=sn+1
  next_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),next_sn)
  set_go_to_atom_molecule(mol_id)
  if (next_res!=-10000 and is_protein_chain_p(mol_id,ch_id)==1):
    set_go_to_atom_chain_residue_atom_name(ch_id,next_res,atom_name)
  elif (next_res!=-10000 and is_nucleotide_chain_p(mol_id,ch_id)==1):
    set_go_to_atom_chain_residue_atom_name(ch_id,next_res,atom_name)

#Go to previous residue in current polymer chain.
def prev_res():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  resn=residue[2]
  atom_name=residue[4]
  sn=get_sn_from_resno(mol_id,ch_id,resn)
  if (sn>=1):
    sn=sn-1
  prev_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),sn)
  set_go_to_atom_molecule(mol_id)
  if (prev_res!=-10000 and is_protein_chain_p(mol_id,ch_id)==1):
    set_go_to_atom_chain_residue_atom_name(ch_id,prev_res,atom_name)
  elif (prev_res!=-10000 and is_nucleotide_chain_p(mol_id,ch_id)==1):
    set_go_to_atom_chain_residue_atom_name(ch_id,prev_res,atom_name)

def mutate_by_entered_code():
  def mutate_single_letter(X):
    entry=str(X).upper()
    residue=_active_residue_or_status()
    if not residue:
      return None
    mol_id=residue[0]
    ch_id=residue[1]
    resno=residue[2]
    ins_code=residue[3]
    resname=residue_name(mol_id,ch_id,resno,ins_code)
    map_id=imol_refinement_map()
    aa_dic={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
    nt_list=['A','C','T','G','U']
    if (resname in aa_dic.values()) and (aa_dic.get(entry,0)!=0):
      mutate(mol_id,ch_id,resno,ins_code,aa_dic.get(entry,0))
    elif (resname in nt_list) and (entry in nt_list):
      mutate_base(mol_id,ch_id,resno,ins_code,entry)
    else:
      info_dialog("Invalid target residue! Must be protein or nucleic acid, and entered code must be single letter.")
  generic_single_entry("New residue? (single letter code)","A","Mutate by single-letter code",mutate_single_letter)
  
def set_map_level_quickly():
  if scroll_wheel_map()!=-1 and map_is_displayed(scroll_wheel_map())!=0:
    current_map_level=get_contour_level_in_sigma(scroll_wheel_map())
    current_map_level="{0:.2f}".format(current_map_level)
    def set_map_level_quickly(X):
      try:
        map_level=float(X)
        map_id=scroll_wheel_map()
        set_contour_level_in_sigma(map_id, map_level)
      except ValueError:
        info_dialog("Has to be a number!") 
    generic_single_entry("New map level in sigma/RMS?",current_map_level,"Set map level",set_map_level_quickly)
  else:
    info_dialog("You need a (scrollable, displayed) map!")
  
def sequence_context():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  resnum=residue[2]
  ins_code=residue[3]
  resname=residue_name(mol_id,ch_id,resnum,ins_code)
  def get_aa_code(resnum):
    if residue_name(mol_id,ch_id,resnum,ins_code):
      aa_code=three_letter_code2single_letter(residue_name(mol_id,ch_id,resnum,ins_code))
      if (len(residue_name(mol_id,ch_id,resnum,ins_code))==1):
        aa_code=residue_name(mol_id,ch_id,resnum,ins_code)
      if (residue_name(mol_id,ch_id,resnum,ins_code)!="ALA") and (aa_code=="A") and (len(residue_name(mol_id,ch_id,resnum,ins_code))!=1):
        aa_code="X"
    else:
      aa_code="-"
    return aa_code
  current_res=get_aa_code(resnum)
  minus_context=""
  for i in range(1,10):
    aa_code=str(get_aa_code(resnum-i))
    minus_context=aa_code + minus_context
  plus_context=""
  for i in range(1,10):
    aa_code=str(get_aa_code(resnum+i))
    plus_context=plus_context + aa_code
  final_string="Residue:  " + resname + "  " + str(resnum) +"  Sequence context: ..." + minus_context + "[" + current_res + "]" + plus_context + "..."
  info_dialog(final_string)
  
#Toggle display of active map
map_disp_flag={0:0}
map_disp_flag_cycle=0
def toggle_map_display():
  global map_disp_flag
  global map_disp_flag_cycle
  if map_disp_flag_cycle==0:
    for map_id in map_molecule_list():
      disp_value=map_is_displayed(map_id)
      map_disp_flag[map_id]=disp_value
      if disp_value==1:
        set_map_displayed(map_id,0) #If any maps are displayed, undisplay them.
    map_disp_flag_cycle=1
  elif map_disp_flag_cycle==1:
    disp_counter=0 
    for map_id in map_molecule_list():
      if map_id not in map_disp_flag: #if the map wasn't present in the previous cycle, assign a disp_value for it
        disp_value=map_is_displayed(map_id)
        map_disp_flag[map_id]=disp_value
      if map_disp_flag[map_id]==1:
        set_map_displayed(map_id,1) #Redisplay any maps that were displayed on the previous cycle.
      disp_counter=disp_counter+map_disp_flag[map_id] #test
    if disp_counter==0: #If no maps were displayed in the prior cycle, display all maps.
      for map_id in map_molecule_list(): 
        set_map_displayed(map_id,1) 
    map_disp_flag_cycle=0
  _ensure_single_displayed_map_is_scrollable()

#Toggle display of modelling toolbar (assumes it is shown by default)
toolbar_toggle_var=0
def toggle_toolbar_display():
  global toolbar_toggle_var
  if toolbar_toggle_var==0:
    hide_modelling_toolbar()
    toolbar_toggle_var=1
  elif toolbar_toggle_var==1:
    show_modelling_toolbar()
    toolbar_toggle_var=0
    
#Colour active segment
def colour_active_segment():
  mol_id=active_residue()[0]
  segments=segment_list(mol_id)
  res_here=active_residue()[2]
  ins_code=active_residue()[3]
  ch_id=active_residue()[1]
  colour_list=[]
  blank_list=[]
  segment_colour=34
  blank_colour=0
  for seg in segments:
    if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      for res in range(res_start,res_end+1):
        res_color_spec=[([ch_id,res,""],segment_colour)]
        colour_list=colour_list+res_color_spec
    else:
      res_start=seg[2]
      res_end=seg[3]
      ch_id_here=seg[1]
      for res in range(res_start,res_end+1):
        blank_color_spec=[([ch_id_here,res,""],blank_colour)]
        blank_list=blank_list+blank_color_spec
  clear_user_defined_atom_colours(mol_id)
  set_user_defined_atom_colour_by_residue_py(mol_id,colour_list)
  set_user_defined_atom_colour_by_residue_py(mol_id,blank_list)
  graphics_to_user_defined_atom_colours_representation(mol_id)


def _all_residue_specs_for_colouring(mol_id):
  residue_specs=[]
  for residue_entry in all_residues_with_serial_numbers(mol_id) or []:
    if not residue_entry or len(residue_entry) < 4:
      continue
    residue_spec=residue_entry[1:]
    chain_id=residue_spec_to_chain_id(residue_spec)
    resno=residue_spec_to_res_no(residue_spec)
    ins_code=residue_spec_to_ins_code(residue_spec)
    if chain_id is False or resno is False or ins_code is False:
      continue
    residue_specs.append([chain_id,resno,ins_code])
  return residue_specs

def _active_polymer_molecule_for_colouring(action_name):
  residue=_active_residue_or_status()
  if not residue:
    return None
  if not _residue_is_polymer(residue[0], residue[1], residue[2], residue[3]):
    info_dialog(action_name+" requires the active residue to be polymer.")
    return None
  return residue[0]

def color_by_rama_native(mol_id):
  try:
    rama_results=all_molecule_ramachandran_score(mol_id)
  except NameError:
    info_dialog("This Coot build does not expose Ramachandran scores.")
    return None
  if not isinstance(rama_results, list) or len(rama_results) < 6:
    info_dialog("Unable to obtain Ramachandran scores.")
    return None
  scored_residues=rama_results[5]
  blank_colour=0
  rama_allowed_colour=27
  rama_outlier_colour=31
  blank_res_list=[]
  rama_colour_list=[]
  for residue_spec in _all_residue_specs_for_colouring(mol_id):
    blank_res_list.append((residue_spec, blank_colour))
  for item in scored_residues:
    if not isinstance(item, list) or len(item) < 3:
      continue
    residue_spec=item[1]
    rama_score=item[2]
    if not isinstance(residue_spec, list):
      continue
    if rama_score < 0.002:
      rama_colour_list.append((residue_spec[1:], rama_outlier_colour))
    elif rama_score < 0.02:
      rama_colour_list.append((residue_spec[1:], rama_allowed_colour))
  try:
    clear_user_defined_atom_colours(mol_id)
    set_user_defined_atom_colour_by_residue_py(mol_id, blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id, rama_colour_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
    info_dialog("Ramachandran coloring:\n\nRed = outlier (<0.2%)\n\nOrange = allowed/disfavored (<2%)")
  except NameError:
    info_dialog("You need a newer Coot build for user-defined residue coloring.")

def color_by_rama_native_for_active_residue():
  mol_id=_active_polymer_molecule_for_colouring("Ramachandran coloring")
  if mol_id is None:
    return None
  return color_by_rama_native(mol_id)

def color_by_density_fit_native(mol_id):
  map_id=imol_refinement_map()
  if map_id==-1:
    info_dialog("You need a refinement map for density-fit coloring.")
    return None
  residue_specs=all_residues_sans_water(mol_id)
  if not residue_specs:
    info_dialog("No residues found for density-fit coloring.")
    return None
  try:
    correlation_results=map_to_model_correlation_per_residue(mol_id, residue_specs, 0, map_id)
  except NameError:
    info_dialog("This Coot build does not expose native per-residue density fit scores.")
    return None
  blank_colour=0
  blank_res_list=[]
  density_colour_list=[]
  for residue_spec in _all_residue_specs_for_colouring(mol_id):
    blank_res_list.append((residue_spec, blank_colour))
  for item in correlation_results:
    if not isinstance(item, list) or len(item) < 2:
      continue
    residue_spec=item[0]
    score=item[1]
    if not isinstance(residue_spec, list):
      continue
    if score < 0.0:
      score=0.0
    if score > 1.0:
      score=1.0
    colour_index=int((1.0-score)*31+2)
    density_colour_list.append((residue_spec[1:], colour_index))
  try:
    clear_user_defined_atom_colours(mol_id)
    set_user_defined_atom_colour_by_residue_py(mol_id, blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id, density_colour_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
    info_dialog("Active molecule colored by model/map correlation, in spectral coloring (blue=CC 1.0, red=CC 0.0)")
  except NameError:
    info_dialog("You need a newer Coot build for user-defined residue coloring.")

def color_by_density_fit_native_for_active_residue():
  mol_id=_active_polymer_molecule_for_colouring("Density-fit coloring")
  if mol_id is None:
    return None
  return color_by_density_fit_native(mol_id)

def color_by_ncs_difference(mol_id):
  if mol_id not in model_molecule_list():
    info_dialog("You need an active model for NCS-difference coloring.")
    return None
  ncs_data=None
  target_chain_id=None
  for chain_id in chain_ids(mol_id):
    try:
      diffs=ncs_chain_differences(mol_id, chain_id)
    except:
      diffs=False
    if diffs:
      ncs_data=diffs
      target_chain_id=chain_id
      break
  if not ncs_data:
    info_dialog("No NCS-difference data were found for the active molecule.")
    return None

  blank_colour=0
  blank_res_list=[]
  ncs_scores={}

  for residue_spec in _all_residue_specs_for_colouring(mol_id):
    blank_res_list.append((residue_spec, blank_colour))

  for i in range(0, len(ncs_data), 3):
    try:
      peer_chain_id=ncs_data[i]
      current_target_chain_id=ncs_data[i+1]
      residue_diffs=ncs_data[i+2]
    except:
      continue
    if not isinstance(residue_diffs, list):
      continue
    for residue_diff in residue_diffs:
      if not isinstance(residue_diff, list) or len(residue_diff) < 3:
        continue
      peer_residue=residue_diff[0]
      target_residue=residue_diff[1]
      mean_diff=residue_diff[2]
      try:
        mean_diff=float(mean_diff)
      except:
        continue
      if isinstance(peer_residue, list) and len(peer_residue) >= 2:
        peer_spec=(peer_chain_id, peer_residue[0], peer_residue[1])
        ncs_scores[peer_spec]=max(ncs_scores.get(peer_spec, 0.0), mean_diff)
      if isinstance(target_residue, list) and len(target_residue) >= 2:
        target_spec=(current_target_chain_id, target_residue[0], target_residue[1])
        ncs_scores[target_spec]=max(ncs_scores.get(target_spec, 0.0), mean_diff)

  if not ncs_scores:
    info_dialog("No NCS-difference values were available for coloring.")
    return None

  ncs_colour_list=[]
  for residue_spec, mean_diff in ncs_scores.items():
    normalized_score=mean_diff/2.0
    if normalized_score < 0.0:
      normalized_score=0.0
    if normalized_score > 1.0:
      normalized_score=1.0
    colour_index=int(normalized_score*31+2)
    ncs_colour_list.append(([residue_spec[0], residue_spec[1], residue_spec[2]], colour_index))

  try:
    clear_user_defined_atom_colours(mol_id)
    set_user_defined_atom_colour_by_residue_py(mol_id, blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id, ncs_colour_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
    info_dialog("Active molecule colored by NCS difference, in spectral coloring (blue=low difference, red=high difference)")
  except NameError:
    info_dialog("You need a newer Coot build for user-defined residue coloring.")

def color_by_ncs_difference_for_active_residue():
  mol_id=_active_polymer_molecule_for_colouring("NCS-difference coloring")
  if mol_id is None:
    return None
  return color_by_ncs_difference(mol_id)

def color_by_clash_score(mol_id):
  if mol_id not in model_molecule_list():
    info_dialog("You need an active model for clash coloring.")
    return None
  try:
    overlap_data=molecule_atom_overlaps(mol_id)
  except:
    overlap_data=False
  if not isinstance(overlap_data, list):
    info_dialog("No clash data were available for the active molecule.")
    return None

  blank_colour=0
  blank_res_list=[]
  max_overlap_by_residue={}

  for residue_spec in _all_residue_specs_for_colouring(mol_id):
    blank_res_list.append((residue_spec, blank_colour))

  def accumulate_overlap_from_atom_spec(atom_spec, overlap_value):
    if not isinstance(atom_spec, list):
      return
    try:
      residue_spec=atom_spec_to_residue_spec(atom_spec)
    except:
      return
    if not isinstance(residue_spec, list) or len(residue_spec) < 3:
      return
    residue_key=(residue_spec_to_chain_id(residue_spec),
                 residue_spec_to_res_no(residue_spec),
                 residue_spec_to_ins_code(residue_spec))
    previous_max=max_overlap_by_residue.get(residue_key, 0.0)
    if overlap_value > previous_max:
      max_overlap_by_residue[residue_key]=overlap_value

  for overlap_item in overlap_data:
    if not isinstance(overlap_item, dict):
      continue
    try:
      overlap_value=float(overlap_item.get('overlap-volume', 0.0))
    except:
      continue
    atom_spec_1=overlap_item.get('atom-1-spec')
    atom_spec_2=overlap_item.get('atom-2-spec')
    accumulate_overlap_from_atom_spec(atom_spec_1, overlap_value)
    accumulate_overlap_from_atom_spec(atom_spec_2, overlap_value)

  if not max_overlap_by_residue:
    info_dialog("No clash data were available for the active molecule.")
    return None

  clash_colour_list=[]
  for residue_spec, max_overlap in max_overlap_by_residue.items():
    normalized_score=max_overlap/2.0
    if normalized_score < 0.0:
      normalized_score=0.0
    if normalized_score > 1.0:
      normalized_score=1.0
    colour_index=int(normalized_score*31+2)
    clash_colour_list.append(([residue_spec[0], residue_spec[1], residue_spec[2]], colour_index))

  try:
    clear_user_defined_atom_colours(mol_id)
    set_user_defined_atom_colour_by_residue_py(mol_id, blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id, clash_colour_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
    info_dialog("Active molecule colored by per-residue maximum clash overlap, in spectral coloring (blue=low clash, red=high clash)")
  except NameError:
    info_dialog("You need a newer Coot build for user-defined residue coloring.")

def color_by_clash_score_for_active_molecule():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  return color_by_clash_score(mol_id)

def show_custom_keybindings_summary():
  info_dialog(
    "Custom keybindings\n\n"
    "G  Toggle global map view\n"
    "A  Refine zone (click 2 atoms)\n"
    "a  Local cylinder refine\n"
    "b  Go to nearest density peak\n"
    "C  Copy active ligand/ion/solvent\n"
    "V  Paste copied ligand/ion/solvent\n"
    "g  Generate smart local extra restraints\n"
    "h  Place helix here\n"
    "J  Jiggle-fit active non-polymer residue\n"
    "m  Measure distance\n"
    "o  Go to next NCS chain\n"
    "O  Go to NCS master chain\n"
    "q  Flip peptide\n"
    "r  Refine triple\n"
    "v  Undo symmetry view\n"
    "w/W  Add water / add water + refine\n"
    "y/Y/T  Add terminal residue / cycle terminus phi / psi\n"
    "z/x  Undo / redo\n"
    "Z  Clear labels and distances\n"
    "[ ]  Cycle representation mode\n"
    "{ }  Cycle symmetry representation mode\n"
    "; and '  Decrease / increase map radius\n"
    "! to (  Set map to 1-9 sigma\n"
    "- and =  Increase / decrease map contour by 0.1 sigma\n"
    "? / ~  Show only active model / map\n"
    "backtick / slash  Toggle map / model display\n"
    "< >  Previous / next residue in chain"
  )
    

#Set refinement map to currently scrollable map
def set_map_to_scrollable_map():
  if map_molecule_list()!=[]:
    set_imol_refinement_map(scroll_wheel_map())
  else:
    info_dialog("You need a map for this to work.")
  
#Toggle display of active model
mol_disp_flag={0:0}
mol_disp_flag_cycle=0
def toggle_mol_display():
  global mol_disp_flag
  global mol_disp_flag_cycle
  if mol_disp_flag_cycle==0:
    for mol_id in model_molecule_list():
      mol_disp_value=mol_is_displayed(mol_id)
      mol_disp_flag[mol_id]=mol_disp_value
      if mol_disp_value==1:
        set_mol_displayed(mol_id,0)
    mol_disp_flag_cycle=1
  elif mol_disp_flag_cycle==1:
    disp_counter=0
    for mol_id in model_molecule_list():
      if mol_id not in mol_disp_flag:
        disp_value=mol_is_displayed(mol_id)
        mol_disp_flag[mol_id]=disp_value
      if mol_disp_flag[mol_id]==1:
        set_mol_displayed(mol_id,1)
      disp_counter=disp_counter+mol_disp_flag[mol_id]
    if disp_counter==0:
      for mol_id in model_molecule_list():
        set_mol_displayed(mol_id,1)
    mol_disp_flag_cycle=0

#Cycle representation mode forward/back
cycle_rep_flag={0:0}
def _cycle_rep_flag_or_default(mol_id):
  flag=cycle_rep_flag.get(mol_id, 0)
  if flag not in [0,1,2,3,4,5]:
    cycle_rep_flag[mol_id]=0
    return 0
  return flag

def cycle_rep_up_active():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  return cycle_rep_up(mol_id, _cycle_rep_flag_or_default(mol_id))

def cycle_rep_up(mol_id,flag):
  global cycle_rep_flag
  cycle_rep_flag[mol_id]=flag
  if cycle_rep_flag[mol_id]==0:
    graphics_to_ca_plus_ligands_representation(mol_id)
    cycle_rep_flag[mol_id]=1
  elif cycle_rep_flag[mol_id]==1:
    graphics_to_ca_plus_ligands_and_sidechains_representation(mol_id)
    cycle_rep_flag[mol_id]=2
  elif cycle_rep_flag[mol_id]==2:
    graphics_to_rainbow_representation(mol_id)
    cycle_rep_flag[mol_id]=3
  elif cycle_rep_flag[mol_id]==3:
    graphics_to_bonds_representation(mol_id)
    cycle_rep_flag[mol_id]=4
  elif cycle_rep_flag[mol_id]==4:
    try:
      graphics_to_user_defined_atom_colours_representation(mol_id)
      cycle_rep_flag[mol_id]=5
    except NameError:
      graphics_to_bonds_representation(mol_id)
      cycle_rep_flag[mol_id]=0
  elif cycle_rep_flag[mol_id]==5:
    try:
      graphics_to_user_defined_atom_colours_all_atoms_representation(mol_id)
      cycle_rep_flag[mol_id]=0
    except NameError:
      graphics_to_bonds_representation(mol_id)
      cycle_rep_flag[mol_id]=0

# def add_partial_water():
#   place_typed_atom_at_pointer("Water")
#   refine_active_residue()
#   mol_id=active_residue()[0]
#   ch_id=active_residue()[1]
#   resno=active_residue()[2]
#   ins_code=active_residue()[3]
#   atom_name=active_residue()[4]
#   alt_conf=""
#   set_atom_attribute(mol_id,ch_id,resno,ins_code,atom_name,alt_conf,"occ",0.5)
# add_key_binding("Add partial water","w",
# lambda: add_partial_water())

# def delete_alt_confs():
#   mol_id=active_residue()[0]
#   for ch_id in chain_ids(mol_id):
#     for resn in (first_residue(mol_id,ch_id),last_residue(mol_id,ch_id)+1):
#      residue_inf=residue_info(mol_id,ch_id,resn,"")      

    
def cycle_rep_down(mol_id,flag):
  global cycle_rep_flag
  cycle_rep_flag[mol_id]=flag
  if cycle_rep_flag[mol_id]==3:
    graphics_to_ca_plus_ligands_and_sidechains_representation(mol_id)
    cycle_rep_flag[mol_id]=2
  elif cycle_rep_flag[mol_id]==2:
    graphics_to_ca_plus_ligands_representation(mol_id)
    cycle_rep_flag[mol_id]=1
  elif cycle_rep_flag[mol_id]==1:
    try:
      graphics_to_user_defined_atom_colours_all_atoms_representation(mol_id)
      cycle_rep_flag[mol_id]=0
    except NameError:
      graphics_to_bonds_representation(mol_id)
      cycle_rep_flag[mol_id]=5
  elif cycle_rep_flag[mol_id]==0:
    try:
      graphics_to_user_defined_atom_colours_representation(mol_id)
      cycle_rep_flag[mol_id]=5
    except NameError:
      graphics_to_bonds_representation(mol_id)
      cycle_rep_flag[mol_id]=5
  elif cycle_rep_flag[mol_id]==5:
    graphics_to_bonds_representation(mol_id)
    cycle_rep_flag[mol_id]=4
  elif cycle_rep_flag[mol_id]==4:
    graphics_to_rainbow_representation(mol_id)
    cycle_rep_flag[mol_id]=3

def cycle_rep_down_active():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  return cycle_rep_down(mol_id, _cycle_rep_flag_or_default(mol_id))


#Refine triple (Paul)
def key_binding_refine_triple():
    active_atom = active_residue()
    if not active_atom:
       print("No active atom")
    else:
       imol = active_atom[0]
       residue_spec = atom_spec_to_residue_spec(active_atom)
       N_terminal_residue_spec = [residue_spec_to_chain_id(residue_spec),
                                  residue_spec_to_res_no(residue_spec)-1,
                                  residue_spec_to_ins_code(residue_spec)]
       C_terminal_residue_spec = [residue_spec_to_chain_id(residue_spec),
                                  residue_spec_to_res_no(residue_spec)+1,
                                  residue_spec_to_ins_code(residue_spec)]
       spec_list = [N_terminal_residue_spec, residue_spec, C_terminal_residue_spec]
       refine_residues(imol, spec_list)
    
#Cycle symmetry represntation mode forward/back
cycle_symm_flag={0:0}
def cycle_symm_up_active():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  return cycle_symm_up(mol_id, cycle_symm_flag.get(mol_id,0))

def cycle_symm_up(mol_id,flag):
  global cycle_symm_flag
  cycle_symm_flag[mol_id]=flag
  if cycle_symm_flag[mol_id]==0:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,1)
    cycle_symm_flag[mol_id]=1
  elif cycle_symm_flag[mol_id]==1:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,0)
    set_symmetry_whole_chain(mol_id,1)
    cycle_symm_flag[mol_id]=2
  elif cycle_symm_flag[mol_id]==2:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,0)
    set_symmetry_whole_chain(mol_id,0)
    cycle_symm_flag[mol_id]=3
  elif cycle_symm_flag[mol_id]==3:
    get_show_symmetry()
    set_show_symmetry_master(0)
    cycle_symm_flag[mol_id]=0
def cycle_symm_down(mol_id,flag):
  global cycle_symm_flag
  cycle_symm_flag[mol_id]=flag
  if cycle_symm_flag[mol_id]==3:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,0)
    set_symmetry_whole_chain(mol_id,1)
    cycle_symm_flag[mol_id]=2
  elif cycle_symm_flag[mol_id]==2:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,1)
    cycle_symm_flag[mol_id]=1
  elif cycle_symm_flag[mol_id]==1:
    get_show_symmetry()
    set_show_symmetry_master(0)
    cycle_symm_flag[mol_id]=0
  elif cycle_symm_flag[mol_id]==0:
    get_show_symmetry()
    set_show_symmetry_master(1)
    symmetry_as_calphas(mol_id,0)
    set_symmetry_whole_chain(mol_id,0)
    cycle_symm_flag[mol_id]=3

def cycle_symm_down_active():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  return cycle_symm_down(mol_id, cycle_symm_flag.get(mol_id,0))
    
def undo_visible():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  set_undo_molecule(mol_id)
  apply_undo()
  
def redo_visible():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  set_undo_molecule(mol_id)
  apply_redo()

#Cycle rotamers for active residue
rotamer_number=0
def cycle_rotamers():
  global rotamer_number
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  res_here=active_residue()[2]
  ins_code=""
  alt_conf=""
  n_rots=n_rotamers(mol_id, ch_id, res_here, ins_code)-1
  turn_off_backup(mol_id)
  update_go_to_atom_from_current_position()
  if rotamer_number>=n_rots:
    rotamer_number=0
    set_residue_to_rotamer_number(mol_id,ch_id,res_here,ins_code,alt_conf,rotamer_number)
  else:
    rotamer_number=rotamer_number+1
    set_residue_to_rotamer_number(mol_id,ch_id,res_here,ins_code,alt_conf,rotamer_number)
  turn_on_backup(mol_id)
  
#Toggle environment distances
def toggle_env_dist():
  if show_environment_distances_state()==1:
    set_show_environment_distances(0)
  else:
    set_show_environment_distances(1)
    
#Check if residue is in a polymer by comparing resname to list
def is_polymer_residue(mol_id,ch_id,sn):
  valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
  resname=resname_from_serial_number(mol_id,ch_id,sn)
  if resname in valid_resnames:
    return 1
  else:
    return 0
    
#Return last residue in polymer
def last_polymer_residue(mol_id,ch_id):
  if ((valid_model_molecule_qm(mol_id)) and (ch_id in chain_ids(mol_id))):
    if not is_solvent_chain_qm(mol_id,ch_id):
      n=chain_n_residues(ch_id,mol_id)-1
      valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
      while resname_from_serial_number(mol_id,"%s"%(ch_id),n) not in valid_resnames:
        n=n-1 
      result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),n)
      return result
  else:
    return -1
    
def first_polymer_residue(mol_id,ch_id):
  if ((valid_model_molecule_qm(mol_id)) and (ch_id in chain_ids(mol_id))):
    if not is_solvent_chain_qm(mol_id,ch_id):
      n=0
      valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
      while resname_from_serial_number(mol_id,"%s"%(ch_id),n) not in valid_resnames:
        n=n+1 
      result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),n)
      return result
  else:
    return -1


#Return last res (polymer or not)
def last_residue(mol_id,ch_id):
          n=chain_n_residues(ch_id,mol_id)-1
          result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),n)
          return result

#Return first res (polymer or not)
def first_residue(mol_id,ch_id):
          result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
          return result
          
#Check if chain is polymer
def is_polymer(mol_id,ch_id):
        a=is_protein_chain_p(mol_id,"%s"%(ch_id))
        b=is_nucleotide_chain_p(mol_id,"%s"%(ch_id))
        if (a==1) or (b==1):
                result=1
        else:
                result=0
        return result
    
#Check if residue is last in polymer
def is_last_polymer_residue_sn(mol_id,ch_id,sn):
  if ((valid_model_molecule_qm(mol_id)) and (ch_id in chain_ids(mol_id))):
    if not is_solvent_chain_qm(mol_id,ch_id):
      valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
      if (resname_from_serial_number(mol_id,"%s"%(ch_id),sn) in valid_resnames) and (resname_from_serial_number(mol_id,"%s"%(ch_id),sn+1) not in valid_resnames):
        return 1
  else:
    return -1

#get serial number from resnum
def get_sn_from_resno(mol_id,ch_id,resno):
  sn=chain_n_residues(ch_id,mol_id)-1
  sn2=0
  resno_out=""
  resno_out_2=""
  if resno>last_residue(mol_id,ch_id) or resno<first_residue(mol_id,ch_id):
    return -1
  elif does_residue_exist_p(mol_id,ch_id,resno,"")==0:
    return -1
  elif resno==first_residue(mol_id,ch_id):
    return 0
  elif resno==chain_n_residues(ch_id,mol_id):
    return chain_n_residues(ch_id,mol_id)
  else:
    while (resno_out!=resno) and (resno_out_2!=resno):
      resno_out=seqnum_from_serial_number(mol_id,ch_id,sn)
      resno_out_2=seqnum_from_serial_number(mol_id,ch_id,sn2)
      sn=sn-1
      sn2=sn+1
    if resno_out_2==resno:
      return int(sn2+1)
    else:
      return int(sn+1)

def get_sn_from_resno_alt(mol_id,ch_id,resno):
  sn=0
  sn_max=chain_n_residues(ch_id,mol_id)-1
  sn_dict={}
  while sn<=sn_max:
    resno_here=seqnum_from_serial_number(mol_id,ch_id,sn)
    sn_dict[resno_here]=sn
    sn=sn+1
  try:
    sn_out=sn_dict[resno]
    return sn_out
  except KeyError:
    return -1

#check if res is at C-term side of break in mid chain
def is_term_type_mc(mol_id,ch_id,resno):
  sn=get_sn_from_resno(mol_id,ch_id,resno)
  if (type(sn) is int) and (is_polymer_residue(mol_id,ch_id,sn)==1) and (is_polymer_residue(mol_id,ch_id,sn+1)==1):
    resn_here=resno
    resn_next=seqnum_from_serial_number(mol_id,ch_id,sn+1)
    diff_next=resn_next-resn_here
    if (diff_next>=2):
      return 1
    else:
      return 0
  else:
    return 0
    

#check if res is at C-term side of break in mid chain
def is_term_type_mc_sn(mol_id,ch_id,sn):
  if sn in range(1,chain_n_residues(ch_id,mol_id)):
    resn_here=seqnum_from_serial_number(mol_id,ch_id,sn)
    if (is_polymer_residue(mol_id,ch_id,sn)==1) and (is_polymer_residue(mol_id,ch_id,sn+1)==1):
      resn_next=seqnum_from_serial_number(mol_id,ch_id,sn+1)
      diff_next=resn_next-resn_here
      if (diff_next>=2):
        return 1
      else:
        return 0
    else:
      return 0
  else:
    return 0
    
#check if res is at N-term side of break in mid chain
def is_term_type_mn(mol_id,ch_id,resno):
  sn=get_sn_from_resno(mol_id,ch_id,resno)
  if type(sn) is int:
    if sn in range(1,chain_n_residues(ch_id,mol_id)) and (is_polymer_residue(mol_id,ch_id,sn)==1):
      resn_here=resno
      resn_prev=seqnum_from_serial_number(mol_id,ch_id,sn-1)
      diff_prev=resn_here-resn_prev
      if (diff_prev>=2):
        return 1
      else:
        return 0
    else:
      return 0
  else:
    return 0
    
def is_term_type_mn_sn(mol_id,ch_id,sn):
  if sn in range(1,chain_n_residues(ch_id,mol_id)):
    resn_here=seqnum_from_serial_number(mol_id,ch_id,sn) 
    if (is_polymer_residue(mol_id,ch_id,sn)==1):
      resn_prev=seqnum_from_serial_number(mol_id,ch_id,sn-1)
      diff_prev=resn_here-resn_prev
      if (diff_prev>=2):
        return 1
      else:
        return 0
    else:
      return 0
  else:
    return 0
  
        
def first_residue_in_seg(mol_id,ch_id,resno):
  sn_here=get_sn_from_resno(mol_id,ch_id,resno)
  sn=sn_here
  while is_term_type_mn_sn(mol_id,ch_id,sn)==0 and sn>0 and seqnum_from_serial_number(mol_id,ch_id,sn-1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn-1)==1:
    sn=sn-1
  res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
  return res_start
  
def last_residue_in_seg(mol_id,ch_id,resno):
  sn_here=get_sn_from_resno(mol_id,ch_id,resno)
  sn=sn_here
  while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
    sn=sn+1
  res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
  return res_end


#Get serial number of active residue
def sn_of_active_res():
  residue=_active_residue_or_status()
  if not residue:
    return None
  sn=0
  resn_to_match=residue[2]
  mol_id=residue[0]
  ch_id=residue[1]
  current_resn=seqnum_from_serial_number(mol_id,"%s"%(ch_id),sn)
  while (current_resn!=resn_to_match):
    sn=sn+1
    current_resn=seqnum_from_serial_number(mol_id,"%s"%(ch_id),sn)
  return sn
    
#Get monomer and delete hydrogens
def get_monomer_no_H(mon):
  get_monomer(mon)
  delete_hydrogens(molecule_number_list()[-1])    
    
#Return list of segments in active mol
def segment_list(mol_id):
  sn=0
  list_out=[]
  for ch_id in chain_ids(mol_id):
    while is_polymer_residue(mol_id,ch_id,sn+1)==1:
      if sn==0:
        res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
        while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
          sn=sn+1
        res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
        list_out.append([mol_id,ch_id,res_start,res_end])
      else:
        while is_term_type_mn_sn(mol_id,ch_id,sn)==0:
          sn=sn+1
        res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
        while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
          sn=sn+1
        res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
        list_out.append([mol_id,ch_id,res_start,res_end])
    sn=0
  return list_out
  
def segment_list_chain(mol_id,ch_id):
  sn=0
  list_out=[]
  while is_polymer_residue(mol_id,ch_id,sn+1)==1:
    if sn==0:
      res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
      while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
        sn=sn+1
      res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
      list_out.append([mol_id,ch_id,res_start,res_end])
    else:
      while is_term_type_mn_sn(mol_id,ch_id,sn)==0:
        sn=sn+1
      res_start=seqnum_from_serial_number(mol_id,ch_id,sn)
      while is_term_type_mc_sn(mol_id,ch_id,sn)==0 and seqnum_from_serial_number(mol_id,ch_id,sn+1)!=-10000 and is_polymer_residue(mol_id,ch_id,sn+1)==1:
        sn=sn+1
      res_end=seqnum_from_serial_number(mol_id,ch_id,sn)
      list_out.append([mol_id,ch_id,res_start,res_end])
  return list_out
  
    
#rigid body refine zone here +/- n+1, then real space refine zone here +/-n
def auto_refine():
  if imol_refinement_map()==-1:
    info_dialog("You must set a refinement map!")
    return None
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  resno=residue[2]
  active_ins_code=residue[3]
  half_window=5
  contact_radius=4.0
  max_gap_to_fill=4
  minimum_fragment_size=3
  main_range_specs=set()
  seed_residue_specs=[]

  if _residue_is_polymer(mol_id, ch_id, resno, active_ins_code):
    fpr=first_polymer_residue(mol_id,ch_id)
    lpr=last_polymer_residue(mol_id,ch_id)
    res_start=max(fpr, resno-half_window)
    res_end=min(lpr, resno+half_window)
    for resn in range(res_start, res_end+1):
      if residue_exists_qm(mol_id,ch_id,resn,""):
        main_range_specs.add((ch_id,resn,""))
        seed_residue_specs.append([ch_id,resn,""])
  else:
    if not residue_exists_qm(mol_id, ch_id, resno, active_ins_code):
      _status_message("No active residue for local cylinder refinement")
      return None
    main_range_specs.add((ch_id,resno,active_ins_code))
    seed_residue_specs.append([ch_id,resno,active_ins_code])

  secondary_residue_specs=set()
  direct_non_polymer_contact_specs=set()
  for spec in seed_residue_specs:
    nearby_residues=residues_near_residue(mol_id, spec, contact_radius) or []
    for nearby_residue in nearby_residues:
      if nearby_residue and len(nearby_residue)==3:
        nearby_spec=(nearby_residue[0], nearby_residue[1], nearby_residue[2])
        if nearby_spec not in main_range_specs:
          secondary_residue_specs.add(nearby_spec)
          if not _residue_is_polymer(mol_id, nearby_spec[0], nearby_spec[1], nearby_spec[2]):
            direct_non_polymer_contact_specs.add(nearby_spec)

  secondary_polymer_specs=set([spec for spec in secondary_residue_specs
                               if _residue_is_polymer(mol_id, spec[0], spec[1], spec[2])])
  secondary_non_polymer_specs=(secondary_residue_specs-secondary_polymer_specs) | direct_non_polymer_contact_specs
  if not _residue_is_polymer(mol_id, ch_id, resno, active_ins_code):
    secondary_polymer_specs=_expand_polymer_contact_windows(mol_id, secondary_polymer_specs, 1)
  secondary_polymer_specs=_fill_short_polymer_gaps(mol_id, secondary_polymer_specs, max_gap_to_fill)
  secondary_polymer_specs=_prune_small_polymer_fragments(mol_id, secondary_polymer_specs, minimum_fragment_size)

  selected_residue_specs=main_range_specs | secondary_polymer_specs | secondary_non_polymer_specs
  residue_list=_sorted_residue_specs(selected_residue_specs)
  refine_residues(mol_id, residue_list)
  _status_message("Local cylinder refinement started")
   

#**** "Custom menu item functions ****
#Deletes active chain
def delete_chain():
  residue=_active_residue_or_status()
  if not residue:
    return None
  active_chain_id=residue[1]
  active_mol_id=residue[0]
  turn_off_backup(active_mol_id)
  while (is_polymer(active_mol_id,active_chain_id)==1) or (
  is_solvent_chain_p(active_mol_id,active_chain_id)!=-1):
    first_res=first_residue(active_mol_id,active_chain_id)
    last_res=last_residue(active_mol_id,active_chain_id)
    delete_residue_range(active_mol_id,active_chain_id,first_res,last_res)
  turn_on_backup(active_mol_id)

#Fits all polymer chains to map
def rigid_fit_all_chains():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  turn_off_backup(mol_id)
  for ch_id in chain_ids(mol_id): #Rigid body refine each chain
    if is_polymer(mol_id,ch_id)==1:
      rigid_body_refine_by_atom_selection(mol_id, "//%s//"%(ch_id))
      accept_regularizement()
  turn_on_backup(mol_id)
  
#Fits active chain to map
def rigid_fit_active_chain():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  rigid_body_refine_by_atom_selection(mol_id, "//%s//"%(ch_id))
  
#Copies active chain
def copy_active_chain():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  new_molecule_by_atom_selection(mol_id, "//%s//"%(ch_id))
  
#Cuts active chain
def cut_active_chain():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  new_molecule_by_atom_selection(mol_id, "//%s//"%(ch_id))
  turn_off_backup(mol_id)
  while (is_polymer(mol_id,ch_id)==1) or (is_solvent_chain_p(mol_id,ch_id)!=-1):
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    delete_residue_range(mol_id,ch_id,first_res,last_res)
  turn_on_backup(mol_id)

#Faster rigid body fit (not relevant for smaller ranges, but much faster for larger ranges)
def fast_rigid_fit(res_start,res_end,ch_id,mol_id):
  if (mol_id in model_molecule_number_list()) and (ch_id in chain_ids(mol_id)):
    turn_off_backup(mol_id)
    ins_code=""
    sn_max=chain_n_residues(ch_id,mol_id)-1
    res_min=seqnum_from_serial_number(mol_id,ch_id,0)
    res_max=seqnum_from_serial_number(mol_id,ch_id,sn_max)
    if res_start>res_end:
      res_start,res_end=res_end,res_start
    while not residue_exists_qm(mol_id,ch_id,res_start,ins_code) and res_start<=res_max:
     res_start=res_start+1
    while not residue_exists_qm(mol_id,ch_id,res_end,ins_code) and res_end>=res_min:
     res_end=res_end-1
    if res_start<=res_end:
      new_molecule_by_atom_selection(mol_id, "//{ch_id}/{res_start}-{res_end}/".format(ch_id=ch_id,res_start=res_start,res_end=res_end))
      mol_id_new=model_molecule_number_list()[-1]
      rigid_body_refine_by_atom_selection(mol_id_new,"/ /")
      accept_regularizement()
      delete_residue_range(mol_id,ch_id,res_start,res_end) #delete copied range from original mol
      merge_molecules([mol_id_new],mol_id) #Merge fit segment back into original mol
      change_chain_id_with_result(mol_id,chain_ids(mol_id)[-1],ch_id,1,res_start,res_end) #Merge chains
      close_molecule(mol_id_new)
    turn_on_backup(mol_id)
  else:
    info_dialog("The specified chain or molecule does not exist!")

#Copy active segment
def copy_active_segment():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  segments=segment_list(mol_id)
  res_here=residue[2]
  ch_id=residue[1]
  for seg in segments:
    if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      new_molecule_by_atom_selection(mol_id, "//{ch_id}/{res_start}-{res_end}/".format(ch_id=ch_id,res_start=res_start,res_end=res_end))

#Cut active segment
def cut_active_segment():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  segments=segment_list(mol_id)
  res_here=residue[2]
  ch_id=residue[1]
  for seg in segments:
    if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      new_molecule_by_atom_selection(mol_id, "//{ch_id}/{res_start}-{res_end}/".format(ch_id=ch_id,res_start=res_start,res_end=res_end))
      delete_residue_range(mol_id,ch_id,res_start,res_end)

#Delete active segment
def delete_active_segment():
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  segments=segment_list(mol_id)
  res_here=residue[2]
  ch_id=residue[1]
  for seg in segments:
    if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      delete_residue_range(mol_id,ch_id,res_start,res_end)


#Jiggle-fits active chain to map
def jiggle_fit_active_chain():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    residue=_active_residue_or_status()
    if not residue:
      return None
    mol_id=residue[0]
    ch_id=residue[1]
    fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)

#Jiggle-fit active chain to B-smoothed map
def jiggle_fit_active_chain_smooth():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    sharpen(imol_refinement_map(),200)
    residue=_active_residue_or_status()
    if not residue:
      sharpen(imol_refinement_map(),0)
      return None
    mol_id=residue[0]
    ch_id=residue[1]
    fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)
    sharpen(imol_refinement_map(),0)
    
#Jiggle-fits active chain to map (more thorough)
def jiggle_fit_active_chain_slow():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    score_string=""
    for scale_fac in range(1,1000,100):
      scale_fac2=scale_fac/10.0
      score=fit_chain_to_map_by_random_jiggle(mol_id,ch_id,100,scale_fac2)
      score_string=score_string+"scale:"+str(scale_fac2)+"score:"+str(score)
    print(score_string)
    
#Fits all polymer chains to map
def jiggle_fit_all_chains():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    mol_id=_active_molecule_or_status()
    if mol_id is None:
      return None
    turn_off_backup(mol_id)
    for ch_id in chain_ids(mol_id): #Rigid body refine each chain
      if is_polymer(mol_id,ch_id)==1:
        fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)
        accept_regularizement()
    turn_on_backup(mol_id)

#Jiggle-fit current molecule to map
def jiggle_fit_active_mol():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    mol_id=_active_molecule_or_status()
    if mol_id is None:
      return None
    fit_molecule_to_map_by_random_jiggle(mol_id,1000,0.1)
    
#Clear labels and distances
def clear_distances_and_labels():
  remove_all_atom_labels()
  try:
    clear_measure_distances()
  except:
    clear_simple_distances()
  
#Delete hydrogens from active molecule
def delete_h_active():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  delete_hydrogens(mol_id)
  
def fit_this_segment():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    mol_id=active_residue()[0]
    segments=segment_list(mol_id)
    res_here=active_residue()[2]
    ch_id=active_residue()[1]
    for seg in segments:
      if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
        res_start=seg[2]
        res_end=seg[3]
        ch_id=seg[1]
        turn_off_backup(mol_id)
        set_refinement_immediate_replacement(1)
        rigid_body_refine_zone(res_start,res_end,ch_id,mol_id)
        accept_regularizement()
        set_refinement_immediate_replacement(0)
        turn_on_backup(mol_id)


#Real space refine for keyboard shortcut
def refine_click():
  def refine_click_a(res1,res2):
    mol_id_1=res1[1]
    mol_id_2=res2[1]
    ch_id_1=res1[2]
    ch_id_2=res2[2]
    resno_1=res1[3]
    resno_2=res2[3]
    alt_conf_1=res1[4]
    alt_conf_2=res2[4]
    if resno_2>=resno_1:
      refine_zone(mol_id_1,ch_id_1,resno_1,resno_2,alt_conf_1)
    else:
      refine_zone(mol_id_1,ch_id_1,resno_2,resno_1,alt_conf_1)
  user_defined_click(2,refine_click_a)
  
def refine_residues_click():
  def refine_residues_click_a(res1,res2):
    mol_id_1=res1[1]
    mol_id_2=res2[1]
    ch_id_1=res1[2]
    ch_id_2=res2[2]
    resno_1=res1[3]
    resno_2=res2[3]
    ins_code_1=res1[6]
    if (mol_id_1==mol_id_2) and (ch_id_1==ch_id_2):
      if resno_1>=resno_2:
        resno_1,resno_2=resno_2,resno_1 #Swap res numbers if wrong way round
      refine_zone(mol_id_1,ch_id_1,resno_1,resno_2,ins_code_1)
    else:
      info_dialog("Sorry, start and end residues must be in same mol and chain!")
  user_defined_click(2,refine_residues_click_a)
  
def refine_residues_range(mol_id,ch_id,resno_1,resno_2,ins_code_1):
  if resno_1>=resno_2:
    resno_1,resno_2=resno_2,resno_1 #Swap res numbers if wrong way round
  res_list=[]
  for resn in range(resno_1,resno_2+1):
    res_triple=[ch_id,resn,ins_code_1] #Does not account for varying ins code. Should be able to fix using residues_matching_criteria(), e.g.:
    #residues_matching_criteria(0, lambda chain_id,resno,ins_code,serial: True), except substitute a test func for true
    #that evaluates to true if resno, mol id and ch_id match, and returns ins_code (third item in output triple)
    res_list.append(res_triple) #append rather than adding, bc making list of lists
  refine_residues(mol_id,res_list)

  
  
#refine range, plus everything in contact with it
def refine_residues_sphere_click():
  def refine_residues_sphere_click_a(res1,res2):
    mol_id_1=res1[1]
    mol_id_2=res2[1]
    ch_id_1=res1[2]
    ch_id_2=res2[2]
    resno_1=res1[3]
    resno_2=res2[3]
    ins_code_1=res1[6]
    ins_code_2=res2[6]
    if (mol_id_1==mol_id_2) and (ch_id_1==ch_id_2) and (imol_refinement_map!=-1):
      if resno_1>=resno_2:
        resno_1,resno_2=resno_2,resno_1 #Swap res numbers if wrong way round
      res_list=[]
      nearby_residues=[]
      for resn in range(resno_1,resno_2+1):
        res_triple=[ch_id_1,resn,ins_code_1] #Does not account for varying ins code. Should be able to fix using residues_matching_criteria(), e.g.:
        #residues_matching_criteria(0, lambda chain_id,resno,ins_code,serial: True), except substitute a test func for true
        #that evaluates to true if resno, mol id and ch_id match, and returns ins_code (third item in output triple)
        res_list.append(res_triple) #append rather than adding, bc making list of lists
      for resn in range(resno_1+1,resno_2-1):
        res_triple=[ch_id_1,resn,ins_code_1] #Does not account for varying ins code. Should be able to fix using residues_matching_criteria(), e.g.:
        new_nearby_residues=residues_near_residue(mol_id_1,res_triple,5)
        nearby_residues.extend(new_nearby_residues)
      if len(res_list)<=2:
        for resn in range(resno_1,resno_2+1):
          res_triple=[ch_id_1,resn,ins_code_1] #Does not account for varying ins code. Should be able to fix using residues_matching_criteria(), e.g.:
          new_nearby_residues=residues_near_residue(mol_id_1,res_triple,5)
          nearby_residues.extend(new_nearby_residues)
      res_list=res_list+nearby_residues
      print(res_list)
      res_list=[list(x) for x in set(tuple(x) for x in res_list)] #Uniquifies list of lists
      print("hopefully unique:",res_list)
      refine_residues(mol_id_1,res_list)
      print(res_list)
    else:
      info_dialog("Sorry, start and end residues must be in same mol and chain! And you need to have a refinement map.")
  user_defined_click(2,refine_residues_sphere_click_a)

#Python version of sphere refine from here:
#https://www.mail-archive.com/coot@jiscmail.ac.uk/msg01463.html
def sphere_refine_active(radius):
    from types import ListType
    active_atom = active_residue()
    if (not active_atom):
        add_status_bar_text("No active residue")
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        centred_residue = active_atom[1:4]
        other_residues = residues_near_residue(imol, centred_residue, radius)
        all_residues = [centred_residue]
        if (type(other_residues) is ListType):
            all_residues += other_residues
        print "imol: %s residues: %s" %(imol, all_residues)
        refine_residues(imol, all_residues)

def stepped_sphere_refine(mol_id,ch_id):
  turn_off_backup(mol_id)
  set_refinement_immediate_replacement(1)
  valid_resnames=['A','C','T','G','U','ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
  if is_polymer(mol_id,ch_id): 
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    if is_protein_chain_p(mol_id,ch_id):
      for res in range(first_res,last_res+1):
        res_exist_flag=does_residue_exist_p(mol_id,ch_id,res,"")
        resname=residue_name(mol_id,ch_id,res,"")
        if (res_exist_flag==1) and (resname in valid_resnames):
          set_go_to_atom_chain_residue_atom_name(ch_id,res," CA ")
          sphere_refine_active(7)
          accept_regularizement()
    elif is_nucleotide_chain_p(mol_id,ch_id):
      for res in range(first_res,last_res+1):
        res_exist_flag=does_residue_exist_p(mol_id,ch_id,res,"")
        resname=residue_name(mol_id,ch_id,res,"")
        if (res_exist_flag==1) and (resname in valid_resnames):
          set_go_to_atom_chain_residue_atom_name(ch_id,res," P ")
          sphere_refine_active(7)
          accept_regularizement()
  set_refinement_immediate_replacement(0)
  turn_on_backup(mol_id)
  info_dialog("Refinement finished - all done!")

  
#Stepped v of cylinder refine:
#get mol from active residue
#for ch_id in mol:
#cylinder_refine four res segments, with two res overlap
#ignore breaks?
#what about waters?


#Copy fragment (click start and end)
def copy_frag_by_click():
  def copy_frag(res1,res2):
    mol_id_1=res1[1]
    mol_id_2=res2[1]
    ch_id_1=res1[2]
    ch_id_2=res2[2]
    resno_1=res1[3]
    resno_2=res2[3]
    if (mol_id_1==mol_id_2) and (ch_id_1==ch_id_2):
      if resno_1>resno_2:
        atom_sel="//%s/%s-%s" %(ch_id_1,resno_2,resno_1)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
      elif resno_2>resno_1:
        atom_sel="//%s/%s-%s" %(ch_id_1,resno_1,resno_2)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
      else:
        atom_sel="//%s/%s" %(ch_id_1,resno_2)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
    else:
      info_dialog("Start and end residues must be in the same molecule and chain!")
  user_defined_click(2,copy_frag)
  
#Cut fragment (click start and end)
def cut_frag_by_click():
  def cut_frag(res1,res2):
    mol_id_1=res1[1]
    mol_id_2=res2[1]
    ch_id_1=res1[2]
    ch_id_2=res2[2]
    resno_1=res1[3]
    resno_2=res2[3]
    if (mol_id_1==mol_id_2) and (ch_id_1==ch_id_2):
      turn_off_backup(mol_id_1)
      if resno_1>resno_2:
        atom_sel="//%s/%s-%s" %(ch_id_1,resno_2,resno_1)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
        delete_residue_range(mol_id_1,ch_id_1,resno_2,resno_1)
      elif resno_2>resno_1:
        atom_sel="//%s/%s-%s" %(ch_id_1,resno_1,resno_2)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
        delete_residue_range(mol_id_1,ch_id_1,resno_1,resno_2)
      else:
        atom_sel="//%s/%s" %(ch_id_1,resno_2)
        new_molecule_by_atom_selection(mol_id_1,atom_sel)
        delete_residue_range(mol_id_1,ch_id_1,resno_2)
      turn_on_backup(mol_id_1)
    else:
      info_dialog("Start and end residues must be in the same molecule and chain!")
  user_defined_click(2,cut_frag)


#Delete sidechain range (click start and end)
def delete_sidechain_range_by_click_a():
  def delete_sidechain_range_by_click_b(res1,res2):
    if (res1[1]!=res2[1]) or (res1[2]!=res2[2]) or (res1[3]==res2[3]):
      info_dialog("Start and end points must be in the same mol and chain!")
    else: 
      if (res1[3] > res2[3]):
        res_start=res2[3]-1
        res_end=res1[3]+1
      else:
        res_start=res1[3]-1
        res_end=res2[3]+1
      mol_id=res1[1]
      ch_id=res1[2]
      turn_off_backup(mol_id)
      delete_sidechain_range(mol_id,ch_id,res_start,res_end)
      turn_on_backup(mol_id)
  user_defined_click(2,delete_sidechain_range_by_click_b)
  
#Mutate range to poly-unk
def mutate_residue_range_by_click_a():
  def mutate_residue_range_by_click_b(res1,res2):
    if (res1[1]!=res2[1]) or (res1[2]!=res2[2]) or (res1[3]==res2[3]):
      info_dialog("Start and end points must be in the same mol and chain!")
    else:
      if (res1[3] > res2[3]):
        res_start=res2[3]
        res_end=res1[3]
        n=res_end-res_start+1
      else:
        res_start=res1[3]
        res_end=res2[3]
        n=res_end-res_start+1
      mol_id=res1[1]
      ch_id=res1[2]
      target_seq=n*"A"
      turn_off_backup(mol_id)
      mutate_residue_range(mol_id,ch_id,res_start,res_end,target_seq)
      for n in range(res_start,res_end+1):
        set_residue_name(mol_id,ch_id,n,"","UNK")
      turn_on_backup(mol_id)
  user_defined_click(2,mutate_residue_range_by_click_b)

#Mutate range to polyala
def mutate_residue_range_by_click_ala_a():
  def mutate_residue_range_by_click_ala_b(res1,res2):
    if (res1[1]!=res2[1]) or (res1[2]!=res2[2]) or (res1[3]==res2[3]):
      info_dialog("Start and end points must be in the same mol and chain!")
    else:
      if (res1[3] > res2[3]):
        res_start=res2[3]
        res_end=res1[3]
        n=res_end-res_start+1
      else:
        res_start=res1[3]
        res_end=res2[3]
        n=res_end-res_start+1
      mol_id=res1[1]
      ch_id=res1[2]
      target_seq=n*"A"
      turn_off_backup(mol_id)
      mutate_residue_range(mol_id,ch_id,res_start,res_end,target_seq)
      for n in range(res_start,res_end+1):
        set_residue_name(mol_id,ch_id,n,"","ALA")
      turn_on_backup(mol_id)
  user_defined_click(2,mutate_residue_range_by_click_ala_b)
  
#Merge two fragments
def merge_fragments():
  def merge_2_fragments(res1,res2):
    mol_daughter=[res2[1]]
    mol_ref=res1[1]
    _ensure_non_polymer_restraints_loaded_from_molecule(res2[1])
    merge_molecules(mol_daughter,mol_ref)
    toggle_display_mol(mol_ref)
    toggle_display_mol(mol_ref)
  user_defined_click(2,merge_2_fragments)
  
#Force addition of residue - useful when
#Coot says "No acceptable position found"
# but density is clear.
def force_add_terminal_residue():
  def force_addition(res1):
    mol_id=res1[1]
    ch_id=res1[2]
    res_no=res1[3]
    res_type="auto"
    add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
    res_type,-57.82,-47)
    if residue_exists_qm(mol_id,ch_id,res_no+1,res1[4]):
      set_b_factor_residue_range(mol_id,ch_id,res_no+1,res_no+1,default_new_atoms_b_factor())
    elif residue_exists_qm(mol_id,ch_id,res_no-1,res1[4]):
      set_b_factor_residue_range(mol_id,ch_id,res_no-1,res_no-1,default_new_atoms_b_factor())
    sort_residues(mol_id)
  user_defined_click(1,force_addition)
  
def force_add_terminal_residue_noclick(mol_id,ch_id,res_no):
  res_type="auto"
  add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
  res_type,-57.82,-47)
  if residue_exists_qm(mol_id,ch_id,res_no+1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no+1,res_no+1,default_new_atoms_b_factor())
  elif residue_exists_qm(mol_id,ch_id,res_no-1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no-1,res_no-1,default_new_atoms_b_factor())
  sort_residues(mol_id)

def force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,res_no,phi,psi):
  res_type="auto"
  add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
  res_type,float(phi),float(psi))
  if residue_exists_qm(mol_id,ch_id,res_no+1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no+1,res_no+1,default_new_atoms_b_factor())
  elif residue_exists_qm(mol_id,ch_id,res_no-1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no-1,res_no-1,default_new_atoms_b_factor())
  sort_residues(mol_id)

def force_add_terminal_residue_noclick_strand(mol_id,ch_id,res_no):
  res_type="auto"
  add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
  res_type,-139,135)
  if residue_exists_qm(mol_id,ch_id,res_no+1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no+1,res_no+1,default_new_atoms_b_factor())
  elif residue_exists_qm(mol_id,ch_id,res_no-1,""):
    set_b_factor_residue_range(mol_id,ch_id,res_no-1,res_no-1,default_new_atoms_b_factor())
  sort_residues(mol_id)

#Paul
def key_binding_terminal_spin():
    active_atom = active_residue()
    if not active_atom:
       print("No active atom")
    else:
       imol = active_atom[0]
       residue_spec = atom_spec_to_residue_spec(active_atom)
       print ('spin_N {} {} {}'.format(imol, residue_spec, 120))
       spin_N_py(imol, residue_spec, 120)

#This works, but should alter to be more flexible. ideally, have one key to tweak phi, one to tweak psi.
#Will need two global cycle variables - residue_phi_cycle and residue_psi_cycle - as well as two variables to keep the current value of phi and psi.
residue_phi_cycle=0
residue_psi_cycle=0
current_phi=-60
current_psi=-50
#should chamnge this to incorporate measurement of starting phi/psi using get_torsion (and maybe setting using set_torsion?)
#get_torsion(0,["A",2393,""," C  ",""], ["A",2394,"", " N  ", ""], ["A", 2394, "", " CA ", ""], ["A", 2394, "", " C  ",""])
#set_torsion(imol, chain_id, res_no, ins_code_alt_conf, atom_name_1,atom_name_2, atom_name_3, atom_name_4)
def cycle_residue_phi():
  global residue_phi_cycle
  global current_phi
  global current_psi
  res_type="auto"
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resn=active_residue()[2]
  atom_name=active_residue()[4]
  set_go_to_atom_molecule(mol_id)
  ins_code=""
  first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
  delta_first=abs(first_in_seg-resn)
  delta_last=abs(last_in_seg-resn)
  set_new_atom_b_fac_to_mean()
  if delta_first<=delta_last:
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
    if (residue_phi_cycle==0):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-180
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=1
    elif (residue_phi_cycle==1):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-150
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=2
    elif (residue_phi_cycle==2):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-120
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=3
    elif (residue_phi_cycle==3):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-90
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=4
    elif (residue_phi_cycle==4):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-60
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=5
    elif (residue_phi_cycle==5):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=-30
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=6
    elif (residue_phi_cycle==6):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=0
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=7
    elif (residue_phi_cycle==7):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=30
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=8
    elif (residue_phi_cycle==8):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=60
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=9
    elif (residue_phi_cycle==9):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=90
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=10
    elif (residue_phi_cycle==10):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=120
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_phi_cycle=0
  else:
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
    if (residue_phi_cycle==0):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-180
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=1
    elif (residue_phi_cycle==1):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-150
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=2
    elif (residue_phi_cycle==2):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-120
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=3
    elif (residue_phi_cycle==3):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-90
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=4
    elif (residue_phi_cycle==4):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-60
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=5
    elif (residue_phi_cycle==5):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=-30
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=6
    elif (residue_phi_cycle==6):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=0
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=7
    elif (residue_phi_cycle==7):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=30
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=8
    elif (residue_phi_cycle==8):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=60
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=9
    elif (residue_phi_cycle==9):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=90
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=10
    elif (residue_phi_cycle==10):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=120
      psi=current_psi
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_phi_cycle=0
  current_phi=phi
def cycle_residue_psi():
  global residue_psi_cycle
  global current_phi
  global current_psi
  res_type="auto"
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resn=active_residue()[2]
  atom_name=active_residue()[4]
  set_go_to_atom_molecule(mol_id)
  ins_code=""
  first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
  delta_first=abs(first_in_seg-resn)
  delta_last=abs(last_in_seg-resn)
  set_new_atom_b_fac_to_mean()
  if delta_first<=delta_last:
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
    if (residue_psi_cycle==0):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-180
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=1
    elif (residue_psi_cycle==1):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-150
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=2
    elif (residue_psi_cycle==2):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-120
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=3
    elif (residue_psi_cycle==3):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-90
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=4
    elif (residue_psi_cycle==4):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-60
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=5
    elif (residue_psi_cycle==5):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=-30
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=6
    elif (residue_psi_cycle==6):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=0
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=7
    elif (residue_psi_cycle==7):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=30
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=8
    elif (residue_psi_cycle==8):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=60
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=9
    elif (residue_psi_cycle==9):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=90
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=10
    elif (residue_psi_cycle==10):
      delete_residue(mol_id,ch_id,first_in_seg,ins_code)
      phi=current_phi
      psi=120
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,first_in_seg+1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      residue_psi_cycle=0
  else:
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
    if (residue_psi_cycle==0):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-180
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=1
    elif (residue_psi_cycle==1):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-150
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=2
    elif (residue_psi_cycle==2):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-120
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=3
    elif (residue_psi_cycle==3):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-90
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=4
    elif (residue_psi_cycle==4):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-60
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=5
    elif (residue_psi_cycle==5):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=-30
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=6
    elif (residue_psi_cycle==6):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=0
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=7
    elif (residue_psi_cycle==7):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=30
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=8
    elif (residue_psi_cycle==8):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=60
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=9
    elif (residue_psi_cycle==9):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=90
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=10
    elif (residue_psi_cycle==10):
      delete_residue(mol_id,ch_id,last_in_seg,ins_code)
      phi=current_phi
      psi=120
      force_add_terminal_residue_noclick_phi_psi(mol_id,ch_id,last_in_seg-1,phi,psi)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      residue_psi_cycle=0
  current_psi=psi


def add_term_shortcut_force():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resn=active_residue()[2]
  atom_name=active_residue()[4]
  set_go_to_atom_molecule(mol_id)
  first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
  delta_first=abs(first_in_seg-resn)
  delta_last=abs(last_in_seg-resn)
  set_new_atom_b_fac_to_mean()
  if delta_first<=delta_last:
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,"CA")
    force_add_terminal_residue_noclick(mol_id,ch_id,first_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg-1,"CA")
  else:
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,"CA")
    force_add_terminal_residue_noclick(mol_id,ch_id,last_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg+1,"CA")



#Grow helix from selected terminus
def grow_helix():
  def grow_helix_post_click(res1):
    def grow_helix_enter_resn(n):
      n_res=_positive_int_from_entry(n,"Number of helix residues")
      if n_res is None:
        return None
      mol_id=res1[1]
      ch_id=res1[2]
      res_no=res1[3]
      _grow_from_selected_terminus(mol_id,ch_id,res_no,n_res,-57.82,-47,1)
    generic_single_entry("How many residues for helix?",
    "10","Grow helix",grow_helix_enter_resn)
  user_defined_click(1,grow_helix_post_click)
  
#Grow strand from selected terminus
def grow_strand():
  def grow_strand_post_click(res1):
    def grow_strand_enter_resn(n):
      n_res=_positive_int_from_entry(n,"Number of strand residues")
      if n_res is None:
        return None
      mol_id=res1[1]
      ch_id=res1[2]
      res_no=res1[3]
      _grow_from_selected_terminus(mol_id,ch_id,res_no,n_res,-139,135,1)
    generic_single_entry("How many residues for strand?",
    "10","Grow strand",grow_strand_enter_resn)
  user_defined_click(1,grow_strand_post_click)
  
#Grow para strand from selected terminus
def grow_parallel_strand():
  def grow_parallel_strand_post_click(res1):
    def grow_parallel_strand_enter_resn(n):
      n_res=_positive_int_from_entry(n,"Number of parallel-strand residues")
      if n_res is None:
        return None
      mol_id=res1[1]
      ch_id=res1[2]
      res_no=res1[3]
      _grow_from_selected_terminus(mol_id,ch_id,res_no,n_res,-119,113,0)
    generic_single_entry("How many residues for parallel strand?",
    "10","Grow parallel strand",grow_parallel_strand_enter_resn)
  user_defined_click(1,grow_parallel_strand_post_click)

#Grow 3-10 helix from selected terminus
def grow_helix_3_10():
  def grow_helix_post_click(res1):
    def grow_helix_enter_resn(n):
      n_res=_positive_int_from_entry(n,"Number of 3-10 helix residues")
      if n_res is None:
        return None
      mol_id=res1[1]
      ch_id=res1[2]
      res_no=res1[3]
      _grow_from_selected_terminus(mol_id,ch_id,res_no,n_res,-49,-26,0)
    generic_single_entry("How many residues for helix?",
    "10","Grow 3-10 helix",grow_helix_enter_resn)
  user_defined_click(1,grow_helix_post_click)

#Renumbers the active chain (from active_residue()). User enters new starting residue number.
def renumber_by_first_res():
  def renum_chain(new_num):
    new_num=int(new_num)
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    offset=new_num-first_res
    renumber_residue_range(mol_id,ch_id,first_res,last_res,int(offset))
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New number for first residue in this chain?",
  str(seqnum_from_serial_number(active_residue()[0],
  "%s"%(active_residue()[1]),0)),"Renumber",renum_chain)
  
#Renumbers the active chain (from active_residue()). User enters new last residue number.
def renumber_by_last_res():
  def renum_chain(new_num):
    new_num=int(new_num)
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    offset=new_num-last_res
    renumber_residue_range(mol_id,ch_id,first_res,last_res,int(offset))
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New number for last residue in this chain?",
  str(seqnum_from_serial_number(active_residue()[0],
  "%s"%(active_residue()[1]),
  (chain_n_residues(active_residue()[1],active_residue()[0])-1))),"Renumber",renum_chain)
  
#Renumbers the active chain (from active_residue()). User enters new active residue number.
def renumber_by_active_res():
  def renum_chain(new_num):
    new_num=int(new_num)
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    current_num=active_residue()[2]
    first_res=first_residue(mol_id,ch_id)
    last_res=last_residue(mol_id,ch_id)
    offset=new_num-current_num
    renumber_residue_range(mol_id,ch_id,first_res,last_res,int(offset))
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New number for this residue?",
  str(active_residue()[2]),"Renumber",renum_chain)

#Renumber from N-term to current residue
def renumber_n_term_segment():
  def renumber_n_term_segment_entry(new_resn):
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    resn=active_residue()[2]
    first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
    offset=int(new_resn)-resn
    renumber_residue_range(mol_id,ch_id,first_res,resn,offset)
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New residue number?",
  str(active_residue()[2]),"Renumber",renumber_n_term_segment_entry)
  
#Renumber from current residue to C-term
def renumber_c_term_segment():
  def renumber_c_term_segment_entry(new_resn):
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    resn=active_residue()[2]
    def last_residue(mol_id,ch_id):
      n=chain_n_residues(ch_id,mol_id)-1
      result=seqnum_from_serial_number(mol_id,"%s"%(ch_id),n)
      return result
    last_res=last_residue(mol_id,ch_id)
    offset=int(new_resn)-resn
    renumber_residue_range(mol_id,ch_id,resn,last_res,offset)
    delete_all_extra_restraints(mol_id)
    set_show_extra_restraints(mol_id,0)
    set_show_extra_restraints(mol_id,1)
  generic_single_entry("New residue number?",
  str(active_residue()[2]),"Renumber",renumber_c_term_segment_entry)
  
#Sharpen or blur map. Slider jerky on large maps.
def sharpen_by_entered_factor():
  if (scroll_wheel_map()==-1):
    info_dialog("You need a map!")
  else:
    def sharpen_by_entry(B):
      B=float(B)
      sharpen(scroll_wheel_map(),B)
    generic_single_entry("New B-factor for map?",
    "-100","Sharpen/Blur",sharpen_by_entry)
    
#Reload map with different hi-res limit
def change_hires_limit():
  if (scroll_wheel_map()==-1):
    info_dialog("You need a map!")
  else:
    def change_hires_by_entry(new_res):
      mol=scroll_wheel_map()
      new_res=float(new_res)
      mtz_file=map_parameters(mol)[0]
      F_col=map_parameters(mol)[1]
      PHI_col=map_parameters(mol)[2]
      make_and_draw_map_with_reso_with_refmac_params(mtz_file,F_col,PHI_col,"",0,0,0,"Fobs:None-specified",
      "SigF:None-specified","RFree:None-specified",0,0,1,1000.0,new_res)
      close_molecule(mol)
    generic_single_entry("New high-res limit for map?",
    "5.0","Change high resolution limit for active map",change_hires_by_entry)
    
def change_hires_limit_copy():
  if (scroll_wheel_map()==-1):
    info_dialog("You need a map!")
  else:
    def change_hires_by_entry(new_res):
      mol=scroll_wheel_map()
      new_res=float(new_res)
      mtz_file=map_parameters(mol)[0]
      F_col=map_parameters(mol)[1]
      PHI_col=map_parameters(mol)[2]
      make_and_draw_map_with_reso_with_refmac_params(mtz_file,F_col,PHI_col,"",0,0,0,"Fobs:None-specified",
      "SigF:None-specified","RFree:None-specified",0,0,1,1000.0,new_res)
    generic_single_entry("New high-res limit for map?",
    "5.0","Make low pass filtered copy of active map",change_hires_by_entry)

#Mutate mets to MSE
def mutate_all_mets_to_mse():
  mol_id=active_residue()[0]
  turn_off_backup(mol_id)
  for ch_id in chain_ids(mol_id):
    for resn in range(0,chain_n_residues(ch_id,mol_id)):
      if (resname_from_serial_number(mol_id,ch_id,resn)=="MET"):
        seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,resn))
        delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"B")
        delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"C")
        mutate(mol_id,ch_id,seqnum,ins_id,"MSE")
        map_id=imol_refinement_map()
        if (map_id!=-1):
          auto_fit_best_rotamer(seqnum,"",ins_id,ch_id,mol_id,map_id,1,0.01)
  turn_on_backup(mol_id)
  
#Mutate MSEs to MET
def mutate_all_mse_to_met():
  mol_id=active_residue()[0]
  turn_off_backup(mol_id)
  for ch_id in chain_ids(mol_id):
    for resn in range(0,chain_n_residues(ch_id,mol_id)):
      if (resname_from_serial_number(mol_id,ch_id,resn)=="MSE"):
        seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,resn))
        delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"B")
        delete_residue_with_full_spec(mol_id,1,ch_id,seqnum,ins_id,"C")
        mutate(mol_id,ch_id,seqnum,ins_id,"MET")
        map_id=imol_refinement_map()
        if (map_id!=-1):
          auto_fit_best_rotamer(seqnum,"",ins_id,ch_id,mol_id,map_id,1,0.01)
  turn_on_backup(mol_id)
  
#Shorten loop by one residue
def shorten_loop():
  active_atom=_active_residue_or_status()
  if not active_atom:
    return None
  mol_id=active_atom[0]
  ch_id=active_atom[1]
  resn=active_atom[2]
  ins_code=active_atom[3]
  delete_residue(mol_id,ch_id,resn,ins_code)
  first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
  renumber_residue_range(mol_id,ch_id,first_res,resn,1)
  delete_all_extra_restraints(mol_id)
  set_show_extra_restraints(mol_id,0)
  set_show_extra_restraints(mol_id,1)
  r1=resn-1
  r2=resn+2
  set_refinement_immediate_replacement(1)
  try:
    refine_zone(mol_id,ch_id,r1,r2,"")
    accept_regularizement()
  finally:
    set_refinement_immediate_replacement(0)

#Lengthen loop by one residue
def lengthen_loop():
  _patch_gap_low_density_average()
  active_atom=active_residue()
  mol_id=active_atom[0]
  ch_id=active_atom[1]
  resn=active_atom[2]
  first_res=seqnum_from_serial_number(mol_id,"%s"%(ch_id),0)
  renumber_residue_range(mol_id,ch_id,first_res,resn,-1)
  delete_all_extra_restraints(mol_id)
  set_show_extra_restraints(mol_id,0)
  set_show_extra_restraints(mol_id,1)
  r1=resn-1
  r2=resn
  set_refinement_immediate_replacement(1)
  fit_gap(mol_id,ch_id,r1,r2,"A",1)
  accept_regularizement()
  set_refinement_immediate_replacement(0)

#Get fractional coordinates of active atom. Useful when inspecting heavy atom sites.
def get_fract_coords():
  a=_active_residue_or_status()
  if not a:
    return None
  atom_info=atom_specs(a[0],a[1],a[2],a[3],a[4],a[5])
  if not atom_info:
    info_dialog("Unable to obtain the active atom coordinates.")
    return None
  x_cart=atom_info[3]
  y_cart=atom_info[4]
  z_cart=atom_info[5]
  mol_id=a[0]
  cell_a=cell(mol_id)[0]
  cell_b=cell(mol_id)[1]
  cell_c=cell(mol_id)[2]
  cell_alpha=math.radians(cell(mol_id)[3])
  cell_beta=math.radians(cell(mol_id)[4])
  cell_gamma=math.radians(cell(mol_id)[5])
  cos_alpha_star=(math.cos(cell_beta)*math.cos(cell_gamma)-math.cos(cell_alpha))/(math.sin(cell_beta)*math.sin(cell_gamma))
  sin_alpha_star=math.sqrt(1-cos_alpha_star**2)
  z_fract=z_cart/(cell_c*math.sin(cell_beta)*sin_alpha_star)
  y_fract=(y_cart-(-1*cell_c*math.sin(cell_beta)*cos_alpha_star)*z_fract)/(cell_b*math.sin(cell_gamma))
  x_fract=(x_cart-(cell_b*math.cos(cell_gamma))*y_fract-(cell_c*math.cos(cell_beta))*z_fract)/cell_a
  x_y_z_string=str("("+str("%.3f" % x_fract)+","+str("%.3f" % y_fract)+","+str("%.3f" % z_fract)+")")
  info_dialog("The fractional coordinates of the active atom are: %s"%(x_y_z_string))
  
#Go to center of scroll wheel map.
def goto_center_of_map():
  if (scroll_wheel_map()==-1):
    info_dialog("You need a map!")
  else:
    mol_id=scroll_wheel_map()
    a=float(cell(mol_id)[0])
    b=float(cell(mol_id)[1])
    c=float(cell(mol_id)[2])
    x=a*0.5
    y=b*0.5
    z=c*0.5
    set_rotation_centre(x,y,z)
    
#Make a button list for inserting common monomers
def pick_common_monomers():
  get_ddm=["DDM", lambda func: get_monomer_no_H("LMT")]
  get_dm=["DM", lambda func: get_monomer_no_H("DMU")]
  get_bog=["OG", lambda func: get_monomer_no_H("HSH")]
  get_ldao=["LDAO", lambda func: get_monomer_no_H("LDA")]
  get_mpg=["Monoolein", lambda func: get_monomer_no_H("MPG")]
  get_glycerol=["Glycerol", lambda func: get_monomer_no_H("GOL")]
  get_eg=["Ethylene glycol", lambda func: get_monomer_no_H("EDO")]
  get_acetate=["Acetate", lambda func: get_monomer_no_H("ACT")]
  get_dmso=["DMSO", lambda func: get_monomer_no_H("DMS")]
  get_tris=["Tris", lambda func: get_monomer_no_H("TAM")]
  get_hepes=["HEPES", lambda func: get_monomer_no_H("EPE")]
  get_mes=["MES", lambda func: get_monomer_no_H("MES")]
  get_cac=["Cacodylate", lambda func: get_monomer_no_H("CAD")]
  get_peg=["PEG", lambda func: get_monomer_no_H("1PE")]
  get_popg=["POPG (lipid)", lambda func: get_monomer_no_H("LHG")]
  get_pe=["PE (lipid)", lambda func: get_monomer_no_H("PEF")]
  get_pc=["PC (lipid)", lambda func: get_monomer_no_H("PLC")]
  get_ps=["PS (lipid)", lambda func: get_monomer_no_H("PSF")]
  get_pip2=["PI(3,4)P2 (lipid)", lambda func: get_monomer_no_H("52N")]
  get_chs=["Cholesterol hemisuccinate", lambda func: get_monomer_no_H("Y01")]
  get_chl=["Cholesterol", lambda func: get_monomer_no_H("CLR")]
  button_list=[get_acetate,get_eg,get_glycerol,get_dmso,get_ddm,get_dm,get_bog,get_ldao,get_mpg,get_tris,get_hepes,get_mes,get_cac,get_peg,get_popg,get_pe,get_pc,get_ps,get_pip2,get_chs,get_chl]
  generic_button_dialog("Common small molecules",button_list)
  
#Switch all models to CA representation
def all_mols_to_ca():
  for mol_id in molecule_number_list():
    graphics_to_ca_plus_ligands_representation(mol_id)
    
#Set b-factor color scaling based on mean B of active mol
def autoscale_b_factor():
  mol_id=_active_molecule_or_status()
  if mol_id is None:
    return None
  mean_b=average_temperature_factor(mol_id)
  if mean_b <= 0.0:
    info_dialog("Mean B-factor is not positive for the active molecule.")
    return None
  scale_fac=50/mean_b
  set_b_factor_bonds_scale_factor(mol_id,scale_fac)
  graphics_to_b_factor_representation(mol_id)
    
#Color molecule by rotamer and missing atom outliers
def color_rotamer_outliers_and_missing_atoms(mol_id):
  missing_atoms_list=[]
  missing_atoms_colour=2
  rotamer_outlier_list=[]
  rotamer_outlier_colour=34
  blank_colour=0
  blank_res_list=[]
  residue_specs=all_residues_sans_water(mol_id) or []
  for residue_spec in residue_specs:
    if not residue_spec or len(residue_spec) < 3:
      continue
    ch_id=residue_spec_to_chain_id(residue_spec)
    resn=residue_spec_to_res_no(residue_spec)
    ins_code=residue_spec_to_ins_code(residue_spec)
    if ch_id is False or resn is False or ins_code is False:
      continue
    blank_res_list.append(([ch_id,resn,ins_code],blank_colour))
    resname=residue_name(mol_id,ch_id,resn,ins_code)
    if resname in ["ALA","GLY","UNK","HOH"]:
      continue
    rot_prob=rotamer_score(mol_id,ch_id,resn,ins_code,"")
    if rot_prob<0.5 and rot_prob>0.0:
      rotamer_outlier_spec=[([ch_id,resn,ins_code],rotamer_outlier_colour)]
      rotamer_outlier_list=rotamer_outlier_list+rotamer_outlier_spec
  for x in missing_atom_info(mol_id) or []:
    if not x or len(x) < 3:
      continue
    missing_atoms_spec=[([x[0],x[1],x[2]],missing_atoms_colour)]
    missing_atoms_list=missing_atoms_list+missing_atoms_spec
  try:
    clear_user_defined_atom_colours(mol_id)
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,rotamer_outlier_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,missing_atoms_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass


#Colors subset of protein residues red, provided by user as string of single-letter ids.
def color_protein_residue_subset():
  def color_from_string(X):
    entry=str(X).upper() #capitalize
    mol_id=active_residue()[0]
    aa_dic={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL','X':'UNK'}
    highlight_colour=34 #red
    blank_colour=0 #light grey
    residue_list=[]
    resname_list=[]
    for char in entry:
      if (aa_dic.get(char,0)!=0): #make a list of 3-letter resnames from user supplied string by checking aa_dic
        resname=[(aa_dic.get(char,0))]
        resname_list=resname_list+resname
    for ch_id in chain_ids(mol_id):
      sn_max=chain_n_residues(ch_id,mol_id)
      for sn in  range(0,sn_max+1):
        resname_here=resname_from_serial_number(mol_id,ch_id,sn)
        for resname in resname_list:
          if resname==resname_here:
            resn=seqnum_from_serial_number(mol_id,ch_id,sn)
            ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
            residue_to_color=[([ch_id,resn,ins_id],highlight_colour)]
            residue_list=residue_list+residue_to_color
        if resname_here not in resname_list:
          resn=seqnum_from_serial_number(mol_id,ch_id,sn)
          ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
          residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
          residue_list=residue_list+residue_to_color
    try:
      set_user_defined_atom_colour_by_residue_py(mol_id,residue_list)
      graphics_to_user_defined_atom_colours_representation(mol_id)
    except NameError:
      info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
      pass
  generic_single_entry("Residues to color? Single letter code (e.g. DE or de will color Asp/Glu)","A","Color entered residue types!",color_from_string)

def color_polars_and_hphobs(mol_id):
  hphob_list=["CYS","ILE","LEU","VAL","TYR","MET","PHE","TRP","ALA"]
  polar_list=["SER","ASN","GLN","HIS","ARG","LYS","GLU","ASP","THR"]
  #based these on Moon&Fleming PNAS 2011 and MacCallum TIBS 2011
  #Gly/Pro colored differently because they are conformationally "special" residues
  polar_colour=5 #light blue
  hphob_colour=28 #orange
  gly_colour=34 #magenta
  pro_colour=15 #green
  blank_colour=0 #light gray
  polar_res_list=[]
  hphob_res_list=[]
  gly_res_list=[]
  pro_res_list=[]
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    for sn in range(0,sn_max+1):
      resname_here=resname_from_serial_number(mol_id,ch_id,sn)
      if resname_here in polar_list:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],polar_colour)]
        polar_res_list=polar_res_list+residue_to_color
      elif resname_here in hphob_list:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],hphob_colour)]
        hphob_res_list=hphob_res_list+residue_to_color
      elif resname_here=="GLY":
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],gly_colour)]
        gly_res_list=gly_res_list+residue_to_color
      elif resname_here=="PRO":
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],pro_colour)]
        pro_res_list=pro_res_list+residue_to_color
      else:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,polar_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,hphob_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,gly_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,pro_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass
  
def color_by_charge(mol_id):
  pos_list=["ARG","LYS","HIS"]
  neg_list=["GLU","ASP"]
  pos_colour=4
  neg_colour=31
  blank_colour=0
  pos_res_list=[]
  neg_res_list=[]
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    for sn in range(0,sn_max+1):
      resname_here=resname_from_serial_number(mol_id,ch_id,sn)
      if resname_here in pos_list:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],pos_colour)]
        pos_res_list=pos_res_list+residue_to_color
      elif resname_here in neg_list:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],neg_colour)]
        neg_res_list=neg_res_list+residue_to_color
      else:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,pos_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,neg_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass

def uncolor_other_chains():
  mol_id=active_residue()[0]
  ch_id_here=active_residue()[1]  
  blank_colour=0
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    if ch_id!=ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass

def color_active_chain():
  mol_id=active_residue()[0]
  ch_id_here=active_residue()[1]  
  blank_colour=0
  chain_colour=22 #yellow
  blank_res_list=[]
  chain_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    if ch_id!=ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
    if ch_id==ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],chain_colour)]
        chain_res_list=chain_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,chain_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass

def color_active_chain_by_num(chain_colour):
  mol_id=active_residue()[0]
  ch_id_here=active_residue()[1]
  blank_colour=0
  blank_res_list=[]
  chain_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    if ch_id!=ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
    if ch_id==ch_id_here:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],chain_colour)]
        chain_res_list=chain_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,chain_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass
    
def color_protein_na(mol_id):
  blank_colour=0
  protein_colour=22 #yellow
  na_colour=31
  protein_res_list=[]
  na_res_list=[]
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    if is_nucleotide_chain_p(mol_id,ch_id):
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],na_colour)]
        na_res_list=na_res_list+residue_to_color
    if is_protein_chain_p(mol_id,ch_id):
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],protein_colour)]
        protein_res_list=protein_res_list+residue_to_color
    else:
      for sn in range(0,sn_max+1):
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,protein_res_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,na_res_list)
    graphics_to_user_defined_atom_colours_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass


def color_waters(mol_id):
  water_colour=31
  blank_colour=0
  water_list=[]
  blank_res_list=[]
  for ch_id in chain_ids(mol_id):
    sn_max=chain_n_residues(ch_id,mol_id)
    for sn in range(0,sn_max+1):
      resname_here=resname_from_serial_number(mol_id,ch_id,sn)
      if resname_here=="HOH":
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],water_colour)]
        water_list=water_list+residue_to_color
      else:
        resn=seqnum_from_serial_number(mol_id,ch_id,sn)
        ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
        residue_to_color=[([ch_id,resn,ins_id],blank_colour)]
        blank_res_list=blank_res_list+residue_to_color
  try:
    set_user_defined_atom_colour_by_residue_py(mol_id,water_list)
    set_user_defined_atom_colour_by_residue_py(mol_id,blank_res_list)
    graphics_to_user_defined_atom_colours_all_atoms_representation(mol_id)
  except NameError:
    info_dialog("You need a newer Coot - custom coloring is only in r6174 and later, sorry.")
    pass

  
#Mutate active chain to entered sequence
default_seq="MAAAA"
def mutate_by_resnum():
  def enter_seq(seq):
    global default_seq
    seq=str(seq).upper()
    seq.replace(" ", "")
    seq_dic={}
    len_seq=len(seq)
    n=0
    nmax=len_seq+1
    aa_dic={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU',
    'Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET',
    'F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
    clean_seq=''
    while (n<len_seq): 
      if ((seq[n].isalpha() and (seq[n] in aa_dic))):
        clean_seq=clean_seq+seq[n]
      n=n+1
    seq=clean_seq
    default_seq=seq
    len_seq=len(seq)
    n=0
    while (n<len_seq):
      value=aa_dic[seq[n]]
      seq_dic[n+1]=value
      n=n+1
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    sn=0
    last_sn=chain_n_residues(ch_id,mol_id)-1
    turn_off_backup(mol_id)
    while (sn<=last_sn):
      res=resname_from_serial_number(mol_id,ch_id,sn)
      seqnum=seqnum_from_serial_number(mol_id,ch_id,sn)
      ins_id=str(insertion_code_from_serial_number(mol_id,ch_id,sn))
      if ((res!=seq_dic.get(seqnum)) and ((res in aa_dic.values()) or (res=="MSE"))):
        if seq_dic.get(seqnum):
          mutate(mol_id,ch_id,seqnum,ins_id,seq_dic.get(seqnum))
          if seq_dic.get(seqnum)!="PRO":
            delete_residue_sidechain(mol_id,ch_id,seqnum,ins_id,0)
      sn=sn+1
    turn_on_backup(mol_id)
  generic_single_entry("Enter raw amino acid sequence (must be complete!)",
  default_seq,"Mutate active chain to match sequence using PDB numbering", enter_seq)
  
#Search PDB by active chain
#Modify example below. Need to get (clean - waters etc removed) sequence of 
#active chain, search PDB using REST (look up RCSB REST for details), parse output line by line (try using first four characters of each line, or everything between numeric and :, 
#load PDBs (get_ebi_pdb pdb_id) and SSM superpose (superpose imol1 imol2 move_imol2_flag;
#superpose-with-chain-selection imol1 imol2 chain_imol1 chain_imol2 chain_used_flag_imol1 chain_used_flag_imol2 
#move_imol2_copy_flag) each on active chain. Would probably be useful to optionally only load the best matching chain
#of each struc; but equally might be good to keep everything, so that ligands, conserved waters, binding partners etc can be seen.
# should probably restrict to top 10, and exclude those that are 90%+identical... or maybe take min_idpct and max_idpct as parameters?
#
# def query_pdb_by_active_chain():
#   import urllib2
#   url = 'http://www.rcsb.org/pdb/rest/search'
#   mol_id=active_residue()[0]
#   ch_id=active_residue()[1]
#   seq=print_sequence_chain(mol_id,ch_id) # print_seq won't work, as does not return seq but pipes to stdout. need to write func.
#   print seq
#   seq.replace("X","")
#   print seq
#   
#   queryText = """
# <?xml version="1.0" encoding="UTF-8"?>
# <orgPdbQuery>
# <queryType>org.pdb.query.simple.SequenceQuery</queryType>
# <sequence>{sequence}</sequence>
# <eCutOff>0.00001</eCutOff>
# <searchTool>blast</searchTool>
# <sequenceIdentityCutoff>50</sequenceIdentityCutoff>
# </orgPdbQuery>
#   """.format(sequence=seq)
#   
#   print "query:\n", queryText
#   print "querying PDB...\n"
#   req = urllib2.Request(url, data=queryText)
#   f = urllib2.urlopen(req)
#   result = f.read()
#   if result:
#     print "Found number of PDB entries:", result.count('\n')
#     print result
#   else:
#     print "Failed to retrieve results" 
  
#Highlight various items with ball-and-stick and lines
def highlight_chain_breaks():
  mol_id=active_residue()[0]
  clear_ball_and_stick(mol_id)
  turn_off_backup(mol_id)
  obj_number=generic_object_with_name("chain_breaks_{mol_id}".format(mol_id=mol_id))
  generic_object_clear(obj_number)
  protein_resnames=['ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
  na_resnames=['A','C','T','G','U']
  missing_segments_list=[]
  for ch_id in chain_ids(mol_id):
    for resn in range(0,chain_n_residues(ch_id,mol_id)):
      resname=resname_from_serial_number(mol_id,ch_id,resn)
      print(resname)
      if (resname in protein_resnames):
        if ((is_term_type_mc_sn(mol_id,ch_id,resn)==1) or (is_term_type_mn_sn(mol_id,ch_id,resn)==1)):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
#         sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
#         make_ball_and_stick(mol_id,sel_string,0,0.5,1
          if (is_term_type_mn_sn(mol_id,ch_id,resn)==1):
            x_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-3]
            y_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-2]
            z_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-1]
            x_mid=(x_mn+x_here)/2
            y_mid=(y_mn+y_here)/2
            z_mid=(z_mn+z_here)/2
            res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
            #formula: sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
            distance=((x_here-x_mn)**2+(y_here-y_mn)**2+(z_here-z_mn)**2)**0.5
            print("distance",distance)
            distance_per_residue=distance/res_missing
            print("per residue distance",distance_per_residue)
            label_string="{ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing)
            if distance_per_residue>3.8:
              label_string="MIND THE GAP! {ch_id}  {res_start}...{res_end} ({res_missing} missing residues for {distance:6.1f} A)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing,distance=distance)
            if res_missing >=50:
              label_string="!!!  {ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing)
            list_entry=[label_string,x_mid,y_mid,z_mid]
            missing_segments_list.append(list_entry)
            if res_missing <=15:
              gap_color="gray"
            elif (res_missing > 15) and (res_missing<50):
              gap_color="orange"
            else:
              gap_color="red"
            if distance_per_residue>3.8:
              gap_color="cyan" 
            line_width=4
            dash_density=3
            to_generic_object_add_dashed_line(obj_number,gap_color,line_width,dash_density,x_here,y_here,z_here,x_mn,y_mn,z_mn)
            #if res_missing>=20:
            #  place_text(str(res_missing),x_mid,y_mid,z_mid,1)
          else: #By definition we always pass through here first (need to hit a term_type_mc befor term_type_mn)
            x_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-3]
            y_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-2]
            z_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-1]
        if (resn==0):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
          make_ball_and_stick(mol_id,sel_string,0,0.5,1)
        if (is_last_polymer_residue_sn(mol_id,ch_id,resn)==1):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
          make_ball_and_stick(mol_id,sel_string,0,0.5,1)
      if (resname in na_resnames):
        if ((is_term_type_mc_sn(mol_id,ch_id,resn)==1) or (is_term_type_mn_sn(mol_id,ch_id,resn)==1)):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
#         sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
#         make_ball_and_stick(mol_id,sel_string,0,0.5,1
          if (is_term_type_mn_sn(mol_id,ch_id,resn)==1):
            x_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-3]
            y_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-2]
            z_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-1]
            x_mid=(x_mn+x_here)/2
            y_mid=(y_mn+y_here)/2
            z_mid=(z_mn+z_here)/2
            res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
            #formula: sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
            distance=((x_here-x_mn)**2+(y_here-y_mn)**2+(z_here-z_mn)**2)**0.5
            distance_per_residue=distance/res_missing
            label_string="{ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing)
            if distance_per_residue>5.9:
              label_string="MIND THE GAP! {ch_id}  {res_start}...{res_end} ({res_missing} missing residues for {distance:6.1f} A)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing,distance=distance)
            if res_missing >=50:
              label_string="!!!  {ch_id}  {res_start}...{res_end} ({res_missing} missing residues)".format(ch_id=ch_id,res_start=seqnum_from_serial_number(mol_id,ch_id,resn-1),res_end=seqnum,res_missing=res_missing)
            list_entry=[label_string,x_mid,y_mid,z_mid]
            missing_segments_list.append(list_entry)
            if res_missing <=15:
              gap_color="gray"
            elif (res_missing > 15) and (res_missing<50):
              gap_color="orange"
            else:
              gap_color="red"
            if distance_per_residue>5.9:
              gap_color="cyan" 
            line_width=4
            dash_density=3
            to_generic_object_add_dashed_line(obj_number,gap_color,line_width,dash_density,x_here,y_here,z_here,x_mn,y_mn,z_mn)
            #if res_missing>=20:
            #  place_text(str(res_missing),x_mid,y_mid,z_mid,1)
          else: #By definition we always pass through here first (need to hit a term_type_mc befor term_type_mn)
            x_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-3]
            y_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-2]
            z_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-1]
        if (resn==0):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/P"
          make_ball_and_stick(mol_id,sel_string,0,0.5,1)
        if (is_last_polymer_residue_sn(mol_id,ch_id,resn)==1):
          seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
          sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/P"
          make_ball_and_stick(mol_id,sel_string,0,0.5,1)
  set_display_generic_object(obj_number,1)
  try:
    attach_generic_object_to_molecule(obj_number, mol_id)
  except NameError:
    info_dialog("attach_generic_object_to_molecule is only present in Coot r6057 and later, sorry.")
    pass
  turn_on_backup(mol_id)
  interesting_things_gui("Missing segments",missing_segments_list)


def highlight_all_chain_breaks():
  protein_resnames=['ALA','UNK','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL']
  na_resnames=['A','C','T','G','U']
  for mol_id in model_molecule_list():
    clear_ball_and_stick(mol_id)
    turn_off_backup(mol_id)
    obj_number=generic_object_with_name("chain_breaks_{mol_id}".format(mol_id=mol_id))
    generic_object_clear(obj_number)
    for ch_id in chain_ids(mol_id):
      for resn in range(0,chain_n_residues(ch_id,mol_id)):
        resname=resname_from_serial_number(mol_id,ch_id,resn)
        print(resname)
        if (resname in protein_resnames):
          if ((is_term_type_mc_sn(mol_id,ch_id,resn)==1) or (is_term_type_mn_sn(mol_id,ch_id,resn)==1)):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
  #         sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
  #         make_ball_and_stick(mol_id,sel_string,0,0.5,1)
            if (is_term_type_mn_sn(mol_id,ch_id,resn)==1):
              x_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-3]
              y_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-2]
              z_mn=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-1]
              x_mid=(x_mn+x_here)/2
              y_mid=(y_mn+y_here)/2
              z_mid=(z_mn+z_here)/2
              res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
              #formula: sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
              distance=((x_here-x_mn)**2+(y_here-y_mn)**2+(z_here-z_mn)**2)**0.5
              distance_per_residue=distance/res_missing
              if res_missing <=15:
                gap_color="gray"
              elif (res_missing > 15) and (res_missing<50):
                gap_color="orange"
              else:
                gap_color="red"
              if distance_per_residue>3.8:
                gap_color="cyan" 
              line_width=4
              dash_density=3
              to_generic_object_add_dashed_line(obj_number,gap_color,line_width,dash_density,x_here,y_here,z_here,x_mn,y_mn,z_mn)
              #res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
              #if res_missing>=20:
              #  place_text(str(res_missing),x_mid,y_mid,z_mid,1)
            else: #By definition we always pass through here first (need to hit a term_type_mc befor term_type_mn)
              x_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-3]
              y_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-2]
              z_here=atom_specs(mol_id,ch_id,seqnum,"","CA","")[-1]
          if (resn==0):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
            make_ball_and_stick(mol_id,sel_string,0,0.5,1)
          if (is_last_polymer_residue_sn(mol_id,ch_id,resn)==1):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/CA"
            make_ball_and_stick(mol_id,sel_string,0,0.5,1)
        if (resname in na_resnames):
          if ((is_term_type_mc_sn(mol_id,ch_id,resn)==1) or (is_term_type_mn_sn(mol_id,ch_id,resn)==1)):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            if (is_term_type_mn_sn(mol_id,ch_id,resn)==1):
              x_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-3]
              y_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-2]
              z_mn=atom_specs(mol_id,ch_id,seqnum,"","P","")[-1]
              x_mid=(x_mn+x_here)/2
              y_mid=(y_mn+y_here)/2
              z_mid=(z_mn+z_here)/2
              res_missing=seqnum-seqnum_from_serial_number(mol_id,ch_id,resn-1)
              #formula: sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
              distance=((x_here-x_mn)**2+(y_here-y_mn)**2+(z_here-z_mn)**2)**0.5
              distance_per_residue=distance/res_missing
              if res_missing <=15:
                gap_color="gray"
              elif (res_missing > 15) and (res_missing<50):
                gap_color="orange"
              else:
                gap_color="red"
              if distance_per_residue>5.9:
                gap_color="cyan" 
              line_width=4
              dash_density=3
              to_generic_object_add_dashed_line(obj_number,gap_color,line_width,dash_density,x_here,y_here,z_here,x_mn,y_mn,z_mn)
              #if res_missing>=20:
              #  place_text(str(res_missing),x_mid,y_mid,z_mid,1)
            else: #By definition we always pass through here first (need to hit a term_type_mc befor term_type_mn)
              x_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-3]
              y_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-2]
              z_here=atom_specs(mol_id,ch_id,seqnum,"","P","")[-1]
          if (resn==0):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/P"
            make_ball_and_stick(mol_id,sel_string,0,0.5,1)
          if (is_last_polymer_residue_sn(mol_id,ch_id,resn)==1):
            seqnum=seqnum_from_serial_number(mol_id,ch_id,resn)
            sel_string="//"+str(ch_id)+"/"+str(seqnum)+"/P"
            make_ball_and_stick(mol_id,sel_string,0,0.5,1)
    set_display_generic_object(obj_number,1) 
    try:
      attach_generic_object_to_molecule(obj_number, mol_id)
    except NameError:
      info_dialog("attach_generic_object_to_molecule is only present in Coot r6057 and later, sorry.")
      pass
    turn_on_backup(mol_id)

def local_hbonds_and_baddies():
  mol_id=active_residue()[0]
  probe_local_sphere(mol_id,5)
  undisplay_list=["wide contact","close contact"]
  display_list=["bad overlap","H-bonds","small overlap"]
  for obj_name in undisplay_list:
    obj_number=generic_object_with_name(obj_name)
    generic_object_clear(obj_number)
  for obj_name in display_list:
    obj_number=generic_object_with_name(obj_name)
    set_display_generic_object(obj_number,1)
    try:
      attach_generic_object_to_molecule(obj_number, mol_id)
    except NameError:
      info_dialog("attach_generic_object_to_molecule is only present in Coot r6057 and later, sorry.")
      pass
      
def global_hbonds_and_baddies():
  mol_id=active_residue()[0]
  probe(mol_id)
  close_molecule(model_molecule_list()[-1])
  undisplay_list=["wide contact","close contact","H-bonds","small overlap"]
  display_list=["bad overlap"]
  for obj_name in undisplay_list:
    obj_number=generic_object_with_name(obj_name)
    generic_object_clear(obj_number)
  for obj_name in display_list:
    obj_number=generic_object_with_name(obj_name)
    set_display_generic_object(obj_number,1)
    try:
      attach_generic_object_to_molecule(obj_number, mol_id)
    except NameError:
      info_dialog("attach_generic_object_to_molecule is only present in Coot r6057 and later, sorry.")
      pass
      
def clear_dots():
  undisplay_list=["wide contact","close contact", "small overlap","bad overlap","H-bonds","H-bond","big-overlap","close-contact","small-overlap","wide-contact","clashes"]
  for obj_name in undisplay_list:
    obj_number=generic_object_with_name(obj_name)
    generic_object_clear(obj_number)

#Make an alkyl chain of entered length and autofit to map if available
def make_alkyl_chain():
  def make_alkyl_chain_length_n(n):
    n_carbons=_positive_int_from_entry(n,"Number of carbons")
    if n_carbons is None:
      return None
    smiles_string=int(n_carbons+1)*"c"
    new_molecule_by_smiles_string("",smiles_string,force_libcheck=True)
    delete_hydrogens(molecule_number_list()[-1])
    mol_id=molecule_number_list()[-1]
    ch_id="A"
    res_no=1
    ins_code=""
    altloc=""
    new_residue_name="CXC"
    set_residue_name(mol_id,ch_id,res_no,ins_code,new_residue_name)
    delete_atom(mol_id,ch_id,res_no,ins_code," C  ",altloc)
    prodrg_ify(mol_id,ch_id,res_no,ins_code)
    close_molecule(molecule_number_list()[-1])
    close_molecule(molecule_number_list()[-1])
    if imol_refinement_map()!=-1:
      fit_chain_to_map_by_random_jiggle(mol_id,ch_id,1000,0.1)
  generic_single_entry("How many carbons do you want in the chain?",
    "10","Make alkyl chain",make_alkyl_chain_length_n)

#Place helix and add helix restraints
def place_helix_with_restraints():
  size_initial=len(model_molecule_number_list())
  if place_helix_here():
    size_after=len(model_molecule_number_list())
    size_diff=size_after-size_initial
    for i in range(1,size_diff+1):
      mol_id=model_molecule_number_list()[-i]
      set_b_factor_molecule(mol_id,default_new_atoms_b_factor())
      graphics_to_rainbow_representation(mol_id)
      ch_id=""
      res1=seqnum_from_serial_number(mol_id,ch_id,0)
      res2=last_polymer_residue(mol_id,ch_id)
      for resn in range(res1,res2-2):
        add_extra_bond_restraint(mol_id,
                                 ch_id, resn    , "", " O  ", "",
                                 ch_id, resn + 3, "", " N  ", "",
                                 3.3, 0.1)
        if (resn + 4 <= res2):
          add_extra_bond_restraint(mol_id,
                                   ch_id, resn    , "", " O  ", "",
                                   ch_id, resn + 4, "", " N  ", "",
                                   2.9, 0.05)

#Make new helix (don't fit)
def place_new_helix():
  def place_new_helix_entry(n):
    n_res=_positive_int_from_entry(n,"Number of helix residues")
    if n_res is None:
      return None
    get_monomer_no_H("ALA")
    mol_id=model_molecule_list()[-1]
    ch_id=chain_ids(mol_id)[0]
    res_no=1
    ins_code=""
    altloc=""
    delete_atom(mol_id,ch_id,1,ins_code," OXT",altloc)
    res_type="auto"
    for i in range(1,n_res):
      add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
      res_type,-57.82,-47)
      res_no=res_no+1
  generic_single_entry("How many residues for helix?",
  "10","Place helix",place_new_helix_entry)
  
#Make new strand (don't fit)
def place_new_strand():
  def place_new_strand_entry(n):
    n_res=_positive_int_from_entry(n,"Number of strand residues")
    if n_res is None:
      return None
    get_monomer_no_H("ALA")
    mol_id=model_molecule_list()[-1]
    ch_id=chain_ids(mol_id)[0]
    res_no=1
    ins_code=""
    altloc=""
    delete_atom(mol_id,ch_id,1,ins_code," OXT",altloc)
    res_type="auto"
    for i in range(1,n_res):
      add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
      res_type,-139,135)
      res_no=res_no+1
  generic_single_entry("How many residues for strand?",
  "10","Place strand",place_new_strand_entry)

#Make new 3-10 helix (don't fit)
def place_new_3_10_helix():
  def place_new_3_10_helix_entry(n):
    n_res=_positive_int_from_entry(n,"Number of 3-10 helix residues")
    if n_res is None:
      return None
    get_monomer_no_H("ALA")
    mol_id=model_molecule_list()[-1]
    ch_id=chain_ids(mol_id)[0]
    res_no=1
    ins_code=""
    altloc=""
    delete_atom(mol_id,ch_id,1,ins_code," OXT",altloc)
    res_type="auto"
    for i in range(1,n_res):
      add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,
      res_type,-49,-26)
      res_no=res_no+1
  generic_single_entry("How many residues for 3-10 helix?",
  "10","Place 3-10 helix",place_new_3_10_helix_entry)

#Merge two chains (throw error msg if they overlap in numbering)
#This needs fixing - says selections overlap in cases where sel1 has one segment on either side from sel2
#Fixed now, I think.
def merge_chains():
  def merge_chains_by_click(res1,res2):
    mol_id_1=res1[1]
    mol_id_2=res2[1]
    ch_id_1=res1[2]
    ch_id_2=res2[2]
    first_res_sel1=first_residue(mol_id_1,ch_id_1)
    last_res_sel1=last_residue(mol_id_1,ch_id_1)
    first_res_sel2=first_residue(mol_id_2,ch_id_2)
    last_res_sel2=last_residue(mol_id_2,ch_id_2)
    size_2=abs(last_res_sel2-first_res_sel2)
    size_1=abs(last_res_sel1-first_res_sel1)
    if (mol_id_1!=mol_id_2) or (ch_id_1==ch_id_2):
      info_dialog("No can do, chains must be in the same mol and have non-overlapping ranges!")
    elif (mol_id_1==mol_id_2) and (ch_id_1!=ch_id_2) and (size_1<=size_2):
      out=change_chain_id_with_result(mol_id_1,ch_id_1,ch_id_2,1,first_res_sel1,last_res_sel1)
      if out[0]==0:
        info_dialog(out[1])
    elif (mol_id_1==mol_id_2) and (ch_id_1!=ch_id_2) and (size_2<size_1):
      out=change_chain_id_with_result(mol_id_1,ch_id_2,ch_id_1,1,first_res_sel2,last_res_sel2)
      if out[0]==0:
        info_dialog(out[1])
    else:
      info_dialog("No can do, chains must be in the same mol and have non-overlapping ranges!")
  user_defined_click(2,merge_chains_by_click)

#Renumber active segment by active residue
def renumber_seg_by_active_res():
  def renum_seg(new_num):
    new_num=int(new_num)
    current_num=active_residue()[2]
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    seg_count=0
    segments=segment_list(mol_id)
    last_res=last_polymer_residue(mol_id,ch_id)
    print(last_res)
    first_res=first_residue(mol_id,ch_id)
    print(first_res)
    new_seg_list=[]
    for seg in segments:
      if ch_id==seg[1]:
        new_seg_list.append(seg) 
    print(new_seg_list)
    for seg in new_seg_list:
      seg_count=seg_count+1
      if (current_num>=seg[2]) and (current_num<=seg[3]):
        res_start=seg[2]
        res_end=seg[3]
        ch_id=seg[1]
        if res_end<last_res:
          seg_next=new_seg_list[seg_count]
        if res_start>first_res:
          seg_prev=new_seg_list[seg_count-2]
        offset=new_num-current_num
        if ((((res_start==first_res) or (res_start+offset)>seg_prev[3])) and (((res_end==last_res) or (res_end+offset)<seg_next[2]))):
          renumber_residue_range(mol_id,ch_id,res_start,res_end,int(offset))
        else:
          info_dialog("No can do, this would result in overlapping sequence numbering!")
        delete_all_extra_restraints(mol_id)
        set_show_extra_restraints(mol_id,0)
        set_show_extra_restraints(mol_id,1)
  generic_single_entry("New number for this residue?",
  str(active_residue()[2]),"Renumber",renum_seg)
  
#Fit all segments to map
def rigid_body_fit_segments():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else:
    mol_id=active_residue()[0]
    segments=segment_list(mol_id)
    for seg in segments:
      res_start=seg[2]
      res_end=seg[3]
      ch_id=seg[1]
      turn_off_backup(mol_id)
      set_refinement_immediate_replacement(1)
      rigid_body_refine_zone(res_start,res_end,ch_id,mol_id)
      accept_regularizement()
      set_refinement_immediate_replacement(0)
      turn_on_backup(mol_id)
      
#Fit current segment
def fit_this_segment():
  if (imol_refinement_map()==-1):
    info_dialog("You must set a refinement map!")
  else: 
    mol_id=active_residue()[0]
    segments=segment_list(mol_id)
    res_here=active_residue()[2]
    ch_id=active_residue()[1]
    for seg in segments:
      if (res_here>=seg[2]) and (res_here<=seg[3]) and (ch_id==seg[1]):
        res_start=seg[2]
        res_end=seg[3]
        ch_id=seg[1]
        turn_off_backup(mol_id)
        set_refinement_immediate_replacement(1)
        rigid_body_refine_zone(res_start,res_end,ch_id,mol_id)
        accept_regularizement()
        set_refinement_immediate_replacement(0)
        turn_on_backup(mol_id)

#Set default b-fac for new atoms to mean B for active mol
def set_new_atom_b_fac_to_mean():
  mol_id=active_residue()[0]
  mean_b=average_temperature_factor(mol_id)
  set_default_temperature_factor_for_new_atoms(mean_b)
  
  
#Shortcut for adding terminal residue (to bind to key)
def add_term_shortcut():
  if imol_refinement_map()!=-1:
    mol_id=active_residue()[0]
    ch_id=active_residue()[1]
    resn=active_residue()[2]
    atom_name=active_residue()[4]
    first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
    last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
    delta_first=abs(first_in_seg-resn)
    delta_last=abs(last_in_seg-resn)
    set_new_atom_b_fac_to_mean()
    if delta_first<=delta_last:
      set_go_to_atom_molecule(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
      add_terminal_residue(mol_id,ch_id,first_in_seg,"auto",1)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg-1,atom_name)
    else:
      set_go_to_atom_molecule(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
      add_terminal_residue(mol_id,ch_id,last_in_seg,"auto",1)
      sort_residues(mol_id)
      set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg+1,atom_name)
  else:
    info_dialog("You must set a refinement map!")
    
def add_term_shortcut_force():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resn=active_residue()[2]
  atom_name=active_residue()[4]
  set_go_to_atom_molecule(mol_id)
  first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
  delta_first=abs(first_in_seg-resn)
  delta_last=abs(last_in_seg-resn)
  set_new_atom_b_fac_to_mean()
  if delta_first<=delta_last:
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
    force_add_terminal_residue_noclick(mol_id,ch_id,first_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg-1,atom_name)
  else:
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
    force_add_terminal_residue_noclick(mol_id,ch_id,last_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg+1,atom_name)
    
def add_term_shortcut_force_strand():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  resn=active_residue()[2]
  atom_name=active_residue()[4]
  set_go_to_atom_molecule(mol_id)
  first_in_seg=first_residue_in_seg(mol_id,ch_id,resn)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,resn)
  delta_first=abs(first_in_seg-resn) 
  delta_last=abs(last_in_seg-resn)
  set_new_atom_b_fac_to_mean()
  if delta_first<=delta_last:
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg,atom_name)
    force_add_terminal_residue_noclick_strand(mol_id,ch_id,first_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,first_in_seg-1,atom_name)
  else:
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg,atom_name)
    force_add_terminal_residue_noclick_strand(mol_id,ch_id,last_in_seg)
    sort_residues(mol_id)
    set_go_to_atom_chain_residue_atom_name(ch_id,last_in_seg+1,atom_name)

def _growth_terminus_state(mol_id, ch_id, res_no):
  first_in_seg=first_residue_in_seg(mol_id,ch_id,res_no)
  last_in_seg=last_residue_in_seg(mol_id,ch_id,res_no)
  delta_first=abs(first_in_seg-res_no)
  delta_last=abs(last_in_seg-res_no)
  if delta_first<=delta_last:
    return [first_in_seg, -1]
  return [last_in_seg, 1]

def _grow_from_selected_terminus(mol_id, ch_id, clicked_res_no, n_res, phi, psi, b_factor_range_qm=0):
  [res_no, step_direction]=_growth_terminus_state(mol_id,ch_id,clicked_res_no)
  res_no_0=res_no

  if step_direction == 1:
    current_term_res_no=last_residue_in_seg(mol_id,ch_id,res_no)
  else:
    current_term_res_no=first_residue_in_seg(mol_id,ch_id,res_no)

  requested_start=current_term_res_no+step_direction
  requested_end=current_term_res_no+(step_direction*n_res)
  requested_min=min(requested_start, requested_end)
  requested_max=max(requested_start, requested_end)

  if step_direction == 1:
    for seg in segment_list_chain(mol_id,ch_id):
      seg_start=seg[2]
      seg_end=seg[3]
      if (seg_start > current_term_res_no) and not ((requested_max < seg_start) or (requested_min > seg_end)):
        info_dialog("Requested growth would overlap the next polymer segment.\n\nPlease renumber or move the later segment first.")
        return 0
  else:
    for seg in segment_list_chain(mol_id,ch_id):
      seg_start=seg[2]
      seg_end=seg[3]
      if (seg_end < current_term_res_no) and not ((requested_max < seg_start) or (requested_min > seg_end)):
        info_dialog("Requested growth would overlap the previous polymer segment.\n\nPlease renumber or move the earlier segment first.")
        return 0

  added_count=0
  for i in range(1,(n_res+1)):
    if step_direction == 1:
      res_no=last_residue_in_seg(mol_id,ch_id,res_no)
    else:
      res_no=first_residue_in_seg(mol_id,ch_id,res_no)
    target_res_no=res_no+step_direction
    if residue_exists_qm(mol_id,ch_id,target_res_no,""):
      if (step_direction == 1) and (not _residue_is_polymer(mol_id, ch_id, target_res_no, "")):
        remaining_growth=(n_res - added_count)
        renumber_offset=remaining_growth + 1
        renumber_residue_range(mol_id,ch_id,target_res_no,last_residue(mol_id,ch_id),renumber_offset)
        sort_residues(mol_id)
        if step_direction == 1:
          res_no=last_residue_in_seg(mol_id,ch_id,res_no)
        else:
          res_no=first_residue_in_seg(mol_id,ch_id,res_no)
        target_res_no=res_no+step_direction
      else:
        _status_message("Stopped growth before overlapping adjacent segment numbering")
        break
    res_type="auto"
    add_terminal_residue_using_phi_psi(mol_id,ch_id,res_no,res_type,phi,psi)
    sort_residues(mol_id)
    if step_direction == 1:
      new_res_no=last_residue_in_seg(mol_id,ch_id,res_no)
    else:
      new_res_no=first_residue_in_seg(mol_id,ch_id,res_no)
    if new_res_no == res_no:
      _status_message("Stopped growth because Coot could not add the next terminal residue")
      break
    res_no=new_res_no
    added_count=added_count+1
  if b_factor_range_qm and added_count>0:
    set_b_factor_residue_range(mol_id,ch_id,min(res_no_0,res_no),max(res_no_0,res_no),default_new_atoms_b_factor())
  return added_count

#Add h-bond restraints to active mol with Prosmart
def run_prosmart_self():
  """
  target is my molecule, ref is the homologous (high-res) model

  extra arg: include_side_chains=False
  """
  imol_target=active_residue()[0]
  dir_stub = "coot-ccp4"
  prosmart_dir="ProSMART_Output"
  make_directory_maybe(dir_stub)
  make_directory_maybe(prosmart_dir)
  target_pdb_file_name = os.path.join(dir_stub,
                                      molecule_name_stub(imol_target, 0).replace(" ", "_") + \
                                      "-prosmart.pdb")
  prosmart_out = os.path.join("ProSMART_Output",
                              molecule_name_stub(imol_target, 0).replace(" ", "_") + \
                              "-prosmart.txt")

  write_pdb_file(imol_target, target_pdb_file_name)
  prosmart_exe = find_exe("prosmart")
  if prosmart_exe:
    l = ["-p1", target_pdb_file_name,
         "-h","-self_restrain","-bond_max","4.0","-bond_override", "2","-o",prosmart_dir]
    popen_command(prosmart_exe,
                  l,
                  [],
                  os.path.join(dir_stub, "prosmart.log"),
                  False)
    if (not os.path.isfile(prosmart_out)):
      print "file not found", prosmart_out
    else:
      print "Reading ProSMART restraints from", prosmart_out
      add_refmac_extra_restraints(imol_target, prosmart_out)
      delete_extra_restraints_worse_than(imol_target,4)
      
#Fix all atoms in active residue
def fix_active_residue():
  ar=active_residue() #ar[0] is mol_id, 1 is ch_id, 2 is resnum, 3 is ins code, 4 is atom name (CA) and 5 is alt conf 
  residue_info_list=residue_info(ar[0],ar[1],ar[2],ar[3])
  atom_list=[]
  for item in residue_info_list: # residue item is list of list of lists. first item in first list of each list in the list is the atom name :-)
    atom_name=item[0][0]
    atom_info=[ar[1],ar[2],ar[3],atom_name,ar[5]] #Make a list of chain id, res number, ins code, atom name and alt conf. Not sure if alt conf and ins code are right way around here...
    atom_list.append(atom_info)
  mark_multiple_atoms_as_fixed(ar[0],atom_list,1)
  
#Save and overwrite active model:
def quicksave_active():
  import shutil
  mol_id=active_residue()[0]
  filename=molecule_name(mol_id)
  filename_bak=filename+".bak"
  shutil.copy(filename,filename_bak)
  save_coordinates(mol_id,filename)
  info_dialog("Saved {filename} (molecule {mol_id}). \n If this was an accident, you can find a backup of the original file here: \n {filename_bak}".format(filename=filename,mol_id=mol_id,filename_bak=filename_bak))

#Refresh model from file (if changed externally for example)
def reload_model():
  mol_id=active_residue()[0]
  filename=molecule_name(mol_id)
  clear_and_update_model_molecule_from_file(mol_id,filename)

#Copy changes from active chain to NCS equivalents.
def copy_ncs_chain_from_active():
  mol_id=active_residue()[0]
  ch_id=active_residue()[1]
  turn_off_backup(mol_id)
  ncs_control_change_ncs_master_to_chain_id(mol_id,ch_id)
  copy_from_ncs_master_to_others(mol_id,ch_id)
  turn_on_backup(mol_id)

def return_seq_as_string(mol_id,ch_id):
  seq=''
  for sn in range(0,chain_n_residues(ch_id,mol_id)):
    resno=seqnum_from_serial_number(mol_id,ch_id,sn)
    ins_code=insertion_code_from_serial_number(mol_id,ch_id,sn)
    residue_name_3_letter=residue_name(mol_id,ch_id,resno,ins_code)
    aa_code=three_letter_code2single_letter(residue_name_3_letter)
    if (len(residue_name(mol_id,ch_id,resno,ins_code))==1):
      aa_code=residue_name_3_letter
    if (residue_name(mol_id,ch_id,resno,ins_code)!="ALA") and (aa_code=="A") and (len(residue_name(mol_id,ch_id,resno,ins_code))!=1):
      aa_code="X"
    seq=seq+aa_code
  return seq

def _compile_sequence_pattern(subseq):
  pattern=[]
  index=0
  while index < len(subseq):
    character=subseq[index]
    if character=="(":
      close_index=subseq.find(")", index+1)
      if close_index==-1:
        raise ValueError("Unmatched '(' in sequence pattern")
      option_text=subseq[index+1:close_index]
      if not option_text:
        raise ValueError("Empty bracketed option in sequence pattern")
      options=[piece.strip().upper() for piece in option_text.split("/") if piece.strip()]
      if not options:
        raise ValueError("Empty bracketed option in sequence pattern")
      allowed=set()
      for option in options:
        if len(option)!=1:
          raise ValueError("Bracketed options must be single-letter residues")
        allowed.add(option)
      pattern.append(("set", allowed))
      index=close_index+1
      continue
    if character=="X":
      pattern.append(("any", None))
    else:
      pattern.append(("set", set([character])))
    index=index+1
  if not pattern:
    raise ValueError("Empty sequence pattern")
  return pattern

def _sequence_pattern_matches_at(seq, start_index, compiled_pattern):
  if start_index+len(compiled_pattern) > len(seq):
    return False
  for offset in range(len(compiled_pattern)):
    mode, allowed=compiled_pattern[offset]
    seq_char=seq[start_index+offset]
    if mode=="any":
      continue
    if seq_char not in allowed:
      return False
  return True

def _sequence_match_context_label(ch_id, seq, sn_start, pattern_length, mol_id):
  sn_end=sn_start+pattern_length-1
  resno_start=seqnum_from_serial_number(mol_id,ch_id,sn_start)
  resno_end=seqnum_from_serial_number(mol_id,ch_id,sn_end)
  left_context_start=max(0, sn_start-3)
  right_context_end=min(len(seq), sn_start+pattern_length+3)
  left_context=seq[left_context_start:sn_start]
  match_context=seq[sn_start:sn_start+pattern_length]
  right_context=seq[sn_start+pattern_length:right_context_end]
  if left_context_start>0:
    left_context="..."+left_context
  if right_context_end<len(seq):
    right_context=right_context+"..."
  return "%s%d-%d: %s*%s*%s" %(ch_id,resno_start,resno_end,left_context,match_context,right_context)

def find_sequence_in_current_chain(subseq):
  subseq=subseq.upper()
  residue=_active_residue_or_status()
  if not residue:
    return None
  mol_id=residue[0]
  ch_id=residue[1]
  seq=return_seq_as_string(mol_id,ch_id)
  try:
    compiled_pattern=_compile_sequence_pattern(subseq)
  except ValueError, e:
    info_dialog(str(e))
    return None
  sn_list=[]
  interesting_list=[]
  pattern_length=len(compiled_pattern)
  index=0
  while index <= len(seq)-pattern_length:
    if _sequence_pattern_matches_at(seq, index, compiled_pattern):
      sn_list.append(index)
      index=index+pattern_length
    else:
      index=index+1
  if len(sn_list)==1:
    sn_start=sn_list[0]
    resno=seqnum_from_serial_number(mol_id,ch_id,sn_start)
    ins_code=insertion_code_from_serial_number(mol_id,ch_id,sn_start)
    alt_conf=residue_alt_confs(mol_id,ch_id,resno,ins_code)[0]
    print("sn_start",sn_start)
    try:
      x=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-3]
      y=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-2]
      z=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-1]
    except TypeError:
      x=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-3]
      y=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-2]
      z=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-1]
    set_rotation_centre(x,y,z)
  elif len(sn_list)==0:
    info_dialog("Sequence not found!")
  elif len(sn_list)>1:
    count=0
    for sn in sn_list:
      count=count+1
      resno=seqnum_from_serial_number(mol_id,ch_id,sn)
      ins_code=insertion_code_from_serial_number(mol_id,ch_id,sn)
      alt_conf=residue_alt_confs(mol_id,ch_id,resno,ins_code)[0]
      print("sn_start",sn)
      try:
        x=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-3]
        y=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-2]
        z=atom_specs(mol_id,ch_id,resno,ins_code,"CA",alt_conf)[-1]
      except TypeError:
        x=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-3]
        y=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-2]
        z=atom_specs(mol_id,ch_id,resno,ins_code,"P",alt_conf)[-1]
      list_entry=[_sequence_match_context_label(ch_id,seq,sn,pattern_length,mol_id),x,y,z]
      interesting_list.append(list_entry)
      print("interesting list",interesting_list)
    print("interesting list",interesting_list)
    interesting_things_gui("Matches to entered sequence",interesting_list)

def find_sequence_with_entry():
  generic_single_entry("Enter sequence fragment to find\n(X=wildcard, (A/S)=either residue)",
  "MAAAA","Find sequence in active chain",find_sequence_in_current_chain)

      
def user_defined_add_arbitrary_length_bond_restraint(bond_length=2.0):
  def make_restr(text_list, continue_qm):
    s = "Now click on 2 atoms to define the additional bond restraint"
    add_status_bar_text(s)
    dist = text_list[0]
    try:
      bl = float(dist)
    except:
      bl = False
      add_status_bar_text("Must define a number for the bond length")
    if bl:
      def make_restr_dist(*args):
        atom_spec_1 = args[0]
        atom_spec_2 = args[1]
        imol = atom_spec_1[1]
        print "BL DEBUG:: imol: %s spec 1: %s and 2: %s" %(imol, atom_spec_1, atom_spec_2)
        add_extra_bond_restraint(imol, atom_spec_1[2], atom_spec_1[3], atom_spec_1[4], atom_spec_1[5], atom_spec_1[6], atom_spec_2[2], atom_spec_2[3], atom_spec_2[4], atom_spec_2[5], atom_spec_2[6], bl, 0.035)
      user_defined_click(2, make_restr_dist)
      if continue_qm:
        user_defined_add_arbitrary_length_bond_restraint(bl)
  def stay_open(*args):
    pass
  #generic_single_entry("Add a User-defined extra distance restraint",
  #                     "2.0",
  #                     "OK...",
  #                     lambda text: make_restr(text))
  generic_multiple_entries_with_check_button(
    [["Add a User-defined extra distance restraint",
      str(bond_length)]],
    ["Stay open?", lambda active_state: stay_open(active_state)],
    "OK...",
    lambda text, stay_open_qm: make_restr(text, stay_open_qm))
      

      
      
            

      
#Test post manipulation background func
# import time, threading
# def background_func():
#   threading.Timer(0.01, post_manip_background).start()
#   if condition:
#     do something
# post_manip_background()

# def post_manipulation_script(imol, mode):
#   print "BL DEBUG:: imol and mode", imol, mode
#   if (mode == DELETED):
#     print "BL DEBUG:: deleted something in mol ", imol
#     return 1
# 
# def post_manipulation_script2(imol, mode):
#   print "BL DEBUG:: imol and mode", imol, mode
#   if (mode == MUTATED):
#     print "BL DEBUG:: moved something in mol ", imol
#     return 1

#****Make and arrange menus****

menu=coot_menubar_menu("Custom")

add_simple_coot_menu_menuitem(menu,
"Custom keybindings...", lambda func: show_custom_keybindings_summary())

#make submenus

submenu_display=gtk.Menu()
menuitem_2=gtk.MenuItem("Display...")
menuitem_2.set_submenu(submenu_display)
menu.append(menuitem_2)
menuitem_2.show()

submenu_colour=gtk.Menu()
menuitem_2b=gtk.MenuItem("Colour...")
menuitem_2b.set_submenu(submenu_colour)
menu.append(menuitem_2b)
menuitem_2b.show()


submenu_fit=gtk.Menu()
menuitem_3=gtk.MenuItem("Fit...")
menuitem_3.set_submenu(submenu_fit)
menu.append(menuitem_3)
menuitem_3.show()


submenu_renumber=gtk.Menu()
menuitem_4=gtk.MenuItem("Renumber...")
menuitem_4.set_submenu(submenu_renumber)
menu.append(menuitem_4)
menuitem_4.show()

submenu_settings=gtk.Menu()
menuitem_5=gtk.MenuItem("Settings...")
menuitem_5.set_submenu(submenu_settings)
menu.append(menuitem_5)
menuitem_5.show()

submenu_build=gtk.Menu()
menuitem_6=gtk.MenuItem("Build...")
menuitem_6.set_submenu(submenu_build)
menu.append(menuitem_6)
menuitem_6.show()

submenu_mutate=gtk.Menu()
menuitem_7=gtk.MenuItem("Mutate...")
menuitem_7.set_submenu(submenu_mutate)
menu.append(menuitem_7)
menuitem_7.show()

submenu_copy=gtk.Menu()
menuitem_8=gtk.MenuItem("Copy...")
menuitem_8.set_submenu(submenu_copy)
menu.append(menuitem_8)
menuitem_8.show()

submenu_delete=gtk.Menu()
menuitem_9=gtk.MenuItem("Delete...")
menuitem_9.set_submenu(submenu_delete)
menu.append(menuitem_9)
menuitem_9.show()

submenu_merge=gtk.Menu()
menuitem_10=gtk.MenuItem("Merge...")
menuitem_10.set_submenu(submenu_merge)
menu.append(menuitem_10)
menuitem_10.show()

submenu_maps=gtk.Menu()
menuitem_11=gtk.MenuItem("Maps...")
menuitem_11.set_submenu(submenu_maps)
menu.append(menuitem_11)
menuitem_11.show()

#**** Populate submenus ****
#"Display..."

add_simple_coot_menu_menuitem(submenu_display, "All Molecules use \"C-alpha\" Symmetry", lambda func: map(lambda imol: valid_model_molecule_qm(imol) and symmetry_as_calphas(imol, 1), molecule_number_list()))

add_simple_coot_menu_menuitem(submenu_display, "Toggle Symmetry", 
lambda func: set_show_symmetry_master(not get_show_symmetry()))

add_simple_coot_menu_menuitem(submenu_display, "Load molecular symmetry copies from file metadata",
lambda func: _load_display_molecular_symmetry_from_metadata())

add_simple_coot_menu_menuitem(submenu_display, "Clear labels and distances", 
lambda func: clear_distances_and_labels())

add_simple_coot_menu_menuitem(submenu_display,
"Switch all mols to CA representation",lambda func: all_mols_to_ca())

add_simple_coot_menu_menuitem(submenu_colour,
"Color active mol by rotamer prob (outliers magenta) and missing atoms (blue)", lambda func: color_rotamer_outliers_and_missing_atoms_for_active_molecule())

add_simple_coot_menu_menuitem(submenu_colour,
"Color active mol by hydrophobics (orange), polars (blue), glys (magenta) and pros (green)", lambda func: color_polars_and_hphobs_for_active_molecule())

add_simple_coot_menu_menuitem(submenu_colour,
"Color active mol by charge (+ve blue, -ve red)", lambda func: color_by_charge_for_active_molecule())

add_simple_coot_menu_menuitem(submenu_colour,
"Uncolor other chains in active mol", lambda func: uncolor_other_chains())

add_simple_coot_menu_menuitem(submenu_colour,
"Color active chain", lambda func: color_active_chain())

add_simple_coot_menu_menuitem(submenu_colour,
"Color active segment", lambda func: colour_active_segment())

add_simple_coot_menu_menuitem(submenu_colour,
"Color by protein/nucleic acid", lambda func: color_protein_na_for_active_molecule())


add_simple_coot_menu_menuitem(submenu_colour,
"Color waters", lambda func: color_waters_for_active_molecule())

add_simple_coot_menu_menuitem(submenu_colour, "Colour entered subset of protein residues for active mol", lambda func: color_protein_residue_subset())

add_simple_coot_menu_menuitem(submenu_colour,
"Color active mol by Ramachandran outliers", lambda func: color_by_rama_native_for_active_residue())

add_simple_coot_menu_menuitem(submenu_colour,
"Color active mol by density fit", lambda func: color_by_density_fit_native_for_active_residue())

add_simple_coot_menu_menuitem(submenu_colour,
"Color active mol by NCS difference", lambda func: color_by_ncs_difference_for_active_residue())

add_simple_coot_menu_menuitem(submenu_colour,
"Color active mol by clash score", lambda func: color_by_clash_score_for_active_molecule())

add_simple_coot_menu_menuitem(submenu_colour, "Highlight chain breaks in active mol", lambda func: highlight_chain_breaks())

add_simple_coot_menu_menuitem(submenu_colour, "Highlight chain breaks in all mols", lambda func: highlight_all_chain_breaks())

add_simple_coot_menu_menuitem(submenu_display, "Show local probe dots (H-bonds, VdW and baddies only)", lambda func: local_hbonds_and_baddies())

add_simple_coot_menu_menuitem(submenu_display, "Show global probe dots (bad overlaps only)", lambda func: global_hbonds_and_baddies())

add_simple_coot_menu_menuitem(submenu_display, "Clear probe dots", lambda func: clear_dots())

add_simple_coot_menu_menuitem(submenu_display, "Find sequence in active chain", lambda func: find_sequence_with_entry())

#"Fit..."
add_simple_coot_menu_menuitem(submenu_fit, "Fit all chains to map", 
lambda func: rigid_fit_all_chains())

add_simple_coot_menu_menuitem(submenu_fit, "Stepped sphere refine active chain",
lambda func: stepped_sphere_refine_active_chain())

add_simple_coot_menu_menuitem(submenu_fit, "Fit current chain to map", 
lambda func: rigid_fit_active_chain())

add_simple_coot_menu_menuitem(submenu_fit, "Fit all segments", lambda func: rigid_body_fit_segments())

add_simple_coot_menu_menuitem(submenu_fit, "Fit this segment", lambda func: fit_this_segment())

add_simple_coot_menu_menuitem(submenu_fit,
"Smart self restrain active mol...",lambda func: prompt_generate_smart_local_extra_restraints())

add_simple_coot_menu_menuitem(submenu_fit,
"Cylinder refine (click start and end of range)",lambda func: refine_residues_sphere_click())

add_simple_coot_menu_menuitem(submenu_fit, "Add distance restraint (click two atoms)",lambda func:  user_defined_add_arbitrary_length_bond_restraint(bond_length=2.0))

#"Renumber..."

add_simple_coot_menu_menuitem(submenu_renumber, "Renumber active chain by first res", 
lambda func: renumber_by_first_res())

add_simple_coot_menu_menuitem(submenu_renumber, 
"Renumber active chain by last res", lambda func: renumber_by_last_res())

add_simple_coot_menu_menuitem(submenu_renumber,
"Renumber active chain by current res", lambda func: renumber_by_active_res())

add_simple_coot_menu_menuitem(submenu_renumber,"Renumber from N-term to active residue",
lambda func: renumber_n_term_segment())

add_simple_coot_menu_menuitem(submenu_renumber,"Renumber from active residue to C-term",
lambda func: renumber_c_term_segment())

add_simple_coot_menu_menuitem(submenu_renumber, "Renumber segment by active res", lambda func: renumber_seg_by_active_res())


#"Settings..."

add_simple_coot_menu_menuitem(submenu_settings,
"Auto-scale B-factor coloring for active mol",lambda func: autoscale_b_factor())

add_simple_coot_menu_menuitem(submenu_settings,
"Set Bfac for new atoms to mean B for active mol",lambda func: set_new_atom_b_fac_to_mean())

add_simple_coot_menu_menuitem(submenu_settings,
"Set proportional editing radius",lambda func: set_proportional_editing_radius())


#"Build..."
add_simple_coot_menu_menuitem(submenu_build,
"Forced addition of terminal residue (click terminus)", 
lambda func: force_add_terminal_residue())

add_simple_coot_menu_menuitem(submenu_build,
"Grow helix (click terminus)",
lambda func: grow_helix())

add_simple_coot_menu_menuitem(submenu_build,
"Grow strand (click terminus)",
lambda func: grow_strand())

add_simple_coot_menu_menuitem(submenu_build,
"Grow parallel strand (click terminus)",
lambda func: grow_parallel_strand())

add_simple_coot_menu_menuitem(submenu_build,
"Grow 3-10 helix (click terminus)",
lambda func: grow_helix_3_10())

add_simple_coot_menu_menuitem(submenu_build, 
"Shorten loop by one residue", lambda func: shorten_loop())

add_simple_coot_menu_menuitem(submenu_build, 
"Lengthen loop by one residue", lambda func: lengthen_loop())

add_simple_coot_menu_menuitem(submenu_build,
"Get fractional coordinates of active atom",lambda func: get_fract_coords()) 

add_simple_coot_menu_menuitem(submenu_build,
"Common monomers",lambda func: pick_common_monomers())

add_simple_coot_menu_menuitem(submenu_build, "Make alkyl chain of length n", lambda func: make_alkyl_chain())

add_simple_coot_menu_menuitem(submenu_build, "Make alpha helix of length n", lambda func: place_new_helix()) 

add_simple_coot_menu_menuitem(submenu_build, "Make 3-10 helix of length n", lambda func: place_new_3_10_helix())

#"Mutate...

add_simple_coot_menu_menuitem(submenu_mutate,
"Mutate range to UNK (click start and end)", lambda func: mutate_residue_range_by_click_a())

add_simple_coot_menu_menuitem(submenu_mutate,
"Mutate range to ALA (click start and end)", lambda func: mutate_residue_range_by_click_ala_a())

add_simple_coot_menu_menuitem(submenu_mutate, "Mutate all Mets to MSE", lambda func: mutate_all_mets_to_mse())

add_simple_coot_menu_menuitem(submenu_mutate, "Mutate all MSEs to Met", lambda func: mutate_all_mse_to_met())

add_simple_coot_menu_menuitem(submenu_mutate,
"Mutate active chain to template sequence (numbering must match sequence!)", lambda func: mutate_by_resnum())

#"Copy..."
add_simple_coot_menu_menuitem(submenu_copy, "Copy current chain", 
lambda func: copy_active_chain())

add_simple_coot_menu_menuitem(submenu_copy, "Cut current chain", 
lambda func: cut_active_chain())

add_simple_coot_menu_menuitem(submenu_copy, "Copy active segment", lambda func: copy_active_segment())

add_simple_coot_menu_menuitem(submenu_copy, "Cut active segment", lambda func: cut_active_segment())

add_simple_coot_menu_menuitem(submenu_copy,
"Copy active ligand/ion/solvent", lambda func: smart_copy_active_non_polymer_residue())

add_simple_coot_menu_menuitem(submenu_copy,
"Paste copied ligand/ion/solvent", lambda func: smart_paste_copied_non_polymer_residue())

add_simple_coot_menu_menuitem(submenu_copy,
"Copy fragment (click start and end)", lambda func: copy_frag_by_click())

add_simple_coot_menu_menuitem(submenu_copy,
"Cut fragment (click start and end)", lambda func: cut_frag_by_click())

add_simple_coot_menu_menuitem(submenu_copy,
"Copy active chain to NCS equivs", lambda func: copy_ncs_chain_from_active())


#"Delete..."
add_simple_coot_menu_menuitem(submenu_delete,
"Delete active chain", lambda func: delete_chain())

add_simple_coot_menu_menuitem(submenu_delete, "Delete active segment", lambda func: delete_active_segment())


add_simple_coot_menu_menuitem(submenu_delete, 
"Delete hydrogens from molecule", lambda func: delete_h_active())

add_simple_coot_menu_menuitem(submenu_delete,
"Delete sidechains in range (click start and end)", lambda func: delete_sidechain_range_by_click_a())  

#"Merge..."

add_simple_coot_menu_menuitem(submenu_merge, 
"Merge two mols (click two; 2nd into 1st)", lambda func: merge_fragments())

add_simple_coot_menu_menuitem(submenu_merge, "Merge chains (click two; 2nd into 1st)", lambda func: merge_chains())


#"Maps..."
add_simple_coot_menu_menuitem(submenu_maps,
"Sharpen (enter B-factor)",lambda func: sharpen_by_entered_factor())

add_simple_coot_menu_menuitem(submenu_maps,
"Change hi-res limit for map",lambda func: change_hires_limit_copy())

add_simple_coot_menu_menuitem(submenu_maps,
"Go to center of scrollable map",lambda func: goto_center_of_map())

add_simple_coot_menu_menuitem(submenu_maps,
"Set refinement map to scrollable map",lambda func: set_map_to_scrollable_map())

add_simple_coot_menu_menuitem(submenu_maps,
"Resample active EM map to 0.5 A/pixel",lambda func: resample_active_map_for_em_half_angstrom())
