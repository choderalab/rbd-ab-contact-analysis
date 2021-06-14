### Author: Ivy Zhang

## Preparing the antibody complex using ISOLDE

Here, we started with a rebuilt version of 7jx3 (contains the RBD, plus the S309, S304, and S2H14 antibodies) that was sent by Tristan Croll. This rebuilt version has a few rotamer adjustments, a handful of peptide flips, and a few weakly resolved waters. It also contains the glycan partially built in (which we removed and rebuilt ourselves, as described below). This was all done in [ISOLDE](https://isolde.cimr.cam.ac.uk/).

### Overview of refinement in ISOLDE

After visual inspection, we discovered a missing loop in the constant domain of S309. In order to build in the missing loop, we used [ISOLDE](https://isolde.cimr.cam.ac.uk/) to copy and paste the residues from another human IgG antibody [4jhw](https://www.rcsb.org/structure/4jhw), which has the same sequence for the chain with the missing loop.

Next, we used  [ISOLDE](https://isolde.cimr.cam.ac.uk/) to build in the glycan, which involved:
1) aligning the FA2G2 (open) structure (`systems/representative_glycan_structures/fa2g2_n_open.pdb`) to the glycan stub present in 7jx3, 
2) adding missing bonds in the sugar and between the sugar and N343, 
3) deleting the glycan stub in 7jx3.

Finally, there were severe clashes introduced when building in the glycan, so we used [ISOLDE](https://isolde.cimr.cam.ac.uk/) to remove these clashes and refine non-ideal rotamers.

### Detailed steps taken in ISOLDE

Below describes the steps taken to refine the full 7jx3 structure (unless otherwise specified, assume the commands are to be run in the ChimeraX command line):

#### Part 1: Build in missing loop in constant domain of S309
1. Load in Tristan's rebuild 7jx3 structure:
`open 01_7jx3_rebuild.pdb`
- Also, load in the crystal density (.mtz file from RCSB)
2. Load in 4jhw:
`open 4jhw`
3. Add hydrogens
`addh`
4. Align 4jhw constant domain (chains H and L) to 7jx3 constant domain (chains A and B)
`match #2 to #1/A:130-150`
5. Delete the two residues on either end (4 residues in total) of the missing loop in 7jx3
```
select #1/A:141-142 #1/A:147-148`
del sel
```
6. Copy 4jhw loop to 7jx3
In ChimeraX command line:
```
sel #2/H:127-134
isolde start
```
In shell (python session):
```python
from chimerax.atomic import selected_atoms
from chimerax.atomic import selected_residues
frag = selected_residues(session)
from chimerax.isolde.atomic.building.merge import merge_fragment
```
In ChimeraX command line:
`select #1/A:140`
In shell (python session):
`anchor_n = selected_residues(session)[0]`
In ChimeraX command line:
`select #1/A:149`
In shell(python session):
```python
anchor_c = selected_residues(session)[0]
m = session.isolde.selected_model
merge_fragment(m, frag, chain_id='A', renumber_from=141, anchor_n=anchor_n, anchor_c=anchor_c, transform=frag.unique_structures[0].position)
```
7. Hide 4jhw
`hide #2`
8. Fix cis bond on GLY148 (chain A) introduced by copy/paste:
`select #1/A:140-149`
Start simulation (click the blue arrow)
Select the cis bond (CTRL + click)
Flip the cis bond using the GUI
9. Remove psuedobond (dashed line in the visualization)
Select the atoms at the two ends of the pseudobond (CTRL + click, then CTRL + SHIFT + click)
```python
pbm = m.pseudobond_group(m.PBG_MISSING_STRUCTURE, create_type=None)
pbm.pseudobonds[pbm.pseudobonds.between_atoms(selected_atoms(session))].delete()
```
10. Identify problems with the structure by checking the different lists in the `Validate` tab. Fix problems and then click `Update list`
11. Save the session and model as a PDB.

#### Part 2: Build in glycan at N343 in the RBD
1. Load the fa2g2 open glycan structure into a new session
`open fa2g2_n_open.pdb` (located in `systems/representative_glycan_structures/`)
2. Rename the glycan residues to their canonical PDB names.
```python
def glycam_to_pdb_residue_names(model):

    from chimerax.atomic import Residue

    non_protein = model.residues[model.residues.polymer_types!=Residue.PT_AMINO]

    from chimerax.isolde.openmm.amberff import glycam

    def is_glycam(rname):

        return (rname[1:] in glycam.glycam_suffix_to_ccd_name.keys() and
            rname[0] in glycam._glycam_prefix.values())

    glycam_name_map = {r: r.name for r in non_protein if is_glycam(r.name)}
    
    import numpy

    anchor_names = list(glycam._anchor_name_map.values())

    
    anchor_names = anchor_names + ['N'+aname for aname in anchor_names] + ['C'+aname for aname in anchor_names]
    
    anchors = m.residues[numpy.in1d(m.residues.names, anchor_names)]
    
    anchor_name_map = {r: r.name for r in anchors}
    
    for r, gname in glycam_name_map.items():
        r.name = glycam.glycam_suffix_to_ccd_name[r.name[1:]]
    glycam_to_pdb_anchor = {val: key for key,val in glycam._anchor_name_map.items()}
    
    for r, gname in anchor_name_map.items():
        r.name = glycam_to_pdb_anchor[r.name[-3:]]
    glycam_name_map.update(anchor_name_map)
    
    return glycam_name_map

m = session.isolde.selected_model
name_map = glycam_to_pdb_residue_names(m)

```
3. Load the structure saved from the Part 1 (and crystal density data)
4. Align the glycan to the first NAG in 7jx3
You'll need to select the atoms from each structure you want to align on:
`align #3/X:2@C1 #3/X:2@O5#3/X:2@C2 #3/X:2@C4 #3/X:2@C5 toAtoms #1.2/R:601@C1 #1.2/R:601@O5 #1.2/R:601@C2 #1.2/R:601@C4 #1.2/R:601@C5`
5. Add missing bonds in the glycan
`sel /X:2@O6 /X:11@C`
Then go to `ISOLDE > model building` and click `Make bond`
`sel /X:4@O3 /X:8@C1`
Then go to `ISOLDE > model building` and click `Make bond`
6. Delete the ROH from the fa2g2/_n/_open structure
```
sel /X:1@O1
del sel
```
7. Combine the models
`select #1 #2`
Go to `ISOLDE > model building` and click `Merge models`
8. Select N343 ND2 (CTRL + CLICK) and NAG C1 (CTRL + SHIFT + CLICK) manually and then go to `ISOLDE > model building` and click `Make bond`
9. Add hydrogens
`addh`
10. Save session and the model as a PDB
11. Select the parts of the structure that have clashes and start the simulation (press the blue arrow) to fix clashes. Use the `Validate` tab to fix problems with rotamers, ramachandran map, etc.
12. Delete hydrogens
`delete H`
13. Rename glycan (and protein) residue names back to match the amber forcefield templates.
```python
from chimerax.isolde.openmm.openmm_interface import find_residue_templates
m = session.isolde.selected_model
ff = session.isolde.forcefield_mgr['amber14']
template_dict = find_residue_templates(m.residues, ff)
for i, template_name in template_dict.items():
    if 'GLYCAM' in template_name:
        m.residues[i].name = template_name.split('_')[1]
    elif '_' not in template_name:
        m.residues[i].name = template_name
```
14. Save the session and the model as a PDB.
