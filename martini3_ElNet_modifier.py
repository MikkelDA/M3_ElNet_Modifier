# %reset -f

'''
This python script modifies the elastic network in a supplied itp file
'''

import os.path
import argparse
import sys
import numpy as np
# import scipy as spy
from scipy.spatial import distance

### Parser
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

### ### Optional arguments
### Output itp file
parser.add_argument("-o", dest = "output_name", default = False,
                    help = "Name of output itp file.")

### Log file
parser.add_argument("-l", dest = "log_name", default = False,
                    help = "Name of output log file.\n"
                    "A log file will not be created unless a name has been given.\n")

### Input pdb file
parser.add_argument('-f', dest = "struc_name", default = False,
                    help = "A pdb or gro file is required to add networks based on pre-existing distances."
                    "Protein must not stretch across the periodic boundary.\n"
                    "NOTE that due to this script handling distances based on the CG structure, the "
                    "calculated distances will differ slightly from those obtained from Martinize2, "
                    "though this primarily has effect at the 5th decimal point, and can therefore "
                    "likely be ignored as being inconsequential.\n")

### Request modifications
parser.add_argument('-m', dest = "requested_modifications", action='extend', type=str, default = [], nargs='+',
                    help = "Modifications will be done in the order that they are given.\n"
                    
                    "RI: Removes internal elastic networks. Examples:\n"
                    "-m RI:245:282\n"
                    "Removes all elastic networks between all residues from residue 245 to residue 282.\n"
                    "-m RI:245:282-E:266:279\n"
                    "Same as above except residue 266 to 279 are exempt from the selection.\n\n"
                    
                    "RE: Removes external elastic networks. Examples:\n"
                    "-m RE:266:279\n"
                    "Removes all elastic networks between any residue from 266 to 279 and all other residues.\n"
                    "-m RE:266:279-E:270:275\n"
                    "Same as above except residue 270 to 275 are exempt from the selection.\n"
                    "-m RE:266:279-E:400:450\n"
                    "Same as first example except residue 400 to 450 are exempt from the non-selected residues.\n\n"
                    
                    "RB: Removes elastic networks between two groups of residues. Examples:\n"
                    "-m RB:30:244,283:500\n"
                    "Removes all elastic network between the two groups of residues.\n"
                    "(group 1: 30 to 244 and group 2: 283 to 500)\n"
                    "-m RB:30:244,283:500-E:350:400\n"
                    "Removes all elastic network between the two groups of residues.\n"
                    "(group 1: 30 to 39 and 41 to 244 and group 2: 283 to 349 and 401 to 500)\n\n"
                    
                    "RA: Removes all elastic networks associated with this selection.\n"
                    "-m RA:30:244\n\n"
                    
                    "AS: Add elastic networks between two groups of residues. Examples:\n"
                    "-m AS:280:433-dis:0.95-fc:700\n"
                    "Adds elastic networks between the two residues.\n"
                    "dis designates the distance for the bond. If no distance is given, then it will be"
                    "calculated based on pdb file.\n"
                    "fc designates the force constant of the bond (default = 700)\n\n"
                    
                    "AG: Add elastic networks for a group of residues. Examples:\n"
                    "Requires pdb file. Distances will be based on said pdb file.\n"
                    "-m AG:50:75-emax:0.85-emin:0.0-fc:700-replace:all-replaceonly:yes-E:300:350\n"
                    "Generates elastic networks for all of the residues in the selection.\n"
                    "emax, emin and fc set the maximum distance, minimum distance and force constant for bonds.\n"
                    "replace:all designates that all existing networks should be replaced instead of skipped.\n"
                    "replaceonly:yes designates that bonds should only be replaced (prevents addition of new bonds).\n")

########################################## NOT IMPLEMENTED - Can be partially done with "-m AG" command
### ### Modify networks
### Modify Fc and/or distance
### parser.add_argument('-ms', dest = "modify_single", action='append', default = [], nargs='+',
###                     help = "Add elastic networks between two groups of residues. Examples:\n"
###                     "-as M:280:433,dis:0.95,fc:700 Changes the Fc and distance networks between the two residues.\n"
###                     "dis designates the distance for the bond. If not distance given, then it will be"
###                     "calculated based on pdb file. fc designates the force constant of the bond\n")

### Modify multiple
### parser.add_argument('-mm', dest = "modify_multiple", action='append', default = [], nargs='+',
###                     help = "Add elastic networks between two groups of residues. Examples:\n"
###                     "-as M:280:433,dis:0.95,fc:700 Changes the Fc and distance networks between the two residues.\n"
###                     "dis designates the distance for the bond. If not distance given, then it will be"
###                     "calculated based on pdb file. fc designates the force constant of the bond\n")
##########################################


### ### Arguments required for script to run
### Input itp file
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', dest = "input_name",
                           help = "Name of input itp file.", required = True)

### ### Handle arguments
### Print parser help if no flags provided:
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

### Parses the arguments (checks if required arguments are present)
args = parser.parse_args()

### Input file naming
input_name = args.input_name

### Requested modifications
requested_modifications = args.requested_modifications

### Output file naming
if args.output_name != False:
    output_name = args.output_name
else:
    output_name = input_name[:-4] + "_modified" + ".itp"

### structure file (pdb or gro)
#if args.struc_name != False:
struc_name = args.struc_name

### Log file naming
log_file = []
log_modifications = []
log_settings = []
log_removed_networks = [] ### Currently not printed to LOG file
if args.log_name != False:
    create_log = True
    log_name = args.log_name
    log_file.extend("Log file: " + log_name + "\n")
else:
    create_log = False

log_file.extend("itp file: " + input_name + "\n")
log_file.extend("Output file: " + output_name + "\n")
if args.struc_name != False:
    log_file.extend("Structure file: " + log_name + "\n")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### Parser Handling done

### Open files
def file_reader(file_name):
    file = open(file_name, "r")
    file_read = [i for i in file]
    file.close()
    return file_read
itp_file = file_reader(input_name)

### Finds the appropriate lines for ElNet modifications
itp_separators_dict = {}
atom_line_found = False
ElNet_line_found = False
print("Processing itp file")
for line_number, line in enumerate(itp_file):
    ### [:-1] is to remove the \n from all lines
    if "[ atoms ]" in line[:-1]:
        itp_separators_dict["atoms_start"] = line_number + 1
        atom_line_found = True
    if atom_line_found == True and line[:-1] == "":
        itp_separators_dict["atoms_end"] = line_number
        atom_line_found = False
    
    if "Rubber band" in line[:-1]:
        itp_separators_dict["ElNet_start"] = line_number + 1
        ElNet_line_found = True
    if ElNet_line_found == True and line[:-1] == "":
        itp_separators_dict["ElNet_end"] = line_number
        ElNet_line_found = False

itp_separators_sorted_list = sorted([(name, line_nr) for name, line_nr in itp_separators_dict.items()], key=lambda x: x[1])

### Creates dictionary of residues and atoms from itp file
residue_dict = {}
list_of_atom_numbers = []

for atom in itp_file[itp_separators_dict["atoms_start"]:itp_separators_dict["atoms_end"]]:
    atom_split = atom.split()
    if "BB" in atom_split:
        residue_dict[str(atom_split[2])] = atom_split[0]
        list_of_atom_numbers.append(atom_split[0])

### Creates list of existing ElNets and ensures that they are properly formatted
ElNets = itp_file[itp_separators_dict["ElNet_start"]:itp_separators_dict["ElNet_end"]]
w0 = max(len(line.split()[0]) for line in ElNets)
w1 = max(len(line.split()[1]) for line in ElNets)
w2 = max(len(line.split()[2]) for line in ElNets)
w3 = max(len(line.split()[3]) for line in ElNets)
w4 = max(len(line.split()[4]) for line in ElNets)
ElNets_split = [line.split() for line in ElNets]
ElNets = ['{i0: >{w0}} {i1: >{w1}} {i2: >{w2}} {i3: >{w3}} {i4: >{w4}}\n'.format(i0=i0, i1=i1, i2=i2, i3=i3, i4=i4, w0=w0, w1=w1, w2=w2, w3=w3, w4=w4)
                    for i0, i1, i2, i3, i4 in ElNets_split]

### Creates list of coordinates for each atom if structure file is present and array of distances between atoms
if struc_name != False:
    struc_file = file_reader(struc_name)
    ### ### Following lines reads structure files and extract the "BB" beads from them
    ## If pdb file
    if struc_name[-3:] == "pdb":
        line = struc_file[0]
        struc_file_split = [(line[22:26].strip(), line[6:11].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip()) for line in struc_file
                           if line[0:6].strip() == "ATOM"
                           and line[6:11].strip() == residue_dict[line[22:26].strip()] ### Bead (atom) in residue dict
                           and line[12:16].strip() == "BB"] ### Double ensure that bead is "BB"
    
    ## If gro file
    if struc_name[-3:] == "gro":
        struc_file_split = [(line[0:5].strip(), line[15:20].strip(), line[20:28].strip(), line[28:36].strip(), line[36:44].strip()) for line in struc_file
                           if line[10:15].strip() == "BB"
                           and line[15:20].strip() == residue_dict[line[0:5].strip()]] ### Bead (atom) in residue dict
    
    coordinates_list = [(float(x) * 0.1, float(y) * 0.1, float(z) * 0.1) for resid, atom, x, y, z in struc_file_split]
    
    
    distances_array = distance.cdist(coordinates_list, coordinates_list, "euclidean")
    distances_dict = {}
    for i in range(len(distances_array)):
        for j in range(len(distances_array)):
            if i != j: 
                distances_dict[(str(list_of_atom_numbers[i]), str(list_of_atom_numbers[j]))] = str(round(distances_array[i, j], 5))

### ### ### ### ### ### ### ### Functions ### ### ### ### ### ### ### ###
### ### ### ### General functions ### ### ### ###
### Removes exemptions from atom list
def exemptions_remover(atom_list, exemptions):
    '''
    Removes atoms from atom_list based on the atoms that are present in exemptions.
    atom_list is a list containing atom numbers.
    exemptions: [(res1, res2)]: A list containing sets of residues that should be excluded from the selection.
    '''
    ### Finds all the atoms that are exempt from the selection
    exemptions_list = []
    for exemption in exemptions:
        start_res, end_res = exemptions
        if int(start_res) > int(end_res):
            start_res, end_res = end_res, start_res
        atom_exempt_list = [residue_dict[str(res)] for res in range(int(start_res), int(end_res) + 1)]
        exemptions_list.extend(atom_exempt_list)

    ### Removes exempt atoms from selection
    for atom_ex in exemptions_list:
        if atom_ex in atom_list:
            atom_list.remove(atom_ex)
    
    return atom_list

### Removes networks from ElNets
def network_remover(remover_list):
    '''
    Removes networks from the ElNets list.
    remover_list is a list containing tuples of the two atoms for which networks should be removed.
    (atom_1, atom_2)
    '''
    log_modifications.append("Will look for " + str(len(remover_list)) + " potential elastic network bonds\n")
    print(log_modifications[-1][:-1])
    
    elastic_remove_counter = 0
    for rubber_band in reversed(ElNets):
        rbs = rubber_band.split()
        if (rbs[0], rbs[1]) in remover_list:
            elastic_remove_counter += 1
            ElNets.remove(rubber_band)
            log_removed_networks.append(rubber_band)
    
    log_modifications.append(str(elastic_remove_counter) + " elastic network bonds were removed\n")
    print(log_modifications[-1])

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### Removals
### ### ### Internal elastic network removal calculator:
def remove_internal(RI_dict):
    '''
    Removes internal elastic networks for a set of residues.
    RI_dict contains the followings keys and values:
    "RI": (res1, res2): A tuple containing two two strings.
    The set of residues for which internal networks between should be removed.
    "E": [(res1,res2)]: A list containing sets of residues that should be excluded from the selection.
    '''
    start_res, end_res = RI_dict["RI"]
    
    ### Finds all the atoms in the removal selection 
    if int(start_res) > int(end_res):
        start_res, end_res = end_res, start_res
    atom_list = [residue_dict[str(res)] for res in range(int(start_res), int(end_res) + 1)]
    
    ### Removes all the atoms that are exempt from the selection
    if RI_dict["E"] != []:
        for exemptions in RI_dict["E"]:
            atom_list = exemptions_remover(atom_list = atom_list, exemptions = exemptions)
    
    ### Removes atoms form ElNets list
    remover_list = [(atom_1, atom_2) for atom_1 in atom_list for atom_2 in atom_list if atom_1 != atom_2]
    network_remover(remover_list)
        
### ### ### External elastic network removal calculator:
def remove_external(RE_dict):
    '''
    Removes external elastic networks for a set of residues.
    RE_dict contains the followings keys and values:
    "RE": (res1, res2): A tuple containing two two strings.
    The set of residues for which external networks between should be removed.
    "E": [(res1,res2)]: A list containing sets of residues that should be excluded from the selection.
    '''
    
    start_res, end_res = RE_dict["RE"]
    
    ### Finds all the atoms in the selection
    if int(start_res) > int(end_res):
        start_res, end_res = end_res, start_res
    atom_list = [residue_dict[str(res)] for res in range(int(start_res), int(end_res) + 1)]
    
    ### Finds all atoms not in selection
    atom_ext_list = [residue_dict[str(res)] for res in list(residue_dict.keys()) if residue_dict[str(res)] not in atom_list]

    ### Removes all the atoms that are exempt from the selection
    if RE_dict["E"] != []:
        for exemptions in RE_dict["E"]:
            atom_list = exemptions_remover(atom_list = atom_list, exemptions = exemptions)
            atom_ext_list = exemptions_remover(atom_list = atom_ext_list, exemptions = exemptions)
    
    remover_list = []
    remover_list.extend([(atom_1, atom_2) for atom_1 in atom_list for atom_2 in atom_ext_list if atom_1 != atom_2])
    remover_list.extend([(atom_1, atom_2) for atom_1 in atom_ext_list for atom_2 in atom_list if atom_1 != atom_2])
    network_remover(remover_list)
        
### ### ### Between groups elastic network removal calculator:
def remove_between(RB_dict):
    '''
    Removes elastic networks between two sets of residues.
    RB_dict contains the followings keys and values:
    "RB": ((res1, res2), (res3, res4)): A tuple containing two tuples each with two strings.
    The two sets of residues for which networks between should be removed.
    "E": [(res1,res2)]: A list containing sets of residues that should be excluded from the selection.
    '''
    
    groups = RB_dict["RB"]
    group_atom_lists = []
        
    for group in groups:
        start_res, end_res = group
        
        ### Finds all the atoms in the selection
        if int(start_res) > int(end_res):
            start_res, end_res = end_res, start_res
        atom_list = [residue_dict[str(res)] for res in range(int(start_res), int(end_res) + 1)]

        ### Removes all the atoms that are exempt from the selection and non-selection
        if RB_dict["E"] != []:
            for exemptions in RB_dict["E"]:
                atom_list = exemptions_remover(atom_list = atom_list, exemptions = exemptions)

        group_atom_lists.append(atom_list)
    
    remover_list = []
    remover_list.extend([(atom_1, atom_2) for atom_1 in group_atom_lists[0] for atom_2 in group_atom_lists[1] if atom_1 != atom_2])
    remover_list.extend([(atom_1, atom_2) for atom_1 in group_atom_lists[1] for atom_2 in group_atom_lists[0] if atom_1 != atom_2])
    network_remover(remover_list)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### Additions
### ### ### Adds a single bond between two specific residues
def add_single(AS_dict):
    '''
    Creates a single elastic network bond between two residues.
    AS_dict contains the followings keys and values:
    "AS": (res1, res2): A tuple containing two strings.
    The two residues for which a network should be created.
    "distance": A string containing either "automatic" or a number.
    "fc": A string containing "default" or a number as a string.
    '''
    res1, res2 = AS_dict["AS"]
    
    atom1, atom2 = residue_dict[res1], residue_dict[res2]
    
    log_modifications.append(" ".join(["Preparing to add a bond between atom", atom1, "and", atom2]) + "\n")
    print(log_modifications[-1][:-1])
    
    ### Distance
    if AS_dict["distance"] == "automatic":
        distance = distances_dict[(atom1, atom2)]
        log_modifications.append("No distance was given. Setting it to " + distance + "\n")
        print(log_modifications[-1][:-1])
    else:
        distance = AS_dict["distance"]

    ### Force constant
    if AS_dict["fc"] == "default":
        force_constant = "700.0"
        log_modifications.append("No force constant given. Setting it to " + force_constant + "\n")
        print(log_modifications[-1][:-1])
    else:
        if "." in AS_dict["fc"]:
            force_constant = AS_dict["fc"]
        else:
            force_constant = AS_dict["fc"] + ".0"

    ### Checking if bond already exists
    log_modifications.append("Checking for existing bond to be replaced.\n")
    print(log_modifications[-1][:-1])
    
    replaced_questionmark = False
    for line, rubber_band in enumerate(ElNets):
        rbs = rubber_band.split()
        if (rbs[0], rbs[1]) == (atom1, atom2):
            ElNets[line] = " ".join([atom1, atom2, "1", distance, force_constant]) + "\n"
            replaced_questionmark = True
            log_modifications.append("Existing bond found. Replacing it.\n")
            print(log_modifications[-1][:-1])
            break
    if replaced_questionmark == False:
        ElNets.append(" ".join([atom1, atom2, "1", distance, force_constant]) + "\n")
        log_modifications.append("No existing bond found. Creating a new one.\n")
        print(log_modifications[-1][:-1])

### ### ### Generates bonds for a range of residues
def add_generated(AG_dict):
    '''
    Creates new elastic networks for all residues in the selection.
    AG_dict contains the followings keys and values:
    "AG": (res1, res2): A tuple containing two strings.
    The set of residues for which new networks should be created.
    "emax": A string containing a number.
    Sets the upper distance limit for generated elastic networks. "Default" = 0.9
    "emin": A string containing a number.
    Sets the lower distance limit for generated elastic networks. "Default" = 0.5
    "fc": A string containing a number.
    Sets the force constant for the generated elastic networks. "Default" = 700.0
    "replace": A string containing "all", "none" or a list of tuples containg either two strings (res1, res2)
    or two tuples with two strings each. One tuple designates residues for which all networks should be replaced
    while two tuples designates sets of residues for which networks between should be replaced.
    "E": [(res1,res2)]: A list containing sets of residues that should be excluded from the selection.
    '''
    start_res, end_res = AG_dict["AG"]
    if int(start_res) > int(end_res):
        start_res, end_res = end_res, start_res
    
    ### ### Handle replacements
    replace_any_list = []
    replace_pair_list = []
    
    ### List of atoms where networks should be replaced
    if AG_dict["replace"] == ["all"]:
        replace_any_list.extend([residue_dict[str(res)] for res in range(int(start_res), int(end_res) + 1)])
    elif AG_dict["replace"] == ["none"]: ### Do nothing if "none"
        pass
    else: # len(AG_dict["replace"]) > 1:
        for replacement in AG_dict["replace"]:
            if len(replacement) == 1:
                res1, res2 = replacement
                if res2 > res1:
                    res1, res2 = res2, res1
                replace_any_list.extend([residue_dict[str(res)] for res in range(int(res1), int(res2) + 1)])
            elif len(replacement) == 2:
                res1, res2, res3, res4 = replacement
                if res2 > res1:
                    res1, res2 = res2, res1
                if res4 > res3:
                    res3, res4 = res4, res3
                replace_pair_list_part1 = [residue_dict[str(res)] for res in range(int(res1), int(res2) + 1)]
                replace_pair_list_part2 = [residue_dict[str(res)] for res in range(int(res3), int(res4) + 1)]
                replace_pair_list.extend([(resA, resB) for resA in replace_pair_list_part1 for resB in replace_pair_list_part2])
                replace_pair_list.extend([(resA, resB) for resA in replace_pair_list_part2 for resB in replace_pair_list_part1])
    
    ### ### Handle standard atom and residue selection processing
    ### Finds all the atoms and residues in the selection
    res_list_selection = [int(resid) for resid in range(int(start_res), int(end_res) + 1)]
    atom_list_selection = [residue_dict[str(res)] for res in range(int(start_res), int(end_res) + 1)]
    
    ### Finds all atoms in the system
    atom_list_all = [residue_dict[str(res)] for res in list(residue_dict.keys())]
    res_list_all = [int(resid) for resid in list(residue_dict.keys())]
    
    ### Removes all the atoms that are exempt from the selection
    if AG_dict["E"] != []:
        for exemptions in AG_dict["E"]:
            atom_list_selection = exemptions_remover(atom_list = atom_list_selection, exemptions = exemptions)
            atom_list_all = exemptions_remover(atom_list = atom_list_all, exemptions = exemptions)
    
    atom_lists_combined = []
    atom_lists_combined.extend([(atom1, atom2) for atom1 in atom_list_selection for atom2 in atom_list_all if atom1 != atom2])
    atom_lists_combined.extend([(atom1, atom2) for atom1 in atom_list_all for atom2 in atom_list_selection if atom1 != atom2])
    
    ### Ensuring that residues are more than 2 residues apart ("i -> i+1 and i -> i+2" rule)
    addition_res_list = []
    addition_res_list.extend([(res1, res2) for res1 in res_list_selection for res2 in res_list_all])
    addition_res_list.extend([(res1, res2) for res1 in res_list_all for res2 in res_list_selection])
    addition_atom_list = sorted(list(set(sorted([(residue_dict[str(res1)], residue_dict[str(res2)])
                                                for res1, res2 in addition_res_list
                                                if (res1 < res2 - 2)],
                                            key=lambda x: (int(x[0]), int(x[1])))
                                        )
                                    ),
                                key=lambda x: (int(x[0]), int(x[1])))
    
    ### ### Handles combining above sections
    ### Creates the list of atom connections taking into account both exemptions, replacements and the no "i -> i+1 and i -> i+2" connections rule
    ### This is somewhat slow
    addition_list = []
    replacement_list = []
    for atom1, atom2 in addition_atom_list: ### Atoms from no "i -> i+1 and i -> i+2" residue connections rule
        if (atom1, atom2) in atom_lists_combined: ### Checks if atom bond is in bonds to be created
            if (atom1, atom2) in [(line.split()[0], line.split()[1]) for line in ElNets]: ### Cecks if bond already exists
                if len(replace_any_list) > 0 and (atom1 in replace_any_list or atom2 in replace_any_list):
                    replacement_list.append((atom1, atom2))
                elif len(replace_pair_list) > 0 and (atom1, atom2) in replace_pair_list:
                    replacement_list.append((atom1, atom2))
            elif AG_dict["replaceonly"] != "yes":
                ### If bond does not already exist then create new one and only replacements not active
                addition_list.append((atom1, atom2))
    
    ### ### Handle distances, pre-existing bonds, upper and lower limits and force constants
    ### Ensuring the force constant has a decimal point for consistency
    if "." in AG_dict["fc"]:
        force_constant = AG_dict["fc"]
    else:
        force_constant = AG_dict["fc"] + ".0"
    
    ### Adding new bonds
    new_additions_counter = 0
    if AG_dict["replaceonly"] != "yes":
        for atom1, atom2 in addition_list:
            distance = distances_dict[(atom1, atom2)]
            if AG_dict["emin"] > distance or distance > AG_dict["emax"]: ### Continue if bond distance outside emin:emax range
                continue
            else:
                ElNets.append(" ".join([atom1, atom2, "1", distance, force_constant]) + "\n")
                new_additions_counter += 1
    
    ### Replacing bonds
    replacement_counter = 0
    for atom1, atom2 in replacement_list:
        distance = distances_dict[(atom1, atom2)]
        if AG_dict["emin"] > distance or distance > AG_dict["emax"]: ### Continue if bond distance outside emin:emax range
            continue
        else:
            for line, rubber_band in enumerate(ElNets):
                rbs = rubber_band.split()
                if (rbs[0], rbs[1]) == (atom1, atom2):
                    ElNets[line] = " ".join([atom1, atom2, "1", distance, force_constant]) + "\n"
                    replacement_counter += 1
                    continue
    
    
    log_modifications.append("Added " + str(new_additions_counter) + " new elastic network bonds\n")
    print(log_modifications[-1][:-1])
    log_modifications.append("Replaced " + str(replacement_counter) + " elastic network bonds\n")
    print(log_modifications[-1][:-1])

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### Processing
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
log_settings.append("Following modifications were requested and will be processed in the shown order:\n")
print(log_settings[-1])
for mod in requested_modifications:
    log_settings.append(mod + "\n")
    log_modifications.append("\nWorking on: " + mod + "\n")
    print(log_modifications[-1][:-1])
    ### Run the respective modification function
    identifier = mod.split(":", maxsplit = 1)[0]
    
    ### Remove INTERNAL elastic networks and Remove ALL elastic networks
    if identifier == "RI" or identifier == "RA": ### Example: "RI:245:282-E:266:279"
        RI_dict = {}
        exemptions = []
        for part in mod.split("-"):
            command, settings = part.split(":", maxsplit = 1)
            if command == "RI" or identifier == "RA":
                res1, res2 = settings.split(":")
                RI_dict["RI"] = (res1, res2)
            if command == "E":
                res1, res2 = settings.split(":")
                exemptions.append((res1, res2))
        RI_dict["E"] = exemptions
        remove_internal(RI_dict)
    
    ### Remove EXTERNAL elastic networks and Remove ALL elastic networks
    if identifier == "RE" or identifier == "RA": ### Example: "RE:463:471-E:350:400"
        RE_dict = {}
        exemptions = []
        for part in mod.split("-"):
            command, settings = part.split(":", maxsplit = 1)
            if command == "RE" or identifier == "RA":
                res1, res2 = settings.split(":")
                RE_dict["RE"] = (res1, res2)
            if command == "E":
                res1, res2 = settings.split(":")
                exemptions.append((res1, res2))
        RE_dict["E"] = exemptions
        remove_external(RE_dict)
    
    ### Remove elastic networks BETWEEN sets of residues
    if identifier == "RB": ### Example: "RB:30:244,283:500-E:350:400"
        RB_dict = {}
        exemptions = []
        for part in mod.split("-"):
            command, settings = part.split(":", maxsplit = 1)
            if command == "RB":
                set1, set2 = settings.split(",")
                res1, res2, res3, res4 = set1.split(":") + set2.split(":")
                RB_dict["RB"] = [(res1, res2), (res3, res4)]
            if command == "E":
                res1, res2 = settings.split(":")
                exemptions.append((res1, res2))
        RB_dict["E"] = exemptions
        remove_between(RB_dict)
    
    ### Add SINGLE elastic network bond
    if identifier == "AS": ### Example: "AS:280:433-dis:0.95-fc:700"
        AS_dict = {}
        AS_dict["distance"] = "automatic"
        AS_dict["fc"] = "default"
        for part in mod.split("-"):
            command, settings = part.split(":", maxsplit = 1)
            if command == "AS":
                res1, res2 = settings.split(":")
                assert abs(int(res1) - int(res2)) > 2, \
                        "Residues are not more than 2 residues apart. Elastic networks are usually not used between residues that are this close."
                AS_dict["AS"] = (res1, res2)
            if command == "dis":
                AS_dict["distance"] = settings
            if command == "fc":
                AS_dict["fc"] = settings
        add_single(AS_dict)
    
    ## Add GENERATED elastic networks for residue selection
    if identifier == "AG": ### Example: "AG:50:75-emax:0.85-emin:0.5-fc:700-replace:60:65-E:300:350"
        AG_dict = {}
        exemptions = []
        AG_dict["emax"] = "0.9" ### Default = 0.9
        AG_dict["emin"] = "0.0" ### Default = 0.0
        AG_dict["fc"] = "700" ### Default = 700
        replacements = []
        AG_dict["replace"] = ["all"] ### "all", "none", selection:(res1:res2), selection between:(res1:res2,res3:res4)
        AG_dict["replaceonly"] = "no" ### "yes", "no": Only replaces, does not create new bonds if "yes"
        for part in mod.split("-"):
            command, settings = part.split(":", maxsplit = 1)
            if command == "AG":
                res1, res2 = settings.split(":")
                AG_dict["AG"] = (res1, res2)
            if command == "E":
                res1, res2 = settings.split(":")
                exemptions.append((res1, res2))
            if command == "emax":
                AG_dict["emax"] = settings
            if command == "emin":
                AG_dict["emin"] = settings
            if command == "fc":
                AG_dict["fc"] = settings
            if command == "replaceonly":
                AG_dict["replaceonly"] = settings
            if command == "replace":
                if len(settings.split(",")) == 1 and type(settings) == tuple:
                    res1, res2 = settings.split(":")
                    replacements.append((res1, res2))
                elif len(settings.split(",")) == 2 and type(settings) == tuple:
                    res1, res2, res3, res4 = settings.split(":")
                    replacements.append(((res1, res2), (res3, res4)))
                else:
                    replacements.append(settings) ### Only here for the assertion below
        if len(replacements) > 0:
            AG_dict["replace"] = replacements
            ### Ensure that "all" and "none" are not combined with residue specific commands
            if len(AG_dict["replace"]) > 1:
                assert "all" in AG_dict["replace"] or "none" in AG_dict["replace"], \
                "You cannot combine 'all' and/or 'none' with residue specific commands for replace"
        AG_dict["E"] = exemptions
        add_generated(AG_dict)

### ### ### Log writer
log_file.extend(log_settings)
log_file.extend(log_modifications)
# log_file.extend("The following lines show the removed elastic networks\n")
# log_file.extend(log_removed_networks)

### ### ### Output file handling
### Formatting and sorting new elastic network system
w0 = max(len(line.split()[0]) for line in ElNets)
w1 = max(len(line.split()[1]) for line in ElNets)
w2 = max(len(line.split()[2]) for line in ElNets)
w3 = max(len(line.split()[3]) for line in ElNets)
w4 = max(len(line.split()[4]) for line in ElNets)
ElNets_split = [line.split() for line in ElNets]
ElNets_formatted = ['{i0: >{w0}} {i1: >{w1}} {i2: >{w2}} {i3: >{w3}} {i4: >{w4}}\n'.format(i0=i0, i1=i1, i2=i2, i3=i3, i4=i4, w0=w0, w1=w1, w2=w2, w3=w3, w4=w4)
                    for i0, i1, i2, i3, i4 in ElNets_split]

ElNets_sorted = sorted(ElNets_formatted, key=lambda x: (x[0], x[1]))
### Inserting new elastic network system file into output itp file
before = itp_file[:itp_separators_dict["ElNet_start"]]
after = itp_file[itp_separators_dict["ElNet_end"]:]
new_itp_file = before + ElNets_sorted + after

### Checks if output file already exists and backs it up
def backupper(output_file_name):
    output_file_split = output_file_name.split("/")
    output_path = ""
    output_name = output_file_split[-1]
    if len(output_file_split) > 1:
        for i in range(len(output_file_split) - 1):
            output_path += output_file_split[i] + "/"
    if os.path.exists(output_file_name):
        print("File " + output_file_name + " already exists. Backing it up")
        number = 1
        while True:
            if os.path.exists(output_path + "#" + output_name + "." + str(number) + "#"):
                number += 1
            else:
                os.rename(output_file_name, output_path + "#" + output_name + "." + str(number) + "#")
                break

print("\n")
### Output itp file
backupper(output_name)
print("Writing output file: " + output_name)
new_file = open(output_name, "w")
for line in new_itp_file:
    new_file.write(line)
new_file.close()

### Output log file
if create_log == True:
    backupper(log_name)
    print("Writing log file: " + output_name)
    new_file = open(log_name, "w")
    for line in log_file:
        new_file.write(line)
    new_file.close()
print("\n")
