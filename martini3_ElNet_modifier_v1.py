# %reset -f

'''
This python script modifies the elastic network in a supplied itp file
Input:
-i [required] Name of input itp file.
-o [optional] Name of output itp file.
-l [optional] Name of output log file containing removed elastic network lines. A log file will not be created unless a name has been given.
------------
At least one of the following commands must be used before the script will run.
For the following commands: All residus mentioned by the commands are included in the selection.
------------
-ri [optional] Remove internal elastic networks. Examples:
    -ri R:245:282" Removes all elastic networks between all residues from residue 245 to residue 282.
    -ri R:245:282-E:266:279 Same as above except residue 266 to 279 are exempt from the selection.

-re [optional] Remove external elastic networks. Examples:
    -re R:245:282 Removes all elastic networks between any residue from 245 to 282 and all other residues.
    -re R:245:282-E:266:279 Same as above except residue 266 to 279 are exempt from the selection.
    -re R:245:282-E:400:450 Same as first example except residue 400 to 450 are exempt from the non-selected residues.

-rb [optional] Remove elastic networks between two groups of residues. Examples:
    -rb R:30:244,R:283:500 Removes all elastic network between the two groups of residues.
    (group 1: 30 to 244 and group 2: 283 to 500)
    -rb R:30:244,R:283:500-E:350:400 Removes all elastic network between the two groups of residues.
    (group 1: 30 to 244 and group 2: 283 to 349 and 401 to 500)
'''

import os.path
import argparse
import sys

### Parser
# parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser = argparse.ArgumentParser()

### Optional arguments
parser.add_argument("-o", dest = "output_name", default = None,
                    help = "Name of output itp file.")
parser.add_argument("-l", dest = "log_name", default = None,
                    help = "Name of output log file containing removed elastic network lines.\n"
                    "A log file will not be created unless a name has been given.\n")

### Remove networks flags
# Remove internal networks
parser.add_argument('-ri', dest = "remove_internal", action='append', default = [], nargs='+',
                    help = "Remove internal elastic networks. Examples:\n"
                    "-ri R:245:282 Removes all elastic networks between all residues from residue 245 to residue 282.\n"
                    "-ri R:245:282-E:266:279 Same as above except residue 266 to 279 are exempt from the selection.\n")
# Remove external networks
parser.add_argument('-re', dest = "remove_external", action='append', default = [], nargs='+',
                    help = "Remove external elastic networks. Examples:\n"
                    "-re R:266:279 Removes all elastic networks between any residue from 266 to 279 and all other residues.\n"
                    "-re R:266:279-E:270:275 Same as above except residue 270 to 275 are exempt from the selection.\n"
                    "-re R:266:279-E:400:450 Same as first example except residue 400 to 450 are exempt from the non-selected residues.\n")
# Remove networks between two groups
parser.add_argument('-rb', dest = "remove_between", action='append', default = [], nargs='+',
                    help = "Remove elastic networks between two groups of residues. Examples:\n"
                    "-rb R:30:244,R:283:500 Removes all elastic network between the two groups of residues.\n"
                    "(group 1: 30 to 244 and group 2: 283 to 500)\n"
                    "-rb R:30:244-E40:50,R:283:500-E:350:400 Removes all elastic network between the two groups of residues.\n"
                    "(group 1: 30 to 39 and 41 to 244 and group 2: 283 to 349 and 401 to 500)\n")

### Arguments required for script to run
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', dest = "input_name", help = "Name of input itp file.", required = True)

### Print parser help if no flags provided:
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

### Parses the arguments (checks if required arguments are present)
args = parser.parse_args()

### Input file naming
input_name = args.input_name
### Networks to be removed:
rem_int_input = args.remove_internal
rem_ext_input = args.remove_external
rem_btw_input = args.remove_between

### Output file naming
if args.output_name != None:
    output_name = args.output_name
else:
    output_name = input_name[:-4] + "_reprotonated" + ".pdb"

### Log file naming
log_file = []
log_file.extend("Output itp file: " + output_name + "\n")
create_log = False
if args.log_name != None:
    create_log = True
    log_name = args.log_name
    log_file.extend("Log file: " + log_name + "\n")

### Checks if changes have been requested
if rem_int_input == [] and rem_ext_input == [] and rem_btw_input == []:
    assert True, "No elastic network changes were requested. Aborting script."

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### Parser Handling done




### ### ### ### ### ### ### ### ### ### Test files

# input_name = "martini3_ElNet_modifier_test_files/AtSUC1_res152_standard.itp"
# output_name = "martini3_ElNet_modifier_test_files/AtSUC1_res152_standard_modified.itp"

# rem_int_input = [[]]
# rem_ext_input = [[]]
# rem_btw_input = [[]]

# # rem_int_input = [["R:245:282-E:266:279"]]
# # rem_ext_input = [["R:245:282"], ["R:463:471"]]
# rem_ext_input = [["R:266:279"]]
# # rem_btw_input = [["R:30:244,R:283:500"]]

# log_file = []
# create_log = True
# log_name = "test.log"
# log_file.extend("Log file: " + log_name + "\n")

### ### ### ### ### ### ### ### ### ### Test files over

print("Following elastic network removals were requested:")
if rem_int_input != [[]]:
    for i in rem_int_input:
        print("Internal: " + str(i))

if rem_ext_input != [[]]:
    for i in rem_ext_input:
        print("External: " + str(i))

if rem_btw_input != [[]]:
    for i in rem_btw_input:
        print("Between: " + str(i))

### Open files
def file_reader(file_name):
    file = open(file_name, "r")
    file_read = [i for i in file]
    file.close()
    return file_read
itp_file = file_reader(input_name)

### Finds the appropriate lines
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

### Creates dictionary of residues and atoms
residue_dict = {}
for atom in itp_file[itp_separators_dict["atoms_start"]:itp_separators_dict["atoms_end"]]:
    atom_split = atom.split()
    if str(atom_split[2]) not in residue_dict.keys() and "BB" in atom_split:
        residue_dict[str(atom_split[2])] = [atom_split[0]]
    elif str(atom_split[2]) in residue_dict.keys() and "BB" in atom_split:
        residue_dict[str(atom_split[2])].append(atom_split[0])
        

### Removes exemptions from atom list
def exemptions_remover(ElNet_list, exemptions):
    ### Finds all the atoms that are exempt from the selection
    exemptions_list = []
    for exemption in exemptions:
        rem_type, start_res, end_res = exemption.split(":")
        atom_exempt_list = [atom for res in range(int(start_res), int(end_res) + 1) for atom in residue_dict[str(res)]]
        exemptions_list.extend(atom_exempt_list)

    ### Removes exempt atoms from selection
    for atom_ex in exemptions_list:
        if atom_ex in ElNet_list:
            ElNet_list.remove(atom_ex)
    
    return ElNet_list



### ### List containing all networks to be removed
print("Creating elastic network removal list")
Main_remover_list = []

### ### Writing requested changes to log file
log_requested_changes = []

### ### ### Internal elastic network calculator:
if rem_int_input != [[]]:
    for internal in rem_int_input:
        exemptions_present = False
        input_split = internal[0].split("-")
        
        log_requested_changes.append("-ri " + internal[0] + "\n")

        remove_list = input_split[0]
        if len(input_split) > 1:
            exemptions = input_split[1:]
            exemptions_present = True

        rem_type, start_res, end_res = remove_list.split(":")

        ### Finds all the atoms in the removal selection 
        atom_list = [str(atom) for res in range(int(start_res), int(end_res) + 1) for atom in residue_dict[str(res)]]

        ### Finds all the atoms that are exempt from the selection
        if exemptions_present == True:
            atom_list = exemptions_remover(ElNet_list = atom_list, exemptions = exemptions)

        Main_remover_list.extend([(atom_1, atom_2) for atom_1 in atom_list for atom_2 in atom_list if atom_1 != atom_2])


### ### ### External elastic network calculator:
if rem_ext_input != [[]]:
    for external in rem_ext_input:
        exemptions_present = False
        input_split = external[0].split("-")
        
        log_requested_changes.append("-re " + external[0] + "\n")
        
        remove_list = input_split[0]
        if len(input_split) > 1:
            exemptions = input_split[1:]
            exemptions_present = True

        rem_type, start_res, end_res = remove_list.split(":")

        ### Finds all the atoms in the selection
        atom_list = [str(atom) for res in range(int(start_res), int(end_res) + 1) for atom in residue_dict[str(res)]]

        ### Finds all atoms not in selection
        atom_ext_list = [str(atom) for res in list(residue_dict.keys()) for atom in residue_dict[str(res)] if atom not in atom_list]
        
        ### Finds all the atoms that are exempt from the selection
        if exemptions_present == True:
            atom_list = exemptions_remover(ElNet_list = atom_list, exemptions = exemptions)
            atom_ext_list = exemptions_remover(ElNet_list = atom_ext_list, exemptions = exemptions)
                
        Main_remover_list.extend([(atom_1, atom_2) for atom_1 in atom_list for atom_2 in atom_ext_list if atom_1 != atom_2])
        Main_remover_list.extend([(atom_1, atom_2) for atom_1 in atom_ext_list for atom_2 in atom_list if atom_1 != atom_2])


### ### ### Between groups elastic network calculator:
if rem_btw_input != [[]]:
    for between in rem_btw_input:
        group_atom_lists = []
        
        log_requested_changes.append("-rb " + between[0] + "\n")
        
        groups = [[group] for group in between[0].split(",")]
        for group in groups:
            exemptions_present = False
            input_split = group[0].split("-")

            remove_list = input_split[0]
            if len(input_split) > 1:
                exemptions = input_split[1:]
                exemptions_present = True

            rem_type, start_res, end_res = remove_list.split(":")

            ### Finds all the atoms in the selection
            atom_list = [str(atom) for res in range(int(start_res), int(end_res) + 1) for atom in residue_dict[str(res)]]

            ### Finds all the atoms that are exempt from the selection and non-selection
            if exemptions_present == True:
                atom_list = exemptions_remover(ElNet_list = atom_list, exemptions = exemptions)

            group_atom_lists.append(atom_list)

        Main_remover_list.extend([(atom_1, atom_2) for atom_1 in group_atom_lists[0] for atom_2 in group_atom_lists[1] if atom_1 != atom_2])
        Main_remover_list.extend([(atom_1, atom_2) for atom_1 in group_atom_lists[1] for atom_2 in group_atom_lists[0] if atom_1 != atom_2])

log_file.extend("The following lines show the requested removals\n")
log_file.extend(log_requested_changes)

print("Removal list calculated")
### Removes duplicates from list
len_MRL_before = len(Main_remover_list)
Main_remover_list = list(dict.fromkeys(Main_remover_list))
len_MRL_after = len(Main_remover_list)

if len_MRL_before - len_MRL_after != 0:
    print("Removed " + str(len_MRL_before - len_MRL_after) + " duplicate entries in the removal list")
print("Will look for " + str(len(Main_remover_list)) + " Potential elastic network bonds")

elastic_remove_counter = 0
log_removed_networks = []
for rubber_band in itp_file[itp_separators_dict["ElNet_start"]:itp_separators_dict["ElNet_end"]]:
    rbs = rubber_band.split()
    if (rbs[0], rbs[1]) in Main_remover_list:
        elastic_remove_counter += 1
        itp_file.remove(rubber_band)
        log_removed_networks.append(rubber_band)
log_file.extend(str(elastic_remove_counter) + " " + "elastic networks were removed")
log_file.extend("The following lines show the removed elastic networks\n")
log_file.extend(log_removed_networks)

print(str(elastic_remove_counter) + " elastic network bonds were removed")

### ### ### Output file handling
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
for line in itp_file:
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
