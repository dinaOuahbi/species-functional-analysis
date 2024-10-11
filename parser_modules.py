#!/usr/bin/python

# Author : Florence Thirion
# Date : 15/09/2019

### Import libraries
import os, sys
import argparse
import re
import datetime
from collections import namedtuple, defaultdict


##################################################
### Create function to parse module defintions ###
##################################################

def load_module_def(module_file, index = [0, 3]):
    "Load the file with the definitions and return a dictionnary based on index"
    ### Initialize dictionary of modules
    module_dict = {}
    with open(module_file) as mf:
        ### Get the definition for each module
        for line in mf:
            mod_name = line.rsplit("\n")[0].rsplit("\t")[index[0]]
            mod_def = line.rsplit("\n")[0].rsplit("\t")[index[1]]
            module_dict[mod_name] = mod_def
    return module_dict

    
def clean_module(mod_def):
    "Remove from the module definition components '-XXX' or ' -- '"
    ### Remove ' -- '
    mod_def = re.sub(r"\s?--\s?", r" ", mod_def).strip()
    ### Remove '-XXX'
    mod_def = re.sub(r"\-\w+", r"", mod_def)
    ### Replace ' ' by '+'
    mod_def = mod_def.replace(" ", "+")
    ### Add parenthesis around the whole mod_def (in case we have such case : KXXX,KYYY
    mod_def = "(" + mod_def + ")"
    ### Remove single parenthesis (do it as long as the module_def is changed, ie, as long as the pattern is found
    ### (the while and the flag are necessary for such patterns that could be found : ((KO1+KO2))
    flag = 0
    while flag == 0:
        mod_def_new = re.sub(r"\(([\w+]+)\)", r"\1", mod_def)
        if mod_def_new == mod_def:
            flag = 1
        else:
            mod_def = mod_def_new
    return mod_def

def replace_submod(mod_def, mod_dict):
    "Replace submodules from mod_def that are in mod_dict"
    mod_def = re.sub(r"\w+", lambda x: mod_dict[x.group()] if x.group() in mod_dict else x.group(), mod_def)
    return mod_def

def find_all_alt(mod_def):
    "Return the list of all alternatives for a module definition"
    ### Define an alternative
    alternative = re.compile(r"(.*)(\(([\w+]+,{1})+[\w+]+\))(.*)")
    list_alt = [mod_def]
    ### Flag to know if some alternatives remain to be solved
    flag = 0
    while flag == 0: ### While alternatives remain
        flag = 1 ### At the beginning of each loop we suppose there are no more alternative
        new_list_alt = [] ### Initialize a new list to keep alternatives that will be solved in the next loop
        for my_def in list_alt:
            my_def = clean_module(my_def)
            ### Look for simple parenthesis (simple parenthesis = alternatives that should be solved)
            res = alternative.match(my_def)
            if res:
                flag = 0 ### Alternatives were found so we suppose other may remain
                new_def = [res.group(1) + k +  res.group(4) for k in res.group(2).strip("()").split(",")]# list of alternatives after resolving simple parenthesis
                new_list_alt += new_def
            else: ### If the def has already been simplified at maximum
                new_list_alt.append(my_def)
        list_alt = new_list_alt
    ### Check there is no more ',' or ')' or '(' in the list of alternatives
    #error_alt = [i for i in list_alt if re.search(r"[,)(]", i)]
    #if len(error_alt) > 0:
    #    print("\n".join(error_alt))
    #    sys.exit(1)
    ### Transform each alternative into set of KOs : module_dict_alt['M000x'] = [set(KO1, KO2), set(KO1, KO3), etc]
    return [set(k.split('+')) for k in set(list_alt)]

def main_parser(module_def_file):
    ### Load file
    module_dict = load_module_def(module_def_file)
    ### Initailize alternative dictionnary
    module_dict_alt = {}
    ### Parse each module
    for my_mod in module_dict:
        ### Which definition
        mod_def = module_dict[my_mod]
        ### Clean defintion
        new_def = clean_module(mod_def)
        ### Remove submodules
        new_def = replace_submod(new_def, module_dict)
        ### Find alternatives
        all_alt = find_all_alt(new_def)
        #print(all_alt)
        #print("\n".join([my_mod + "\t" + " ".join(x) for x in all_alt]))
        module_dict_alt[my_mod] = all_alt
    return module_dict_alt



def main():
    ##########################
    ### Parse command line ###
    ##########################

    ### Create parser
    parser = argparse.ArgumentParser()

    ### Create options
    parser.add_argument("kegg_def", help = "Path to the KEGG module definition file.")

    ### Parse command line
    args = parser.parse_args()

    ### Compute alternatives
    module_dict_alt = main_parser(args.kegg_def)
    ### Print results
    for my_mod in module_dict_alt:
        print("\n".join([my_mod + "\t" + " ".join(x) for x in module_dict_alt[my_mod]]))


if __name__ == "__main__":
    main()
