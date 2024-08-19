"""
Script to perform macro substitutions to the source code for bns_nurates
and export the resulting files for use with Weakrates2 thorn

Usage: python export_thc.py /path/to/bns_nurates /path/to/destination
"""

import re
import os
import shutil
import numpy as np
import subprocess
import argparse

# Function to define macro for KOKKOS_INLINE_FUNCTION as normal inline in bns_nurates.hpp
def replace_macro_bns_nurates_hpp(input_filepath, output_filepath):
    replacement_line = '#define KOKKOS_INLINE_FUNCTION inline\n'
    
    with open(input_filepath, 'r') as file:
        input_lines = file.readlines()

    idx = 0
    with open(output_filepath, 'w') as file:
        inside_block = False
        for line in input_lines:
            if(idx < 16):
                if re.match(r'#ifndef KOKKOS_INLINE_FUNCTION', line):
                    inside_block = True
                    file.write(replacement_line)
                    continue
                elif re.match(r'#endif', line):
                    if inside_block:
                        inside_block = False
                    continue 
                elif inside_block:
                    continue
                
            file.write(line)
            idx = idx + 1

# Function to define macro for tgamma with tgamma in math.h in kernel_pair.hpp
def replace_macro_kernel_pair_hpp(input_filepath, output_filepath):
    replacement_line = '#define func_tgamma(x) tgamma(x)\n'
    
    with open(input_filepath, 'r') as file:
        input_lines = file.readlines()

    idx = 0
    with open(output_filepath, 'w') as file:
        inside_block = False
        for line in input_lines:
            if(idx < 15):
                if re.match(r'#ifdef KOKKOS_FLAG', line):
                    inside_block = True
                    file.write(replacement_line)
                    continue
                elif re.match(r'#endif', line):
                    if inside_block and idx > 0:
                        inside_block = False
                    continue 
                elif inside_block:
                    continue
            
            file.write(line)
            idx = idx + 1

# Compute the diff between a file in the source and destination
def diff_files(source_filepath, destination_filepath):

    filediff = subprocess.run(['git', 'diff', '--no-index', source_filepath, destination_filepath], capture_output=True, text=True)

    return filediff.stdout

# function to copy hpp and cpp files into a destination directory fixing the macros
# Needs full path of top level bns_nurates directory and the destination directory
def refactor_and_copy(source_dir, destination_dir):
    
    files_include = os.listdir(source_dir + '/include/')
    files_include_full = np.array([os.path.join(source_dir + '/include/', file) for file in files_include])
    files_src_integration = os.listdir(source_dir + '/src/integration/')
    files_src_integration_full = np.array([os.path.join(source_dir + '/src/integration/', file) for file in files_src_integration])
    files = np.concatenate((files_include, files_src_integration))
    files_full = np.concatenate((files_include_full, files_src_integration_full))
    diff_stdout = []

    for i in range(0,len(files)):
    
        file = files[i]
        file_full = files_full[i]
        destination_file = destination_dir + '/' + file
    
        print(f"Copying {file_full} to {destination_file}")
        
        if(file == "bns_nurates.hpp"):
            replace_macro_bns_nurates_hpp(file_full, destination_file)
        elif(file == "kernel_pair.hpp"):
            replace_macro_kernel_pair_hpp(file_full, destination_file)
        else:    
            shutil.copy(file_full, destination_file)

        diff_stdout.append(diff_files(file_full, destination_file))

    print("=============")
    print("Changes made:")
    print("=============")
    for diff in diff_stdout:
        print(diff)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare include and src files from bns_nurates for THC")
    parser.add_argument('source_dir', type=str, help="Full path for the top level bns_nurates directory")
    parser.add_argument('destination_dir', type=str, help="Full path of the destination directory")

    args = parser.parse_args()

    source_dir = args.source_dir
    destination_dir = args.destination_dir

    if not os.path.isdir(source_dir):
        print(f"Error: The source directory {source_dir} does not exist.")
        exit(1)
    
    if not os.path.isdir(destination_dir):
        os.makedirs(destination_dir)
        print(f"Created destination directory: {destination_dir}")

    refactor_and_copy(source_dir, destination_dir)