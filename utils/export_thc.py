"""
Script to perform macro substitutions to the source code for bns_nurates
and export the resulting files for use with Weakrates2 thorn

Usage: python export_thc.py /path/to/bns_nurates /path/to/destination
"""

import os
import argparse
import subprocess

# Replacement rules go here (only change this)
replacement_rules = [("BS_PRINTF", "printf"),("KOKKOS_INLINE_FUNCTION","inline")]

# Text block to replace
text_block = """#ifdef KOKKOS_INLINE_FUNCTION
#include <Kokkos_Core.hpp>
#define BS_PRINTF Kokkos::printf
#else
#define KOKKOS_INLINE_FUNCTION inline
#define BS_PRINTF printf
#endif
"""

# Compute the diff between a file in the source and destination
def diff_files(source_filepath, destination_filepath):

    filediff = subprocess.run(['git', 'diff', '--no-index', source_filepath, destination_filepath], capture_output=True, text=True)

    return filediff.stdout

# Perform text subsitutions on file content
def replace_file_content(file_content, replacements):

    modified_content = file_content
    for find_string, replace_string in replacements:
        modified_content = modified_content.replace(find_string, replace_string)

    return modified_content

# Remove Kokkos block from bns_nurates.hpp
def remove_bns_nurates_text_block(file_content):
    modified_content = file_content.replace(text_block,"")
    
    return modified_content
    
# Perform text substitions on all files in a folder
def refactor_and_copy(source_dir, destination_dir, replacements):
    
    source_include_dir = os.path.join(source_dir, "include")
    if not os.path.isdir(source_include_dir):
        print(f"Source directory {source_include_dir} does not exist!")
        print("Aborting.")
        return False

    try:
        files = os.listdir(source_include_dir)
    except OSError as e:
        print(f"Error getting filenames in {source_include_dir}: {e}")
        print("Aborting.")
        return False

    diff_stdout = []
    for filename in files:
        source_filepath = os.path.join(source_include_dir, filename)
        destination_filepath = os.path.join(destination_dir, filename)

        if not os.path.isfile(source_filepath):
            printf(f"Error {source_filepath} is not a file!")
            printf("Aborting.")
            return False

        with open(source_filepath, 'r') as f:
            content = f.read()
            
            if(filename == "bns_nurates.hpp"):
                content = remove_bns_nurates_text_block(content)
                
            modified_content = replace_file_content(content, replacements)

            with open(destination_filepath, 'w') as out:
                out.write(modified_content)

        diff_stdout.append(diff_files(source_filepath, destination_filepath))

    print("=============")
    print("Changes made:")
    print("=============")
    for diff in diff_stdout:
        print(diff)

    return True
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare include files from bns_nurates for THC installation")
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

    result = refactor_and_copy(source_dir, destination_dir, replacement_rules)
    print(result)
