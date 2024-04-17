######################################
# Generate Python bindings with SWIG #
# for the bns_nurates C library      #
######################################

To generate the Python bindings for bns_nurates with SWIG:

1) Launch the setup script by running 'python3 setup.py install' in the current directory.
   This should generate a module that can be imported in a python environment with 'import
   bns_nurates', if the PYHTONPATH is correctly set.

2) Test the bindings by running the 'test_bns_nurates.py' script in the current directory.
