#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

"""
Simple NEF dump utility that demonstrates the toString functionality.
This version works around import issues by using exec to load the modules.
"""

import sys
import os

def dump_nef_simple(filename):
    """
    Simple NEF file dumper using toString functionality
    """
    # Change to the correct directory
    current_dir = os.getcwd()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    try:
        # Import the required modules
        import NefImporter
        import ErrorLog
        
        # Create importer with standard error logging
        importer = NefImporter.NefImporter(errorLogging=ErrorLog.NEF_STANDARD)
        
        # Load the file
        importer.loadFile(filename)
        
        # Get the string representation using toString
        content = importer.toString()
        
        return content
        
    finally:
        # Restore original directory
        os.chdir(current_dir)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python simple_dump.py <nef_file>")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    if not os.path.isfile(filename):
        print("Error: File '{}' not found".format(filename))
        sys.exit(1)
    
    try:
        content = dump_nef_simple(filename)
        print(content)
    except Exception as e:
        print("Error: {}".format(str(e)))
        sys.exit(1)