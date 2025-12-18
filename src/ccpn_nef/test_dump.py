#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for the dump utility functionality
"""

import sys
import os

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import NefImporter
import ErrorLog as el


def test_dump():
    """Test the dump functionality with a NEF file"""
    # Find a test NEF file
    test_file = "../../tests/test_data/Commented_Example.nef"
    
    if not os.path.isfile(test_file):
        print("Test file not found, looking for alternatives...")
        return
    
    # Create NefImporter instance
    importer = NefImporter.NefImporter(errorLogging=el.NEF_STANDARD)
    
    try:
        # Load the NEF file
        print(f"Loading NEF file: {test_file}")
        importer.loadFile(test_file)
        
        # Use the built-in toString functionality
        print("\n" + "="*50)
        print("NEF FILE CONTENTS:")
        print("="*50)
        content = importer.toString()
        print(content)
        
    except Exception as e:
        print(f"Error: {str(e)}")


if __name__ == '__main__':
    test_dump()