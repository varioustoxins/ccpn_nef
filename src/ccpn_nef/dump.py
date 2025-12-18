#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
"""
NEF file dump utility

This utility uses the NefImporter to parse a NEF file and dumps its contents
as a string displaying the nested dictionary structure using the library's
built-in toString functionality.

Usage:
    python dump.py <nef_file>
    or
    from ccpn_nef.dump import dump_nef_file
    dump_nef_file('path/to/file.nef')
"""

import sys
import os

# Handle both module import and standalone execution
try:
    from . import NefImporter
    from . import ErrorLog as el
except ImportError:
    # If running as standalone script, try direct import
    try:
        import NefImporter
        import ErrorLog as el
    except ImportError:
        # Last resort - add current directory to path and import
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import NefImporter
        import ErrorLog as el


def dump_nef_file(filename, error_logging=el.NEF_STANDARD):
    """
    Parse a NEF file and dump its contents using the library's toString functionality.
    
    :param filename: Path to the NEF file to parse
    :param error_logging: Error logging mode (default: NEF_STANDARD)
    :return: String representation of the NEF file contents
    """
    # Create NefImporter instance
    importer = NefImporter.NefImporter(errorLogging=error_logging)
    
    # Load the NEF file
    try:
        importer.loadFile(filename)
        
        # Use the built-in toString functionality to get string representation
        content = importer.toString()
        
        return content
        
    except Exception as e:
        error_msg = "Error loading NEF file '{}': {}".format(filename, str(e))
        if error_logging != el.NEF_SILENT:
            print(error_msg, file=sys.stderr)
        raise


def dump_nef_text(text, error_logging=el.NEF_STANDARD):
    """
    Parse NEF-formatted text and dump its contents using the library's toString functionality.
    
    :param text: NEF-formatted text string
    :param error_logging: Error logging mode (default: NEF_STANDARD)
    :return: String representation of the NEF contents
    """
    # Create NefImporter instance
    importer = NefImporter.NefImporter(errorLogging=error_logging)
    
    # Load the NEF text
    try:
        importer.loadText(text)
        
        # Use the built-in toString functionality to get string representation
        content = importer.toString()
        
        return content
        
    except Exception as e:
        error_msg = "Error loading NEF text: {}".format(str(e))
        if error_logging != el.NEF_SILENT:
            print(error_msg, file=sys.stderr)
        raise


def main():
    """Command-line interface for the dump utility."""
    if len(sys.argv) != 2:
        print("Usage: python dump.py <nef_file>", file=sys.stderr)
        sys.exit(1)
    
    filename = sys.argv[1]
    
    # Check if file exists
    if not os.path.isfile(filename):
        print("Error: File '{}' not found".format(filename), file=sys.stderr)
        sys.exit(1)
    
    try:
        # Dump the NEF file contents
        content = dump_nef_file(filename)
        print(content)
        
    except Exception as e:
        print("Error: {}".format(str(e)), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()