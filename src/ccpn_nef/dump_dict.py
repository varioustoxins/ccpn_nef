#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

"""
NEF dictionary structure dump utility

This utility uses the NefImporter to parse a NEF file and dumps its contents
as a nested dictionary structure, showing the actual Python object hierarchy
rather than the NEF format output.
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


def dump_nef_structure(filename, max_depth=3, show_data=True, error_logging=el.NEF_STANDARD):
    """
    Parse a NEF file and dump its nested dictionary structure.
    
    :param filename: Path to the NEF file to parse
    :param max_depth: Maximum depth to traverse (default: 3)
    :param show_data: Whether to show actual data values (default: True)
    :param error_logging: Error logging mode (default: NEF_STANDARD)
    :return: String representation of the NEF dictionary structure
    """
    # Create NefImporter instance
    importer = NefImporter.NefImporter(errorLogging=error_logging)
    
    try:
        # Load the NEF file
        importer.loadFile(filename)
        
        # Get the internal dictionary structure
        nef_dict = importer._nefDict
        
        # Build structure representation
        output = []
        output.append("NEF Dictionary Structure")
        output.append("=" * 50)
        output.append("DataBlock: {} (type: {})".format(nef_dict.name, type(nef_dict).__name__))
        output.append("")
        
        # Traverse the structure
        _traverse_structure(nef_dict, output, indent="", depth=0, max_depth=max_depth, show_data=show_data)
        
        return "\n".join(output)
        
    except Exception as e:
        error_msg = "Error loading NEF file '{}': {}".format(filename, str(e))
        if error_logging != el.NEF_SILENT:
            print(error_msg, file=sys.stderr)
        raise


def _traverse_structure(obj, output, indent="", depth=0, max_depth=3, show_data=True):
    """Recursively traverse and display the NEF structure."""
    
    if depth >= max_depth:
        output.append("{}[max depth reached]".format(indent))
        return
    
    if hasattr(obj, 'keys'):  # Dict-like object
        for key in obj.keys():
            value = obj[key]
            type_name = type(value).__name__
            
            if hasattr(value, 'keys'):  # Nested dict-like
                output.append("{}{}: {} ({})".format(indent, key, getattr(value, 'name', ''), type_name))
                _traverse_structure(value, output, indent + "  ", depth + 1, max_depth, show_data)
            elif hasattr(value, '__len__') and not isinstance(value, str):  # List-like
                output.append("{}{}: [{}] ({})".format(indent, key, len(value), type_name))
                if show_data and len(value) > 0 and depth < max_depth - 1:
                    # Show first few items
                    for i, item in enumerate(value[:3]):
                        output.append("{}  [{}]: {} ({})".format(indent, i, _format_value(item), type(item).__name__))
                    if len(value) > 3:
                        output.append("{}  ... and {} more items".format(indent, len(value) - 3))
            else:  # Simple value
                if show_data:
                    formatted_value = _format_value(value)
                    output.append("{}{}: {} ({})".format(indent, key, formatted_value, type_name))
                else:
                    output.append("{}{}: ({})".format(indent, key, type_name))


def _format_value(value):
    """Format a value for display, truncating if necessary."""
    str_val = str(value)
    if len(str_val) > 100:
        return str_val[:100] + "..."
    return str_val


def main():
    """Command-line interface for the dictionary dump utility."""
    if len(sys.argv) < 2:
        print("Usage: python dump_dict.py <nef_file> [max_depth] [show_data]", file=sys.stderr)
        print("  max_depth: Maximum traversal depth (default: 3)", file=sys.stderr)
        print("  show_data: Show data values - true/false (default: true)", file=sys.stderr)
        sys.exit(1)
    
    filename = sys.argv[1]
    max_depth = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    show_data = sys.argv[3].lower() != 'false' if len(sys.argv) > 3 else True
    
    # Check if file exists
    if not os.path.isfile(filename):
        print("Error: File '{}' not found".format(filename), file=sys.stderr)
        sys.exit(1)
    
    try:
        # Dump the NEF dictionary structure
        content = dump_nef_structure(filename, max_depth=max_depth, show_data=show_data)
        print(content)
        
    except Exception as e:
        print("Error: {}".format(str(e)), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()