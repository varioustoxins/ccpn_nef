#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

"""
NEF dictionary dump utility - treats NEF structures as dictionaries

This utility uses the NefImporter to parse a NEF file and dumps its contents
by iterating over the structures as dictionaries, showing the true hierarchical
structure as nested OrderedDicts with Loop objects containing columns and data.
"""

import sys
import os
try:
    from tabulate import tabulate
except ImportError:
    tabulate = None

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


def dump_nef_as_dict(filename, max_depth=None, show_loop_data=True, max_loop_rows=None, error_logging=el.NEF_STANDARD):
    """
    Parse a NEF file and dump its contents as constructor-style representation.
    
    :param filename: Path to the NEF file to parse
    :param max_depth: Maximum depth to traverse (default: None = unlimited)
    :param show_loop_data: Whether to show actual loop data (default: True)
    :param max_loop_rows: Maximum number of loop rows to display (default: None = all rows)
    :param error_logging: Error logging mode (default: NEF_STANDARD)
    :return: String representation of the NEF as constructor-style objects
    """
    # Create NefImporter instance
    importer = NefImporter.NefImporter(errorLogging=error_logging)
    
    try:
        # Load the NEF file
        importer.loadFile(filename)
        
        # Get the internal dictionary structure
        nef_dict = importer._nefDict
        
        # Build constructor-style representation
        output = []
        output.append("# source={}".format(filename))
        output.append("")
        
        # Generate constructor-style output
        output.append("parse = {}".format(_format_as_constructor(nef_dict, show_loop_data, max_loop_rows, max_depth, depth=0)))
        
        return "\n".join(output)
        
    except Exception as e:
        error_msg = "Error loading NEF file '{}': {}".format(filename, str(e))
        if error_logging != el.NEF_SILENT:
            print(error_msg, file=sys.stderr)
        raise


def _format_as_constructor(obj, show_loop_data=True, max_loop_rows=None, max_depth=None, depth=0):
    """Format NEF structures as constructor-style representation."""
    
    if max_depth is not None and depth >= max_depth:
        return "[max depth reached]"
    
    base_indent = "  " * depth
    param_indent = "  " * (depth + 1)
    child_indent = "  " * (depth + 2)
    
    # Check if it's a dictionary-like object
    if hasattr(obj, 'keys'):
        class_name = type(obj).__name__
        
        # Get constructor parameters
        constructor_params = []
        
        # Add name parameter if available
        if hasattr(obj, 'name'):
            constructor_params.append("name='{}'".format(getattr(obj, 'name', '')))
        
        # Build children dictionary
        children_items = []
        for key in obj.keys():
            value = obj[key]
            
            # Check if it's a Loop object
            if hasattr(value, 'columns') and hasattr(value, 'data'):
                if show_loop_data:
                    loop_repr = _format_loop_as_constructor(value, max_loop_rows, depth + 2)
                    children_items.append("{}'{}': {}".format(child_indent, key, loop_repr))
                else:
                    children_items.append("{}'{}': Loop(name='{}', columns={}, rows={})".format(
                        child_indent, key, getattr(value, 'name', ''), len(getattr(value, 'columns', [])), len(getattr(value, 'data', []))))
            
            # Check if it's another dictionary-like (SaveFrame)
            elif hasattr(value, 'keys'):
                child_repr = _format_as_constructor(value, show_loop_data, max_loop_rows, max_depth, depth + 2)
                children_items.append("{}'{}': {}".format(child_indent, key, child_repr))
            
            # Simple value
            else:
                value_repr = _format_constructor_value(value)
                children_items.append("{}'{}': {} # {}".format(child_indent, key, value_repr, type(value).__name__))
        
        # Build constructor call with proper formatting
        if children_items:
            children_str = "{{\n{}\n{}}}".format("\n".join(children_items), param_indent)
            
            if len(constructor_params) > 0:
                param_str = "\n{}{}, \n{}children={}".format(
                    param_indent, constructor_params[0], param_indent, children_str)
                return "{}({}\n{})".format(class_name, param_str, base_indent)
            else:
                return "{}(\n{}children={}\n{})".format(class_name, param_indent, children_str, base_indent)
        else:
            # No children, simple constructor
            if constructor_params:
                return "{}({})".format(class_name, ", ".join(constructor_params))
            else:
                return "{}()".format(class_name)
    
    else:
        # Not dictionary-like
        return "{} # {}".format(_format_constructor_value(obj), type(obj).__name__)


def _format_loop_as_constructor(loop_obj, max_loop_rows=None, depth=0):
    """Format Loop objects as constructor-style representation."""
    
    base_indent = "  " * depth
    param_indent = "  " * (depth + 1)
    data_indent = "  " * (depth + 2)
    
    constructor_params = []
    
    # Add name parameter
    if hasattr(loop_obj, 'name'):
        constructor_params.append("name='{}'".format(getattr(loop_obj, 'name', '')))
    
    # Add columns - format nicely if long
    if hasattr(loop_obj, 'columns'):
        columns_repr = repr(loop_obj.columns)
        if len(columns_repr) > 60:
            # Multi-line columns if too long
            col_items = ["'{}'".format(col) for col in loop_obj.columns]
            columns_str = "(\n{}{}\n{})".format(
                param_indent + "  ", 
                (",\n" + param_indent + "  ").join(col_items),
                param_indent
            )
            constructor_params.append("columns={}".format(columns_str))
        else:
            constructor_params.append("columns={}".format(columns_repr))
    
    # Add data with proper indentation
    if hasattr(loop_obj, 'data'):
        data_len = len(loop_obj.data)
        if data_len > 0:
            rows_to_show = data_len if max_loop_rows is None else min(max_loop_rows, data_len)
            
            # Check if we should use tabulate formatting for large datasets
            use_tabulate = (tabulate is not None and 
                           data_len >= 5 and  # At least 5 rows
                           len(loop_obj.columns) > 6)  # More than 6 columns
            
            
            if use_tabulate and rows_to_show > 0:
                data_str = _format_loop_data_as_table(loop_obj, rows_to_show, max_loop_rows, data_len, param_indent, data_indent)
            
            else:
                # Use original dict-style formatting for smaller datasets
                data_items = []
                for i in range(rows_to_show):
                    row = loop_obj.data[i]
                    if hasattr(row, 'keys'):  # OrderedDict-like row
                        row_items = []
                        for col_key, col_val in row.items():
                            row_items.append("'{}': {}".format(col_key, _format_constructor_value(col_val)))
                        data_items.append("{}{{{} }}".format(data_indent, ", ".join(row_items)))
                    else:
                        data_items.append("{}{}".format(data_indent, _format_constructor_value(row)))
                
                if max_loop_rows is not None and data_len > max_loop_rows:
                    data_items.append("{}# ... ({} more rows)".format(data_indent, data_len - max_loop_rows))
                
                data_str = "[\n{}\n{}]".format("\n".join(data_items), param_indent)
            
            constructor_params.append("data={}".format(data_str))
        else:
            constructor_params.append("data=[]")
    
    # Format constructor with proper multi-line indentation
    if len(constructor_params) > 1:
        param_lines = []
        for param in constructor_params:
            param_lines.append("{}{}".format(param_indent, param))
        
        return "Loop(\n{}\n{})".format(",\n".join(param_lines), base_indent)
    else:
        return "Loop({})".format(", ".join(constructor_params))


def _format_loop_data_as_table(loop_obj, rows_to_show, max_loop_rows, data_len, param_indent, data_indent):
    """Format loop data with column alignment while maintaining valid Python syntax."""
    
    # First pass: collect all key-value strings and calculate max widths
    max_widths = {}
    all_rows_data = []
    
    for i in range(rows_to_show):
        row = loop_obj.data[i]
        if hasattr(row, 'keys'):  # OrderedDict-like row
            row_data = {}
            for col_key in loop_obj.columns:
                col_val = row.get(col_key, None)
                formatted_val = _format_constructor_value(col_val)
                key_val_str = "'{}': {}".format(col_key, formatted_val)
                row_data[col_key] = key_val_str
                max_widths[col_key] = max(max_widths.get(col_key, 0), len(key_val_str))
            all_rows_data.append(row_data)
    
    # Second pass: format with proper alignment
    data_items = []
    for row_data in all_rows_data:
        dict_parts = []
        for col_key in loop_obj.columns:
            key_val_str = row_data.get(col_key, '')
            # Pad to align nicely (except for the last column)
            if col_key == loop_obj.columns[-1]:
                # Don't pad the last column
                dict_parts.append(key_val_str)
            else:
                padded_str = key_val_str.ljust(max_widths[col_key])
                dict_parts.append(padded_str)
        
        data_items.append("{}{{ {} }}".format(data_indent, ", ".join(dict_parts)))
    
    if max_loop_rows is not None and data_len > max_loop_rows:
        data_items.append("{}# ... ({} more rows)".format(data_indent, data_len - max_loop_rows))
    
    return "[\n{}\n{}]".format("\n".join(data_items), param_indent)


def _format_constructor_value(value):
    """Format a value for constructor-style display."""
    if value is None:
        return "None"
    elif isinstance(value, str):
        if len(value) > 80:
            return "'{}'...".format(value[:77])
        else:
            return "'{}'".format(value)
    else:
        str_val = str(value)
        if len(str_val) > 80:
            return "{}...".format(str_val[:77])
        return str_val




def main():
    """Command-line interface for the dictionary dump utility."""
    if len(sys.argv) < 2:
        print("Usage: python dump_as_dict.py <nef_file> [max_depth] [show_loop_data] [max_loop_rows]", file=sys.stderr)
        print("  max_depth: Maximum traversal depth (default: unlimited)", file=sys.stderr)
        print("  show_loop_data: Show loop data - true/false (default: true)", file=sys.stderr)
        print("  max_loop_rows: Max loop rows to show (default: all rows)", file=sys.stderr)
        print("", file=sys.stderr)
        print("Examples:", file=sys.stderr)
        print("  python dump_as_dict.py file.nef              # Show everything", file=sys.stderr)
        print("  python dump_as_dict.py file.nef 3            # Limit depth to 3", file=sys.stderr)
        print("  python dump_as_dict.py file.nef 3 false      # Limit depth, hide loop data", file=sys.stderr)
        print("  python dump_as_dict.py file.nef 3 true 10    # Limit depth, show max 10 rows", file=sys.stderr)
        sys.exit(1)
    
    filename = sys.argv[1]
    max_depth = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != 'none' else None
    show_loop_data = sys.argv[3].lower() != 'false' if len(sys.argv) > 3 else True
    max_loop_rows = int(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] != 'none' else None
    
    # Check if file exists
    if not os.path.isfile(filename):
        print("Error: File '{}' not found".format(filename), file=sys.stderr)
        sys.exit(1)
    
    try:
        # Dump the NEF as dictionaries
        content = dump_nef_as_dict(filename, max_depth=max_depth, 
                                  show_loop_data=show_loop_data, 
                                  max_loop_rows=max_loop_rows)
        print(content)
        
    except Exception as e:
        print("Error: {}".format(str(e)), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()