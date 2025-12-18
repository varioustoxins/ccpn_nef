#!/usr/bin/env python3
"""
Test script for tabulate formatting functionality
"""

from tabulate import tabulate

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


def test_tabulate_formatting():
    """Test the tabulate formatting approach with dummy loop data."""
    
    # Create dummy loop-like data
    columns = ['index', 'chain_code', 'sequence_code', 'residue_name', 'linking', 'residue_variant', 'cis_peptide']
    
    data = [
        {'index': 1, 'chain_code': 'A', 'sequence_code': '13', 'residue_name': 'ALA', 'linking': 'start', 'residue_variant': None, 'cis_peptide': None},
        {'index': 2, 'chain_code': 'A', 'sequence_code': '14', 'residue_name': 'TRP', 'linking': 'middle', 'residue_variant': None, 'cis_peptide': None},
        {'index': 3, 'chain_code': 'A', 'sequence_code': '15', 'residue_name': 'GLY', 'linking': 'middle', 'residue_variant': None, 'cis_peptide': None},
        {'index': 4, 'chain_code': 'A', 'sequence_code': '16', 'residue_name': 'ASN', 'linking': 'middle', 'residue_variant': None, 'cis_peptide': None},
        {'index': 5, 'chain_code': 'A', 'sequence_code': '17', 'residue_name': 'VAL', 'linking': 'middle', 'residue_variant': None, 'cis_peptide': None}
    ]
    
    data_indent = "            "
    param_indent = "          "
    
    print("=== Test 1: Original non-tabulated formatting ===")
    for row in data:
        row_items = []
        for col_key, col_val in row.items():
            row_items.append("'{}': {}".format(col_key, _format_constructor_value(col_val)))
        print("{}{{ {} }}".format(data_indent, ", ".join(row_items)))
    
    print("\n=== Test 2: Tabulated approach - key-value pairs as columns ===")
    
    # Approach 1: Each key-value pair as a column
    table_data = []
    for row in data:
        formatted_row = []
        for col_key in columns:
            col_val = row.get(col_key, None)
            formatted_val = _format_constructor_value(col_val)
            formatted_row.append("'{}': {}".format(col_key, formatted_val))
        table_data.append(formatted_row)
    
    # Use tabulate to align the key-value pairs
    table_str = tabulate(table_data, tablefmt='plain', stralign='left')
    table_lines = table_str.split('\n')
    
    for line in table_lines:
        if line.strip():
            # Split by multiple spaces and rejoin with commas
            parts = [part.strip() for part in line.split('  ') if part.strip()]
            dict_content = ", ".join(parts)
            print("{}{{ {} }}".format(data_indent, dict_content))
    
    print("\n=== Test 3: Simpler approach - align whole dict strings ===")
    
    # Approach 2: Build complete dict strings and align them
    dict_strings = []
    for row in data:
        row_items = []
        for col_key, col_val in row.items():
            row_items.append("'{}': {}".format(col_key, _format_constructor_value(col_val)))
        dict_str = "{{ {} }}".format(", ".join(row_items))
        dict_strings.append([dict_str])
    
    # Use tabulate just for consistent spacing
    table_str = tabulate(dict_strings, tablefmt='plain')
    table_lines = table_str.split('\n')
    
    for line in table_lines:
        if line.strip():
            print("{}{}".format(data_indent, line.strip()))
    
    print("\n=== Test 4: Column-wise alignment ===")
    
    # Approach 3: Align values within each column position
    # Create table where each column is just the value (for alignment)
    value_table = []
    for row in data:
        value_row = []
        for col_key in columns:
            col_val = row.get(col_key, None)
            formatted_val = _format_constructor_value(col_val)
            value_row.append(formatted_val)
        value_table.append(value_row)
    
    # Get aligned values
    aligned_table = tabulate(value_table, tablefmt='plain', stralign='left')
    aligned_lines = aligned_table.split('\n')
    
    # Reconstruct dicts with aligned spacing
    for line in aligned_lines:
        if line.strip():
            values = line.split()
            if len(values) == len(columns):
                dict_pairs = []
                for i, col_name in enumerate(columns):
                    # Pad values to maintain alignment
                    dict_pairs.append("'{}': {}".format(col_name, values[i]))
                print("{}{{ {} }}".format(data_indent, ", ".join(dict_pairs)))


def test_wide_data():
    """Test with wider data that would benefit from tabulation."""
    
    print("\n" + "="*80)
    print("=== Test with Wide Data (many columns) ===")
    
    columns = ['index', 'restraint_id', 'chain_code_1', 'sequence_code_1', 'residue_name_1', 'atom_name_1', 
               'chain_code_2', 'sequence_code_2', 'residue_name_2', 'atom_name_2', 'weight', 'target_value']
    
    data = [
        {'index': 1, 'restraint_id': 1, 'chain_code_1': 'A', 'sequence_code_1': '17', 'residue_name_1': 'VAL', 'atom_name_1': 'H', 'chain_code_2': 'A', 'sequence_code_2': '21', 'residue_name_2': 'ALA', 'atom_name_2': 'HB%', 'weight': 1, 'target_value': 3.7},
        {'index': 2, 'restraint_id': 1, 'chain_code_1': 'A', 'sequence_code_1': '17', 'residue_name_1': 'VAL', 'atom_name_1': 'H', 'chain_code_2': 'A', 'sequence_code_2': '22', 'residue_name_2': 'THR', 'atom_name_2': 'HG2%', 'weight': 1, 'target_value': 3.7},
        {'index': 3, 'restraint_id': 1, 'chain_code_1': 'A', 'sequence_code_1': '18', 'residue_name_1': 'LEU', 'atom_name_1': 'H', 'chain_code_2': 'A', 'sequence_code_2': '21', 'residue_name_2': 'ALA', 'atom_name_2': 'HB%', 'weight': 1, 'target_value': 3.7}
    ]
    
    data_indent = "            "
    
    print("\n--- Without tabulation ---")
    for row in data:
        row_items = []
        for col_key, col_val in row.items():
            row_items.append("'{}': {}".format(col_key, _format_constructor_value(col_val)))
        print("{}{{ {} }}".format(data_indent, ", ".join(row_items)))
    
    print("\n--- With column-wise tabulation ---")
    
    # Create table for alignment
    value_table = []
    max_widths = {}
    
    # First pass: collect all values and calculate max widths
    for row in data:
        for col_key in columns:
            col_val = row.get(col_key, None)
            formatted_val = _format_constructor_value(col_val)
            key_val_str = "'{}': {}".format(col_key, formatted_val)
            max_widths[col_key] = max(max_widths.get(col_key, 0), len(key_val_str))
    
    # Second pass: format with proper spacing
    for row in data:
        dict_parts = []
        for col_key in columns:
            col_val = row.get(col_key, None)
            formatted_val = _format_constructor_value(col_val)
            key_val_str = "'{}': {}".format(col_key, formatted_val)
            # Pad to align nicely
            padded_str = key_val_str.ljust(max_widths[col_key])
            dict_parts.append(padded_str)
        
        print("{}{{ {} }}".format(data_indent, ", ".join(dict_parts)))


if __name__ == "__main__":
    test_tabulate_formatting()
    test_wide_data()