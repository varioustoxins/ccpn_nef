# NEF File Processing Modules Summary

This directory contains a suite of Python modules for handling, parsing, comparing, and validating NEF (NMR Exchange Format) files. Below is a summary of the main components:

## Main Modules

- **CompareNef**: Routines to compare the contents of two NEF files. See `CompareNef.py` for details.
- **GenericStarParser**: A parser for STAR-type files, used for NEF and related formats. See `GenericStarParser.py` for details.
- **Specification**: Code for handling the NEF specification.
- **StarIo**: Handles input/output for NEF and NmrStar formats. See `StarIo.py` for details.
- **StarTokeniser**: Tokenizer for STAR files.
- **NefImporter**: Routines for reading NEF files and examining their contents.

## Additional Components

- **Validation Dictionaries**: Files such as `mmcif_nef.dic` and `mmcif_nef_v1_1.dic` for NEF file validation.
- **Test Data**: Example NEF files and batch test folders in `testdata/`.
- **Documentation**: `README`, `README.md`, `Overview.md`, and related documents.
- **Testing Scripts**: Located in the `testing/` directory (e.g., `Test_Compare_files.py`).
- **NEF Specification Files**: Found in `NEF/specification/`.
- **Python Environment**: Provided in the `venv/` folder.

## Usage

- Parse NEF files, extract and manipulate their contents using `NefImporter` and related modules.
- Compare NEF files or directories of NEF files using `CompareNef` routines.
- Validate NEF files against specification dictionaries.
- Extend or test functionality using provided scripts and test data.

## Note

All required modules and files appear present for full NEF file processing functionality. If you encounter missing imports or errors, check for the presence of the referenced modules and dictionary files.
