# Change Log
Updates to BAD_Mutations will be documented here. Versioning follows the
form X.Y, where X is the major version, and Y is the minor version. Major
version increments involve alterations to the way the user interacts with the
program, or fundamentally change how the program behaves. Examples are
switching of supported file formats or addition of new subcommands. Minor
version increments are small feature additions or bug fixes. Examples are
modifications to existing subcommands or alterations to output file formats.
Modifications to documentation do not alter version numbers.

## 0.2 - 2016-01-05
### Modified
- Sequences with internal stop codons are skipped from PASTA alignment
- Colons(:) are replaced by underscores(_) in alignment and prediction
- MSA and phylogenetic trees are cleaned for HyPhy input prior to being copied to the output directory

## 0.1 - 2015-12-22
### Added
- User manuals in plain text and GitHub-Flavored Markdown
- Brief section on parallel processing to the manual
- This CHANGELOG to document changes

### Removed
- Old workflows
- Old package source files
