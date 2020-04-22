# Contributing to BAD_Mutations
Thank you for your interest in contributing to BAD_Mutations. We welcome suggestions and tips for improvement. You can contribute to the project by either reporting issues or by submitting changes to our code.

## Reporting Issues
Please use the GitHub issues tracker for reporting problems. If you file a new issue, please use one of the templates that we provided. Depending on the nature of the issue, we may ask you to provide a small data set so that we can test our patches.

## Submitting Changes
Please use the GitHub pull request functionality to submit changes to BAD_Mutations. If you do, please follow our stylistic guidelines described below, and include both a description and a rationale for the changes you made.

## Our Stylistic Guidelines
### Python
For Python, we mostly adhere to standard Python style as defined in PEP8. We use four spaces for indentation. Function names are written in `small_case` with underscores between words, class/module names are written in `CamelCase`, and constants are written in `ALL_CAPS`. Lines are wrapped at 80 characters.

### Bash
For bash, we use the `-e`, `-u`, `-o pipefail` options. We double-quote all variables, and place them inside curly braces (e.g., `${VARIABLE}`). Our variable names are `ALL_CAPS`. We do not set a line length limit on bash like we do with Python.
