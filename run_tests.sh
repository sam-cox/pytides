#!/bin/bash -e

function run_pep8_style_checks {
    flake8 .
}


function run_python_unit_tests {
    nosetests
}

# run_pep8_style_checks  # disabled until PEP8 changes are merged
run_python_unit_tests
