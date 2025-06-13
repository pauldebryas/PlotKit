#!/usr/bin/env bash

action() {
    local shell_is_zsh="$( [ -z "${ZSH_VERSION}" ] && echo "false" || echo "true" )"
    local this_file="$( ${shell_is_zsh} && echo "${(%):-%x}" || echo "${BASH_SOURCE[0]}" )"
    local this_dir="$( cd "$( dirname "${this_file}" )" && pwd )"
    local parent_dir="$( dirname "${this_dir}" )"  

    export PYTHONPATH="${parent_dir}:${PYTHONPATH}"
    export RUN_PATH="${parent_dir}"

    source /afs/cern.ch/work/p/pdebryas/miniconda3/bin/activate MLenv
}
action
