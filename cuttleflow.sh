#!/bin/bash

usage() {
    echo "Usage: $0 --descriptionFolder path_to_workflow --output path_to_output"
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --descriptionFolder) path_to_workflow="$2"; shift ;;
        --output) path_to_output="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

if [ -z "$path_to_workflow" ] || [ -z "$path_to_output" ]; then
    usage
fi

python source/cuttleflow.py --descriptionFolder "$path_to_workflow" --output "$path_to_output"
