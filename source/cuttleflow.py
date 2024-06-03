import os
import glob
import json
from infra import Infra
import input
import infra
from DAW import DAW
from to_nextflow import to_nextflow
from annotation import AnnotationDB
import argparse
import subprocess


def main(path_to_workflow, path_to_output):
    current_directory = subprocess.check_output(['pwd'], universal_newlines=True).strip() + "/source/"
    annotation_files = []
    annotation_path = current_directory + 'annotation_files'
    for i in os.listdir(annotation_path):
        if (i.endswith('.json')):
            full_path = '%s/%s' % (annotation_path, i)
            annotation_files.append(full_path)

    annotDB = AnnotationDB(annotation_files)

    with open(path_to_workflow+'/INFRA.json') as jsonfile:
        infra_description = json.load(jsonfile)

    with open(path_to_workflow+'/DAW.json') as jsonfile:
        daw_description = json.load(jsonfile)

    with open(path_to_workflow+'/SPLIT_MERGE_TASKS.json') as jsonfile:
        split_merge_tasks = json.load(jsonfile)

    with open(path_to_workflow+'/INPUT.json') as jsonfile:
        input_description = json.load(jsonfile)

    my_infra = Infra(infra_description)
    my_DAW = DAW(daw_description, input_description, my_infra, annotDB)
    my_DAW = my_DAW.rewrite(annotationdb=annotDB, input_description=input_description, split_merge_annotation=split_merge_tasks)

    to_nextflow(my_DAW, path_to_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--descriptionFolder', required=True, help='Path to the workflow description files')
    parser.add_argument('--output', required=True, help='Path to the output directory where the DAW is written')

    args = parser.parse_args()

    main(args.descriptionFolder, args.output)
