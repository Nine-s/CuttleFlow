from task import Task
import numpy as np
import math
import re
import pandas as pd
import statistics
from sklearn.preprocessing import PolynomialFeatures

class ToolException(Exception):
    pass


#returns a degree of parallelism from which splitting is reducing runtime
def min_beneficial_split(alignment_task, annotation_database, median_input_size, ram, cpu):
    align_time_no_split = annotation_database.predict_runtime(alignment_task.tool, ram, cpu, median_input_size)
    split_time = annotation_database.predict_runtime("split-merge", ram, cpu, median_input_size)
    
    #max alignment runtime after splitting to reduce runtime overall
    max_align_time_split = align_time_no_split - split_time
    
    #search for split param from where predicted alignment time is beneficial considering the newly added split time
    for i in range(1,100):
        input_size_chunked = median_input_size/i
        align_time_chunked = annotation_database.predict_runtime(alignment_task.tool, ram, cpu, input_size_chunked)
        if(align_time_chunked < max_align_time_split):
            return i 
            
    return 0
    
    
    #return ReLU of min_split, when min_split is negative, splitting increases runtime


def find_split_input(DAW, channeled_inputs, split_input_type):
    for c_input in channeled_inputs:
        finished = False
        no_parents = False
        task = next(task for task in DAW.tasks if c_input.split(".")[0]==task.module_name)
        print(task.tool)
        print(task.operation)
        print(task.inputs)
        print(task.input_description)
        for i in task.inputs:
            print(i)
            print(i.input_type)
            print(i.paths)
            if i.input_type==split_input_type:
                return c_input
            child_task_requirements = task.require_input_from
            #return c_input
        while (finished == False) & (no_parents == False):
            parent_tasks = [task for task in DAW.tasks if task.module_name in [req.split(".")[0] for req in child_task_requirements]]
            print(parent_tasks)
            print("F")
            if parent_tasks == []:
                no_parents = True
                #check for all parent tasks whether one of the has samples as input
            for index, p in enumerate(parent_tasks):     
                for i in p.inputs:
                    if i.input_type==split_input_type:
                        return child_tasks_requirements[index]
                #otherwise: set child = parents to prepare next loop iteration
            child_tasks_requirements = []
            for p in parent_tasks:
                child_tasks_requirements.extend(p.require_input_from)
          
"""
input parameters should be:
x task operation to be split 
x use annotation to predict runtime yes/no?
--> paths to merge and split tasks (if needed)
"""

def split(DAW, annotation_database, input_description, split_operation, predict_runtime, input_type_split, i_split_task, i_merge_task=None):
    
    if input_type_split == "sample":
        if ( len(DAW.infra.list_nodes) < len(DAW.input.input_samples) ):
            print(str(len(DAW.infra.list_nodes)) + " < " + str(len(DAW.input.input_samples)))
            print('Less nodes than input samples available')
            return DAW
        
    to_split_tasks = [task for task in DAW.tasks if task.operation == split_operation]
    if to_split_tasks == []:
        print('No task with operation "' + split_operation + '" was found.')
        return DAW
 

    for to_split_task in to_split_tasks:
        try:
            #print(alignment_task.module_name)
            try:
                annotation_split = next(annot for annot in annotation_database.annotation_db if annot.toolname.lower() == to_split_task.tool.lower())
            except StopIteration:
                raise ToolException("No annotation for tool " + str(to_split_task.tool) + " was found.")
            
            if annotation_split.is_splittable == "False": #if tool does not support splitting, return DAW (no changes)
                raise ToolException("Tool " + str(to_split_task.tool) + " is not splittable according to the annotation database.")
            if predict_runtime == "True":
                median_input_size = statistics.median(DAW.input.size_of_samples)
                ram = DAW.infra.RAM
                cpus = [node.cpu for node in DAW.infra.list_nodes]
                for i, element in enumerate(cpus):
                    cpus[i] = int(element.replace("m",""))
                cpu = statistics.median(cpus)
            
                min_split = min_beneficial_split(to_split_task, annotation_database, median_input_size, ram, cpu)
                if((DAW.infra.number_nodes < min_split) | (min_split < 1)): #splitting not beneficial
                    raise ToolException("Not enough nodes available or splitting not beneficial at all for " + str(to_split_task.tool) + ".")
                else:
                    split_number = DAW.infra.number_nodes
                    DAW.wf_level_params.append(("split", split_number))
                    #cores = [int(node.cores) for node in DAW.infra.list_nodes] 
            
            #find splittable tasks, first is align
            first_split_task = to_split_task
            last_split_task = to_split_task
            task_splittable = True
        
            while task_splittable == True:
                output_last_split_task = re.compile(last_split_task.module_name + ".out_channel.*")
                next_tasks = [task for task in DAW.tasks if [requirement for requirement in task.require_input_from if output_last_split_task.match(requirement)] != []]
                if next_tasks != []:    
                    for task in next_tasks:
                        print(task.tool)
                        annotation_next_task = [annotation for annotation in annotation_database.annotation_db if annotation.toolname == task.tool]
                        if annotation_next_task != []:
                            if annotation_next_task[0].is_splittable == "True":
                                last_split_task = task 
                                task_splittable = True
                                continue              
                            else:
                                task_splittable = False 
                        else: 
                            task_splittable = False
                else: #no next task found
                    task_splittable = False       

            output_last_split_task = last_split_task.module_name + ".out_channel." + last_split_task.outputs[0]
            child_tasks = [task for task in DAW.tasks if output_last_split_task in task.require_input_from]
            
            try:
                input_to_split = next(input for input in first_split_task.inputs if input.input_type == input_type_split)
                input_to_split_channel = first_split_task.inputs_from_DAW[first_split_task.inputs_from_DAW.index(input_to_split.name)]
                
            except StopIteration:
            #samples not direct input of align task: search along DAG to find task that provides preprocessed input
                channeled_inputs = first_split_task.require_input_from
                print(channeled_inputs)
                split_input_channel = find_split_input(DAW, channeled_inputs, input_type_split)
                input_to_split = split_input_channel
                input_to_split_channel = first_split_task.inputs_from_DAW[first_split_task.inputs_from_DAW.index(input_to_split)]
            full_module_path = first_split_task.module_path
        
            last_slash_index = full_module_path.rfind("/")
            module_path = full_module_path[:last_slash_index]


            split_task = Task("split", "fastqsplit", [input_to_split_channel], ["split_reads"], [], "split", ("FASTQSPLIT_" + to_split_task.tool.upper()),  module_path + "/FASTQSPLIT.nf", input_description, {"include_from": "FASTQSPLIT"}) 
            split_task_output = split_task.module_name + ".out_channel." + split_task.outputs[0]
            first_split_task.change_input(split_task_output, input_to_split)
            
            add_s_inputs = {}
            if "channel_operators" in i_split_task:
                add_s_inputs["channel_operators"] = i_split_task["channel_operators"]
            print(add_s_inputs)
            print(i_split_task)
            add_s_inputs["include_from"] = i_split_task["module_name"]
            split_task = Task(i_split_task["name"], i_split_task["tool"], [input_to_split], i_split_task["outputs"],\
                             i_split_task["parameters"], i_split_task["operation"], \
                            (i_split_task["module_name"] + "_" + to_split_task.tool.upper()), i_split_task["module_path"], \
                            input_description, add_s_inputs)
            print(split_task)
            for task in DAW.tasks:
                print(task.name)
            DAW.insert_tasks(split_task) 
            print("###########################")
            print(to_split_task.tool.upper())
            print("###########################")
            if(i_merge_task!=None):  
                add_m_inputs = {}
                if "channel_operators" in i_merge_task:
                    add_m_inputs["channel_operators"] = i_merge_task["channel_operators"]
                add_m_inputs["include_from"] = i_merge_task["module_name"]
                merge_task = Task(i_merge_task["name"], i_merge_task["tool"], [output_last_split_task], i_merge_task["outputs"],\
                             i_merge_task["parameters"], i_merge_task["operation"], \
                            (i_merge_task["module_name"] + "_" + to_split_task.tool.upper()), i_merge_task["module_path"], \
                            input_description, add_m_inputs)
                merge_task_output = merge_task.module_name + ".out_channel." + merge_task.outputs[0]
                for child_task in child_tasks:
                    child_task.change_input(merge_task_output, output_last_split_task)
                DAW.insert_tasks(merge_task)
        except ToolException:
            continue
        return DAW
        
    
