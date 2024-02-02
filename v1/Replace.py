from task import Task
import pandas as pd
from sklearn.preprocessing import PolynomialFeatures

def find_alternative_tool(annotation_db, tool_to_replace):
    alt_tool_list = []
    # print(annotation_db)
    # print(tool_to_replace)
    for tool in annotation_db:
        #find tool with matching operations
        #print(tool.operation)
        #print(tool_to_replace.operation)
        if (tool.operation == tool_to_replace.operation):
        #    print("OK")
            #find tool that runs with the same input/out as original 
            if ( input_output_matches(tool, tool_to_replace) ):
                alt_tool_list.append(tool)
        # else:
        #     print("NOT OK")
    return alt_tool_list

def find_index_tool(annotationdb, final_tool):
    for tool in annotationdb:
        #find tool with matching operations
        # print(tool.operation)
        # print(tool.toolname)
        if (tool.operation == "index"):
            if ( final_tool.toolname.casefold() == tool.toolname.casefold() ):
                return tool
        #print("ERROR: no tool")
            



def input_output_matches(tool, tool_to_replace):
    
    input_to_replace = tool_to_replace.mendatory_input_list
    input_tool_from_database = tool.mendatory_input_list

    output_to_replace = tool_to_replace.output_list
    output_tool_from_database = tool.output_list
    
    if ( set(input_to_replace).issubset(set(input_tool_from_database)) ):
        if ( set(output_to_replace).issubset(set(output_tool_from_database)) ):
            return True
    else:
        return False


def is_tool_runnable( tool, RAM, reference_size, model ):
    min_RAM = model.predict([[reference_size]])
    if ( RAM > min_RAM[0] + 1 ):
        return True
    else:
        return False

def choose_best_tool(list_alt_tools, annot, input_of_daw):
    dataset_size = 18.5
    ram = 256
    cpu = 2400.0000
    #TODO: those values should not be fixed, right?
    list_predicted_runtimes = []

    for tool in list_alt_tools:
        try:
            list_predicted_runtimes.append(annot.predict_runtime(tool.toolname, ram, cpu, dataset_size))
        except ValueError:
            print("No runtime model found for " + tool.toolname + " on most similar cluster.")
            list_predicted_runtimes.append(float("inf"))
    
    min_number = min(list_predicted_runtimes)
    min_index = list_predicted_runtimes.index(min_number)
    best_tool = list_alt_tools[min_index]
    return best_tool

def create_new_task(task, new_tool, input_description, operation_):
    name = new_tool.toolname + "_" + operation_
    tool = new_tool.toolname
    outputs = task.outputs
    parameters = [] #TODO
    module_path = new_tool.module_path
    module_name = new_tool.module_name#(module_path.split("/")[-1].split(".")[0]).upper()
    inputs_from_DAW = task.inputs_from_DAW
    require_input_from_list = task.require_input_from #TODO

    new_task = Task(name=name, tool=tool, inputs_from_DAW=inputs_from_DAW, outputs=outputs, parameters=parameters, operation=operation_, module_name=module_name, module_path=module_path, input_description=input_description)
    #print(self.tasks)
    #new_task.my_print()
    return new_task

def task_needing_input_change(task, old_task):
    print(task.name)
    old_input = old_task.module_name
    is_old_input_used_by_task = False
    for require in task.require_input_from: 
        if (old_input in require):
            #print(old_input)
            #print(require)
            is_old_input_used_by_task = True
    #print(old_input)
    #print(task.require_input_from)
    #print([mtask for mtask in task.require_input_from])
    if ( is_old_input_used_by_task ):
        #print("##################################### TRUE")
        return True
    # else:
    #     return False
    
    # print("$$$$$$$")
    # print(task.name)
    # print("$$$$$$$")
    
def change_inputs_task(task, old_task, new_task):
    for i in range(len( task.require_input_from)):
        if(old_task.module_name in task.require_input_from[i]):
            task.require_input_from[i] = task.require_input_from[i].replace(old_task.module_name, new_task.module_name) 
    for i in range(len( task.inputs_task )):
        if(old_task.module_name in task.inputs_task[i]):
            task.inputs_task[i] = task.inputs_task[i].replace(old_task.module_name, new_task.module_name)
    #i = task.inputs.index(old_task.outputs)
    #task.inputs[i:i + 1] = new_task.output_list
    return task

def replace_tool(daw, annotations, input_description, input_of_daw):
    for i in range(len(daw.tasks)):
        ### check if task is "align"
        replaced = False
        task = daw.tasks[i]
        if (task.operation == "align"):
            # print("++++++")
            # print(task.tool)
            annotation_tools_list = [tool for tool in annotations.annotation_db]# if tool.operation=="align"] 
            
            # print("$$$$$$$$")
            # print(annotation_tools_list)
            # print("$$$$$$$$$$")

            # find the tool of the DAW task in the annot DB
            tool_to_replace = next((tool_in_annot for tool_in_annot in annotation_tools_list if (tool_in_annot.toolname.casefold() == task.tool.casefold() and tool_in_annot.operation == "align")), None)
            if(tool_to_replace == None):
                print("Tool to replace not found in DB")
                continue

            # find the tools in the db that match the requirements
            alternative_tools_list = find_alternative_tool(annotations.annotation_db, tool_to_replace)
            
            if ( len(alternative_tools_list) < 1 ):
                print("No alternative found")
                continue

            else:    
                RAM = daw.infra.RAM 
                reference_size = input_of_daw.size_of_reference_genome_max
                # print("*******")
                # print(alternative_tools_list)
                # alternative_tools_list = [tool_alt for tool_alt in alternative_tools_list if (is_tool_runnable(tool_alt, RAM, reference_size, tool_alt.RAM_requirements_model))] 
                # print(alternative_tools_list)
                # print("*******")
                if (len(alternative_tools_list) > 1):
                    replaced = True
                    final_tool = choose_best_tool(alternative_tools_list, annotations, input_of_daw) #TODO
                    # print("######")
                    # print(final_tool.toolname)
                    # print("######")
                    new_align_task = create_new_task(task, final_tool, input_description, "align")
                    daw.tasks[i] = new_align_task
                    #new_align_task.my_print()
                    #daw.tasks[i-2].my_print()
                    for z in range(len(daw.tasks)):
                        if (task_needing_input_change(daw.tasks[z], old_task=task)):
                            daw.tasks[z] = change_inputs_task(daw.tasks[z], old_task=task, new_task=final_tool)
                            #daw.tasks[z].my_print()
                    for j in range(len(daw.tasks)):
                        if (tool_to_replace.toolname.casefold() in daw.tasks[j].tool.casefold()):
                            if (daw.tasks[j].operation == "index"):
                                old_task_index = daw.tasks[j]
                                index_tool = find_index_tool(annotations.annotation_db, final_tool)
                                new_task = create_new_task(old_task_index, index_tool, input_description, "index")
                                daw.tasks[j] = new_task
                                # print("$$$$$$$$$$$$$$$$$$$$$")
                                # new_task.my_print()
                                # print("$$$$$$$$$$$$$$$$$$$$$")
                                ##old_task_index.my_print()
                                for y in range(len(daw.tasks)):
                                    if(task_needing_input_change(daw.tasks[y], old_task=old_task_index)):
                                        daw.tasks[y] = change_inputs_task(daw.tasks[y], old_task_index, new_task=index_tool)
                                        # print("##### OK ######")
                                        # daw.tasks[y].my_print()
                                        # print("##### OK ######")
                            #daw = change_inputs_DAW(daw, old_task=old_task_index, new_task=index_tool)
                else:
                    continue
    # for z in range(len(daw.tasks)):
    #     daw.tasks[z].my_print()
    #print(daw)
    print("##### :) ######")
    return daw