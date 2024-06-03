import pandas as pd
from itertools import chain
import os

class to_nextflow:
    DAW = None
    PATH_TO_DAW = "generated_workflow/"

    def write_workflow(self, input_tasks_list):
        start = """\nworkflow{
        read_pairs_ch = Channel
            .fromPath( params.csv_input )
            .splitCsv(header: true, sep: ',')
            .map {row -> tuple(row.sample, [row.path_r1, row.path_r2], row.condition)}
            .view()
        """
        end = "\n}"
        core = self.write_core_workflow(input_tasks_list)
        return start + core + end

    def write_core_workflow(self, input_tasks_list):
        core = "\n"
        tasks = self.DAW.tasks
        for i in range(len(tasks)):
            priority_index = self.DAW.tasks_priority.index(i)
            tmp_task = tasks[priority_index]
            tmp = str(tmp_task.module_name) + "(" 
            for index, my_input in enumerate(tmp_task.inputs_task):
                if (my_input == "samples"):
                        tmp += "read_pairs_ch"            
                elif((".out" not in my_input) and (my_input != "reads")):
                    tmp += "params." + my_input
                    
                else: 
                    tmp_input = my_input
                    tmp_input = tmp_input.replace("out_channel", "out")
                    tmp += tmp_input
                if tmp_task.channel_operators != None:
                    tmp += tmp_task.channel_operators[index]
                tmp += ", "
            tmp = tmp[:-2] + ")\n"
            core += tmp
        return core
    
    def generate_include_modules(self):
        include_string = ""
        include_dictionnary = {}
        for task in self.DAW.tasks:
            if (task.module_path in include_dictionnary):
                include_dictionnary[task.module_path].append((task.module_name, task.include_from))
            else:
                include_dictionnary[task.module_path] = []
                include_dictionnary[task.module_path].append((task.module_name,task.include_from))
        for path in include_dictionnary:
            module_names_string = ""
            for module_name in include_dictionnary[path]:
                if(module_name[1] == ""):
                    module_names_string += module_name[0] + " ; "
                else: 
                    module_names_string += module_name[1] + " as " + module_name[0] + "; "
            include = "include { " + module_names_string[:-2] + " } from '" + path + "'\n"
            include_string = include_string + include
        return include_string    


    def create_config_file(self):
        with open("./resources/default_nextflow_config.config") as f:
            base_config = f.read()

        params_string = "\n\nparams {\n"
        params_string += "\tstrand = " + "'" + self.DAW.input.first_strand + "'\n"
        params_string += "\toutdir = 'results'\n"
        params_string += "\tcsv_input = './input.csv'\n"
        for ref in self.DAW.input.input_references:
            params_string += "\t"+ str(ref.name) + " = '"+ str(ref.paths[0]) + "'\n"

        for i in range(len(self.DAW.tasks)):
            task_inputs = self.DAW.tasks[i].inputs
            for j in range(len(task_inputs)):
                _input = task_inputs[j]
        for additional_param in self.DAW.wf_level_params:
            param, value = additional_param
            params_string += "\t" + param + " = " + str(value) + "\n"
        for input_param in self.DAW.input.input_parameters:
            value = self.DAW.input.input_parameters[input_param]
            params_string += "\t" + input_param + " = " + str(value) + "\n"
        params_string += "\tthreads = " + str(self.DAW.infra.threads) + "\n"
        params_string += "\tbasedir = '" + os.popen('pwd').read().strip() + "/generated_workflow'\n"
        params_string += "}\n"
        with open(self.PATH_TO_DAW +'/nextflow.config', 'w') as config_file:
            config_file.write(base_config + params_string)
            config_file.close()
    
    def write_input_csv(self, DAW):
        mandatory_columns = ['sample', 'path_r1', 'path_r2']
        additional_columns = list(set(chain.from_iterable(list(sample.additional_columns.keys()) for sample in DAW.input.input_samples)))
        columns = mandatory_columns + additional_columns
        input_tasks_list = []
        df = pd.DataFrame(columns=columns)
        for i in range(len(DAW.input.input_samples)):
            inputDAW = DAW.input.input_samples[i]
            mandatory_values = [inputDAW.name, inputDAW.paths[0], inputDAW.paths[1]]
            additional_values = []
            for additional_element in additional_columns:
                if(additional_element in inputDAW.additional_columns):
                    additional_values.append(inputDAW.additional_columns[additional_element])
                else:
                    additional_values.append("")
            df.loc[i] = mandatory_values + additional_values
            input_tasks_list.append(inputDAW.name)
        df.to_csv(self.PATH_TO_DAW +"/input.csv", index=False)
        return input_tasks_list

    
    def __init__(self, DAW,PATH_TO_DAW):
        self.DAW = DAW
        self.PATH_TO_DAW = PATH_TO_DAW
        os.makedirs(PATH_TO_DAW, exist_ok=True)
        self.create_config_file()
        input_tasks_list = self.write_input_csv(DAW)
        nextflow_header = "nextflow.enable.dsl = 2\n"
        include = self.generate_include_modules()
        header = nextflow_header + "\n\n" + include
        workflow = self.write_workflow(input_tasks_list)
        with open(self.PATH_TO_DAW +"/main.nf", "w") as f:
            f.write("")
            f.write(header + workflow)
            f.close()
        print("Workflow written in the folder: "+ self.PATH_TO_DAW)

