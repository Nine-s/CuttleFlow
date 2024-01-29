from task import Task
from input import Input_of_DAW
from infra import Infra
from Replace import replace_tool
from Compress import compress_before_file_transfer
from Split import split
#import Split
#import Compress

class DAW:  

    tasks:list = []
    tasks_priority:list = []
    is_scatter_gather:bool = False
    input:Input_of_DAW = None
    infra:Infra = None
    wf_level_params:list = []
    
    @staticmethod
    def build_dict_from_dependencies(tasks_list):
        dependencies_dict = {}
        #get names of tasks
        task_name_list = []
        task_modules_name = []
        for task in tasks_list:
            task_name_list.append(task.name)
            task_modules_name.append(task.module_name)
        for module_name in task_modules_name:
            dependencies_dict[module_name] = []
        for task in tasks_list:    
            for input in task.require_input_from:   
                if (".out_channel." in input):
                    #print(input)
                    dependency_name = input.split(".out")[0]
                    #print(dependency_name)
                    #modify dependency name so it is align instead of STAR_ALIGN
                    #XXXX
                    if(dependency_name in task_modules_name):
                        dependencies_dict[dependency_name].append(task.module_name)
                    else:
                        raise Exception("error: " + dependency_name + " is not declared in " + str(task_modules_name))# str(input)) 
        return dependencies_dict

    @staticmethod
    def build_DAG_from_dict(dependencies_dict):
        #https://stackoverflow.com/questions/42195291/topological-sort-kahns-algorithm-trouble
        # Find number of incoming edges for each vertex
        in_degree = {}
        for x, neighbors in dependencies_dict.items():
            in_degree.setdefault(x, 0)
            for n in neighbors:
                in_degree[n] = in_degree.get(n, 0) + 1
        # Iterate over edges to find vertices with no incoming edges
        empty = {v for v, count in in_degree.items() if count == 0}
        result = []
        while empty:
            # Take random vertex from empty set
            v = empty.pop()
            result.append(v)
            # Remove edges originating from it, if vertex not present
            # in adjacency list use empty list as neighbors
            for neighbor in dependencies_dict.get(v, []):
                in_degree[neighbor] -= 1
                # If neighbor has no more incoming edges add it to empty set
                if in_degree[neighbor] == 0:
                    empty.add(neighbor)
        if len(result) != len(in_degree):
            return None # Not DAG
        else:
            return result

    def define_tasks_priority(self):
        tasks_dict = self.build_dict_from_dependencies(self.tasks)
        ordered_tasks_list = self.build_DAG_from_dict(tasks_dict)
        tasks_priority = []
        for mtask in self.tasks:
            task_name = mtask.module_name
            #print(task_name)
            index = ordered_tasks_list.index(task_name)
            tasks_priority.append(index)
        return tasks_priority
        
    def insert_tasks(self, new_tasks):
        self.tasks.append(new_tasks)
        self.tasks_priority = self.define_tasks_priority()
    
    # Creates a DAW object from the description 


    def __init__    (self, DAW_description, input_description, infra, annotDB):
        # create task objects
        self.infra = infra
        self.input = Input_of_DAW(input_description)
        tasks_list = []
        for i in range(len(DAW_description["tasks"])):
            json_task = DAW_description["tasks"][i]
            additional_inputs = {}
            if "channel_operators" in json_task:
                additional_inputs["channel_operators"] = json_task["channel_operators"]
            if "include_from" in json_task:
                additional_inputs["include_from"] = json_task["include_from"]
            new_task = Task(json_task["name"], json_task["toolname"], json_task["inputs"], json_task["outputs"], 
                    json_task["parameters"], json_task["operation"], json_task["module_name"], json_task["module_path"], 
                    input_description, additional_inputs)
            
            tasks_list.append(new_task)
        self.tasks = tasks_list
        # define their priority
        self.tasks_priority = self.define_tasks_priority() 
        #self.rewrite(annotDB, input_description)
    
    def rewrite(self, annotationdb, input_description):
        new_daw = replace_tool(self, annotationdb, input_description, self.input)
        #for task in new_daw.tasks:
         #   print(task.name)
        new_daw = split(self, annotationdb, input_description)
        return new_daw

    def my_print(self):
        print("DAW")
        print(str(self.name))
        print(str(self.tool))
        print(str(self.inputs))
        print(str(self.outputs))
        print(str(self.parameters))
        print(str(self.operation))
        print(self.require_input_from)
