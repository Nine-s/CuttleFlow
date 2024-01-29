from input import Input

class Task:  

    name:str = "task_default"
    tool:str = "none"
    operation:str = "none"
    inputs:list = []
    outputs:list = []
    parameters:list = []
    require_input_from:list = [] 
    module_name:str = ""
    module_path:str = ""
    inputs_from_DAW = ""
    channel_operators = None
    include_from:str = ""

    def create_input(self, inputs_from_DAW, input_description):
        list_inputs = []
        for input in inputs_from_DAW:
            is_input_described = False
            if(".out" in input):
                continue #TODO: check if the name before "out" exists 
            for sample in input_description["samples"]:
                # if (sample["name"] == input):
                mInput = Input(input, "sample", [sample["path_r1"], sample["path_r2"]], sample["strand"], "", sample["uncompressed_size"])
                list_inputs.append(mInput)
                is_input_described = True
            if(is_input_described):
                continue
            for reference in input_description["references"]:
                if (reference["name"] == input):
                    mInput = Input(input, "reference", [reference["path"]], "", reference["reference_type"], reference["uncompressed_size"])
                    list_inputs.append(mInput)   
                    is_input_described = True
            for parameter in input_description["parameters"]:
                if(parameter["name"] == input): 
                    is_input_described = True
            if(input == "csv_input"):
                is_input_described = True
            if not is_input_described:
                    raise Exception( 'The input "' + input + '" of the task '+ self.name + ' was not described in the input description file')
        return list_inputs

    def change_input(self, new_input, old_input):
        if isinstance(old_input, Input): 
            old_input_position = self.inputs.index(old_input)
            self.inputs.pop(old_input_position)
        elif isinstance(old_input, str) and ".out" in old_input:
            self.require_input_from.remove(old_input)
            old_input_position = self.inputs_task.index(old_input)
        self.inputs_task[old_input_position] = new_input
        if isinstance(new_input, str) and ".out" in new_input and new_input not in self.require_input_from:
            self.require_input_from.append(new_input)
    
    def __init__(self, name, tool, inputs_from_DAW, outputs, parameters, operation, module_name, module_path, input_description, additional_inputs={}):#, require_input_from): 
        self.name:str = name
        self.tool:str = tool
        self.outputs:list = outputs
        self.parameters:list = parameters
        self.operation:str = operation
        self.module_name:str = module_name
        self.module_path:str = module_path
        self.inputs_from_DAW = inputs_from_DAW
        self.input_description = input_description
        if "channel_operators" in additional_inputs:
            self.channel_operators = additional_inputs["channel_operators"]            
        if "include_from" in additional_inputs:
            self.include_from = additional_inputs["include_from"]
        inputs = self.create_input(inputs_from_DAW, input_description)
        self.inputs = inputs
        self.inputs_task:list = inputs_from_DAW
        require_input_from_list = []
        for input in self.inputs_task:        
            if (".out_channel." in input):
                require_input_from_list.append(input)
        self.require_input_from:list = require_input_from_list #build_task_dependencies(require_input_from)

    def my_print(self):
        print("\n#+~ TASK #+~ ")
        print("name:    " + str(self.name))
        print("tool:    " + str(self.tool))
        print("inputs:  " + str([minputs.name for minputs in self.inputs]))
        print("outputs: " + str(self.outputs))
        print("params:  " + str(self.parameters))
        print("operation:   " + str(self.operation))
        print("require input from:  " + str([str(minputs) for minputs in self.require_input_from]))
        print("inputs_task: " + str([minputs for minputs in self.inputs_task]))
        print("  #+~  \n\n")
