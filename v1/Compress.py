import DAW #from DAW import add task to daw
from statistics import median
from task import Task


def create_compress_tasks(compression_type, task, input_description):
    
    ## Compress task
    operation = "compression"
    name = compression_type + "_" + operation
    if (compression_type == "naf"):
        tool = "naf" 
        outputs = "naf"
    elif (compression_type == "cram"):
        tool = "samtools" 
        outputs = "cram"
    module_path = "./nextflow_modules/" + operation + "_" + compression_type + ".nf"
    module_name = module_path.split("/")[-1].split(".")[0]
    parameters = []
    operation = operation
    inputs = [(task.name + ".out_channel" + out) for out in task.outputs]
    inputs_from_DAW = task.inputs_from_DAW
    

    new_task_compress = Task(name=name, tool=tool, inputs_from_DAW=inputs_from_DAW, outputs=outputs, parameters=parameters, operation=operation, module_name=module_name, module_path=module_path, input_description=input_description)
    new_task_compress.my_print()

    ## decompress task
    operation = "decompression"
    name = compression_type + "_" + operation
    if (compression_type == "naf"):
        tool = "naf" 
        outputs = "naf"
    elif (compression_type == "cram"):
        tool = "samtools" 
        outputs = "bam"#
    module_path = "./nextflow_modules/" + operation + "_" + compression_type + ".nf"
    module_name = module_path.split("/")[-1].split(".")[0]
    inputs = [new_task_compress.name + ".out_channel" + new_task_compress.outputs]
    parameters = []
    operation = operation
    inputs_from_DAW = task.inputs_from_DAW
    input_description = task.input_description

    new_task_decompress= Task(name=name, tool=tool, inputs_from_DAW=inputs_from_DAW, outputs=outputs, parameters=parameters, operation=operation, module_name=module_name, module_path=module_path, input_description=input_description)
    return [new_task_compress, new_task_decompress]


def is_output_fastq(task):
    if ('fastq' in task.outputs):
        return True
    else:
        return False
    
def is_output_sam(task):
    if ('sam' in task.outputs):
        return True
    else:
        return False

def compress_before_file_transfer(daw, input_description):
    for task in daw.tasks:
        print(task.outputs)
        if ( median(daw.input.size_of_samples) < 1):
            continue
        elif( is_output_fastq(task)):
            print("### fastq output found")    
            new_tasks = create_compress_tasks("naf", task, input_description)
            new_task_compress  = new_tasks[0]
            new_task_decompress = new_tasks[1]
        # elif (is_output_sam(task)):
        #     print("### sam output found")
        #     new_tasks = create_compress_tasks("cram", task, input_description)
        #     new_task_compress  = new_tasks[0]
        #     new_task_decompress = new_tasks[1]
        # TODO: change input of child tasks to output of decompress task
        #child_tasks = [task for task in daw.tasks if output_last_split_task in task.require_input_from]
        #for child_task in child_tasks:
        #    child_task.change_input(merge_task_output, output_last_split_task)
        # TODO: add decompress task
            daw.insert_tasks(new_task_compress)
            daw.insert_tasks(new_task_decompress)
    return daw