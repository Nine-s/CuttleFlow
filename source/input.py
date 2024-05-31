import os
import pathlib 
import csv


    
class Input:

    uncompressed_size:int = 1 #GB
    name:str = "sample"
    input_type:str = "reads"
    paths = []
    size:float = 1
    paired:bool = True
    strand:str = "forward"
    ref_type:str = ""
    additional_columns:dict = None


    def __init__(self, name, input_type, paths, strand, ref_type, uncompressed_size, additional_columns=None):
        
        self.name = name
        #print(name) 
        possible_input_types = ["sample", "reference"]
        self.uncompressed_size = uncompressed_size
        if not (input_type in possible_input_types): raise Exception( 'The variable "input_type" should be either "sample" or "reference"')
        else: self.input_type = input_type
        #TODO after testing: number of lines and uncompressed size from os.command
        # for path in paths:
        #     if not (os.path.exists(path)): raise Exception( 'The input file provided ' + path + ' does not exist') 
        self.paths = paths
        
        possible_strand_types = ["forward", "reverse", "unstranded"]

        possible_ref_types = ["genome", "annotation", "transcriptome", "others"]
        

        if (input_type=="sample"):
            if (strand not in possible_strand_types): 
                raise Exception( 'The variable "strand" should be "forward", "reverse" or "unstranded", not '+strand)
            else:
                self.strand = strand

        elif (input_type=="reference"):
            
            if not (ref_type in possible_ref_types): raise Exception( 'The variable "ref_type" should be "genome", "annotation" or "others"')
            self.ref_type = ref_type
            strand = ""
        self.additional_columns = additional_columns 
        
class Input_of_DAW:  

    number_of_samples:int = 0
    size_of_samples:list = []
    input_samples:list = []
    input_references:list = []
    input_parameters:dict = {}
    first_strand = "foward"
    size_of_reference_genome_max:int = 0 #TODO: consider metagenomics: add a condition to take the max number?

    def __init__(self, input_description):
        _input_samples = []
        _input_references = []
        _size_of_samples = []
        _size_of_reference_genome_max = -1
        
        for sample in input_description["samples"]:
            standard_columns = ["name","type","path_r1","path_r2","strand","uncompressed_size"]
            additional_columns = {description: value for description, value in sample.items() if description not in standard_columns}
            mInput = Input(sample["name"], "sample", [sample["path_r1"], sample["path_r2"]], sample["strand"], "", sample["uncompressed_size"], additional_columns)
            _input_samples.append(mInput) 
            self.first_strand = mInput.strand    
            _size_of_samples.append(float(sample["uncompressed_size"]))
                
        for reference in input_description["references"]:
            mInput = Input(reference["name"], "reference", [reference["path"]], "", reference["reference_type"], reference["uncompressed_size"])
            _input_references.append(mInput)
            if (reference["reference_type"] == "genome"):
                if (_size_of_reference_genome_max == -1):
                    _size_of_reference_genome_max = float(reference["uncompressed_size"])
                else:
                    _size_of_reference_genome_max = float(max(_size_of_reference_genome_max, float(reference["uncompressed_size"])))
        self.number_of_samples = len(_input_samples)         
        self.input_references = _input_references
        self.size_of_samples = _size_of_samples
        self.input_samples = _input_samples     
        self.size_of_reference_genome_max = _size_of_reference_genome_max
        for param in input_description["parameters"]:
            value = isinstance(param["value"], str) * "'" + str(param["value"]) + isinstance(param["value"], str) * "'"
            self.input_parameters[param["name"]] = value     
       
