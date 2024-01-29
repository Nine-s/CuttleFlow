
import os
import json
import numpy as np
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
from scipy.optimize import nnls


class AnnotationDB:  

    annotation_db:list = []
    runtime_estimation_models:object = None


    
    @staticmethod
    def create_alignment_runtime_estimation_model_alt(runtime_measured):

        df_runtime = pd.read_csv(runtime_measured, delimiter=",")

        aligners = (df_runtime.columns[5:])
        #print(aligners)

        columns_to_normalize = ['dataset_size', "RAM", "CPUMHz", "ref_genome_size"]
        X = df_runtime[columns_to_normalize]

        scaler = StandardScaler()
        X_normalized = scaler.fit_transform(X)

        df_normalized = pd.DataFrame(X_normalized, columns=columns_to_normalize)
        df_runtime[columns_to_normalize] = df_normalized
        #print(df_runtime)

        list_models = {}
        for aligner in aligners:
            df = pd.DataFrame({'dataset_size': df_runtime["dataset_size"], 'runtime': df_runtime[aligner], 
                               "RAM": df_runtime["RAM"], "CPU": df_runtime["CPUMHz"], "ref_size": df_runtime["ref_genome_size"]})

            X = df[['dataset_size', "RAM", "CPU", "ref_size"]]
            y = df['runtime']

            poly = PolynomialFeatures(degree=2)
            X_poly = poly.fit_transform(X)
            model = LinearRegression()
            model.fit(X_poly, y)

            list_models[aligner] = model
            
            test = False
            if(aligner == "Salmon" and test is True):
                print("SALMON")
                #test 2 with FONDA cluster: D1
                new_data = pd.DataFrame({'dataset_size': [(3.5)],
                                    'RAM': [(251)], 
                                    'CPUMHz': [(3400.0000)], 
                                    'ref_genome_size': [(0.137)] })
                normalized_new_data = scaler.transform(new_data)

                new_data_point_poly = poly.transform(normalized_new_data)
                predicted_runtime = model.predict(new_data_point_poly)
                print("$$$$ test 1")
                print("real runtime: 8.3"  )
                print("predicted runtime: "+ str(predicted_runtime))

                #test 2 with FONDA cluster: D2
                new_data = pd.DataFrame({'dataset_size': [(13)],
                                    'RAM': [(251)], 
                                    'CPUMHz': [(3400.0000)], 
                                    'ref_genome_size': [(0.137)] })
                normalized_new_data = scaler.transform(new_data)

                new_data_point_poly = poly.transform(normalized_new_data)
                predicted_runtime = model.predict(new_data_point_poly)
                print("$$$$ test 2")
                print("real runtime: 8.3"  )
                print("predicted runtime: "+ str(predicted_runtime))

                #test 2 with FONDA cluster: D3
                new_data = pd.DataFrame({'dataset_size': [(48)],
                                    'RAM': [(251)], 
                                    'CPUMHz': [(3400.0000)], 
                                    'ref_genome_size': [(0.137)] })
                normalized_new_data = scaler.transform(new_data)

                new_data_point_poly = poly.transform(normalized_new_data)
                predicted_runtime = model.predict(new_data_point_poly)
                print("$$$$ test 3")
                print("real runtime: 81"  )
                print("predicted runtime: "+ str(predicted_runtime))

            if(aligner == "HISAT2" and test is True):
                print("HISAT2")
                #test 2 with FONDA cluster: D1
                new_data = pd.DataFrame({'dataset_size': [(3.5)],
                                    'RAM': [(251)], 
                                    'CPUMHz': [(3400.0000)], 
                                    'ref_genome_size': [(0.137)] })
                normalized_new_data = scaler.transform(new_data)

                new_data_point_poly = poly.transform(normalized_new_data)
                predicted_runtime = model.predict(new_data_point_poly)
                print("$$$$ test 1")
                print("real runtime: 292"  )
                print("predicted runtime: "+ str(predicted_runtime))

                #test 2 with FONDA cluster: D2
                new_data = pd.DataFrame({'dataset_size': [(13)],
                                    'RAM': [(251)], 
                                    'CPUMHz': [(3400.0000)], 
                                    'ref_genome_size': [(0.137)] })
                normalized_new_data = scaler.transform(new_data)

                new_data_point_poly = poly.transform(normalized_new_data)
                predicted_runtime = model.predict(new_data_point_poly)
                print("$$$$ test 2")
                print("real runtime: 1178.1"  )
                print("predicted runtime: "+ str(predicted_runtime))

                #test 2 with FONDA cluster: D3
                new_data = pd.DataFrame({'dataset_size': [(48)],
                                    'RAM': [(251)], 
                                    'CPUMHz': [(3400.0000)], 
                                    'ref_genome_size': [(0.137)] })
                normalized_new_data = scaler.transform(new_data)

                new_data_point_poly = poly.transform(normalized_new_data)
                predicted_runtime = model.predict(new_data_point_poly)
                print("$$$$ test 3")
                print("real runtime: 3916"  )
                print("predicted runtime: "+ str(predicted_runtime))

        #    dataset_size = df[["dataset_size"]]
        #     plt.scatter(dataset_size, y, label=aligner)
        #     plt.plot(dataset_size, model.predict(X_poly))
        #     plt.legend()
        # plt.title('Aligner: '+ aligner)
        # plt.show()

        return [list_models, scaler]

    @staticmethod
    def create_split_runtime_estimation_model(runtime_measured):

        df_runtime = pd.read_csv(runtime_measured, delimiter=",")
        columns_to_normalize = ['dataset_size', "RAM", "CPUMHz", "ref_genome_size"]
        X = df_runtime[columns_to_normalize]

        scaler = StandardScaler()
        X_normalized = scaler.fit_transform(X)

        df_normalized = pd.DataFrame(X_normalized, columns=columns_to_normalize)
        df_runtime[columns_to_normalize] = df_normalized
        #print(df_runtime)

        X = df_runtime[['dataset_size', "RAM", "CPUMHz", "ref_genome_size"]]
        y = df_runtime['duration']

        poly = PolynomialFeatures(degree=2)
        X_poly = poly.fit_transform(X)
        model = LinearRegression()
        model.fit(X_poly, y)
        return model, scaler


    def __init__ (self, annotation_files_list):
        
        path_runtimes_align = "./annotation_files/runtime_aligners_with_CPU_RAM.csv"
        if( os.path.isfile(path_runtimes_align) == False):
            raise Exception("File missing: "+path_runtimes_align) 
        else:
            with open(path_runtimes_align) as mfile:
                runtime_model = self.create_alignment_runtime_estimation_model_alt(mfile)          
                self.runtime_estimation_models = runtime_model[0]
                self.sandardscaler = runtime_model[1]
        path_runtimes_split = "./annotation_files/runtime_split_merge.csv"
        if( os.path.isfile(path_runtimes_split) == False):
            raise Exception("File missing: "+path_runtimes_split) 
        else:
            with open(path_runtimes_split) as mfile:
                split_estimators = self.create_split_runtime_estimation_model(mfile)   
        #for infra in split_estimators:
        self.runtime_estimation_models["split_merge"], self.split_merge_scaler = split_estimators
        annotation_db = []
        for file_path in annotation_files_list:
            with open(file_path) as json_file:
                #print(file_path)                

                tool_annotated = ToolAnnotation(json.load(json_file))
            annotation_db.append(tool_annotated)
        self.annotation_db = annotation_db
        

class ToolAnnotation:

    toolname:str = ""
    operation:list = []
    domain_specific_features:list = []
    is_splittable:bool = False
    mendatory_input_list:list = []
    output_list:list = []
    #optional_inputs_list:list = []
    RAM_requirements_model:object = None

    @staticmethod
    def create_resource_requirements_RAM (reference_sizes, ram_used):
        df = pd.DataFrame({'reference_sizes': reference_sizes, 'ram_used': ram_used })
        df.head()
        X = df.iloc[:,:-1].values # feature matrix: reference_sizes
        y = df.iloc[:,1].values # response vector: ram_used

        model = LinearRegression()
        model.fit(X, y)

        ###plot
        # plt.scatter(X, y,color='g')
        # plt.plot(X, model.predict(X),color='k')
        # plt.show()

        return model
    
    
    def __init__ (self, tool_description):
    
        self.toolname = tool_description["toolname"]
        self.operation = tool_description["operation"]
        self.domain_specific_features = tool_description["domain_specific_features"]
        self.is_splittable = tool_description["is_splittable"]
        self.mendatory_input_list = tool_description["mendatory_input_list"]
        #self.optional_inputs_list = tool_description["optional_inputs_list"]
        self.output_list = tool_description["output_list"]
        self.module_path =  tool_description["module_path"]
        self.module_name =  tool_description["module_name"]
        reference_sizes = []
        ram_used = []
        for resource_requirements in tool_description["resource_requirements_RAM"]:
            RAM_require = resource_requirements["RAM"].split("GB")[0][:-1]
            ref_size = resource_requirements["reference_size"][:-1]
            reference_sizes.append(float(ref_size))
            ram_used.append(float(RAM_require)) #in GB
        self.RAM_requirements_model = self.create_resource_requirements_RAM(reference_sizes, ram_used)
