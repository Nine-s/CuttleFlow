
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
from scipy.spatial import distance

class AnnotationDB:  

    annotation_db:list = []
    runtime_estimation_models:object = None
    dataset_scaler = []
    CPU_scaler = []
    RAM_scaler = []

    @staticmethod
    def create_alignment_runtime_estimation_model_alt(runtime_measured):

        df_runtime = pd.read_csv(runtime_measured, delimiter=",")

        aligners = (df_runtime.columns[6:])
        #print(aligners)
        infrastructures = set(df_runtime.infrastructure.tolist())


        columns_to_normalize = ['dataset_size']
        dataset_sizes = df_runtime[columns_to_normalize]
        print(dataset_sizes)
        ram_values = df_runtime[["RAM"]]
        cpu_values = df_runtime[["CPUMHz"]]
        
        ram_scaler = StandardScaler()
        RAM_normalized = ram_scaler.fit_transform(ram_values)
        df_RAM = pd.DataFrame(RAM_normalized, columns=["RAM"])
        df_runtime.loc[:,["RAM"]] = df_RAM

        cpu_scaler = StandardScaler()
        CPU_normalized = cpu_scaler.fit_transform(cpu_values)
        df_CPU = pd.DataFrame(CPU_normalized, columns=["CPUMHz"])
        df_runtime.loc[:,["CPUMHz"]] = df_CPU

        dataset_size_scaler = StandardScaler()
        dataset_size_scaler.fit(dataset_sizes)
        list_models = {}

        for infra in infrastructures:
            df_infra = df_runtime[df_runtime.infrastructure==infra]
            list_models[infra] = {}
            list_models[infra]["RAM"] = np.mean(df_infra["RAM"].tolist())
            list_models[infra]["CPU"] = np.mean(df_infra["CPUMHz"].tolist())

            X_normalized = pd.Series(dataset_size_scaler.transform(df_infra[columns_to_normalize]).flatten())
            #X_normalized.reset_index(inplace=True)
            df_infra.reset_index(inplace=True)
            df_normalized = pd.DataFrame(X_normalized, columns=columns_to_normalize)
            df_infra.loc[:, columns_to_normalize] = df_normalized
            
            for aligner in aligners:
                
                df = pd.DataFrame({'dataset_size': df_infra["dataset_size"], 'runtime': df_infra[aligner], 
                                "RAM": df_infra["RAM"], "CPU": df_infra["CPUMHz"], "ref_size": df_infra["ref_genome_size"]})
                df = df[df['runtime'].notnull()]
                if not df.empty:
                    X = df[['dataset_size']]
                    y = df['runtime']

                    poly = PolynomialFeatures(degree=1)
                    X_poly = poly.fit_transform(X)
                    model = LinearRegression()
                    model.fit(X_poly, y)

                    list_models[infra][aligner] = model

        return [list_models, dataset_size_scaler, ram_scaler, cpu_scaler]

    def create_split_runtime_estimation_model(self, runtime_measured):

        df_runtime = pd.read_csv(runtime_measured, delimiter=",")

        for infra in self.runtime_estimation_models:
            print(infra)
            df_infra = df_runtime[df_runtime.infrastructure == infra]
            if df_infra.empty:
                print("No annotation data to fit a runtime estimation model for split-merge-processes for infrastructure " + infra)
            else:
                column_to_normalize = ['dataset_size']
                X_normalized = pd.Series(self.dataset_scaler.transform(df_infra[column_to_normalize]).flatten())
                df_infra.reset_index(inplace=True)
                df_normalized = pd.DataFrame(X_normalized, columns=column_to_normalize)
                df_infra.loc[:, column_to_normalize] = df_normalized
                #print(df_runtime)
                X = df_infra[['dataset_size']]
                y = df_infra['duration']

                poly = PolynomialFeatures(degree=1)
                X_poly = poly.fit_transform(X)
                model = LinearRegression()
                model.fit(X_poly, y)
                print(model.coef_)
                print("______")
                self.runtime_estimation_models[infra]["split-merge"] = model
        return None

    def find_most_similar_infrastructure(self, RAM, CPU):
        RAM_norm = self.RAM_scaler.transform(np.array([RAM]).reshape(-1, 1))
        CPU_norm = self.CPU_scaler.transform(np.array([CPU]).reshape(-1, 1))

        min_dist = float("inf")
        most_similar_infrastructure = ""
        for infra in self.runtime_estimation_models:
            if(isinstance(self.runtime_estimation_models[infra], dict)):
                infra_CPU = self.runtime_estimation_models[infra]["CPU"]
                infra_RAM = self.runtime_estimation_models[infra]["RAM"]
                infra_distance = distance.euclidean([RAM_norm, CPU_norm], [infra_RAM, infra_CPU])
                if(infra_distance < min_dist):
                    min_dist = infra_distance
                    most_similar_infrastructure = infra
        
        return most_similar_infrastructure

    def predict_runtime(self, task_name, RAM, CPU, dataset_size):
        poly = PolynomialFeatures(degree=1)

        most_similar_infrastructure = self.find_most_similar_infrastructure(RAM, CPU)
        align_scaler = self.dataset_scaler
        infra_runtime_models = self.runtime_estimation_models[most_similar_infrastructure]
        runtime_models = {k.lower(): infra_runtime_models[k] for k in infra_runtime_models}
        task_runtime_model = runtime_models[task_name.lower()] 
        scaled_dataset_size = align_scaler.transform(np.array(dataset_size).reshape(-1, 1))
        model_input = poly.fit_transform(scaled_dataset_size)
        estimated_task_time = task_runtime_model.predict(model_input)

        return estimated_task_time


        

    def __init__ (self, annotation_files_list):
        
        path_runtimes_align = "./annotation_files/runtime_aligners_with_CPU_RAM.csv"
        if( os.path.isfile(path_runtimes_align) == False):
            raise Exception("File missing: "+path_runtimes_align) 
        else:
            with open(path_runtimes_align) as mfile:
                runtime_model = self.create_alignment_runtime_estimation_model_alt(mfile)          
                self.runtime_estimation_models = runtime_model[0]
                self.dataset_scaler = runtime_model[1]
                self.RAM_scaler = runtime_model[2]
                self.CPU_scaler = runtime_model[3]
        path_runtimes_split = "./annotation_files/runtime_split_merge.csv"
        if( os.path.isfile(path_runtimes_split) == False):
            raise Exception("File missing: "+path_runtimes_split) 
        else:
            with open(path_runtimes_split) as mfile:
                self.create_split_runtime_estimation_model(mfile)   

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
        if self.operation == "align":
            self.domain_specific_features = tool_description["domain_specific_features"]
            reference_sizes = []
            ram_used = []
            for resource_requirements in tool_description["resource_requirements_RAM"]:
                RAM_require = resource_requirements["RAM"].split("GB")[0][:-1]
                ref_size = resource_requirements["reference_size"][:-1]
                reference_sizes.append(float(ref_size))
                ram_used.append(float(RAM_require)) #in GB
            self.RAM_requirements_model = self.create_resource_requirements_RAM(reference_sizes, ram_used)

        
        self.is_splittable = tool_description["is_splittable"]
        
        self.mendatory_input_list = tool_description["mendatory_input_list"]
        #self.optional_inputs_list = tool_description["optional_inputs_list"]
        self.output_list = tool_description["output_list"]
        self.module_path =  tool_description["module_path"]
        self.module_name =  tool_description["module_name"]
        

"""
 #TODO: test cases must be rewritten
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
"""
