import pandas as pd
import numpy as np

def data_parsing(DATA_PATH , GRAPH_FILE ,TARGET , INDEX_COL) :
    
    META_DATA_PATH = [f'{DATA_PATH}/datMeta_{mod}.csv' for mod in GRAPH_FILE[:-4].split('_')[1:-1]]

    meta = pd.Series(dtype=str)
    for path in META_DATA_PATH : 
        meta_tmp = pd.read_csv(path , index_col=0)
        
        if INDEX_COL == '' :
            pass
        else :
            meta_tmp = meta_tmp.set_index(INDEX_COL)
            
        meta = pd.concat([meta , meta_tmp[TARGET]])

    meta = meta[~meta.index.duplicated(keep='first')] # Remove duplicated entries
    meta.index = [str(i) for i in meta.index] # Ensures the patient ids are strings

    TRAIN_DATA_PATH = [f'{DATA_PATH}/datExpr_{mod}.csv' for mod in GRAPH_FILE[:-4].split('_')[1:-1]] # Looks for all expr file names
    datModalities = {}
    for path in TRAIN_DATA_PATH : 
        print('Importing \t %s \n' % path)
        
        dattmp = np.genfromtxt(path , delimiter=',' , dtype = str)
        if len(set(meta.index.astype(str)) & set(np.core.defchararray.strip(dattmp[1: , 0], '"'))) > 0 :
            dattmp = pd.DataFrame(dattmp[1: , 1:] , columns=np.core.defchararray.strip(dattmp[0 ,1:], '"') , index = [int(i.strip('"')) for i in dattmp[1: , 0]])
        else : 
            dattmp = pd.DataFrame(dattmp[1: , 1:] , columns=[int(i.strip('"')) for i in dattmp[0,1:]] , index = np.core.defchararray.strip(dattmp[1: , 0], '"'))
            dattmp = dattmp.T

        dattmp.columns = [str(i) for i in dattmp.columns] # Ensures the patient ids are strings
        dattmp.name = path.split('.')[0].split('_')[-1] #Assumes there is no '.' in file name as per specified naming convention. Can lead to erros down stream. Files should be modality_datEXpr.csv e.g. mRNA_datExpr.csv
        datModalities[dattmp.name] = dattmp

    return datModalities , meta
