import pandas as pd
from scipy.stats import hypergeom
import statsmodels.api as sm
import statsmodels
import re
import numpy as np

q005 = pd.read_csv('0.05.csv')
q001 = pd.read_csv('0.01.csv')
q0001 = pd.read_csv('0.001.csv')
alias_data = pd.read_csv('gdac_entrez.csv')


def analyze(dbType, geneSetType, geneSet, qValueCutoff, inputCancerLevel):

    geneSet = re.sub(r"\s", "", geneSet)
    geneSet = geneSet.split(',')

    if(geneSetType == 'str'):
        alias_data2 = alias_data.set_index(alias_data['symbol'])
        intersec = list(set(alias_data['symbol'].to_list()) & set(geneSet))

    if(geneSetType == 'int'):
        alias_data2 = alias_data.set_index(alias_data['entrez'])
        geneSet_list = [int(v) for v in geneSet]
        intersec = list(set(alias_data['entrez'].to_list()) & set(geneSet_list))
        
    geneSet_list = alias_data2.loc[intersec, 'entrez'].to_list()

    # STEP 0.fdr matrix
    lenGL = len(geneSet_list)-1
    q005_selec = q005.iloc[:, lenGL]
    q001_selec = q001.iloc[:, lenGL]
    q0001_selec = q0001.iloc[:, lenGL]
    new_fdr = pd.DataFrame()
    new_fdr['q005'] = q005_selec
    new_fdr['q001'] = q001_selec
    new_fdr['q0001'] = q0001_selec


    # STEP.1 pathway DB
    pathwayDB = pd.read_csv(f"{dbType}.csv")
    pathwayDB = pathwayDB.dropna()
    pathwayDB = pathwayDB.fillna(0)
    if(inputCancerLevel == '1'):
        pathwayDB = pathwayDB[pathwayDB['cancer_level'] == 1]
    elif (inputCancerLevel == '1&2'):
        pathwayDB = pathwayDB[(pathwayDB['cancer_level'] == 1) | (pathwayDB['cancer_level'] == 2)]
    elif (inputCancerLevel == '1&2&3'):
        pathwayDB = pathwayDB[(pathwayDB['cancer_level'] == 1) | (pathwayDB['cancer_level'] == 2) | (pathwayDB['cancer_level'] == 3)]
    

    # STEP.2 hypergeom
    pvalue = []
    pvalue_string = []
    qvalue = []
    pathway = []
    source = []
    size = []
    overlap = []
    cancerLevel = []
    q_005 = []
    q_001 = []
    q_0001 = []
    _genes_ = []

    for id in range(len(pathwayDB)):
        selected_genes = pathwayDB.iloc[id, 2]
        selected_genes = str(selected_genes).split(';')  # string으로 바꿔야 가능
        selected_list = [int(s) for s in selected_genes]
        intersection = list(set(selected_list) & set(geneSet_list))
        selected_pvalue = (hypergeom.sf(
            len(intersection)-1, 27214-len(geneSet_list), len(geneSet_list), len(selected_list)))
        pvalue.append(selected_pvalue)
        selected_pvalue_string = f"{selected_pvalue:.2E}"
        pvalue_string.append(selected_pvalue_string)
        selected_pathway = pathwayDB.iloc[id, 0]
        pathway.append(selected_pathway)
        size.append(len(selected_list))
        overlap.append(len(intersection))
        selected_cancerLevel = pathwayDB.iloc[id, 3]
        cancerLevel.append(selected_cancerLevel)

        # size가 len(selected_list)
        gsize = len(selected_list)-1
        fdr005 = new_fdr.iloc[gsize, 0]
        fdr001 = new_fdr.iloc[gsize, 1]
        fdr0001 = new_fdr.iloc[gsize, 2]
        
        q_005.append(fdr005)
        q_001.append(fdr001)
        q_0001.append(fdr0001)
        
        # source add
        selected_source = pathwayDB.iloc[id, 1]
        source.append(selected_source)
        
        # gene add
        intersection_str = ','.join([str(i) for i in intersection])
        _genes_.append(intersection_str)

    qvalue = statsmodels.stats.multitest.fdrcorrection(pvalue, alpha=0.05, method='indep', is_sorted=False)[1]
    
    # STEP.3 result table
    result_table = pd.DataFrame()
    result_table['pvalue'] = pvalue
    result_table['qval'] = qvalue
    
    result_table['q_005_val'] = q_005
    result_table['q_001_val'] = q_001
    result_table['q_0001_val'] = q_0001

    result_table['pathway'] = pathway
    result_table['source'] = source
    result_table['size'] = size
    result_table['overlap'] = overlap
    result_table['cancerLevel'] = cancerLevel
    result_table['_genes_'] = _genes_

    result_table = result_table[result_table['overlap'] > 1]
    result_table = result_table[result_table['pvalue'] < result_table[qValueCutoff]]
    result_table = result_table[result_table['qval'] < 0.2]
    
    result_json = result_table.transpose().to_json()

    # return result
    return {
        "data": result_json,
    }
