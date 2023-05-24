import pandas as pd
from scipy.stats import hypergeom
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

    print(geneSet_list)
    print(len(geneSet_list),'gensenlist')
    print(len(geneSet),'total')

    # STEP 0.fdr matrix
    lenGL = len(geneSet_list)-1
    q005_selec = q005.iloc[:, lenGL]
    q001_selec = q001.iloc[:, lenGL]
    q0001_selec = q0001.iloc[:, lenGL]
    new_fdr = pd.DataFrame()
    new_fdr['q005'] = q005_selec
    new_fdr['q001'] = q001_selec
    new_fdr['q0001'] = q0001_selec
    # new_fdr = new_fdr.fillna(1)

    # STEP.1 pathway DB
    pathwayDB = pd.read_csv(f"{dbType}.csv")
    pathwayDB = pathwayDB.dropna()
    pathwayDB = pathwayDB.fillna(0)

    # STEP.2 gene symbol or entrez ID
    pathwayDB_sort = pathwayDB.iloc[:, [1, 2, 3, 7]]
    selected_genes = pathwayDB_sort.iloc[:, 1][0].split(';')

    # STEP.3 hypergeom
    pvalue = []
    pvalue_string = []
    qvalue = []
    qvalue_level = []
    pathway = []
    source = []
    size = []
    overlap = []
    cancerLevel = []

    for id in range(len(pathwayDB)):
        selected_genes = pathwayDB_sort.iloc[id, 2]
        selected_genes = str(selected_genes).split(';')  # string으로 바꿔야 가능
        selected_list = [int(s) for s in selected_genes]
        intersection = list(set(selected_list) & set(geneSet_list))
        selected_pvalue = (hypergeom.sf(
            len(intersection)-1, 27214, len(geneSet_list), len(selected_list)))
        pvalue.append(selected_pvalue)
        selected_pvalue_string = f"{selected_pvalue:.2E}"
        pvalue_string.append(selected_pvalue_string)
        selected_pathway = pathwayDB_sort.iloc[id, 0]
        pathway.append(selected_pathway)
        size.append(len(selected_list))
        overlap.append(len(intersection))
        selected_cancerLevel = pathwayDB_sort.iloc[id, 3]
        cancerLevel.append(selected_cancerLevel)

        # size가 len(selected_list)
        gsize = len(selected_list)-1
        fdr005 = new_fdr.iloc[gsize, 0]
        fdr001 = new_fdr.iloc[gsize, 1]
        fdr0001 = new_fdr.iloc[gsize, 2]

        if selected_pvalue < fdr0001:
            selected_qval = fdr0001
            selected_qval_level = 0.001
        elif selected_pvalue < fdr001:
            selected_qval = fdr001
            selected_qval_level = 0.01
        elif selected_pvalue < fdr005:
            selected_qval = fdr005
            selected_qval_level = 0.05
        else:
            selected_qval = 1
            selected_qval_level = 1
        qvalue.append(round(selected_qval, 4))
        qvalue_level.append(selected_qval_level)

        # source add
        selected_source = pathwayDB_sort.iloc[id, 1]
        source.append(selected_source)

    # STEP.4 result table
    result_table = pd.DataFrame()
    result_table['pvalue'] = pvalue
    result_table['pvalue_string'] = pvalue_string
    result_table['qvalue'] = qvalue
    result_table['qvalue_level'] = qvalue_level
    result_table['pathway'] = pathway
    result_table['source'] = source
    result_table['size'] = size
    result_table['overlap'] = overlap
    result_table['cancerLevel'] = cancerLevel

    result_table = result_table[result_table['overlap'] > 1]
    result_table = result_table[result_table['pvalue'] < 0.05]
    result_table = result_table[result_table['qvalue_level'] <= qValueCutoff]
    
    countTB = result_table['cancerLevel'].value_counts()
    countTB_df = pd.DataFrame(countTB)
    countTB_df['level'] = countTB_df.index
    countTB_lv1 = countTB_df[countTB_df['level'] == 1]
    countTB_lv2 = countTB_df[countTB_df['level'] == 2]
    countTB_lv3 = countTB_df[countTB_df['level'] == 3]
    countTB_lv4 = countTB_df[countTB_df['level'] == 4]

    result_table = result_table[result_table['cancerLevel']
                                == inputCancerLevel]
    print(result_table)
    # result_table = result_table.sort_values(by=['pvalue'])
    # print(result_table.head(10))
    result_json = result_table.transpose().to_json()

    # return result
    return {
        "data": result_json,
        "countDB": {
            # "lv1" : int(countTB_lv1['cancerLevel'].values[0]),    
            # "lv2" : int(countTB_lv2['cancerLevel'].values[0]),    
            # "lv3" : int(countTB_lv3['cancerLevel'].values[0]),
            # "lv4" : int(countTB_lv4['cancerLevel'].values[0])
            "lv1" : 10,
            "lv2" : 10,
            "lv3" : 10,
            "lv4" : 10,
        },
        "matchInfo" : {
            "total" : len(geneSet),
            "intersec" : len(geneSet_list),
            "percent" : round((len(geneSet_list)/len(geneSet)) * 100)
        }
    }
