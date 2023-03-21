import pandas as pd
from scipy.stats import hypergeom
import re


q005 = pd.read_csv('0.05.csv')
q001 = pd.read_csv('0.01.csv')
q0001 = pd.read_csv('0.001.csv')
# STEP.1 pathway DB
pathwayDB = pd.read_csv("CPDB.csv")
pathwayDB = pathwayDB.dropna()
pathwayDB = pathwayDB.fillna(0)


def analyze(geneSetType, geneSet, qValueCutoff, inputCancerLevel):

    geneSet = re.sub(r"\s", "", geneSet)
    geneSet_list = geneSet.split(',')
    geneSet_list = [v for v in geneSet_list if v]

    if(geneSetType == 'int'):
        geneSet_list = [int(s) for s in geneSet_list]

    # STEP 1.fdr matrix
    lenGL = len(geneSet_list)-1
    q005_selec = q005.iloc[:, lenGL]
    q001_selec = q001.iloc[:, lenGL]
    q0001_selec = q0001.iloc[:, lenGL]
    new_fdr = pd.DataFrame()
    new_fdr['q005'] = q005_selec
    new_fdr['q001'] = q001_selec
    new_fdr['q0001'] = q0001_selec
    new_fdr = new_fdr.fillna(1)

    # STEP.2 gene symbol or entrez ID
    if(geneSetType == 'int'):
        pathwayDB_sort = pathwayDB.iloc[:, [1, 2, 3, 7]]
    else:
        pathwayDB_sort = pathwayDB.iloc[:, [1, 2, 4, 7]]

    selected_genes = pathwayDB_sort.iloc[:, 1][0].split(';')

    # STEP.3 hypergeom
    pvalue = []
    qvalue = []
    pathway = []
    source = []
    size = []
    overlap = []
    cancerLevel = []

    if(geneSetType == 'int'):
        for id in range(len(pathwayDB)):
            selected_genes = pathwayDB_sort.iloc[id, 2]
            selected_genes = str(selected_genes).split(';')  # string으로 바꿔야 가능
            selected_list = [int(s) for s in selected_genes]
            intersection = list(set(selected_list) & set(geneSet_list))

            selected_pvalue = (round(hypergeom.sf(len(intersection)-1, 20501,
                                                  len(geneSet_list), len(selected_list)), 4))
            pvalue.append(selected_pvalue)
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
                selected_qval = 0.001
            elif selected_pvalue < fdr001:
                selected_qval = 0.01
            elif selected_pvalue < fdr005:
                selected_qval = 0.05
            else:
                selected_qval = 1
            qvalue.append(selected_qval)

            # source add
            selected_source = pathwayDB_sort.iloc[id, 1]
            source.append(selected_source)

    if(geneSetType == 'str'):
        for id in range(len(pathwayDB)):
            selected_genes = pathwayDB_sort.iloc[id, 2]
            selected_list = str(selected_genes).split(';')  # string으로 바꿔야 가능
            intersection = list(set(selected_list) & set(geneSet_list))
            selected_pvalue = (round(hypergeom.sf(len(intersection)-1, 20501,
                                                  len(geneSet_list), len(selected_list)), 4))

            pvalue.append(selected_pvalue)
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
                selected_qval = 0.001
            elif selected_pvalue < fdr001:
                selected_qval = 0.01
            elif selected_pvalue < fdr005:
                selected_qval = 0.05
            else:
                selected_qval = 1
            qvalue.append(selected_qval)

            # source add
            selected_source = pathwayDB_sort.iloc[id, 1]
            source.append(selected_source)

    # STEP.4 result table
    result_table = pd.DataFrame()
    result_table['pvalue'] = pvalue
    result_table['qvalue'] = qvalue
    result_table['pathway'] = pathway
    result_table['source'] = source
    result_table['size'] = size
    result_table['overlap'] = overlap
    result_table['cancerLevel'] = cancerLevel

    result_table = result_table[result_table['overlap'] > 4]
    result_table = result_table[result_table['pvalue'] < 0.05]
    result_table = result_table[result_table['qvalue'] <= qValueCutoff]
    result_table = result_table[result_table['cancerLevel']
                                == inputCancerLevel]
    # result_table = result_table.sort_values(by=['pvalue'])
    # print(result_table.head(10))
    result_json = result_table.transpose().to_json()

    # return result
    return {
        "data": result_json,
    }
