import pandas as pd
from scipy.stats import hypergeom
import json
import re

q3000 = pd.read_csv('3000.csv')
ccg = pd.read_csv('ccg.csv')
alias_data = pd.read_csv('gdac_entrez.csv')


def analyze(geneSetType, geneSet):
    
    geneSet = re.sub(r"\s", "", geneSet)
    geneSet = geneSet.split(',')

    if(geneSetType == 'str'):
        alias_data2 = alias_data.set_index(alias_data['symbol'])
        intersec = list(set(alias_data['symbol'].to_list()) & set(geneSet))

    if(geneSetType == 'int'):
        alias_data2 = alias_data.set_index(alias_data['entrez'])
        geneSet_list = [int(v) for v in geneSet]
        intersec = list(set(alias_data['entrez'].to_list()) & set(geneSet_list))
    
    alias_data3 = alias_data.set_index(alias_data['entrez'])
    geneSet_list = alias_data2.loc[intersec, 'entrez'].to_list()
    
    print(geneSet_list)
    print(len(geneSet),'total')
    print(len(geneSet_list),'intersec')

    ccgList = ccg.iloc[:, 1].tolist()
    intersection = list(set(geneSet_list) & set(ccgList))
    pvalue = hypergeom.sf(len(intersection)-1, 27214,
                          len(geneSet_list), len(ccgList))
    selected_pvalue_string = f"{pvalue:.2E}"
    size = len(geneSet_list)
    overlap = len(intersection)
    intersection_str = ';'.join(str(e) for e in intersection)
    entrezId = intersection_str
    symbol_str = alias_data3.loc[intersection, 'symbol'].to_list()
    symbol = ';'.join(str(e) for e in symbol_str)
    
    # size가 len(geneSet_list)
    gsize = len(geneSet_list)-5
    fdr005 = q3000.iloc[gsize, 1]
    fdr001 = q3000.iloc[gsize, 2]
    fdr0001 = q3000.iloc[gsize, 3]
    if pvalue < fdr0001:
        cancerLevel = 'Level 1 - Cancer pathway (highest confidence)'
    elif pvalue < fdr001:
        cancerLevel = 'Level 2 - Cancer pathway (high confidence)'
    elif pvalue < fdr005:
        cancerLevel = 'Level 3 - Cancer pathway (moderate confidence)'
    else:
        cancerLevel = 'Level 4 - Pathway with weak cancer association'

    # result table
    result_table = {
        'pvalue': pvalue,
        'pvalue_string': selected_pvalue_string,
        'size': size,
        'overlap': overlap,
        'entrezId': entrezId,
        'symbol': symbol,
        'cancerLevel': cancerLevel,
    }

    result_json = json.dumps(result_table)
    # return result
    return {
        "data": result_json,
        "matchInfo" : {
            "total" : len(geneSet),
            "intersec" : len(geneSet_list),
            "percent" : round((len(geneSet_list)/len(geneSet)) * 100)
        }
    }
