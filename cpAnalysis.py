import pandas as pd
from scipy.stats import hypergeom
import json
import re

q3000 = pd.read_csv('3000.csv')
ccg = pd.read_csv('ccg.csv')
alias_data = pd.read_csv('gdac_entrez.csv')


def analyze(geneSetType, geneSet):

    geneSet = re.sub(r"\s", "", geneSet)
    geneSet_list = geneSet.split(',')
    print(geneSet_list)

    if(geneSetType == 'int'):
        geneSet_list = [int(s) for s in geneSet_list]
        alias_data2 = alias_data.set_index(alias_data['entrez'])

    if(geneSetType == 'str'):
        alias_data2 = alias_data.set_index(alias_data['symbol'])

    if(geneSetType == 'int'):
        ccgList = ccg.iloc[:, 1].tolist()
        intersection = list(set(geneSet_list) & set(ccgList))
        pvalue = (round(hypergeom.sf(len(intersection)-1, 20501,
                                     len(geneSet_list), len(ccgList)), 9))
        size = len(geneSet_list)
        overlap = len(intersection)
        intersection_str = ';'.join(str(e) for e in intersection)
        entrezId = intersection_str
        symbol_str = alias_data2.loc[intersection, 'symbol'].to_list()
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

    if(geneSetType == 'str'):
        ccgList = ccg.iloc[:, 0].tolist()
        intersection = list(set(geneSet_list) & set(ccgList))
        pvalue = (round(hypergeom.sf(len(intersection)-1, 20501,
                                     len(geneSet_list), len(ccgList)), 10))
        size = len(geneSet_list)
        overlap = len(intersection)
        intersection_str = ';'.join(str(e) for e in intersection)
        symbol = intersection_str
        entrezId_str = alias_data2.loc[intersection, 'entrez'].to_list()
        entrezId = ';'.join(str(e) for e in entrezId_str)

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
        'size': size,
        'overlap': overlap,
        'entrezId': entrezId,
        'symbol': symbol,
        'cancerLevel': cancerLevel,
    }

    result_json = json.dumps(result_table)
    print(result_json)
    # return result
    return {
        "data": result_json,
    }
