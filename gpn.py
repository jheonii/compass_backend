import pandas as pd

#Symbol_entrez
symbol_entrez_df = pd.read_csv('gdac_entrez.csv')
symbol_entrez_df = symbol_entrez_df.astype(str)
mapping_dict = symbol_entrez_df.set_index('entrez')['symbol'].to_dict()


def getGenePathwayTable(df):
  df['overlapGenes'] = df['overlapGenes'].str.split(',')
  newdf = df.explode('overlapGenes')
  newdf['overlapSymbol'] = newdf['overlapGenes'].apply(lambda x: ','.join(mapping_dict[i] for i in x.split(',')))
  newdf2 = newdf.loc[:,['pathway', 'overlapGenes', 'overlapSymbol']]
  newdf3=newdf2.drop_duplicates()
  
  df_grouped_entrez = newdf3.groupby('overlapGenes')['pathway'].apply(list).reset_index(name='pathways')
  df_grouped_entrez['pathwayN'] = df_grouped_entrez['pathways'].str.len()
  df_grouped_entrez = df_grouped_entrez.sort_values('pathwayN')
  df_grouped_entrez.loc[:,'pathways'] = [','.join(s) for s in df_grouped_entrez.loc[:,'pathways'].values]
  df_grouped_entrez = df_grouped_entrez.groupby('pathwayN')['overlapGenes'].apply(', '.join).reset_index(name='entrez')
  
  df_grouped_symbol = newdf3.groupby('overlapSymbol')['pathway'].apply(list).reset_index(name='pathways')
  df_grouped_symbol['pathwayN'] = df_grouped_symbol['pathways'].str.len()
  df_grouped_symbol = df_grouped_symbol.sort_values('pathwayN')
  df_grouped_symbol.loc[:,'pathways'] = [','.join(s) for s in df_grouped_symbol.loc[:,'pathways'].values]
  df_grouped_symbol = df_grouped_symbol.groupby('pathwayN')['overlapSymbol'].apply(', '.join).reset_index(name='symbol')
  
  df_result = pd.merge(df_grouped_symbol, df_grouped_entrez, on='pathwayN')
  return {
    'df_result' :df_result,
    'newdf3' : newdf3
    }



def getGenePathwayNetwork(df, pathwayN):  
  result = getGenePathwayTable(df)
  gene_pathway_table = result['df_result']
  filtered_gpt = gene_pathway_table[gene_pathway_table['pathwayN'] == pathwayN]
  entrez_list = gene_pathway_table.loc[:,'entrez'].values
  nested_list = [s.split(', ') for s in entrez_list]
  flattened_list = [item for sublist in nested_list for item in sublist]
  
  newdf3=result['newdf3']
  gpd_df = newdf3[newdf3['overlapGenes'].isin(flattened_list)]
  gpd_df2 = gpd_df.loc[:,['pathway', 'overlapSymbol']]
  gpd_df2.columns = ['pathway', 'gene']
  gpd_df2 = gpd_df2.reset_index(drop=True)
  return {
    'gpd_df2' : gpd_df2,
    'gpt_list' : gene_pathway_table['pathwayN']
  }

