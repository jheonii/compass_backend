from typing import Union
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional
from typing import List
import oraAnalysis
import oraAnalysis_download
import cpAnalysis
import gpn
import pandas as pd


app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# 1. ORA analysis
class ORAInput(BaseModel):
    dbType: Optional[str]
    geneSetType: Optional[str]
    geneSet: Optional[str]
    qValueCutoff: Optional[str]
    cancerLevel: Optional[str]


@app.post("/ora")
async def analyzeORA(input: ORAInput):
    dicted_input = dict(input)
    dbType = dicted_input['dbType']
    geneSetType = dicted_input['geneSetType']
    geneSet = dicted_input['geneSet']
    qValueCutoff = dicted_input['qValueCutoff']
    cancerLevel = dicted_input['cancerLevel']
    data = oraAnalysis.analyze(
        dbType, geneSetType, geneSet, qValueCutoff, cancerLevel)
    return data


# 1-download ORA analysis
class ORAInputDownload(BaseModel):
    dbType: Optional[str]
    geneSetType: Optional[str]
    geneSet: Optional[str]
    qValueCutoff: Optional[str]
    cancerLevel: Optional[str]

@app.post("/ora_download")
async def downlaodORA(input: ORAInputDownload):
    dicted_input = dict(input)
    dbType = dicted_input['dbType']
    geneSetType = dicted_input['geneSetType']
    geneSet = dicted_input['geneSet']
    qValueCutoff = dicted_input['qValueCutoff']
    cancerLevel = dicted_input['cancerLevel']
    
    data = oraAnalysis_download.analyze(
        dbType,
        geneSetType, 
        geneSet, 
        qValueCutoff,
        cancerLevel)
    return data

# 2. Cancer prioritization


class CPInput(BaseModel):
    geneSetType: Optional[str]
    geneSet: Optional[str]


@app.post("/cp")
async def analyzeORA(input: CPInput):
    dicted_input = dict(input)
    geneSetType = dicted_input['geneSetType']
    geneSet = dicted_input['geneSet']
    data = cpAnalysis.analyze(geneSetType, geneSet)
    return data


# 3. Gene Pathway Network

class getORAInput(BaseModel):
    pathway : str
    overlapGenes : str
    

@app.post('/network', status_code=200)
def getGenePathwayNetworkRes(_input: List[getORAInput], pathwayN: int):
    _dict = [s.dict() for s in _input]
    df = pd.DataFrame(_dict)
    result = gpn.getGenePathwayNetwork(df, pathwayN)
    gene_pathway_network = result['gpd_df2']
    gene_pathway_table_list= result['gpt_list']
    return {
        "gene_pathway_network" : gene_pathway_network.transpose().to_json(),
        "gene_pathway_table_list" : gene_pathway_table_list.to_json()
    }