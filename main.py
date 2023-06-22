from typing import Union
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional
import oraAnalysis
import oraAnalysis_download
import cpAnalysis


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
    cancerLevel: Optional[int]


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
    dbType: str | None = None
    geneSetType: str | None = None
    geneSet: str | None = None
    qValueCutoff: str | None = None

@app.post("/ora_download")
async def downlaodORA(input: ORAInputDownload):
    data = oraAnalysis_download.analyze(
        input.dbType,
        input.geneSetType, 
        input.geneSet, 
        input.qValueCutoff)
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


@app.get("/test")
async def testAnalyze():
    return 'data'
