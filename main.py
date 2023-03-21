from typing import Union
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional
import oraAnalysis
import oraAnalysis2
import cpAnalysis
import test

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# 1-(1). ORA analysis - genesetType === int
class ORAInput(BaseModel):
    geneSetType: Optional[str]
    geneSet: Optional[str]
    qValueCutoff: Optional[float]
    cancerLevel: Optional[int]


@app.post("/ora")
async def analyzeORA(input: ORAInput):
    dicted_input = dict(input)
    geneSetType = dicted_input['geneSetType']
    geneSet = dicted_input['geneSet']
    qValueCutoff = dicted_input['qValueCutoff']
    cancerLevel = dicted_input['cancerLevel']
    data = oraAnalysis.analyze(geneSetType, geneSet, qValueCutoff, cancerLevel)
    return data

# 1-(2). ORA analysis - genesetType === str


class ORAInput2(BaseModel):
    geneSetType: Optional[str]
    geneSet: Optional[str]
    qValueCutoff: Optional[float]
    cancerLevel: Optional[int]


@app.post("/ora2")
async def analyzeORA2(input: ORAInput2):
    dicted_input = dict(input)
    geneSetType = dicted_input['geneSetType']
    geneSet = dicted_input['geneSet']
    qValueCutoff = dicted_input['qValueCutoff']
    cancerLevel = dicted_input['cancerLevel']
    data = oraAnalysis2.analyze(
        geneSetType, geneSet, qValueCutoff, cancerLevel)
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
    data = test.analyze()
    return data
