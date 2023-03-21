import pandas as pd


def analyze():

    GSEC = pd.read_csv('GSEC.csv')
    print(GSEC.head(10))

    return {
        "data": 'test'
    }
