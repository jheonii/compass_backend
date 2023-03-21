import pandas as pd


def analyze():

    GSEC = pd.read_csv('GSEC.csv')
    result_data = GSEC.head(1000)
    print(GSEC.head(10))
    result_json = result_data.transpose().to_json()
    return {
        "data": result_json
    }
