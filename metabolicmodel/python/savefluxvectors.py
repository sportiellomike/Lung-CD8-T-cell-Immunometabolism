# this script is to turn the flux vector fromfbasolution to json
# https://escher.readthedocs.io/en/latest/getting_started.html

# importing pandas and json modules
import pandas as pd
import json

# reading csv file for DP
data = pd.read_csv("./fluxvectors/vectortojsontabledp.csv")
# dropping null value columns to avoid errors
data.dropna(inplace = True)
# converting to dict
fluxvectordp=dict(zip(data.rxns,data.fluxvector))
# save a single flux vector as JSON
flux_dictionary = fluxvectordp
with open('./fluxvectors/fluxvectordp.json', 'w') as f:
    json.dump(flux_dictionary, f)

#now do the same for DN
data = pd.read_csv("./fluxvectors/vectortojsontabledn.csv")
data.dropna(inplace = True)
fluxvectordp=dict(zip(data.rxns,data.fluxvector))
flux_dictionary = fluxvectordp
with open('./fluxvectors/fluxvectordn.json', 'w') as f:
    json.dump(flux_dictionary, f)