#if you have not yet installed cobra, uncomment the below line and
#run it in terminal (not in python console)
#pip install cobra
#do the rest of this script in python console. I used pycharm. 
#IMPORTANT: you might have to slightly change the pathing below to 
#make sure you are actually accessing the correct files.

#import libraries
import cobra
import pandas as pd
import sys
#convert models from sbml (which end in '.xml') to jsons (useable by escher)
exec(open('./python/convertsbmlstojsons.py').read())
#convert flux vector csvs to jsons
exec(open('./pythonsavefluxvectors.py').read())

#versions
for module in sys.modules:
    try:
        print(module,sys.modules[module].__version__)
    except:
        try:
            if  type(modules[module].version) is str:
                print(module,sys.modules[module].version)
            else:
                print(module,sys.modules[module].version())
        except:
            try:
                print(module,sys.modules[module].VERSION)
            except:
                pass
