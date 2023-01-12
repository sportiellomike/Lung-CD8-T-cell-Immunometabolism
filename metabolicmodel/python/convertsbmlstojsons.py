#configure cobra
cobra_config=cobra.Configuration()
#path in the models saved from matlab
imm1865='./models/iMM1865.xml'
tissuemodeldp='./models/tissuemodeldp.xml'
tissuemodeldn='./models/tissuemodeldn.xml'

#read in these models
imm1865=cobra.io.read_sbml_model(imm1865)
tissuemodeldp=cobra.io.read_sbml_model(tissuemodeldp)
tissuemodeldn=cobra.io.read_sbml_model(tissuemodeldn)

#save them as jsons. Note: this may take several minutes
cobra.io.save_json_model(tissuemodeldp,'./models/tissuemodeldp.json')
cobra.io.save_json_model(tissuemodeldn,'./models/tissuemodeldn.json')
cobra.io.save_json_model(imm1865,'./models/imm1865.json')