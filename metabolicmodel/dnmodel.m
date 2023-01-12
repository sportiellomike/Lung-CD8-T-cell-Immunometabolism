% Read in gene expression data. this is reading in the basemean, lfc, stat, pvalue,padj,ensemblnames,fdr, genesymbol, and entrezid
dntpmtable=xlsread('./necessaryfiles/dn.genesymboltpm.xlsx')

% Convert gene ID lists to cell strings instead of numbers
%this is taking just the entrezid column
dngenes = cellstr(num2str(dntpmtable(:,1)));
dngenes = strtrim(dngenes);
%dngenes = strtrim(unique(dngenes));
modelgenes=model.genes
modelgenes = cellstr((modelgenes));
modelgenes = strtrim(unique(modelgenes));

%get tpm
dntpm = dntpmtable(:,3);

% build expression data and map to reactions
expressionData=struct()
expressionData.gene=cellstr(dngenes)
expressionData.value=dntpm
model.genes=cellstr(modelgenes)
[expressionRxns parsedGPR]=mapExpressionToReactions(model,expressionData)
%expressionRxns(isnan(expressionRxns))=-1;

%get mean expression value of 1,865 metabolism genes
test=expressionData
test.modelgenes=model.genes
intersect=ismember(test.gene,test.modelgenes)
test.intersect=intersect
index1=find(intersect == 1)

dntpmsubset=dntpm(index1)
dngenessubset=dngenes(index1)

meanmetabtpmdn=mean(dntpmsubset)
medianmetabtpmdn=median(dntpmsubset)
stdmetabtpmdn=std(dntpmsubset)
halfstdmetabtpmdn=stdmetabtpmdn*0.5

lowerthreshold=medianmetabtpmdn
upperthreshold=meanmetabtpmdn

options=struct()
options.expressionRxns = expressionRxns
options.solver = 'GIMME'
options.threshold=medianmetabtpmdn
options.runtime= 288000

%use GIMME
tissuemodeldn=createTissueSpecificModel(model,options)
FBAsolutiontissuedn = optimizeCbModel(tissuemodeldn,'max')
%make sure to uncomment the line below this. Very time consuming so it is commented out.
writeCbModel(tissuemodeldn, 'sbml','./models/tissuemodeldn.xml')

%look at table of fluxes for each reaction
tabledn=table(tissuemodeldn.rxns,tissuemodeldn.rxnNames,tissuemodeldn.subSystems,tissuemodeldn.rxnKEGGID,tissuemodeldn.rxnReactomeID,FBAsolutiontissuedn.v)
tabledn.Properties.VariableNames = {'rxns','rxnnames','subsystems','keggid','reactomeid','fbasolutiontissuev'}
writetable(tabledn,'tissuemodeldntable.txt','Delimiter','\t');

%make table for python script to turn table into jsons GIMME
vectortojsontabledn=table(tissuemodeldn.rxns, FBAsolutiontissuedn.v)
vectortojsontabledn.Properties.VariableNames = {'rxns','fluxvector'}
writetable(vectortojsontabledn,'./fluxvectors/vectortojsontabledn.csv','Delimiter','comma');

%%%%%%%%%%%Proceed to next script%%%%%%%%%%%