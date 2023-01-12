% Read in gene expression data. this is reading in the basemean, lfc, stat, pvalue,padj,ensemblnames,fdr, genesymbol, and entrezid
dptpmtable=xlsread('./necessaryfiles/dp.genesymboltpm.xlsx')

% Convert gene ID lists to cell strings instead of numbers
%this is taking just the entrezid column
dpgenes = cellstr(num2str(dptpmtable(:,1)));
dpgenes = strtrim(dpgenes);
%dpgenes = strtrim(unique(dpgenes));
modelgenes=model.genes
modelgenes = cellstr((modelgenes));
modelgenes = strtrim(unique(modelgenes));

%get tpm
dptpm = dptpmtable(:,3);

% build expression data and map to reactions
expressionData=struct()
expressionData.gene=cellstr(dpgenes)
expressionData.value=dptpm
model.genes=cellstr(modelgenes)
[expressionRxns parsedGPR]=mapExpressionToReactions(model,expressionData)
%expressionRxns(isnan(expressionRxns))=-1;

%get mean expression value of 1,865 metabolism genes
test=expressionData
test.modelgenes=model.genes
intersect=ismember(test.gene,test.modelgenes)
test.intersect=intersect
index1=find(intersect == 1)

dptpmsubset=dptpm(index1)
dpgenessubset=dpgenes(index1)

meanmetabtpmdp=mean(dptpmsubset)
medianmetabtpmdp=median(dptpmsubset)
stdmetabtpmdp=std(dptpmsubset)
halfstdmetabtpmdp=stdmetabtpmdp*0.5

lowerthreshold=medianmetabtpmdp
upperthreshold=meanmetabtpmdp

options=struct()
options.expressionRxns = expressionRxns
options.solver = 'GIMME'
options.threshold=medianmetabtpmdp
options.runtime= 288000

%use GIMME
tissuemodeldp=createTissueSpecificModel(model,options)
FBAsolutiontissuedp = optimizeCbModel(tissuemodeldp,'max')

%make sure to uncomment the line below this. Very time consuming so it is commented out.
writeCbModel(tissuemodeldp, 'sbml','./models/tissuemodeldp.xml')

%look at table of fluxes for each reaction
tabledp=table(tissuemodeldp.rxns,tissuemodeldp.rxnNames,tissuemodeldp.subSystems,tissuemodeldp.rxnKEGGID,tissuemodeldp.rxnReactomeID,FBAsolutiontissuedp.v)
tabledp.Properties.VariableNames = {'rxns','rxnnames','subsystems','keggid','reactomeid','fbasolutiontissuev'}
writetable(tabledp,'tissuemodeldptable.txt','Delimiter','\t');

%make table for python script to turn table into jsons GIMME
vectortojsontabledp=table(tissuemodeldp.rxns, FBAsolutiontissuedp.v)
vectortojsontabledp.Properties.VariableNames = {'rxns','fluxvector'}
writetable(vectortojsontabledp,'./fluxvectors/vectortojsontabledp.csv','Delimiter','comma');

%%%%%%%%%%%Proceed to next script%%%%%%%%%%%