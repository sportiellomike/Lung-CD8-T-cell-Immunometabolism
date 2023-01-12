% Read in gene expression data. this is reading in the basemean, lfc, stat, pvalue,padj,ensemblnames,fdr, genesymbol, and entrezid
cd49asptpmtable=xlsread('./necessaryfiles/cd49asp.genesymboltpm.xlsx')

% Convert gene ID lists to cell strings instead of numbers
%this is taking just the entrezid column
cd49aspgenes = cellstr(num2str(cd49asptpmtable(:,1)));
cd49aspgenes = strtrim(cd49aspgenes);
%cd49aspgenes = strtrim(unique(cd49aspgenes));
modelgenes=model.genes
modelgenes = cellstr((modelgenes));
modelgenes = strtrim(unique(modelgenes));

%get tpm
cd49asptpm = cd49asptpmtable(:,3);

% build expression data and map to reactions
expressionData=struct()
expressionData.gene=cellstr(cd49aspgenes)
expressionData.value=cd49asptpm
model.genes=cellstr(modelgenes)
[expressionRxns parsedGPR]=mapExpressionToReactions(model,expressionData)
%expressionRxns(isnan(expressionRxns))=-1;

%get mean expression value of 1,865 metabolism genes
test=expressionData
test.modelgenes=model.genes
intersect=ismember(test.gene,test.modelgenes)
test.intersect=intersect
index1=find(intersect == 1)

cd49asptpmsubset=cd49asptpm(index1)
cd49aspgenessubset=cd49aspgenes(index1)

meanmetabtpmcd49asp=mean(cd49asptpmsubset)
medianmetabtpmcd49asp=median(cd49asptpmsubset)
stdmetabtpmcd49asp=std(cd49asptpmsubset)
halfstdmetabtpmcd49asp=stdmetabtpmcd49asp*0.5

lowerthreshold=medianmetabtpmcd49asp
upperthreshold=meanmetabtpmcd49asp

options=struct()
options.expressionRxns = expressionRxns
options.solver = 'GIMME'
options.threshold=medianmetabtpmcd49asp
options.runtime= 288000

%use GIMME
tissuemodelcd49asp=createTissueSpecificModel(model,options)
FBAsolutiontissuecd49asp = optimizeCbModel(tissuemodelcd49asp,'max')

%make sure to uncomment the line below this. Very time consuming so it is commented out.
writeCbModel(tissuemodelcd49asp, 'sbml','./models/tissuemodelcd49asp.xml')

%look at table of fluxes for each reaction
tablecd49asp=table(tissuemodelcd49asp.rxns,tissuemodelcd49asp.rxnNames,tissuemodelcd49asp.subSystems,tissuemodelcd49asp.rxnKEGGID,tissuemodelcd49asp.rxnReactomeID,FBAsolutiontissuecd49asp.v)
tablecd49asp.Properties.VariableNames = {'rxns','rxnnames','subsystems','keggid','reactomeid','fbasolutiontissuev'}
writetable(tablecd49asp,'tissuemodelcd49asptable.txt','Delimiter','\t');

%make table for python script to turn table into jsons GIMME
vectortojsontablecd49asp=table(tissuemodelcd49asp.rxns, FBAsolutiontissuecd49asp.v)
vectortojsontablecd49asp.Properties.VariableNames = {'rxns','fluxvector'}
writetable(vectortojsontablecd49asp,'./fluxvectors/vectortojsontablecd49asp.csv','Delimiter','comma');

%%%%%%%%%%%Proceed to next script%%%%%%%%%%%