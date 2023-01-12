% Read in gene expression data. this is reading in the basemean, lfc, stat, pvalue,padj,ensemblnames,fdr, genesymbol, and entrezid
cd103sptpmtable=xlsread('./necessaryfiles/cd103sp.genesymboltpm.xlsx')

% Convert gene ID lists to cell strings instead of numbers
%this is taking just the entrezid column
cd103spgenes = cellstr(num2str(cd103sptpmtable(:,1)));
cd103spgenes = strtrim(cd103spgenes);
%cd103spgenes = strtrim(unique(cd103spgenes));
modelgenes=model.genes
modelgenes = cellstr((modelgenes));
modelgenes = strtrim(unique(modelgenes));

%get tpm
cd103sptpm = cd103sptpmtable(:,3);

% build expression data and map to reactions
expressionData=struct()
expressionData.gene=cellstr(cd103spgenes)
expressionData.value=cd103sptpm
model.genes=cellstr(modelgenes)
[expressionRxns parsedGPR]=mapExpressionToReactions(model,expressionData)
%expressionRxns(isnan(expressionRxns))=-1;

%get mean expression value of 1,865 metabolism genes
test=expressionData
test.modelgenes=model.genes
intersect=ismember(test.gene,test.modelgenes)
test.intersect=intersect
index1=find(intersect == 1)

cd103sptpmsubset=cd103sptpm(index1)
cd103spgenessubset=cd103spgenes(index1)

meanmetabtpmcd103sp=mean(cd103sptpmsubset)
medianmetabtpmcd103sp=median(cd103sptpmsubset)
stdmetabtpmcd103sp=std(cd103sptpmsubset)
halfstdmetabtpmcd103sp=stdmetabtpmcd103sp*0.5

lowerthreshold=medianmetabtpmcd103sp
upperthreshold=meanmetabtpmcd103sp

options=struct()
options.expressionRxns = expressionRxns
options.solver = 'GIMME'
options.threshold=medianmetabtpmcd103sp
options.runtime= 288000

%use GIMME
tissuemodelcd103sp=createTissueSpecificModel(model,options)
FBAsolutiontissuecd103sp = optimizeCbModel(tissuemodelcd103sp,'max')

%make sure to uncomment the line below this. Very time consuming so it is commented out.
writeCbModel(tissuemodelcd103sp, 'sbml','./models/tissuemodelcd103sp.xml')

%look at table of fluxes for each reaction
tablecd103sp=table(tissuemodelcd103sp.rxns,tissuemodelcd103sp.rxnNames,tissuemodelcd103sp.subSystems,tissuemodelcd103sp.rxnKEGGID,tissuemodelcd103sp.rxnReactomeID,FBAsolutiontissuecd103sp.v)
tablecd103sp.Properties.VariableNames = {'rxns','rxnnames','subsystems','keggid','reactomeid','fbasolutiontissuev'}
writetable(tablecd103sp,'tissuemodelcd103sptable.txt','Delimiter','\t');

%make table for python script to turn table into jsons GIMME
vectortojsontablecd103sp=table(tissuemodelcd103sp.rxns, FBAsolutiontissuecd103sp.v)
vectortojsontablecd103sp.Properties.VariableNames = {'rxns','fluxvector'}
writetable(vectortojsontablecd103sp,'./fluxvectors/vectortojsontablecd103sp.csv','Delimiter','comma');

%%%%%%%%%%%Proceed to next script%%%%%%%%%%%