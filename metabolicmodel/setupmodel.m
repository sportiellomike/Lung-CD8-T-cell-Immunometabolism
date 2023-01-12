%initCobraToolbox
%read this in, but I later saved it and loaded it so the original readCbModel is commented out. When you run it for the first time, you will need to uncomment it.
modelori=readCbModel('necessaryfiles/models/iMM1865.xml')
model=modelori

%read in this table produced by 'changereactionbounds.m' and 'checkintersectofboundsandmodelorirxns.R'
newbounds=readtable('necessaryfiles/exchangeboundsfromR.txt')

%grabrxns
listofrxnstochange = table2cell(newbounds(:,1));
listofrxnstochange = strtrim(listofrxnstochange);

%get lowerbounds
lowerbounds = table2array(newbounds(:,2));
%get upperbounds
upperbounds = table2array(newbounds(:,3));

%change reaction bounds
model = changeRxnBounds(model,listofrxnstochange,lowerbounds,'l');
model = changeRxnBounds(model,listofrxnstochange,upperbounds,'u');

FBAsolution = optimizeCbModel(model,'max')

%%%%%%%%%%%Proceed to next script%%%%%%%%%%%