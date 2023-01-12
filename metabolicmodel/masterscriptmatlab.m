%Read the damn ReadMe and do what it says

%create necessary directories
status = mkdir('models')
status = mkdir('fluxvectors')
status = mkdir('fluxvectorsums')

%%%matlab stuff%%%
%run setup
setupmodel

%this script will load model for CD49a+CD103+ (Double Positive=DP) T cells
dpmodel

%this script will load model for CD49a-CD103- (Double Positive=DN) T cells
dnmodel

%this script will load model for CD49a+CD103- (CD49a Single Positive=CD49a SP) T cells
cd49aspmodel

%this script will load model for CD49a-CD103+ (CD103 Single Positive=CD103 SP) T cells
cd103spmodel

%%%Python Stuff
%run the masterscriptpy.py to get the flux vector jsons and map jsons for use with escher

%%%R stuff
%run masterscriptR.R to recreate data column plot in figure 2B