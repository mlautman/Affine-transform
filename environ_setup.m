%% addpaths
addpath(genpath('./maxflow/'));
addpath(genpath('./solution_GC_2D/'));
addpath(genpath('./data/'));
addpath(genpath('./mrfcode/'));

%% Mex compilation
cd maxflow/
make 
cd ..

