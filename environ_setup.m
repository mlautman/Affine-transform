%% addpaths
if isunix()
    addpath(genpath('./maxflow/'));
    addpath(genpath('./solution_GC_2D/'));
    addpath(genpath('./data/'));
    addpath(genpath('./mrfcode/'));

    cd maxflow/
    make 
    cd ..
else 
    addpath(genpath('.\maxflow\'));
    addpath(genpath('.\solution_GC_2D\'));
    addpath(genpath('.\data\'));
    addpath(genpath('.\mrfcode\'));
end    

