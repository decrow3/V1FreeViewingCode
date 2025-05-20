Exp=load('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode/b20220809re2.mat');


%%
Exp0=Exp;


%% MTfoveal
Exp=Exp0;
fovbool=((Exp.osp.clusterDepths>5000));
V1fovsub=Exp.osp.cids(fovbool);
Exp.osp.clusterAmps=Exp.osp.clusterAmps(fovbool);
Exp.osp.clusterDepths=Exp.osp.clusterDepths(fovbool);
Exp.osp.firingRates=Exp.osp.firingRates(fovbool);

Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1fovsub));
Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1fovsub));

Exp.osp.cids = unique(Exp.osp.clu); 

Exp.FileTag = 'brie/b20220809_MTshallow.mat';
save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220809_MTshallow.mat', '-v7.3', '-struct', 'Exp')

 %% V1calc shallow
% 
% Exp=Exp0;
% calcbool=((Exp.osp.clusterDepths>3250)&(Exp.osp.clusterDepths<5300));
% V1calcsub=Exp.osp.cids(calcbool);
% Exp.osp.clusterAmps=Exp.osp.clusterAmps(calcbool);
% Exp.osp.clusterDepths=Exp.osp.clusterDepths(calcbool);
% Exp.osp.firingRates=Exp.osp.firingRates(calcbool);
% 
% Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1calcsub));
% Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1calcsub));
% 
% Exp.osp.cids = unique(Exp.osp.clu); 
% 
% Exp.spikeIds = Exp.osp.clu;
% Exp.spikeTimes = Exp.osp.st;
% Exp.FileTag = 'brie/b20220625_probe2_calc1.mat';
% save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220625_probe2_calc1.mat', '-v7.3', '-struct', 'Exp')
% 
% 
% %% V1calc
% 
% Exp=Exp0;
% calcbool=((Exp.osp.clusterDepths>1500)&(Exp.osp.clusterDepths< 3250));
% V1calcsub=Exp.osp.cids(calcbool);
% Exp.osp.clusterAmps=Exp.osp.clusterAmps(calcbool);
% Exp.osp.clusterDepths=Exp.osp.clusterDepths(calcbool);
% Exp.osp.firingRates=Exp.osp.firingRates(calcbool);
% 
% Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1calcsub));
% Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1calcsub));
% 
% 
% Exp.osp.cids = unique(Exp.osp.clu); 
% 
% Exp.spikeIds = Exp.osp.clu;
% Exp.spikeTimes = Exp.osp.st;
% Exp.FileTag = 'brie/b20220625_probe2_calc2.mat';
% save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220625_probe2_calc2.mat', '-v7.3', '-struct', 'Exp')
% 
% 
% 
% %% V?fovealdeep
% Exp=Exp0;
% fovbool=((Exp.osp.clusterDepths<1500));
% V1fovsub=Exp.osp.cids(fovbool);
% Exp.osp.clusterAmps=Exp.osp.clusterAmps(fovbool);
% Exp.osp.clusterDepths=Exp.osp.clusterDepths(fovbool);
% Exp.osp.firingRates=Exp.osp.firingRates(fovbool);
% 
% Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1fovsub));
% Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1fovsub));
% 
% Exp.osp.cids = unique(Exp.osp.clu); 
% 
% Exp.spikeIds = Exp.osp.clu;
% Exp.spikeTimes = Exp.osp.st;
% Exp.FileTag = 'brie/b20220625_probe2_fov2.mat';
% save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220625_probe2_fov2.mat', '-v7.3', '-struct', 'Exp')
