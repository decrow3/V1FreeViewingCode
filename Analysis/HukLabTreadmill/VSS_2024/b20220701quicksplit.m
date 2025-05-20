HUKDATA = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');
Exp.FileTag = 'brie/b20220701_probe2.mat';
Exp=load(fullfile(HUKDATA, Exp.FileTag));


V1fovsub=Exp.osp.cids(Exp.osp.clusterDepths<1000);
subset=V1fovsub;
nsub = numel(subset);

%%
Exp0=Exp;

%% THESE depths seem reversed??? Could be we are 'too deep'
%% V1foveal
Exp=Exp0;
fovbool=((Exp.osp.clusterDepths<1000));
V1fovsub=Exp.osp.cids(fovbool);
Exp.osp.clusterAmps=Exp.osp.clusterAmps(fovbool);
Exp.osp.clusterDepths=Exp.osp.clusterDepths(fovbool);
Exp.osp.firingRates=Exp.osp.firingRates(fovbool);

Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1fovsub));
Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1fovsub));

Exp.osp.cids = unique(Exp.osp.clu); 

Exp.FileTag = 'brie/b20220701_probe2_fov2.mat';
save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220701_probe2_fov2.mat', '-v7.3', '-struct', 'Exp')
%% V1calc

Exp=Exp0;
calcbool=((Exp.osp.clusterDepths>2000)&(Exp.osp.clusterDepths<4000));
V1calcsub=Exp.osp.cids(calcbool);
Exp.osp.clusterAmps=Exp.osp.clusterAmps(calcbool);
Exp.osp.clusterDepths=Exp.osp.clusterDepths(calcbool);
Exp.osp.firingRates=Exp.osp.firingRates(calcbool);

Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1calcsub));
Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1calcsub));


Exp.osp.cids = unique(Exp.osp.clu); 

Exp.spikeIds = Exp.osp.clu;
Exp.spikeTimes = Exp.osp.st;
Exp.FileTag = 'brie/b20220701_probe2_calc2.mat';
save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220701_probe2_calc2.mat', '-v7.3', '-struct', 'Exp')


%% V1calc

Exp=Exp0;
calcbool=((Exp.osp.clusterDepths>4000)&(Exp.osp.clusterDepths<5500));
V1calcsub=Exp.osp.cids(calcbool);
Exp.osp.clusterAmps=Exp.osp.clusterAmps(calcbool);
Exp.osp.clusterDepths=Exp.osp.clusterDepths(calcbool);
Exp.osp.firingRates=Exp.osp.firingRates(calcbool);

Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1calcsub));
Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1calcsub));

Exp.osp.cids = unique(Exp.osp.clu); 

Exp.spikeIds = Exp.osp.clu;
Exp.spikeTimes = Exp.osp.st;
Exp.FileTag = 'brie/b20220701_probe2_calc1.mat';
save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220701_probe2_calc1.mat', '-v7.3', '-struct', 'Exp')


%% V1foveal
Exp=Exp0;
fovbool=((Exp.osp.clusterDepths>5500));
V1fovsub=Exp.osp.cids(fovbool);
Exp.osp.clusterAmps=Exp.osp.clusterAmps(fovbool);
Exp.osp.clusterDepths=Exp.osp.clusterDepths(fovbool);
Exp.osp.firingRates=Exp.osp.firingRates(fovbool);

Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1fovsub));
Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1fovsub));

Exp.osp.cids = unique(Exp.osp.clu); 

Exp.FileTag = 'brie/b20220701_probe2_fov1.mat';
save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220701_probe2_fov1.mat', '-v7.3', '-struct', 'Exp')