HUKDATA = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');
Exp.FileTag = 'brie/b20220625_probe1.mat';
Exp=load(fullfile(HUKDATA, Exp.FileTag));



%%
Exp0=Exp;


%% V1foveal
Exp=Exp0;
fovbool=((Exp.osp.clusterDepths>5110));
V1fovsub=Exp.osp.cids(fovbool);
Exp.osp.clusterAmps=Exp.osp.clusterAmps(fovbool);
Exp.osp.clusterDepths=Exp.osp.clusterDepths(fovbool);
Exp.osp.firingRates=Exp.osp.firingRates(fovbool);

Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1fovsub));
Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1fovsub));

Exp.osp.cids = unique(Exp.osp.clu); 

Exp.FileTag = 'brie/b20220625_probe1_fov.mat';
save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220625_probe1_fov.mat', '-v7.3', '-struct', 'Exp')

%% V1calc shallow

Exp=Exp0;
calcbool=((Exp.osp.clusterDepths>4417)&(Exp.osp.clusterDepths<5110));
V1calcsub=Exp.osp.cids(calcbool);
Exp.osp.clusterAmps=Exp.osp.clusterAmps(calcbool);
Exp.osp.clusterDepths=Exp.osp.clusterDepths(calcbool);
Exp.osp.firingRates=Exp.osp.firingRates(calcbool);

Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1calcsub));
Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1calcsub));

Exp.osp.cids = unique(Exp.osp.clu); 

Exp.spikeIds = Exp.osp.clu;
Exp.spikeTimes = Exp.osp.st;
Exp.FileTag = 'brie/b20220625_probe1_calc1.mat';
save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220625_probe1_calc1.mat', '-v7.3', '-struct', 'Exp')


%% V1calc

Exp=Exp0;
calcbool=((Exp.osp.clusterDepths>430)&(Exp.osp.clusterDepths< 4417));
V1calcsub=Exp.osp.cids(calcbool);
Exp.osp.clusterAmps=Exp.osp.clusterAmps(calcbool);
Exp.osp.clusterDepths=Exp.osp.clusterDepths(calcbool);
Exp.osp.firingRates=Exp.osp.firingRates(calcbool);

Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1calcsub));
Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1calcsub));


Exp.osp.cids = unique(Exp.osp.clu); 

Exp.spikeIds = Exp.osp.clu;
Exp.spikeTimes = Exp.osp.st;
Exp.FileTag = 'brie/b20220625_probe1_calc2.mat';
save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220625_probe1_calc2.mat', '-v7.3', '-struct', 'Exp')



%% V?fovealdeep
Exp=Exp0;
fovbool=((Exp.osp.clusterDepths<430));
V1fovsub=Exp.osp.cids(fovbool);
Exp.osp.clusterAmps=Exp.osp.clusterAmps(fovbool);
Exp.osp.clusterDepths=Exp.osp.clusterDepths(fovbool);
Exp.osp.firingRates=Exp.osp.firingRates(fovbool);

Exp.osp.st  = Exp.osp.st(ismember(Exp.osp.clu,V1fovsub));
Exp.osp.clu = Exp.osp.clu(ismember(Exp.osp.clu,V1fovsub));

Exp.osp.cids = unique(Exp.osp.clu); 

Exp.spikeIds = Exp.osp.clu;
Exp.spikeTimes = Exp.osp.st;
Exp.FileTag = 'brie/b20220625_probe1_fov2.mat';
save('/media/huklab/Data/NPX/HuklabTreadmill/VSS/brie_20220625_probe1_fov2.mat', '-v7.3', '-struct', 'Exp')
