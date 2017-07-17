% BatchMapping.m
myCluster = parcluster('local');

if getenv('ENVIRONMENT')
   myCluster.JobStorageLocation = getenv('TMPDIR'); 
end

parpool(myCluster,3);

Animals = [84932,85362];
Dates = [20170717];

parfor ii=1:length(Animals)
    MapRetinotopy(Animals(ii),Dates);
end

delete(gcp);
