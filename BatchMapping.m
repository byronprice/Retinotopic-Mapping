% BatchMapping.m
myCluster = parcluster('local');

if getenv('ENVIRONMENT')
   myCluster.JobStorageLocation = getenv('TMPDIR'); 
end

parpool(myCluster,4);

Animals = [34271,34272,43650,43652,43653,62500,62501,62502];
Dates = [20170607,20170607];

parfor ii=1:legnth(Animals)
    MapRetinotopy(Animals(ii),Dates(ii));
end

delete(gcp);