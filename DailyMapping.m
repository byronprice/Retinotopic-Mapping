%DailyMapping.m
%   Used along with a shell script to automatically map retinotopy for
%    all the animals from the previous day
%
%INPUT: Will run through all of the files in the CloudStation folder
%        ByronExp/RetinoExp and see if any new files exist, if they do, it
%        will run MapRetinotopy.m on those
%OUTPUT: 
%
% Created: 2016/08/02, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/08/08
%  By: Byron Price 

cd('~/CloudStation/ByronExp/Retino');

fileStart = sprintf('RetinoData*.plx');

files = dir(fileStart);
numFiles = size(files,1);

for ii=1:numFiles
    index = regexp(files(ii).name,'_');
    Date = str2double(files(ii).name(index-8:index-1));
    AnimalName = str2double(files(ii).name(index+1:end-4));
    
    fileCheck = sprintf('RetinoMapBayes*_%d',AnimalName);
    checkFiles = dir(fileCheck);
    
    if isempty(checkFiles)==1
        MapRetinotopy(AnimalName,Date);
    end
end

fprintf('Ran cron job\n');
