%DailyMapping.m
%   Used along with a shell script to automatically run all the animals
%    from the current day, or all those whose data has not been analyzed
%
%INPUT: Will run through all of the files in the CloudStation folder
%        ByronExp/RetinoExp and see if any new files exist, if they do, it
%        will run MapRetinotopy.m on those
%OUTPUT: 
%
% Created: 2016/08/02, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/08/02
%  By: Byron Price 

cd('~/CloudStation/ByronExp/RetinoExp');

today = datetime('today','Format','yyyy-MM-dd');
today = char(today); today = strrep(today,'-','');
today = str2double(today);

fileStart = 'RetinoData*.plx';

fileList = dir(fileStart);
numFiles = size(fileList,1);

for ii=1:numFiles
    AnimalName = str2double(fileList(ii).name(end-8:end-4));
    Date = str2double(fileList(ii).name(end-17:end-10));
    if Date == today
        [~,~,~] = MapRetinotopy(AnimalName,Date);
    end
end
