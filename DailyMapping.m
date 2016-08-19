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
% Updated: 2016/08/18
%  By: Byron Price 

cd('~/CloudStation/ByronExp/Retino');

yest = datetime('yesterday','Format','yyyy-MM-dd');
yest = char(yest); yest = strrep(yest,'-','');
yest = str2double(yest);

fileStart = sprintf('RetinoData*%d*.plx',yest);

fileList = dir(fileStart);
numFiles = size(fileList,1);

datelen = 8;
idlen = 5;
for ii=1:numFiles
    index = regexp(fileList(ii).name,'_');
    Date = str2double(fileList(ii).name(index-datelen:index-1));
    AnimalName = str2double(fileList(ii).name(index+1:index+idlen));
    [~,~,~] = MapRetinotopy(AnimalName,Date);
end
