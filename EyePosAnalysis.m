function [] = EyePosAnalysis(AnimalName,Date)
%EyePosition.m
%  Display a series of flashing gabors to determine where the mouse is
%   looking
%  Each circle will occupy an ~ 2.5-degree radius of visual space
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%        Date - date of experiment
%
% OUTPUT: a file with results named EyePosResultDate_AnimalName.mat and
%          some figures
% Created: 2018/08/22 at 24 Cummington, Boston
%  Byron Price
% Updated: 2018/08/22
%  By: Byron Price

eyeTrackFile = sprintf('EyePosition_%d-%-MLP.mat',Date,AnimalName);
load(eyeTrackFile,'pupilDiameter','pupilRotation','blink','Fs');
eyeFs = Fs;

ephysDate = sprintf('%d',Date);
ephysDate = [ephysDate(1:4),'-',ephysDate(5:6),'-',ephysDate(7:8)];
ephysFile = dir(sprintf('CompiledData*%s*_%d.mat',ephysDate,AnimalName));

load(ephysFile(1).name,'auxData','events','eventTimes','lowpassData',...
    'lowpassTimes','lpFs');

stimFile = sprintf('EyePosStim%d_%d.mat',Date,AnimalName);
load(stimFile,'centerVals','stimStrobeNums','DistToScreen','mmPerPixel',...
    'blocks','directions','numStimuli');

directionDegrees = linspace(0,360,directions+1);
directionDegrees = directionDegrees(1:end-1);

vidInds = find(auxData(:,2)==1);

avgTime = 2;
posResult = zeros(directions,blocks,eyeFs*avgTime,2);
diamResult = zeros(directions,blocks,eyeFs*avgTime);

try
    cm = magma(100);
catch
    cm = 'jet';
end
for ii=1:directions
    times = eventTimes(events==stimStrobeNums(ii,1));
    count = 0;
    figure;
    for kk=1:2:length(times)
        count = count+1;
        [~,ind] = min(abs(times(kk)-lowpassTimes));
        [~,frameInd] = min(abs(ind-vidInds));
        diamResult(ii,count,:) = pupilDiameter(frameInd-1:frameInd+eyeFs*avgTime);
        posResult(ii,count,:,:) = pupilRotation(frameInd-1:frameInd+eyeFs*avgTime,:);
        
        colorscale = linspace(0,avgTime,eyeFs*avgTime);
        scatter(squeeze(posResult(ii,count,:,1)),squeeze(posResult(ii,count,:,2)),...
            25,colorscale,'filled');colormap(cm);hold on;
    end
    title(sprintf('Stimulus-Triggered Eye Position, %d Degrees',directionDegrees(ii)));
    ylabel('Vertical Eye Position (pixels)');
    xlabel('Horizontal Eye Position (pixels)');
end

saveFile = sprintf('EyePosResult%d_%d.mat',Date,AnimalName);
save(saveFile,'diamResult','posResult','directions','blocks','numStimuli',...
    'centerVals','DistToScreen','mmPerPixel','directionDegrees');

end