function [h] = MakeRetinoCallMaps(AnimalName,Date)
%MakeRetinoMaps.m
%  Plot the data from an LFP retinotopic mapping experiment given the animal's
%   name and the date of the experiment (see RetinotopyCallaway.m for the
%   stimulus)
%
%Created: 2017/10/25, 24 Cummington Mall, Boston
%Updated: 2017/10/25
% By: Byron Price

fileName = sprintf('RetinoCallResults%d_%d.mat',Date,AnimalName);
load(fileName);

DistToScreen = DistToScreen*10;

numChans = size(Response,1);

% for ii=1:numChans
%     h(ii) = figure;
% end
figure();
for ii=1:numChans
   xPos = Results.ScreenPos{ii,1};
   yPos = Results.ScreenPos{ii,2};
   bHorz = Results.b{ii,1};
   bVert = Results.b{ii,2};
   
   horzDesign = [ones(length(xPos),1),(xPos-xPos(1))./(xPos(end)-xPos(1)),((xPos-xPos(1))./(xPos(end)-xPos(1))).^2];
   vertDesign = [ones(length(yPos),1),(yPos-yPos(1))./(yPos(end)-yPos(1)),((yPos-yPos(1))./(yPos(end)-yPos(1))).^2];
   
   muHorz = exp(horzDesign*bHorz);
   muVert =exp(vertDesign*bVert);
   
   muHorz = repmat(muHorz',[dsStimLen(2),1]);
   muVert = repmat(muVert,[1,dsStimLen(1)]);
   
   x = linspace(xPos(1),xPos(end),dsStimLen(1));
   y = linspace(yPos(1),yPos(end),dsStimLen(2));
   centerX = (xPos(end)-xPos(1))/2; centerY = (yPos(end)-yPos(1))/2;
   
   xForDeg = (x-centerX)*mmPerPixel;yForDeg = (y-centerY)*mmPerPixel;
   
   xDeg = 2*atand(xForDeg/(2*DistToScreen));
   yDeg = 2*atand(yForDeg/(2*DistToScreen));
   
   finalIm = muHorz.*muVert;
   subplot(numChans,1,ii);axis([min(xDeg) max(xDeg) min(yDeg) max(yDeg)]);
   imagesc(xDeg,yDeg,finalIm);
   set(gca,'YDir','normal');colormap('jet');
   title(sprintf('LFP Retinotopy: Chan %d, Animal %d',ii,AnimalName));
   xlabel('Degrees Azimuth');
   ylabel('Degrees Elevation');
end


end