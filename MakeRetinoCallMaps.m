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
load(fileName)

numChans = size(Response,1);

for ii=1:numChans
    h(ii) = figure;
end

for ii=1:numChans
   xPos = Results.ScreenPos{ii,1};
   yPos = Results.ScreenPos{ii,2};
   bHorz = Results.b{ii,1};
   bVert = Results.b{ii,2};
   
   horzDesign = [ones(length(xPos),1),xPos./w_pixels,(xPos./w_pixels).^2];
   vertDesign = [ones(length(xPos),1),yPos./h_pixels,(yPos./h_pixels).^2];
   
   muHorz = exp(horzDesign*bHorz);
   muVert =exp(vertDesign*bVert);
   
   muHorz = repmat(muHorz',[dsStimLen(1),1]);
   muVert = repmat(muVert,[1,dsStimLen(2)]);
   
   finalIm = muHorz.*muVert;
   figure(h(ii));imagesc(linspace(0,w_pixels,dsStimLen(1)),linspace(0,h_pixels,dsStimLen(2)),finalIm);
   set(gca,'YDir','normal');colormap('jet');
   title(sprintf('LFP Retinotopy: Chan %d, Animal %d',ii,AnimalName));
   xlabel('Horizontal Screen Position (pixels)');
   ylabel('Vertical Screen Position (pixels)');
end


end