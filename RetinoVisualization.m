function [ ] = RetinoVisualization(AnimalName,Date)
%RetinoVisualization.m
%  Look at the data from a retinotopic mapping experiment in two 
%   ways.
%INPUT: AnimalName - animal's unique identifier as a number
%       Date - date of the experiment as a number

%OUTPUT: Plots
%
% Created: 2017/03/10, 24 Cummington Mall
%  Byron Price
% Updated: 2017/03/10
% By: Byron Price

fileName = sprintf('RetinoMap%d_%d.mat',Date,AnimalName);
load(fileName);

numChans = size(dimReduceData,1);
stimLen = size(vepResponse,3);

xAxis = linspace(0,w_pixels,12);
yAxis = linspace(0,h_pixels,8);
xDiff = mean(diff(xAxis));
yDiff = mean(diff(yAxis));

xconv = stimLen/xDiff;
yconv = 800/yDiff; % for height of the stimulus

positivity = 120:300;
negativity = 50:120;
for ii=1:numChans
%     data = max(squeeze(vepResponse(ii,:,positivity)),[],2)-min(squeeze(vepResponse(ii,:,negativity)),[],2);
    data = dimReduceData{ii}(:,3);
    xPos = dimReduceData{ii}(:,1);
    yPos = dimReduceData{ii}(:,2);
    
    figure();
    scatter(xPos,yPos,[],data);colormap('jet');
    title(sprintf('%d %d',Date,AnimalName));
    
%     figure();
%     for jj=1:length(xAxis)-1
%         for kk=1:length(yAxis)-1
%             xRange = [xAxis(jj),xAxis(jj+1)];
%             yRange = [yAxis(kk),yAxis(kk+1)];
%             
%             xOne = find(xPos > xRange(1));
%             xTwo = find(xPos <= xRange(2));
%             xInds = intersect(xOne,xTwo);
%             yOne = find(yPos > yRange(1));
%             yTwo = find(yPos <=yRange(2));
%             yInds = intersect(yOne,yTwo);
%             
%             commonInds = intersect(xInds,yInds);
%             
%             if isempty(commonInds) ~= 1
%                 meanVep = mean(squeeze(vepResponse(ii,commonInds,:)),1);
%                 
%                 plot((1:1:stimLen)./xconv+mean(xRange)-0.5*xDiff,...
%                     (meanVep)./yconv+mean(yRange),'k','LineWidth',2);hold on;
%             end
%         end
%     end
    
end

end
