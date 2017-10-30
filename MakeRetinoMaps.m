function [] = MakeRetinoMaps(AnimalName,Date)
%MakeRetinoMaps.m
%  Plot the data from a retinotopic mapping experiment given the animal's
%   name and the date of the experiment
%
%Created: 2017/02/22, 10 Lawrence Street, Cambridge
%  Byron Price
%Updated: 2017/06/14
% By: Byron Price

fileName = sprintf('RetinoMapBayes%d_%d.mat',Date,AnimalName);
load(fileName)

DistToScreen = 250;

numChans = size(posteriorSample,1);
x = 1:1:w_pixels;y = 1:1:h_pixels;
[X,Y] = meshgrid(x,y);

centerX = w_pixels/2;centerY = h_pixels/4;
xForDeg = (x-centerX)*mmPerPixel;yForDeg = (y-centerY)*mmPerPixel;
xDeg = 2*atand(xForDeg/(2*DistToScreen));
yDeg = 2*atand(yForDeg/(2*DistToScreen));

% for ii=1:numChans
%     h(ii) = figure;
% end
numSamples = size(posteriorSample,3);
figure();
for ii=1:numChans
    subplot(numChans,1,ii);axis([min(xDeg) max(xDeg) min(yDeg) max(yDeg)]);
    title(sprintf('LFP Retinotopy: Chan %d, Animal %d',ii,AnimalName));
    xlabel('Degrees Azimuth');ylabel('Degrees Elevation');
    hold on;
    
    finalIm = zeros(length(y),length(x));
    samples = squeeze(posteriorSample(ii,:,:));
    N = 1000;
    for ll=1:N
        index = random('Discrete Uniform',numSamples);
        parameterVec = samples(:,index);
        b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(6)];
        distX = X-parameterVec(2);distY = Y-parameterVec(3);
        finalIm = finalIm+b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-...
                     (distY.^2)./(2*b(3)*b(3)))+b(4);
    end
    finalIm = finalIm./N;
    imagesc(xDeg,yDeg,finalIm);set(gca,'YDir','normal');w=colorbar;
    caxis([b(4) b(4)+150]);
%     caxis([310 470]);
    ylabel(w,'Mean VEP Magnitude (\muV)');colormap('jet');hold off;
    
end

end

