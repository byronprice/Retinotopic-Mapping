function [h] = MakeRetinoMaps(AnimalName,Date)
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

numChans = size(posteriorSample,1);
x = 1:10:w_pixels;y = 1:10:h_pixels;

for ii=1:numChans
    h(ii) = figure;
end
numSamples = size(posteriorSample,3);
for ii=1:numChans
    figure(h(ii));axis([0 w_pixels 0 h_pixels]);
    title(sprintf('LFP Retinotopy: Chan %d, Animal %d',ii,AnimalName));
    xlabel('Horizontal Screen Position (pixels)');ylabel('Vertical Screen Position (pixels)');
    hold on;
    
    finalIm = zeros(length(x),length(y));
    samples = squeeze(posteriorSample(ii,:,:));
    N = 1000;
    for ll=1:N
        index = random('Discrete Uniform',numSamples);
        parameterVec = samples(:,index);
        b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(6)];
        for jj=1:length(x)
            for kk=1:length(y)
                distX = x(jj)-parameterVec(2);
                distY = y(kk)-parameterVec(3);
                
                finalIm(jj,kk) = finalIm(jj,kk)+b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-...
                    (distY.^2)./(2*b(3)*b(3)))+b(4);
            end
        end
    end
    finalIm = finalIm./N;
    imagesc(x,y,finalIm');set(gca,'YDir','normal');w=colorbar;
    caxis([b(4) b(4)+150]);
    ylabel(w,'Log Mean VEP Magnitude (\muV)');colormap('jet');hold off;
    
end

end

