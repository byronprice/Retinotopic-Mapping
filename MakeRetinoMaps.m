function [h] = MakeRetinoMaps(AnimalName,Date)
%MakeRetinoMaps.m
%  Plot the data from a retinotopic mapping experiment given the animal's
%   name and the date of the experiment
%
%Created: 2017/02/22, 10 Lawrence Street, Cambridge
%  Byron Price
%Updated: 2017/02/22
% By: Byron Price

fileName = sprintf('RetinoMap%d_%d.mat',Date,AnimalName);
load(fileName)

numChans = size(dimReduceData,1);
x = 1:w_pixels;y = 1:h_pixels;

for ii=1:numChans
    h(ii) = figure;
end

for ii=1:numChans
    figure(h(ii));axis([0 w_pixels 0 h_pixels]);
    title(sprintf('VEP Retinotopy: Chan %d, Animal %d',ii,AnimalName));
    xlabel('Horizontal Screen Position (pixels)');ylabel('Vertical Screen Position (pixels)');
    hold on;
    
    finalIm = zeros(length(x),length(y));
    parameterVec = finalParameters(ii,:);
    for jj=1:length(x)
        for kk=1:length(y)
            distX = x(jj)-parameterVec(2);
            distY = y(kk)-parameterVec(3);
            b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(7)];
            finalIm(jj,kk) = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4);
        end
    end
    imagesc(x,y,finalIm');set(gca,'YDir','normal');w=colorbar;
    caxis([parameterVec(7) parameterVec(7)+500]);
    ylabel(w,'Mean VEP Magnitude (\muV)');colormap('jet');hold off;
    
end

end

