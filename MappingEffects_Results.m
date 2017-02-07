function [ ] = MappingEffects_Results(Animals)
%MappingEffects_Results.m
%
%  Use with MappingEffects.m stimuli and MappingEffects_Analysis.m
%   Once MappingEffects_Analysis.m has calculated receptive fields and
%    gotten SRP data, this code combines information across animals.
%
%INPUT: Animals - numeric array with animal IDs
%
%Created: 2017/02/07
%  Byron Price
%Updated: 2017/02/07
% By: Byron Price

numAnimals = length(Animals);
% filename = sprintf('MappingEffectsResults_%d.mat',Animals(1));
% 
% load(filename);
% [numDays,numChans,numParameters] = size(dailyParameters);
w_pixels = 2500; h_pixels = 1400;
% clear dailyParameters parameterCI srpVEP fisherInformation ConditionNumber;

modelSlope_retinoSRP = [];
modelSlope_greySRP = [];
Y = [];Design = [];

linespecs = {'or','+b'};
figure();hold on;
for ii=1:numAnimals
    filename = sprintf('MappingEffectsResults_%d.mat',Animals(ii));
    load(filename);
    
    [numDays,numChans,numParameters] = size(dailyParameters);
    if ConditionNumber == 1
        for jj=1:numChans
            xdata = squeeze(dailyParameters(:,jj,2));
            ydata = squeeze(dailyParameters(:,jj,3));
            magData = squeeze(dailyParameters(:,jj,1)+dailyParameters(:,jj,7));
        
            xErr = squeeze(parameterCI(:,jj,2))./1.96;
            yErr = squeeze(parameterCI(:,jj,3))./1.96;
            magErr = sqrt((squeeze(parameterCI(:,jj,1))./1.96).^2+(squeeze(parameterCI(:,jj,7))./1.96).^2);
            
            subplot(numAnimals,2,1+(ii-1)*2);errorbar(xdata,ydata,yErr,yErr,xErr,xErr,linespecs{ii});
            axis([0 w_pixels 0 h_pixels]);
            title('Center of Mass, Retinotopy-Grey');
            xlabel('Horizontal Screen Position');ylabel('Vertical Position');
            subplot(numAnimals,2,2+(ii-1)*2);errorbar(1:numDays,magData,magErr,linespecs{ii});
            axis([0 numDays+1 -500 -100]);
            title('Single-Trial VEP Magnitude');xlabel('Experimental Day');ylabel('Negativity (\muVolts');
        end
    elseif ConditionNumber == 2 
        for jj=1:numChans
            xdata = squeeze(dailyParameters(:,jj,2));
            ydata = squeeze(dailyParameters(:,jj,3));
            magData = squeeze(dailyParameters(:,jj,1)+dailyParameters(:,jj,7));
        
            xErr = squeeze(parameterCI(:,jj,2))./1.96;
            yErr = squeeze(parameterCI(:,jj,3))./1.96;
            magErr = sqrt((squeeze(parameterCI(:,jj,1))./1.96).^2+(squeeze(parameterCI(:,jj,7))./1.96).^2);
            
            subplot(numAnimals,2,1+(ii-1)*2);errorbar(xdata,ydata,yErr,yErr,xErr,xErr,linespecs{ii});
            axis([0 w_pixels 0 h_pixels]);
            title('Center of Mass, Retinotopy-SRP');
            xlabel('Horizontal Screen Position');ylabel('Vertical Position');
            subplot(numAnimals,2,2+(ii-1)*2);errorbar(1:numDays,magData,magErr,linespecs{ii});
            axis([0 numDays+1 -500 -100]);
            title('Single-Trial VEP Magnitude');xlabel('Experimental Day');ylabel('Negativity (\muVolts');
        end
        
        dayVec = [0:numDays-1]';
        for jj=1:numChans
            data = squeeze(srpVEP(:,jj,:));data = -min(data,[],2);
            [b,dev,stats] = glmfit(data,dayVec,'normal');
            modelSlope_retinoSRP = [modelSlope_retinoSRP,b(2)];
            
            Y = [Y;data];Design = [Design;dayVec,ones(numDays,1)];
        end
    elseif ConditionNumber == 3
        dayVec = [0:numDays-1]';
        for jj=1:numChans
            data = squeeze(srpVEP(:,jj,:));data = -min(data,[],2);
            [b,dev,stats] = glmfit(data,dayVec,'normal');
            modelSlope_greySRP = [modelSlope_greySRP,b(2)];
            
            Y = [Y;data];Design = [Design;dayVec,zeros(numDays,1)];
        end
    end
end
hold off;

nbins = 10;
figure();h1 = histogram(modelSlope_retinoSRP,nbins);hold on
h2 = histogram(modelSlope_greySRP,nbins);
h1.EdgeAlpha = 0.5;h2.EdgeAlpha = 0.5;
title('Histogram of Linear Regression Slope Coefficients for SRP');
xlabel('Slope Coefficient');ylabel('Count');legend('Retinotopy-SRP','Grey-SRP');

numModels = 4;
AIC = zeros(numModels,1);modelParams = 1:numModels;
modelFit = struct('b',cell(numModels,1),'se',cell(numModels,1));
[b,dev,stats] = glmfit(Y,ones(length(Y),1),'normal');
AIC(1) = dev+2*modelParams(1);
modelFit.b{1} = b;modelFit.se{1} = stats.se;

[b,dev,stats] = glmfit(Y,Design(:,1),'normal');
AIC(2) = dev+2*modelParams(2);
modelFit.b{2} = b;modelFit.se{2} = stats.se;

[b,dev,stats] = glmfit(Y,Design,'normal');
AIC(3) = dev+2*modelParams(3);
modelFit.b{3} = b;modelFit.se{3} = stats.se;

[b,dev,stats] = glmfit(Y,[Design,Design(:,1).*Design(:,2)],'normal');
AIC(4) = dev+2*modelParams(4);
modelFit.b{4} = b;modelFit.se{4} = stats.se;

figure();scatter(modelParams,AIC,'*r');title('Model Comparison with AIC');
xlabel('Number of Model Parameters');ylabel('AIC');
end

