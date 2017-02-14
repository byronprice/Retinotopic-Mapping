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
w_pixels = 2560; h_pixels = 1440;
% clear dailyParameters parameterCI srpVEP fisherInformation ConditionNumber;

Y = [];Design = [];

linespecs = {'or','*b'};
h(1) = figure();h(2) = figure();h(3) = figure();
for ii=1:numAnimals
    filename = sprintf('MappingEffectsResults_%d.mat',Animals(ii));
    load(filename);
    
    numDays = numFiles;
    
    if ConditionNumber == 1
        [numChans,numParameters] = size(dailyParameters{1});
        parameters = zeros(numDays,numChans,numParameters);
        confIntervals = zeros(numDays,numChans,numParameters);
        for jj=1:numDays
            parameters(jj,:,:) = dailyParameters{jj};
            confIntervals(jj,:,:) = parameterCI{jj};
        end
        
        for jj=1:numChans
            xdata = squeeze(parameters(:,jj,2));
            ydata = squeeze(parameters(:,jj,3));
            sigmaXdata = squeeze(parameters(:,jj,4));
            sigmaYdata = squeeze(parameters(:,jj,5));
            riseData = squeeze(parameters(:,jj,1));
            baselineData = squeeze(parameters(:,jj,7));
        
            xErr = squeeze(confIntervals(:,jj,2));
            yErr = squeeze(confIntervals(:,jj,3));
            sigmaXerr = squeeze(confIntervals(:,jj,4));
            sigmaYerr = squeeze(confIntervals(:,jj,5));
            riseErr = squeeze(confIntervals(:,jj,1));
            baselineErr = squeeze(confIntervals(:,jj,7));
            
            figure(h(1));
            subplot(numAnimals,4,1+(ii-1)*4);hold on;
%             errorbar(xdata,ydata,yErr,yErr,xErr,xErr,linespecs{jj});
            scatter(xdata,ydata);
            axis([0 w_pixels 0 h_pixels]);
            title('CoM, Retinotopy-Grey');
            xlabel('Horizontal Screen Position');ylabel('Vertical Position');
            
            subplot(numAnimals,4,2+(ii-1)*4);hold on;
%             errorbar(sigmaXdata,sigmaYdata,...
%                 sigmaYerr,sigmaYerr,sigmaXerr,sigmaXerr,linespecs{jj});
            scatter(sigmaXdata,sigmaYdata);
            axis([100 600 100 600]);
            title('Region Size');
            xlabel('Parameter Sigma_x');ylabel('Parameter Sigma_y');
            
            subplot(numAnimals,4,3+(ii-1)*4);hold on;
%             errorbar(1:numDays,riseData,riseErr,linespecs{jj});
            scatter(1:numDays,riseData);
            axis([0 numDays+1 0 500]);
            title('Negativity Rise at CoM');xlabel('Experimental Day');
            ylabel('VEP Negativity (\muV)');
            
            subplot(numAnimals,4,4+(ii-1)*4);hold on;
%             errorbar(1:numDays,baselineData,baselineErr,linespecs{jj});
            scatter(1:numDays,baselineData);
            axis([0 numDays+1 0 500]);
            title('Baseline Negativity');xlabel('Experimental Day');
            ylabel('VEP Negativity (\muV)');
        end
        hold off;
    elseif ConditionNumber == 2 
        [numChans,numParameters] = size(dailyParameters{1});
        srp_reps = size(srpNegativity{1},2);
        parameters = zeros(numDays,numChans,numParameters);
        confIntervals = zeros(numDays,numChans,numParameters);
        negativity = zeros(numDays,numChans,srp_reps);
        for jj=1:numDays
            parameters(jj,:,:) = dailyParameters{jj};
            confIntervals(jj,:,:) = parameterCI{jj};
            negativity(jj,:,:) = srpNegativity{jj};
        end
        
        for jj=1:numChans
             xdata = squeeze(parameters(:,jj,2));
            ydata = squeeze(parameters(:,jj,3));
            sigmaXdata = squeeze(parameters(:,jj,4));
            sigmaYdata = squeeze(parameters(:,jj,5));
            riseData = squeeze(parameters(:,jj,1));
            baselineData = squeeze(parameters(:,jj,7));
        
            xErr = squeeze(confIntervals(:,jj,2));
            yErr = squeeze(confIntervals(:,jj,3));
            sigmaXerr = squeeze(confIntervals(:,jj,4));
            sigmaYerr = squeeze(confIntervals(:,jj,5));
            riseErr = squeeze(confIntervals(:,jj,1));
            baselineErr = squeeze(confIntervals(:,jj,7));
            
            figure(h(1));
            subplot(numAnimals,4,1+(ii-1)*4);hold on;
%             errorbar(xdata,ydata,yErr,yErr,xErr,xErr,linespecs{jj});
            scatter(xdata,ydata);
            axis([0 w_pixels 0 h_pixels]);
            title('CoM, Retinotopy-SRP');
            xlabel('Horizontal Screen Position');ylabel('Vertical Position');
            
            subplot(numAnimals,4,2+(ii-1)*4);hold on;
%             errorbar(sigmaXdata,sigmaYdata,...
%                 sigmaYerr,sigmaYerr,sigmaXerr,sigmaXerr,linespecs{jj});
            scatter(sigmaXdata,sigmaYdata);
            axis([100 600 100 600]);
            title('Region Size');
            xlabel('Parameter Sigma_x');ylabel('Parameter Sigma_y');
            
            subplot(numAnimals,4,3+(ii-1)*4);hold on;
%             errorbar(1:numDays,riseData,riseErr,linespecs{jj});
            scatter(1:numDays,riseData);
            axis([0 numDays+1 0 500]);
            title('Negativity Rise at CoM');xlabel('Experimental Day');
            ylabel('VEP Negativity (\muV)');
            
            subplot(numAnimals,4,4+(ii-1)*4);hold on;
%             errorbar(1:numDays,baselineData,baselineErr,linespecs{jj});
            scatter(1:numDays,baselineData);
            axis([0 numDays+1 0 500]);
            title('Baseline Negativity');xlabel('Experimental Day');
            ylabel('VEP Negativity (\muV)');
        end
        hold off;
        
        figure(h(2));
         for jj=1:numChans
             data = zeros(numDays,srp_reps);dayVec = zeros(numDays,srp_reps);
             meanVals = zeros(numDays,1);
            for kk=1:numDays
                for ll=1:srp_reps
                    data(kk,ll) = abs(negativity(kk,jj,ll));
                    dayVec(kk,ll) = kk;
                end
                meanVals(kk) = mean(data(kk,:));
            end
                randJitter = randn([numDays*srp_reps,1]).*0.05;
                subplot(numAnimals,2,jj+(ii-1)*2);hold on;
                scatter(dayVec(:)+randJitter,data(:),linespecs{jj});hold on;
                scatter(1:numDays,meanVals,50,'^k','filled');
                axis([0 numDays+1 0 800]);
                title('SRP Negativity, Retinotopy-SRP');xlabel('Experimental Day');ylabel('Negativity (\muV)');
                Y = [Y;data(:)];Design = [Design;dayVec(:),ones(srp_reps*numDays,1)];
         end
        hold off;
    elseif ConditionNumber == 3
        [numChans,srp_reps] = size(srpNegativity{1});
        negativity = zeros(numDays,numChans,srp_reps);
        for jj=1:numDays
            negativity(jj,:,:) = srpNegativity{jj};
        end
        
        figure(h(2));
         for jj=1:numChans
             data = zeros(numDays,srp_reps);dayVec = zeros(numDays,srp_reps);
             meanVals = zeros(numDays,1);
            for kk=1:numDays
                for ll=1:srp_reps
                    data(kk,ll) = abs(negativity(kk,jj,ll));
                    dayVec(kk,ll) = kk;
                end
                meanVals(kk) = mean(data(kk,:));
            end
                randJitter = randn([numDays*srp_reps,1]).*0.05;
                subplot(numAnimals,2,jj+(ii-1)*2);hold on;
                scatter(dayVec(:)+randJitter,data(:),linespecs{jj});hold on;
                scatter(1:numDays,meanVals,50,'^k','filled');
                axis([0 numDays+1 0 800]);
                title('SRP Negativity, Grey-SRP');xlabel('Experimental Day');ylabel('Negativity (\muV)');
                Y = [Y;data(:)];Design = [Design;dayVec(:),zeros(srp_reps*numDays,1)];
         end
        hold off;
    end
end

numModels = 4;
AIC = zeros(numModels,1);modelParams = 1:numModels;
modelFit = struct('b',cell(numModels,1),'se',cell(numModels,1),...
    'p',cell(numModels,1),'dev',zeros(numModels,1));

[b,dev,stats] = glmfit(ones(length(Y),1),Y,'normal','link','log','constant','off');
AIC(1) = dev+2*modelParams(1);
modelFit(1).b = b;modelFit(1).se = stats.se;
modelFit(1).p = stats.p;
modelFit(1).dev = dev;
clear b dev stats;

[b,dev,stats] = glmfit(Design(:,1),Y,'normal','link','log');
AIC(2) = dev+2*modelParams(2);
modelFit(2).b = b;modelFit(2).se = stats.se;
modelFit(2).p = stats.p;
modelFit(2).dev = dev;
clear b dev stats;

[b,dev,stats] = glmfit(Design,Y,'normal','link','log');
AIC(3) = dev+2*modelParams(3);
modelFit(3).b = b;modelFit(3).se = stats.se;
modelFit(3).p = stats.p;
modelFit(3).dev = dev;
clear b dev stats;

[b,dev,stats] = glmfit([Design,Design(:,1).*Design(:,2)],Y,'normal','link','log');
AIC(4) = dev+2*modelParams(4);
modelFit(4).b = b;modelFit(4).se = stats.se;
modelFit(4).p = stats.p;
modelFit(4).dev = dev;
clear b dev stats;

figure(h(3));plot(modelParams,AIC,'b','LineWidth',2);title('Model Comparison with AIC');
xlabel('Number of Model Parameters');ylabel('AIC');

chiSquareTest = zeros(numModels-1,1);
for ii=2:numModels
    statistic = (modelFit(ii-1).dev-modelFit(ii).dev);
    k = modelParams(ii)-modelParams(ii-1);
    chiSquareTest(ii-1) = 1-chi2cdf(statistic,k);
end
chiSquareTest

for ii=1:numModels
    display(modelFit(ii).b);
    display(modelFit(ii).se);
end
end

