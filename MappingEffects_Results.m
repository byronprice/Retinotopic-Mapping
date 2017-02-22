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
% clear dailyParameters parameterCI srpVEP fisherInformation ConditionNumber;

Y = [];Design = [];

conditionNames = {'Retino-Grey','Retino-SRP','Grey-SRP'};
linespecs = {'or','*b','xg'};
h(1) = figure();h(2) = figure();h(3) = figure();h(4) = figure();
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
            
            figure(h(3));
            errorbar(xdata-xdata(1),ydata-ydata(1),yErr,yErr,xErr,xErr,'o');hold on;
            title('Center of Mass Relative to Day One');xlabel('Horizontal Distance (pixels)');
            ylabel('Vertical Distance (pixels)');
            axis([-w_pixels/2 w_pixels/2 -h_pixels/2 h_pixels/2]);
            
            figure(h(1));
%             subplot(numAnimals,5,1+(ii-1)*5);hold on;
%             errorbar(xdata,ydata,yErr,yErr,xErr,xErr,linespecs{jj});
% %             scatter(xdata,ydata);
%             axis([0 w_pixels 0 h_pixels]);
%             title('CoM, Retinotopy-Grey');
%             xlabel('Horizontal Screen Position');ylabel('Vertical Position');
            
            subplot(2,2,1);hold on;
            errorbar(1:numDays,sigmaXdata,sigmaXerr,linespecs{ConditionNumber});
%             scatter(sigmaXdata,sigmaYdata);
            axis([0 numDays+1 0 800]);
            title('Region Size');
            xlabel('Experimental Day');ylabel('Parameter Sigma_x');
            
            subplot(2,2,2);hold on;
            errorbar(1:numDays,sigmaYdata,sigmaYerr,linespecs{ConditionNumber});
%             scatter(sigmaXdata,sigmaYdata);
            axis([0 numDays+1 0 800]);
            title('Region Size');
            xlabel('Experimental Day');ylabel('Parameter Sigma_y');
            
            subplot(2,2,3);hold on;
            errorbar(1:numDays,riseData,riseErr,linespecs{ConditionNumber});
%             scatter(1:numDays,riseData);
            axis([0 numDays+1 0 500]);
            title('VEP Rise at CoM');xlabel('Experimental Day');
            ylabel('Positivity-Negativity (\muV)');
            
            subplot(2,2,4);hold on;
            errorbar(1:numDays,baselineData,baselineErr,linespecs{ConditionNumber});
%             scatter(1:numDays,baselineData);
            axis([0 numDays+1 100 700]);
            title('Baseline');xlabel('Experimental Day');
            ylabel('Positivity-Negativity (\muV)');
        end
        hold off;
    elseif ConditionNumber == 2 
        [numChans,numParameters] = size(dailyParameters{1});
        srp_reps = size(srpSize{1},2);
        parameters = zeros(numDays,numChans,numParameters);
        confIntervals = zeros(numDays,numChans,numParameters);
        vepSize = zeros(numDays,numChans,srp_reps);
        for jj=1:numDays
            parameters(jj,:,:) = dailyParameters{jj};
            confIntervals(jj,:,:) = parameterCI{jj};
            vepSize(jj,:,:) = srpSize{jj};
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
            
            figure(h(3));
            errorbar(xdata-xdata(1),ydata-ydata(1),yErr,yErr,xErr,xErr,'o');hold on;
            title('Center of Mass Relative to Day One');xlabel('Horizontal Distance (pixels)');
            ylabel('Vertical Distance (pixels)');
            axis([-w_pixels/2 w_pixels/2 -h_pixels/2 h_pixels/2]);
           
            figure(h(1));
%             subplot(numAnimals,5,1+(ii-1)*5);hold on;
%             errorbar(xdata,ydata,yErr,yErr,xErr,xErr,linespecs{jj});
% %             scatter(xdata,ydata);
%             axis([0 w_pixels 0 h_pixels]);
%             title('CoM, Retinotopy-SRP');
%             xlabel('Horizontal Screen Position');ylabel('Vertical Position');
            
            subplot(2,2,1);hold on;
            errorbar(1:numDays,sigmaXdata,sigmaXerr,linespecs{ConditionNumber});
%             scatter(sigmaXdata,sigmaYdata);
            axis([0 numDays+1 0 800]);
            title('Region Size');
            xlabel('Experimental Day');ylabel('Parameter Sigma_x');
            
            subplot(2,2,2);hold on;
            errorbar(1:numDays,sigmaYdata,sigmaYerr,linespecs{ConditionNumber});
%             scatter(sigmaXdata,sigmaYdata);
            axis([0 numDays+1 0 800]);
            title('Region Size');
            xlabel('Experimental Day');ylabel('Parameter Sigma_y');
            
            subplot(2,2,3);hold on;
            errorbar(1:numDays,riseData,riseErr,linespecs{ConditionNumber});
%             scatter(1:numDays,riseData);
            axis([0 numDays+1 0 500]);
            title('VEP Rise at CoM');xlabel('Experimental Day');
            ylabel('Positivity-Negativity (\muV)');
            
            subplot(2,2,4);hold on;
            errorbar(1:numDays,baselineData,baselineErr,linespecs{ConditionNumber});
%             scatter(1:numDays,baselineData);
            axis([0 numDays+1 100 700]);
            title('Baseline');xlabel('Experimental Day');
            ylabel('Positivity-Negativity (\muV)');
        end
        hold off;
        
        figure(h(2));
         for jj=1:numChans
             data = zeros(numDays,srp_reps);dayVec = zeros(numDays,srp_reps);
             meanVals = zeros(numDays,1);
             stdVals = zeros(numDays,1);
            for kk=1:numDays
                for ll=1:srp_reps
                    data(kk,ll) = vepSize(kk,jj,ll);
                    dayVec(kk,ll) = kk;
                end
                meanVals(kk) = mean(data(kk,:));
                stdVals(kk) = std(data(kk,:));
            end
                randJitter = randn([numDays*srp_reps,1]).*0.05;
                subplot(1,2,1);hold on;
%                 scatter(dayVec(:)+randJitter,data(:),linespecs{jj});hold on;
%                 scatter(1:numDays,meanVals,50,'^k','filled');
                scatter(1:numDays,meanVals,linespecs{ConditionNumber});
                axis([0 numDays+1 400 900]);
                title('SRP');xlabel('Experimental Day');ylabel('Log(Pos-Neg) (\muV)');
                Y = [Y;data(:)];Design = [Design;dayVec(:),ones(srp_reps*numDays,1)];
                
                subplot(1,2,2);hold on;
                scatter(meanVals,stdVals,linespecs{ConditionNumber});axis([400 900 0 500]);
                title('Mean vs STD VEP Size');
                xlabel('Mean Log(Pos-Neg)');ylabel('STD Log(Pos-Neg)');
         end
        hold off;
    elseif ConditionNumber == 3
        [numChans,srp_reps] = size(srpSize{1});
        vepSize = zeros(numDays,numChans,srp_reps);
        for jj=1:numDays
            vepSize(jj,:,:) = srpSize{jj};
        end
        
        figure(h(2));
         for jj=1:numChans
             data = zeros(numDays,srp_reps);dayVec = zeros(numDays,srp_reps);
             meanVals = zeros(numDays,1);stdVals = zeros(numDays,1);
            for kk=1:numDays
                for ll=1:srp_reps
                    data(kk,ll) = vepSize(kk,jj,ll);
                    dayVec(kk,ll) = kk;
                end
                meanVals(kk) = mean(data(kk,:));
                stdVals(kk) = std(data(kk,:));
            end
                randJitter = randn([numDays*srp_reps,1]).*0.05;
                subplot(1,2,1);hold on;
%                 scatter(dayVec(:)+randJitter,data(:),linespecs{jj});hold on;
%                 scatter(1:numDays,meanVals,50,'^k','filled');
                scatter(1:numDays,meanVals,linespecs{ConditionNumber});
                axis([0 numDays+1 400 900]);
                title('SRP');xlabel('Experimental Day');ylabel('Log(Pos-Neg) (\muV)');
                Y = [Y;data(:)];Design = [Design;dayVec(:),zeros(srp_reps*numDays,1)];
                
                subplot(1,2,2);hold on;
                scatter(meanVals,stdVals,linespecs{ConditionNumber});axis([400 900 0 500]);
                title('Mean vs STD VEP Size');
                xlabel('Mean Log(Pos-Neg)');ylabel('STD Log(Pos-Neg)');
         end
        hold off;
    end
end

numModels = 5;
AIC = zeros(numModels,1);modelParams = [1,2,3,4,3];
modelFit = struct('b',cell(numModels,1),'se',cell(numModels,1),...
    'p',cell(numModels,1),'dev',zeros(numModels,1),'resid',cell(numModels,1),...
    'dfe',zeros(numModels,1));

% constant model, VEP ~ b0
[b,dev,stats] = glmfit(ones(length(Y),1),Y,'normal','link','identity','constant','off');
AIC(1) = dev+2*modelParams(1);
modelFit(1).b = b;modelFit(1).se = stats.se;
modelFit(1).p = stats.p;
modelFit(1).dev = dev;
modelFit(1).resid = stats.resid;
modelFit(1).dfe = stats.dfe;
% figure();scatter(Design(:,1),modelFit(1).resid);
clear b dev stats;

% linear model, VEP ~ b0+b1*day
[b,dev,stats] = glmfit(Design(:,1),Y,'normal','link','identity');
AIC(2) = dev+2*modelParams(2);
modelFit(2).b = b;modelFit(2).se = stats.se;
modelFit(2).p = stats.p;
modelFit(2).dev = dev;
modelFit(2).resid = stats.resid;
modelFit(2).dfe = stats.dfe;
% figure();scatter(Design(:,1),modelFit(2).resid);
clear b dev stats;

% linear model with two intercepts, VEP ~ b0 + b1*day + b2*I(retino-srp)
[b,dev,stats] = glmfit(Design,Y,'normal','link','identity');
AIC(3) = dev+2*modelParams(3);
modelFit(3).b = b;modelFit(3).se = stats.se;
modelFit(3).p = stats.p;
modelFit(3).dev = dev;
modelFit(3).resid = stats.resid;
modelFit(3).dfe = stats.dfe;
% figure();scatter(Design(:,1),modelFit(3).resid);
clear b dev stats;

% linear model with two intercepts, two slopes, VEP ~ b0 + b1*day +
% b2*I(retino-srp) + b3*I(retino-srp)*day
[b,dev,stats] = glmfit([Design,Design(:,1).*Design(:,2)],Y,'normal','link','identity');
AIC(4) = dev+2*modelParams(4);
modelFit(4).b = b;modelFit(4).se = stats.se;
modelFit(4).p = stats.p;
modelFit(4).dev = dev;
modelFit(4).resid = stats.resid;
modelFit(4).dfe = stats.dfe;
% figure();scatter(Design(:,1),modelFit(4).resid);
clear b dev stats;

% linear model with one intercept, two slopes, VEP ~ b0 + b1*day +
% b2*I(retino-srp)*day
[b,dev,stats] = glmfit([Design(:,1),Design(:,1).*Design(:,2)],Y,'normal','link','identity');
AIC(5) = dev+2*modelParams(5);
modelFit(5).b = b;modelFit(5).se = stats.se;
modelFit(5).p = stats.p;
modelFit(5).dev = dev;
modelFit(5).resid = stats.resid;
modelFit(5).dfe = stats.dfe;
% figure();scatter(Design(:,1),modelFit(5).resid);
clear b dev stats;

figure(h(4));plot(modelParams,AIC,'b','LineWidth',2);title('Model Comparison with AIC');
xlabel('Number of Model Parameters');ylabel('AIC');

% chiSquareTest = zeros(numModels-1,1);
% for ii=2:numModels
%     statistic = (modelFit(ii-1).dev-modelFit(ii).dev);
%     k = modelParams(ii)-modelParams(ii-1);
%     chiSquareTest(ii-1) = 1-chi2cdf(statistic,k);
% end
% chiSquareTest

display('F tests: ');

display('Significant increase across days: ');
SS1 = sum((modelFit(1).resid).^2);
SS2 = sum((modelFit(2).resid).^2);
df1 = modelFit(1).dfe;df2 = modelFit(2).dfe;

F = ((SS1-SS2)/(df1-df2))/(SS2/df2);
p = 1-fcdf(F,df1-df2,df2)

display('Retino-SRP and Grey-SRP have different day one and same slope: ');
SS1 = sum((modelFit(2).resid).^2);
SS2 = sum((modelFit(3).resid).^2);
df1 = modelFit(2).dfe;df2 = modelFit(3).dfe;

F = ((SS1-SS2)/(df1-df2))/(SS2/df2);
p = 1-fcdf(F,df1-df2,df2)

display('Retino-SRP and Grey-SRP have same day one and different slope: ');
SS1 = sum((modelFit(2).resid).^2);
SS2 = sum((modelFit(5).resid).^2);
df1 = modelFit(2).dfe;df2 = modelFit(5).dfe;

F = ((SS1-SS2)/(df1-df2))/(SS2/df2);
p = 1-fcdf(F,df1-df2,df2)

display('Retino-SRP and Grey-SRP have different day one and different slope: ');
SS1 = sum((modelFit(2).resid).^2);
SS2 = sum((modelFit(4).resid).^2);
df1 = modelFit(2).dfe;df2 = modelFit(4).dfe;

F = ((SS1-SS2)/(df1-df2))/(SS2/df2);
p = 1-fcdf(F,df1-df2,df2)

% display('Chi-square tests: ');
% 
% display('Significant increase across days: ');
% statistic = (modelFit(1).dev-modelFit(2).dev);
% k = modelParams(2)-modelParams(1);
% p = 1-chi2cdf(statistic,k)
% 
% display('Retino-SRP and Grey-SRP have different day one VEP magnitude: ');
% statistic = (modelFit(2).dev-modelFit(3).dev);
% k = modelParams(3)-modelParams(2);
% p = 1-chi2cdf(statistic,k)
% 
% display('Retino-SRP and Grey-SRP have same day one and different slope: ');
% statistic = (modelFit(2).dev-modelFit(5).dev);
% k = modelParams(5)-modelParams(2);
% p = 1-chi2cdf(statistic,k)
% 
% display('Retino-SRP and Grey-SRP have different day one and different slope: ');
% statistic = (modelFit(2).dev-modelFit(4).dev);
% k = modelParams(4)-modelParams(2);
% p = 1-chi2cdf(statistic,k)
end

