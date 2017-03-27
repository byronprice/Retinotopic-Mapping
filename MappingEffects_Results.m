function [ ] = MappingEffects_Results(Animals,Channels)
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
%Updated: 2017/03/27
% By: Byron Price

numAnimals = length(Animals);
% filename = sprintf('MappingEffectsResults_%d.mat',Animals(1));
% 
% load(filename);
% [numDays,numChans,numParameters] = size(dailyParameters);
% clear dailyParameters parameterCI srpVEP fisherInformation ConditionNumber;

Y_srp = [];Design_srp = [];
Y_retino = [];Design_retino = [];

conditionNames = {'Retino-Grey','Retino-SRP','Grey-SRP'};
linespecs = {'or','*b','xg'};
h(1) = figure();h(2) = figure();h(3) = figure();h(4) = figure();
h(5) = figure();

allRelativeDistToCenterMass_1 = [];
allRelativeDistToCenterMass_2 = [];
for ii=1:numAnimals
    filename = sprintf('MappingEffectsResults_%d.mat',Animals(ii));
    load(filename);
    
    currentChannels = Channels{ii};
    numDays = numFiles;
    
    if ConditionNumber == 1
        numChans = length(currentChannels);
        [~,numParameters] = size(dailyParameters{1});
        parameters = zeros(numDays,numChans,numParameters);
        confIntervals = zeros(numDays,numChans,numParameters);
        
        for jj=1:numDays
            parameters(jj,:,:) = dailyParameters{jj};
            confIntervals(jj,:,:) = parameterCI{jj};
            for kk=currentChannels
                tempNum = length(mapData{jj}{kk}(:,3));
                Y_retino = [Y_retino;mapData{jj}{kk}(:,3)];
                Design_retino = [Design_retino;mapData{jj}{kk}(:,1:2),(Animals(ii)+kk)*ones(tempNum,1)];
            end
        end
        
        for jj=currentChannels
            xdata = atand(squeeze(parameters(:,jj,2)).*pix_to_degree{1});
            ydata = atand(squeeze(parameters(:,jj,3)).*pix_to_degree{1});
        
            xErr = atand(squeeze(confIntervals(:,jj,2)).*pix_to_degree{1});
            yErr = atand(squeeze(confIntervals(:,jj,3)).*pix_to_degree{1});
            
            figure(h(1));
            xDist = xdata-xdata(1);yDist = ydata-ydata(1);
            for zz=1:length(xDist)
                allRelativeDistToCenterMass_1 = [allRelativeDistToCenterMass_1,sqrt(xDist(zz)^2+yDist(zz)^2)];
            end
            errorbar(xDist,yDist,yErr,yErr,xErr,xErr,'o','LineWidth',2);hold on;
            title(sprintf('Center of Mass Relative to Day 1: %s',conditionNames{ConditionNumber}));
            xlabel('Horizontal dva');
            ylabel('Vertical dva');
            axis([-w_pixels/2 w_pixels/2 -h_pixels/2 h_pixels/2]);
        end
        hold off;
    elseif ConditionNumber == 2 
        numChans = length(currentChannels);
        [~,numParameters] = size(dailyParameters{1});
        srp_reps = size(srpSize{1},2);
        stimLen = size(srpVEP{1},2);
        parameters = zeros(numDays,numChans,numParameters);
        confIntervals = zeros(numDays,numChans,numParameters);
        vepSize = zeros(numDays,numChans,srp_reps);
        meanVEP = zeros(numDays,numChans,stimLen);
        for jj=1:numDays
            parameters(jj,:,:) = dailyParameters{jj};
            confIntervals(jj,:,:) = parameterCI{jj};
            vepSize(jj,:,:) = srpSize{jj};
            meanVEP(jj,:,:) = srpVEP{jj};
            for kk=currentChannels
                tempNum = length(mapData{jj}{kk}(:,3));
                Y_retino = [Y_retino;mapData{jj}{kk}(:,3)];
                Design_retino = [Design_retino;mapData{jj}{kk}(:,1:2),(Animals(ii)+kk)*ones(tempNum,1)];
            end
        end
        
        for jj=1:currentChannels
            xdata = atand(squeeze(parameters(:,jj,2)).*pix_to_degree{1});
            ydata = atand(squeeze(parameters(:,jj,3)).*pix_to_degree{1});
        
            xErr = atand(squeeze(confIntervals(:,jj,2)).*pix_to_degree{1});
            yErr = atand(squeeze(confIntervals(:,jj,3)).*pix_to_degree{1});
            
            figure(h(2));
            xDist = xdata-xdata(1);yDist = ydata-ydata(1);
            for zz=1:length(xDist)
                allRelativeDistToCenterMass_2 = [allRelativeDistToCenterMass_2,sqrt(xDist(zz)^2+yDist(zz)^2)];
            end
            errorbar(xDist,yDist,yErr,yErr,xErr,xErr,'o','LineWidth',2);hold on;
            title(sprintf('Center of Mass Relative to Day 1: %s',conditionNames{ConditionNumber}));
            xlabel('Horizontal dva');
            ylabel('Vertical dva');
            axis([-w_pixels/2 w_pixels/2 -h_pixels/2 h_pixels/2]);
        end
        hold off;
        
         for jj=currentChannels
             data = zeros(numDays,srp_reps);dayVec = zeros(numDays,srp_reps);
            for kk=1:numDays
                for ll=1:srp_reps
                    data(kk,ll) = vepSize(kk,jj,ll);
                    dayVec(kk,ll) = kk;
                end
            end
            Y_srp = [Y_srp;data(:)];Design_srp = [Design_srp;dayVec(:),ones(srp_reps*numDays,1)];
         end

    elseif ConditionNumber == 3
        numChans = length(currentChannels);
        [~,srp_reps] = size(srpSize{1});
        stimLen = size(srpVEP{1},2);
        vepSize = zeros(numDays,numChans,srp_reps);
        meanVEP = zeros(numDays,numChans,stimLen);
        for jj=1:numDays
            vepSize(jj,:,:) = srpSize{jj};
            meanVEP(jj,:,:) = srpVEP{jj};
        end
        
         for jj=currentChannels
            data = zeros(numDays,srp_reps);dayVec = zeros(numDays,srp_reps);
            for kk=1:numDays
                for ll=1:srp_reps
                    data(kk,ll) = vepSize(kk,jj,ll);
                    dayVec(kk,ll) = kk;
                end
            end
            
            Y_srp = [Y_srp;data(:)];Design_srp = [Design_srp;dayVec(:),zeros(srp_reps*numDays,1)];
                
         end
    end
end

figure(h(3));
day1 = find(Design_srp(:,1)==1);
day2 = find(Design_srp(:,1)==2);
day3 = find(Design_srp(:,1)==3);
day4 = find(Design_srp(:,1)==4);
day5 = find(Design_srp(:,1)==5);

cond2 = find(Design_srp(:,2)==1);
cond3 = find(Design_srp(:,2)==0);

barMean = [mean(Y_srp(day1 & cond2)),mean(Y_srp(day1 & cond3));...
    mean(Y_srp(day2 & cond2)),mean(Y_srp(day2 & cond3));...
    mean(Y_srp(day3 & cond2)),mean(Y_srp(day3 & cond3));...
    mean(Y_srp(day4 & cond2)),mean(Y_srp(day4 & cond3));...
    mean(Y_srp(day5 & cond2)),mean(Y_srp(day5 & cond3))];

barError = [std(Y_srp(day1 & cond2))/sqrt(length(Y_srp(day1 & cond2))),...
    std(Y_srp(day1 & cond3))/sqrt(length(Y_srp(day1 & cond3)));...
    std(Y_srp(day2 & cond2))/sqrt(length(Y_srp(day2 & cond2))),...
    std(Y_srp(day2 & cond3))/sqrt(length(Y_srp(day2 & cond3)));...
    std(Y_srp(day3 & cond2))/sqrt(length(Y_srp(day3 & cond2))),...
    std(Y_srp(day3 & cond3))/sqrt(length(Y_srp(day3 & cond3)));...
    std(Y_srp(day4 & cond2))/sqrt(length(Y_srp(day4 & cond2))),...
    std(Y_srp(day4 & cond3))/sqrt(length(Y_srp(day4 & cond3)));...
    std(Y_srp(day5 & cond2))/sqrt(length(Y_srp(day5 & cond2))),...
    std(Y_srp(day5 & cond3))/sqrt(length(Y_srp(day5 & cond3)))];

superbar(1:5,barMean','E',barError');title('SRP Protocol Mean VEP Magnitude');
xlabel('Experimental Day');ylabel('Magnitude (\muV)');
legend(conditionNames{2},conditionNames{3});

figure(h(4));
histogram(allRelativeDistToCenterMass_1);hold on;
histogram(allRelativeDistToCenterMass_2);
legend(conditionNames{1},conditionNames{2});
title('Center of Mass (CoM) Relative To Day 1');
xlabel('Distance to Day 1 CoM (dva)');
ylabel('Count');

savefig(h,'MappingEffects_FinalResults.fig');

save('MappingEffects_FinalResults.mat','Y_srp','Design_srp','Y_retino','Design_retino');
end

