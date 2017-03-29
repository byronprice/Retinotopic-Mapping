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
    
    numDays = numFiles;
    
    if ConditionNumber == 1
        [numChans,numParameters] = size(dailyParameters{1});
        parameters = zeros(numDays,numChans,numParameters);
        confIntervals = zeros(numDays,numChans,numParameters);
        
        for jj=1:numDays
            parameters(jj,:,:) = dailyParameters{jj};
            confIntervals(jj,:,:) = parameterCI{jj};
            for kk=1:numChans
                tempNum = length(mapData{jj}{kk}(:,3));
                Y_retino = [Y_retino;squeeze(mapData{jj}{kk}(:,3))];
                tempMapXPos = squeeze(mapData{jj}{kk}(:,1));
                tempMapYPos = squeeze(mapData{jj}{kk}(:,2));
                tempMapXPos = tempMapXPos-squeeze(parameters(jj,kk,2));
                tempMapYPos = tempMapYPos-squeeze(parameters(jj,kk,3));
                Design_retino = [Design_retino;tempMapXPos,tempMapYPos,(Animals(ii)+kk)*ones(tempNum,1),jj*ones(tempNum,1)];
            end
        end
        
        for jj=1:numChans
            xdata = atand(squeeze(parameters(:,jj,2)).*pix_to_degree{1});
            ydata = atand(squeeze(parameters(:,jj,3)).*pix_to_degree{1});
        
            xErr = atand(squeeze(confIntervals(:,jj,2)).*pix_to_degree{1});
            yErr = atand(squeeze(confIntervals(:,jj,3)).*pix_to_degree{1});
            
            figure(h(1));
            xDist = xdata-xdata(1);yDist = ydata-ydata(1);
            for zz=1:length(xDist)
                allRelativeDistToCenterMass_1 = [allRelativeDistToCenterMass_1,sqrt(xDist(zz)^2+yDist(zz)^2)];
            end
            errorbar(xDist,yDist,yErr,yErr,xErr,xErr,linespecs{ConditionNumber},'LineWidth',2);hold on;
            title(sprintf('Center of Mass Relative to Day 1: %s',conditionNames{ConditionNumber}));
            xlabel('Horizontal dva');
            ylabel('Vertical dva');
            axis([-w_pixels/2 w_pixels/2 -h_pixels/2 h_pixels/2]);
        end
        hold off;
    elseif ConditionNumber == 2 
        [numChans,numParameters] = size(dailyParameters{1});
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
            for kk=1:numChans
                tempNum = length(mapData{jj}{kk}(:,3));
                Y_retino = [Y_retino;mapData{jj}{kk}(:,3)];
                tempMapXPos = squeeze(mapData{jj}{kk}(:,1));
                tempMapYPos = squeeze(mapData{jj}{kk}(:,2));
                tempMapXPos = tempMapXPos-squeeze(parameters(jj,kk,2));
                tempMapYPos = tempMapYPos-squeeze(parameters(jj,kk,3));
                Design_retino = [Design_retino;tempMapXpos,tempMapYPos,(Animals(ii)+kk)*ones(tempNum,1),jj*ones(tempNum,1)];
            end
        end
        
        for jj=1:numChans
            xdata = atand(squeeze(parameters(:,jj,2)).*pix_to_degree{1});
            ydata = atand(squeeze(parameters(:,jj,3)).*pix_to_degree{1});
        
            xErr = atand(squeeze(confIntervals(:,jj,2)).*pix_to_degree{1});
            yErr = atand(squeeze(confIntervals(:,jj,3)).*pix_to_degree{1});
            
            figure(h(2));
            xDist = xdata-xdata(1);yDist = ydata-ydata(1);
            for zz=1:length(xDist)
                allRelativeDistToCenterMass_2 = [allRelativeDistToCenterMass_2,sqrt(xDist(zz)^2+yDist(zz)^2)];
            end
            errorbar(xDist,yDist,yErr,yErr,xErr,xErr,linespecs{ConditionNumber},'LineWidth',2);hold on;
            title(sprintf('Center of Mass Relative to Day 1: %s',conditionNames{ConditionNumber}));
            xlabel('Horizontal dva');
            ylabel('Vertical dva');
            axis([-w_pixels/2 w_pixels/2 -h_pixels/2 h_pixels/2]);
        end
        hold off;
        
         for jj=1:numChans
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
day1 = Design_srp(:,1)==1;
day2 = Design_srp(:,1)==2;
day3 = Design_srp(:,1)==3;
day4 = Design_srp(:,1)==4;
day5 = Design_srp(:,1)==5;

cond2 = Design_srp(:,2)==1;
cond3 = Design_srp(:,2)==0;

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

figure(h(5));numDays = 5;
allData = zeros(numDays,3);
alpha = 0.05;
for ii=1:numDays
    dayData = Y_srp(Design_srp(:,1)==ii);
    allQuantiles = quantile(dayData,[alpha/2,0.5,1-alpha/2]);
    allQuantiles(2) = mean(dayData);
    allData(ii,:) = allQuantiles;
end
subplot(2,1,1);boundedline(1:numDays,allData(:,2),[abs(allData(:,1)-allData(:,2)),abs(allData(:,3)-allData(:,2))]);
title('Combined SRP Data Across Days');ylabel('Max-Min Statistic Magnitude (\muV)');
xlabel('Experimental Day');

cond2Data = zeros(numDays,3);
cond3Data = zeros(numDays,3);
cond2y = Y_srp(Design_srp(:,2)==1);
cond3y = Y_srp(Design_srp(:,2)==0);
cond2Design = Design_srp(Design_srp(:,2)==1,:);
cond3Design = Design_srp(Design_srp(:,2)==0,:);
for ii=1:numDays
    cond2temp = cond2y(cond2Design(:,1)==ii);
    cond3temp = cond3y(cond3Design(:,1)==ii);
    cond2quantiles = quantile(cond2temp,[alpha/2,0.5,1-alpha/2]);
    cond3quantiles = quantile(cond3temp,[alpha/2,0.5,1-alpha/2]);
    cond2quantiles(2) = mean(cond2temp);
    cond3quantiles(2) = mean(cond3temp);
    cond2Data(ii,:) = cond2quantiles;
    cond3Data(ii,:) = cond3quantiles;
end
subplot(2,1,2);boundedline(1:numDays,cond2Data(:,2),[abs(cond2Data(:,1)-cond2Data(:,2)),abs(cond2Data(:,3)-cond2Data(:,2))]);hold on;
boundedline(1:numDays,cond3Data(:,2),[abs(cond3Data(:,1)-cond3Data(:,2)),abs(cond3Data(:,3)-cond3Data(:,2))]);
title('Separated SRP Data');ylabel('Max-Min Statistic Magnitude (\muV)');
xlabel('Experimental Day');legend(conditionNames{2},conditionNames{3});

savefig(h,'MappingEffects_FinalResults.fig');

save('MappingEffects_FinalResults.mat','Y_srp','Design_srp','Y_retino','Design_retino');
end

