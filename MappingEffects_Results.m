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
h(1) = figure();h(3) = figure();h(4) = figure();
h(5) = figure();

allRelativeDistToCenterMass_1 = [];
allRelativeDistToCenterMass_2 = [];
allStd_1 = [];
allStd_2 = [];
allHeight_1 = [];
allHeight_2 = [];
allParams = [];
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
                if (squeeze(parameters(jj,kk,1))-squeeze(confIntervals(jj,kk,1))) > 0.05
                    tempNum = length(mapData{jj}{kk}(:,3));
                    
                    tempMapXPos = squeeze(mapData{jj}{kk}(:,1));
                    tempMapYPos = squeeze(mapData{jj}{kk}(:,2));
                    tempMapXPos = tempMapXPos-squeeze(parameters(jj,kk,2));
                    tempMapYPos = tempMapYPos-squeeze(parameters(jj,kk,3));
                    
                
                    Design_retino = [Design_retino;tempMapXPos,tempMapYPos,zeros(tempNum,1),jj*ones(tempNum,1)];
                    Y_retino = [Y_retino;squeeze(mapData{jj}{kk}(:,3))];
                    temp = squeeze(parameters(jj,kk,:));
                    tempmle = mle(mapData{jj}{kk}(:,3),'distribution','loglogistic');
                    temp(6) = tempmle(1);temp(7) = tempmle(2);
                    allParams = [allParams,temp];
                end
            end
        end
        figure(h(1));hold on;
        for jj=1:numChans
            xdata = atand(squeeze(parameters(:,jj,2)).*pix_to_degree{1});
            ydata = atand(squeeze(parameters(:,jj,3)).*pix_to_degree{1});
        
            xstd = atand(4.*squeeze(parameters(:,jj,4)).*pix_to_degree{1});
            xstdErr = atand(4.*squeeze(confIntervals(:,jj,4)).*pix_to_degree{1});
            ystd = atand(4.*squeeze(parameters(:,jj,5)).*pix_to_degree{1});
            ystdErr = atand(4.*squeeze(confIntervals(:,jj,5).*pix_to_degree{1}));
            
            height = exp(squeeze(parameters(:,jj,1))+squeeze(parameters(:,jj,6)));
            heightErr = exp(squeeze(parameters(:,jj,1))+squeeze(parameters(:,jj,6))+...
                squeeze(confIntervals(:,jj,1))+squeeze(confIntervals(:,jj,6)));
            
            xErr = atand(squeeze(confIntervals(:,jj,2)).*pix_to_degree{1});
            yErr = atand(squeeze(confIntervals(:,jj,3)).*pix_to_degree{1});
            
            
            xDist = xdata-xdata(1);yDist = ydata-ydata(1);
            for zz=2:length(xDist)
                if (parameters(zz,jj,1)-confIntervals(zz,jj,1)) > 0 && sum(squeeze(confIntervals(zz,jj,:))) < 1500
                    allRelativeDistToCenterMass_1 = [allRelativeDistToCenterMass_1;xDist(zz),yDist(zz)];
                else
                    xDist(zz) = -1000;
                    yDist(zz) = - 1000;
                end
            end
            for zz=1:length(xDist)
                if (parameters(zz,jj,1)-confIntervals(zz,jj,1)) > 0 && sum(squeeze(confIntervals(zz,jj,:))) < 1500
                    allStd_1 = [allStd_1;xstd(zz),xstdErr(zz),ystd(zz),ystdErr(zz),zz];
                    allHeight_1 = [allHeight_1;height(zz),heightErr(zz),zz];
                end
            end
            errorbar(xDist,yDist,yErr,yErr,xErr,xErr,linespecs{ConditionNumber},'LineWidth',2);hold on;
            title(sprintf('Center of Mass Relative to Day 1'));
            xlabel('Horizontal dva');
            ylabel('Vertical dva');
            axis([-20 20 -20 20]);
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
                if (squeeze(parameters(jj,kk,1))-squeeze(confIntervals(jj,kk,1))) > 0.05 
                    tempNum = length(mapData{jj}{kk}(:,3));
                    
                    tempMapXPos = squeeze(mapData{jj}{kk}(:,1));
                    tempMapYPos = squeeze(mapData{jj}{kk}(:,2));
                    tempMapXPos = tempMapXPos-squeeze(parameters(jj,kk,2));
                    tempMapYPos = tempMapYPos-squeeze(parameters(jj,kk,3));
                    

                    Design_retino = [Design_retino;tempMapXPos,tempMapYPos,ones(tempNum,1),jj*ones(tempNum,1)];
                    Y_retino = [Y_retino;mapData{jj}{kk}(:,3)];
                    temp = squeeze(parameters(jj,kk,:));
                    tempmle = mle(mapData{jj}{kk}(:,3),'distribution','loglogistic');
                    temp(6) = tempmle(1);temp(7) = tempmle(2);
                    allParams = [allParams,temp];
                end
            end
        end
        
        figure(h(1));hold on;
        for jj=1:numChans
            xdata = atand(squeeze(parameters(:,jj,2)).*pix_to_degree{1});
            ydata = atand(squeeze(parameters(:,jj,3)).*pix_to_degree{1});
            
            xstd = atand(4.*squeeze(parameters(:,jj,4)).*pix_to_degree{1});
            xstdErr = atand(4.*squeeze(confIntervals(:,jj,4)).*pix_to_degree{1});
            ystd = atand(4.*squeeze(parameters(:,jj,5)).*pix_to_degree{1});
            ystdErr = atand(4.*squeeze(confIntervals(:,jj,5).*pix_to_degree{1}));
            
            height = exp(squeeze(parameters(:,jj,1))+squeeze(parameters(:,jj,6)));
            heightErr = exp(squeeze(parameters(:,jj,1))+squeeze(parameters(:,jj,6))+...
                squeeze(confIntervals(:,jj,1))+squeeze(confIntervals(:,jj,6)));
        
            xErr = atand(squeeze(confIntervals(:,jj,2)).*pix_to_degree{1});
            yErr = atand(squeeze(confIntervals(:,jj,3)).*pix_to_degree{1});
            
            
            xDist = xdata-xdata(1);yDist = ydata-ydata(1);
            for zz=2:length(xDist)
                if (parameters(zz,jj,1)-confIntervals(zz,jj,1)) > 0 && sum(squeeze(confIntervals(zz,jj,:))) < 1500
                    allRelativeDistToCenterMass_2 = [allRelativeDistToCenterMass_2;xDist(zz),yDist(zz)];
                else
                    xDist(zz) = -1000;
                    yDist(zz) = -1000;
                end
            end
            for zz=1:length(xDist)
                if (parameters(zz,jj,1)-confIntervals(zz,jj,1)) > 0 && sum(squeeze(confIntervals(zz,jj,:))) < 1500
                    allStd_2 = [allStd_2;xstd(zz),xstdErr(zz),ystd(zz),ystdErr(zz),zz];
                    allHeight_2 = [allHeight_2;height(zz),heightErr(zz),zz];
                end
            end
            errorbar(xDist,yDist,yErr,yErr,xErr,xErr,linespecs{ConditionNumber},'LineWidth',2);hold on;
            title(sprintf('Center of Mass Relative to Day 1'));
            xlabel('Horizontal dva');
            ylabel('Vertical dva');
            axis([-20 20 -20 20]);
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
        [numChans,srp_reps] = size(srpSize{1});
        stimLen = size(srpVEP{1},2);
        vepSize = zeros(numDays,numChans,srp_reps);
        meanVEP = zeros(numDays,numChans,stimLen);
        for jj=1:numDays
            vepSize(jj,:,:) = srpSize{jj};
            meanVEP(jj,:,:) = srpVEP{jj};
        end
        
         for jj=1:numChans
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
figure(h(1));
legend(conditionNames{1},conditionNames{2});

% 
% figure(h(3));
% day1 = Design_srp(:,1)==1;
% day2 = Design_srp(:,1)==2;
% day3 = Design_srp(:,1)==3;
% day4 = Design_srp(:,1)==4;
% day5 = Design_srp(:,1)==5;
% 
% cond2 = Design_srp(:,2)==1;
% cond3 = Design_srp(:,2)==0;
% 
% barMean = [mean(Y_srp(day1 & cond2)),mean(Y_srp(day1 & cond3));...
%     mean(Y_srp(day2 & cond2)),mean(Y_srp(day2 & cond3));...
%     mean(Y_srp(day3 & cond2)),mean(Y_srp(day3 & cond3));...
%     mean(Y_srp(day4 & cond2)),mean(Y_srp(day4 & cond3));...
%     mean(Y_srp(day5 & cond2)),mean(Y_srp(day5 & cond3))];
% 
% barError = [std(Y_srp(day1 & cond2))/sqrt(length(Y_srp(day1 & cond2))),...
%     std(Y_srp(day1 & cond3))/sqrt(length(Y_srp(day1 & cond3)));...
%     std(Y_srp(day2 & cond2))/sqrt(length(Y_srp(day2 & cond2))),...
%     std(Y_srp(day2 & cond3))/sqrt(length(Y_srp(day2 & cond3)));...
%     std(Y_srp(day3 & cond2))/sqrt(length(Y_srp(day3 & cond2))),...
%     std(Y_srp(day3 & cond3))/sqrt(length(Y_srp(day3 & cond3)));...
%     std(Y_srp(day4 & cond2))/sqrt(length(Y_srp(day4 & cond2))),...
%     std(Y_srp(day4 & cond3))/sqrt(length(Y_srp(day4 & cond3)));...
%     std(Y_srp(day5 & cond2))/sqrt(length(Y_srp(day5 & cond2))),...
%     std(Y_srp(day5 & cond3))/sqrt(length(Y_srp(day5 & cond3)))];
% 
% days = [1:5;1:5];
% superbar(days,barMean','E',barError')
% title('SRP Protocol Mean VEP Magnitude');
% xlabel('Experimental Day');ylabel('Magnitude (\muV)');
% legend(conditionNames{2},conditionNames{3});
figure(h(3));
subplot(2,1,1);
histogram(allStd_1(:,1),10:5:70,'FaceColor','r');hold on;
histogram(allStd_2(:,1),10:5:70,'FaceColor','b');
legend(conditionNames{1},conditionNames{2});
title('Approximate Horizontal Map Size');
xlabel('4*sigma_x (dva)');
ylabel('Count');

subplot(2,1,2);
histogram(allStd_1(:,2),10:5:70,'FaceColor','r');hold on;
histogram(allStd_2(:,2),10:5:70,'FaceColor','b');
legend(conditionNames{1},conditionNames{2});
title('Approximate Vertical Map Size');
xlabel('4*sigma_y (dva)');
ylabel('Count');


figure(h(4));
subplot(2,1,1);
histogram(allRelativeDistToCenterMass_1(:,1),-10:1:10,'FaceColor','r');hold on;
histogram(allRelativeDistToCenterMass_2(:,1),-10:1:10,'FaceColor','b');
legend(conditionNames{1},conditionNames{2});
title('Horizontal Center of Mass (hCoM) Relative To Day 1');
xlabel('Distance to Day 1 hCoM (dva)');
ylabel('Count');

subplot(2,1,2);
histogram(allRelativeDistToCenterMass_1(:,2),-10:1:10,'FaceColor','r');hold on;
histogram(allRelativeDistToCenterMass_2(:,2),-10:1:10,'FaceColor','b');
legend(conditionNames{1},conditionNames{2});
title('Vertical Center of Mass (vCoM) Relative To Day 1');
xlabel('Distance to Day 1 vCoM (dva)');
ylabel('Count');
% 
figure(h(5));numDays = 5;
allData = zeros(numDays,3);
alpha = 0.05;

for ii=1:numDays
    dayData = Y_srp(Design_srp(:,1)==ii);
    meanVals = zeros(1000,1);
    for jj=1:3000
        indeces = random('Discrete Uniform',length(dayData),[length(dayData),1]);
        meanVals(jj) = mean(dayData(indeces));
    end
    allQuantiles = quantile(meanVals,[alpha/2,0.5,1-alpha/2]);
    allQuantiles(2) = mean(dayData);
    allData(ii,:) = allQuantiles;
end
subplot(2,1,1);boundedline(1:numDays,allData(:,2),[abs(allData(:,1)-allData(:,2)),abs(allData(:,3)-allData(:,2))]);
title('SRP Protocol (Combined Data)');ylabel('Max-Min Statistic Magnitude (\muV)');
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
    mean2 = zeros(1000,1);
    mean3 = zeros(1000,1);
    for jj=1:3000
       ind2 = random('Discrete Uniform',length(cond2temp),[length(cond2temp),1]);
       ind3 = random('Discrete Uniform',length(cond3temp),[length(cond3temp),1]);
        mean2(jj) = mean(cond2temp(ind2));
        mean3(jj) = mean(cond3temp(ind3));
    end
    cond2quantiles = quantile(mean2,[alpha/2,0.5,1-alpha/2]);
    cond3quantiles = quantile(mean3,[alpha/2,0.5,1-alpha/2]);
    cond2quantiles(2) = mean(cond2temp);
    cond3quantiles(2) = mean(cond3temp);
    cond2Data(ii,:) = cond2quantiles;
    cond3Data(ii,:) = cond3quantiles;
end

subplot(2,1,2);boundedline(1:numDays,cond2Data(:,2),[abs(cond2Data(:,1)-cond2Data(:,2)),abs(cond2Data(:,3)-cond2Data(:,2))],linespecs{2},...
       1:numDays,cond3Data(:,2),[abs(cond3Data(:,1)-cond3Data(:,2)),abs(cond3Data(:,3)-cond3Data(:,2))],linespecs{3});
title('SRP Protocol (Separated Data)');ylabel('Max-Min Statistic Magnitude (\muV)');
xlabel('Experimental Day');
legend(conditionNames{2},conditionNames{3},'location','northwest');


figure();
errorbar(allHeight_1(:,3),allHeight_1(:,1),abs(allHeight_1(:,2)-allHeight_1(:,1)),linespecs{1},'Linewidth',2);hold on;
errorbar(allHeight_2(:,3),allHeight_2(:,1),abs(allHeight_2(:,2)-allHeight_2(:,1)),linespecs{2},'Linewidth',2);
legend(conditionNames{1},conditionNames{2},'location','northwest');
title('Height Parameter Across Days');
% savefig(h,'MappingEffects_FinalResults.fig');

figure();
errorbar(allStd_1(:,end),allStd_1(:,1),allStd_1(:,2),linespecs{1},'Linewidth',2);hold on;
errorbar(allStd_2(:,end),allStd_2(:,3),allStd_2(:,4),linespecs{2},'LineWidth',2);
title('Map Size (Horizontal) Across Days');

save('MappingEffects_FinalResults.mat','Y_srp','Design_srp','Y_retino','Design_retino');
end

% xBins = linspace(-1500,1500,50);
% x = discretize(Design_retino(:,1),xBins);figure(1);
% for ii=1:length(unique(x))
% data = Y_retino(x==ii);xPos = mean(Design_retino(x==ii,1));
% figure(1);errorbar(xPos,mean(data),2*std(data)/sqrt(length(data)));hold on;
% end
% yBins = linspace(-1000,1000,25);
% y = discretize(Design_retino(:,2),yBins);figure(2);
% for ii=1:length(unique(y))
% data = Y_retino(y==ii);yPos = mean(Design_retino(y==ii,2));
% figure(2);errorbar(yPos,mean(data),2*std(data)/sqrt(length(data)));hold on;
% end
