function [] = OptimizeVEPstat()
% OptimizeVEPstat.m
%   ROC curves to determine which metric best discriminates a VEP
%    from a not-VEP

load('VEPS_For_Comparison_With_Grey_Screen.mat');

% TEST 1 - simple (max - min) within a window statistic
% posMin = 120:10:200;
% posMax = 250:10:350;
% negMin = 50:10:120;
% negMax = 100:10:180;

% fullAccuracy = zeros(length(posMin),length(posMax),length(negMin),length(negMax));
% for aa=1:length(posMin)
%     for bb=1:length(posMax)
%         for cc=1:length(negMin)
%             for dd=1:length(negMax)

                positivityTime = 180:270;
                negativityTime = 120:160;

% train on training data first
                vepStat = max(trainingVEPs(:,positivityTime),[],2)-min(trainingVEPs(:,negativityTime),[],2);
%                 vepLabels = trainingLabels(vepStat>0);vepStat = vepStat(vepStat>0);
                vepLabels = trainingLabels;
                
                if isempty(vepStat) == 0
                    threshold = 100:1500;
                    truePositives = zeros(length(threshold),1);
                    falsePositives = zeros(length(threshold),1);
                    accuracy = zeros(length(threshold),1);
                    for ii=1:length(threshold)
                        %                     labels = zeros(length(vepStat),1);
                        %                     for jj=1:length(vepStat)
                        %                         if vepStat(jj) <= threshold(ii)
                        %                             labels(jj) = 0;
                        %                         else
                        %                             labels(jj) = 1;
                        %                         end
                        %                     end
                        labels = vepStat>threshold(ii);
                        combined = labels+vepLabels;
                        truePositives(ii) = sum(combined==2)/sum(vepLabels==1);
                        
                        trueNegs = vepLabels==0;
                        labelledPos = labels==1;
                        tempCombined = trueNegs+labelledPos;
                        falsePositives(ii) = sum(tempCombined==2)/sum(trueNegs);
                        accuracy(ii) = (sum(combined==0)+sum(combined==2))/length(combined);
                    end
                    figure(1);subplot(2,1,1);plot(falsePositives,truePositives,'LineWidth',2);
                    title(sprintf('ROC Curve: Max[%d,%d]-Min[%d,%d]',positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end)));
                    figure(1);subplot(2,1,2);plot(threshold,accuracy,'LineWidth',2);axis([threshold(1) threshold(end) 0 1]);
                    title(sprintf('Accuracy: Max[%d,%d]-Min[%d,%d]',positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end)));
%                 fullAccuracy(aa,bb,cc,dd) = max(accuracy);
                else
%                     fullAccuracy(aa,bb,cc,dd) = 0.5;
                end
%             end
%         end
%     end
% end
% [maxAccuracy,I] = max(fullAccuracy(:));
% [I1,I2,I3,I4] = ind2sub(size(fullAccuracy),I);
% maximum accuracy of 0.6796 achieved with positivity times from 180:270
%   and negativity times from 120:150, a threshold of 298, and a false
%   positive rate of 0.30

[~,index] = max(accuracy);
testThreshold = threshold(index);

% test on the test data
testStat = max(testVEPs(:,positivityTime),[],2)-min(testVEPs(:,negativityTime),[],2);
labels = testStat>testThreshold;

combined = labels+testLabels;
truePos = sum(combined==2)/sum(testLabels==1);

trueNegs = testLabels==0;
labelledPos = labels==1;
tempCombined = trueNegs+labelledPos;

falsePos = sum(tempCombined==2)/sum(trueNegs);
accuracy = (sum(combined==0)+sum(combined==2))/length(combined);

fprintf('Max[%d,%d]-Min[%d,%d] Statistic Results\n\n',positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end));
fprintf('True Positive Rate: %3.3f\n',truePos);
fprintf('False Positive Rate: %3.3f\n',falsePos);
fprintf('Accuracy: %3.3f\n',accuracy);


% TEST 2 -  (median(before) + max - min) within windows statistic
posMax = 200:10:350;
negMin = 60:10:150;
negMax = 100:10:200;

fullAccuracy = zeros(length(posMax),length(negMin),length(negMax));
for bb=1:length(posMax)
    for cc=1:length(negMin)
        for dd=1:length(negMax)
            
            positivityTime = negMax(dd):posMax(bb);
            negativityTime = negMin(cc):negMax(dd);
            medianTime = 1:negMin-50;
            
            %train on training data first
            vepStat = median(trainingVEPs(:,medianTime),2)+...
                max(trainingVEPs(:,positivityTime),[],2)-min(trainingVEPs(:,negativityTime),[],2);
            vepLabels = trainingLabels;
            
            if isempty(vepStat) == 0
                threshold = 50:1000;
                truePositives = zeros(length(threshold),1);
                falsePositives = zeros(length(threshold),1);
                accuracy = zeros(length(threshold),1);
                for ii=1:length(threshold)
                    labels = vepStat>threshold(ii);
                    combined = labels+vepLabels;
                    truePositives(ii) = sum(combined==2)/sum(vepLabels==1);
                    
                    trueNegs = vepLabels==0;
                    labelledPos = labels==1;
                    tempCombined = trueNegs+labelledPos;
                    falsePositives(ii) = sum(tempCombined==2)/sum(trueNegs);
                    accuracy(ii) = (sum(combined==0)+sum(combined==2))/length(combined);
                end
%                 figure(2);subplot(2,1,1);plot(falsePositives,truePositives,'LineWidth',2);
%                 title(sprintf('ROC Curve: Max[%d,%d]-Min[%d,%d]',positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end)));
%                 figure(2);subplot(2,1,2);plot(threshold,accuracy,'LineWidth',2);axis([threshold(1) threshold(end) 0 1]);
%                 title(sprintf('Accuracy: Max[%d,%d]-Min[%d,%d]',positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end)));
                fullAccuracy(bb,cc,dd) = max(accuracy);
            else
                fullAccuracy(bb,cc,dd) = 0.5;
            end
        end
    end
end
[maxAccuracy,I] = max(fullAccuracy(:));
[I1,I2,I3] = ind2sub(size(fullAccuracy),I);

save('OptimizeResults.mat','maxAccuracy','I1','I2','I3')
% maximum accuracy of 0.6796 achieved with positivity times from 180:270
%   and negativity times from 120:150, a threshold of 298, and a false
%   positive rate of 0.30

% [~,index] = max(accuracy);
% testThreshold = threshold(index);
% 
% % test on the test data
% testStat = max(testVEPs(:,positivityTime),[],2)-min(testVEPs(:,negativityTime),[],2);
% labels = testStat>testThreshold;
% 
% combined = labels+testLabels;
% truePos = sum(combined==2)/sum(testLabels==1);
% 
% trueNegs = testLabels==0;
% labelledPos = labels==1;
% tempCombined = trueNegs+labelledPos;
% 
% falsePos = sum(tempCombined==2)/sum(trueNegs);
% accuracy = (sum(combined==0)+sum(combined==2))/length(combined);
% 
% fprintf('Max[%d,%d]-Min[%d,%d] Statistic Results\n\n',positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end));
% fprintf('True Positive Rate: %3.3f\n',truePos);
% fprintf('False Positive Rate: %3.3f\n',falsePos);
% fprintf('Accuracy: %3.3f\n',accuracy);

% downSampledVeps = trainingVEPs(:,1:3:end);
% tempVepMu = mean(downSampledVEPs(trainingLabels==1,:),1);
% tempVepCov = cov(downSampledVEPs(trainingLabels==1,:));
% tempGreyMu = mean(downSampledVEPs(trainingLabels==0,:),1);
% tempGreyCov = cov(downSampledVEPs(trainingLabels==0,:));
% 
% threshold = 0.01:0.01:0.99;
% truePositives = zeros(length(threshold),1);
% falsePositives = zeros(length(threshold),1);
% accuracy = zeros(length(threshold),1);
% for ii=1:length(threshold)
%     
%     
% end

end