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
                    maxMinAccuracy = zeros(length(threshold),1);
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
                        maxMinAccuracy(ii) = (sum(combined==0)+sum(combined==2))/length(combined);
                    end
                    figure(1);subplot(2,1,1);plot(falsePositives,truePositives,'LineWidth',2);
                    title(sprintf('ROC Curve: Max[%d,%d]-Min[%d,%d]',positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end)));
                    figure(1);subplot(2,1,2);plot(threshold,maxMinAccuracy,'LineWidth',2);axis([threshold(1) threshold(end) 0 1]);
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

[~,index] = max(maxMinAccuracy);
maxMinThreshold = threshold(index);

% test on the test data
testStat = max(testVEPs(:,positivityTime),[],2)-min(testVEPs(:,negativityTime),[],2);
labels = testStat>maxMinThreshold;

combined = labels+testLabels;
truePos = sum(combined==2)/sum(testLabels==1);

trueNegs = testLabels==0;
labelledPos = labels==1;
tempCombined = trueNegs+labelledPos;

falsePos = sum(tempCombined==2)/sum(trueNegs);
maxMinAccuracy = (sum(combined==0)+sum(combined==2))/length(combined);

fprintf('Max[%d,%d]-Min[%d,%d] Statistic Results\n\n',positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end));
fprintf('True Positive Rate: %3.3f\n',truePos);
fprintf('False Positive Rate: %3.3f\n',falsePos);
fprintf('Accuracy: %3.3f\n\n',maxMinAccuracy);


% TEST 2 -  (median(before) + max - min) within windows statistic
% posMin = 150:10:250;
% posMax = 200:10:350;
% negMin = 60:10:150;
% negMax = 100:10:200;

% fullAccuracy = zeros(length(posMin),length(posMax),length(negMin),length(negMax));
% for aa=1:length(posMin)
%     for bb=1:length(posMax)
%         for cc=1:length(negMin)
%             for dd=1:length(negMax)
                
                positivityTime = 150:270;
                negativityTime = 120:150;
                medianTime = 1:(120-50);
                
%                 if isempty(positivityTime) == 0 && isempty(negativityTime) == 0
                    %train on training data first
                    vepStat = median(trainingVEPs(:,medianTime),2)+...
                        max(trainingVEPs(:,positivityTime),[],2)-min(trainingVEPs(:,negativityTime),[],2);
                    vepLabels = trainingLabels;
                    
                    threshold = 50:1000;
                    truePositives = zeros(length(threshold),1);
                    falsePositives = zeros(length(threshold),1);
                    medianMaxMinAccuracy = zeros(length(threshold),1);
                    for ii=1:length(threshold)
                        labels = vepStat>threshold(ii);
                        combined = labels+vepLabels;
                        truePositives(ii) = sum(combined==2)/sum(vepLabels==1);
                        
                        trueNegs = vepLabels==0;
                        labelledPos = labels==1;
                        tempCombined = trueNegs+labelledPos;
                        falsePositives(ii) = sum(tempCombined==2)/sum(trueNegs);
                        medianMaxMinAccuracy(ii) = (sum(combined==0)+sum(combined==2))/length(combined);
                    end
                    figure(2);subplot(2,1,1);plot(falsePositives,truePositives,'LineWidth',2);
                    title(sprintf('ROC Curve: Median[1,%d]+Max[%d,%d]-Min[%d,%d]',negativityTime(1)-50,positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end)));
                    figure(2);subplot(2,1,2);plot(threshold,medianMaxMinAccuracy,'LineWidth',2);axis([threshold(1) threshold(end) 0 1]);
                    title(sprintf('Accuracy: Median[1,%d]+Max[%d,%d]-Min[%d,%d]',negativityTime(1),positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end)));
%                     fullAccuracy(aa,bb,cc,dd) = max(mvnAccuracy);
%                 else
%                     fullAccuracy(aa,bb,cc,dd) = 0.5;
%                 end
%             end
%         end
%     end
% end
% [maxAccuracy,I] = max(fullAccuracy(:));
% [I1,I2,I3,I4] = ind2sub(size(fullAccuracy),I);


% maximum accuracy of 0.665 achieved with positivity times from 150:270
%   and negativity times from 130:150, a threshold of 298, and a false
%   positive rate of 0.30

[~,index] = max(medianMaxMinAccuracy);
medianMaxMinThreshold = threshold(index);

% test on the test data
testStat = median(testVEPs(:,medianTime),2)+...
                        max(testVEPs(:,positivityTime),[],2)-min(testVEPs(:,negativityTime),[],2);
labels = testStat>medianMaxMinThreshold;

combined = labels+testLabels;
truePos = sum(combined==2)/sum(testLabels==1);

trueNegs = testLabels==0;
labelledPos = labels==1;
tempCombined = trueNegs+labelledPos;

falsePos = sum(tempCombined==2)/sum(trueNegs);
medianMaxMinAccuracy = (sum(combined==0)+sum(combined==2))/length(combined);

fprintf('Median[1,%d]+Max[%d,%d]-Min[%d,%d] Statistic Results\n\n',negativityTime(1)-50,positivityTime(1),positivityTime(end),negativityTime(1),negativityTime(end));
fprintf('True Positive Rate: %3.3f\n',truePos);
fprintf('False Positive Rate: %3.3f\n',falsePos);
fprintf('Accuracy: %3.3f\n\n',medianMaxMinAccuracy);

% TEST 3 - Multivariate normal distribution for Bayesian hypothesis testing
%   fit a multivariate normal to the training data for the VEPs and a
%   second distribution to the training data for the grey screen ... then
%   calculate the likelihood ratio: P(D | D is a grey response)
%    / P(D | D is a VEP)  ... try a variety of threshold across possible
%   likelihood ratios and test with ROC, accuracy, etc.

downSampledVEPs = trainingVEPs(:,1:3:end);
tempVepMu = mean(downSampledVEPs(trainingLabels==1,:),1)';
tempVepCov = cov(downSampledVEPs(trainingLabels==1,:));
tempGreyMu = mean(downSampledVEPs(trainingLabels==0,:),1)';
tempGreyCov = cov(downSampledVEPs(trainingLabels==0,:));

threshold = -5:1:5;
truePositives = zeros(length(threshold),1);
falsePositives = zeros(length(threshold),1);
mvnAccuracy = zeros(length(threshold),1);
for ii=1:length(threshold)
    statistic = zeros(length(testVEPs),1);
    for jj=1:length(testVEPs)
        statistic(jj) = log(mvnpdf(testVEPs(jj,1:3:end)',tempVepMu,tempVepCov)/mvnpdf(testVEPs(jj,1:3:end)',tempGreyMu,tempGreyCov));
    end
    labels = statistic>threshold(ii);
    combined = labels+testLabels;
    truePositives(ii) = sum(combined==2)/sum(testLabels==1);
    
    trueNegs = testLabels==0;
    labelledPos = labels==1;
    tempCombined = trueNegs+labelledPos;
    falsePositives(ii) = sum(tempCombined==2)/sum(trueNegs);
    mvnAccuracy(ii) = (sum(combined==0)+sum(combined==2))/length(combined);
end

figure(3);subplot(2,1,1);plot(falsePositives,truePositives,'LineWidth',2);
title(sprintf('ROC Curve: Multivariate Normal'));
figure(3);subplot(2,1,2);plot(threshold,mvnAccuracy,'LineWidth',2);axis([threshold(1) threshold(end) 0 1]);

[mvnAccuracy,index] = max(mvnAccuracy);
mvnThreshold = threshold(index);

fprintf('Multivariate Normal Results\n');
fprintf('True Positive Rate: %3.3f\n',truePositives(index));
fprintf('False Positive Rate: %3.3f\n',falsePositives(index));
fprintf('Accuracy: %3.3f\n\n',mvnAccuracy);


% TEST 4 - train a perceptron to output "grey" or "VEP" on the training
% data and then test it's accuracy on the test data

% CREATE THE NETWORK WITH RANDOMIZED WEIGHTS AND BIASES
% [N,stimLen] = size(trainingVEPs);
% 
% numOutputs = 2;
% numHidden1 = 200;
% numHidden2 = 50;
% numHidden3 = 10;
% myNet = Network([stimLen,numHidden1,numHidden2,numHidden3,numOutputs]); % from a function
% % in this directory, builds a 4-layer network
% 
% % [Output,~] = Feedforward(trainingVEPs(1,:)',myNet);
% % Output{end}
% % trainingLabels(1)
% % get desired output matrix
% DesireOutput = zeros(numOutputs,N);
% 
% for ii=1:N
%     if trainingLabels(ii) == 1
%         DesireOutput(1,ii) = 1;
%         DesireOutput(2,ii) = 0;
%     elseif trainingLabels(ii) == 0
%         DesireOutput(1,ii) = 0;
%         DesireOutput(2,ii) = 1;
%     end
% end
% 
% % STOCHASTIC GRADIENT DESCENT
% batchSize = 10; % make mini batches and run the algorithm
% % on those "runs" times
% runs = 1e6;

etaRange = [1e-4,1e-3,1e-2,1e-1,1e0,0.5e1,1e1];
lambdaRange = [1e-4,1e-3,1e0,0.5e1,1e1,1e2];
perceptronAccuracy = zeros(length(etaRange),length(lambdaRange));
for bigEta = 1:length(etaRange)
    for bigLambda = 1:length(lambdaRange)
        downSampledVEPs = trainingVEPs(:,1:2:end);
        [N,stimLen] = size(downSampledVEPs);
        downSampledVEPs = downSampledVEPs./2000;
        downSampledVEPs(downSampledVEPs>1) = 1;
        
        numOutputs = 2;
        numHidden1 = 150;
        numHidden2 = 50;
        numHidden3 = 10;
        myNet = Network([stimLen,numHidden1,numHidden2,numHidden3,numOutputs]); % from a function
        % in this directory, builds a 4-layer network
        
        % [Output,~] = Feedforward(trainingVEPs(1,:)',myNet);
        % Output{end}
        % trainingLabels(1)
        % get desired output matrix
        DesireOutput = zeros(numOutputs,N);
        
        for ii=1:N
            if trainingLabels(ii) == 1
                DesireOutput(1,ii) = 1;
                DesireOutput(2,ii) = 0;
            elseif trainingLabels(ii) == 0
                DesireOutput(1,ii) = 0;
                DesireOutput(2,ii) = 1;
            end
        end
        
        % STOCHASTIC GRADIENT DESCENT
        batchSize = 10; % make mini batches and run the algorithm
        % on those "runs" times
        runs = 1e5;

        numCalcs = size(myNet.Weights,2);
        dCostdWeight = cell(1,numCalcs);
        dCostdBias = cell(1,numCalcs);
        
        for ii=1:runs
            indeces = ceil(rand([batchSize,1]).*(N-1));
            for jj=1:numCalcs
                layer1 = size(myNet.Weights{jj},1);
                layer2 = size(myNet.Weights{jj},2);
                dCostdWeight{jj} = zeros(layer1,layer2);
                dCostdBias{jj} = zeros(layer2,1);
            end
            for jj=1:batchSize
                index = indeces(jj);
                [costweight,costbias] = BackProp(downSampledVEPs(index,:)',myNet,...
                    DesireOutput(:,index));
                for kk=1:numCalcs
                    dCostdWeight{kk} = dCostdWeight{kk}+costweight{kk};
                    dCostdBias{kk} = dCostdBias{kk}+costbias{kk};
                end
            end
            [myNet] = GradientDescent(myNet,dCostdWeight,dCostdBias,batchSize,etaRange(bigEta),N,lambdaRange(bigLambda));

            
            clear indeces;% dCostdWeight dCostdBias;
            check = isnan(myNet.Weights{1});
            if sum(sum(check)) > 0
               display(myRange(bigEta));
               display(myRange(bigLambda));
               break; 
            end
        end
        [Output,~] = Feedforward(downSampledVEPs(13456,:)',myNet);
         display(Output{end})
         display(DesireOutput(:,13456))
        
        % COMPARE ON TEST DATA
        downSampledTestVEPs = testVEPs(:,1:2:end);
        downSampledTestVEPs = downSampledTestVEPs./2000;
        downSampledTestVEPs(downSampledTestVEPs>1) = 1;
        [N,~] = size(downSampledTestVEPs);
        numOutputs = 2;
        
        DesireOutput = zeros(numOutputs,N);
        
        for ii=1:N
            if testLabels(ii) == 1
                DesireOutput(1,ii) = 1;
                DesireOutput(2,ii) = 0;
            elseif testLabels(ii) == 0
                DesireOutput(1,ii) = 0;
                DesireOutput(2,ii) = 1;
            end
        end
        
        classifiedVals = zeros(N,1);
        count = 0;
        for ii=1:N
            [Output,~] = Feedforward(downSampledTestVEPs(ii,:)',myNet);
            [~,realVal] = max(DesireOutput(:,ii));
            [~,netVal] = max(Output{end});
            classifiedVals(ii) = netVal;
            if realVal == netVal
                count = count+1;
            end
        end
        perceptronAccuracy(bigEta,bigLambda) = count/N;
        fprintf('Acurracy: %3.2f\n',count/N);
    end
end
save('OptimizeResults.mat','maxMinAccuracy','maxMinThreshold','medianMaxMinAccuracy','medianMaxMinThreshold',...
    'mvnAccuracy','mvnThreshold','perceptronAccuracy','myNet');
end