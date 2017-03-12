% RetinotopyAnalysis.m
%  Use data from 2016 summer retinotopic mapping experiments with new
%   retinotopic region models (Gaussian distribution and binomial
%   distribution) ... try to figure which most accurately portrays the 
%   data
%  fileStart = 'RetinoMap*.mat';
function [] = RetinotopyAnalysis(fileStart)
xaxis = 1:2560; yaxis = 1:1440;

files = dir(fileStart);
numFiles = length(files);

numModels = 9;
numParameters = [2,6,7,8,8,10,11,10,9];
AIC = zeros(numModels,1);
for zz=1:numFiles
   load(files(zz).name);
   temp = files(zz).name;
   index = regexp(temp,'_');
   date = temp(10:index-1);
   name = temp(index+1:end-4);
   

       load(sprintf('RetinoStim%s_%s.mat',date,name));
       [numChans,numStimuli,stimLen] = size(vepResponse);    
       
       display(strcat('Running file-',files(zz).name));
       negParams = [50,120];
       posParams = [120,stimLen];
       
       for ii=1:numChans
           check = sum(ninetyfiveErrors(ii,:));
           
           if check < 2000
               %                dataFile = [dataFile;files(zz).name];
               gaussResponse = zeros(numStimuli,1);
               
%                xPos = unique(centerVals(:,1));yPos = unique(centerVals(:,2));
%                VEPs = zeros(length(xPos),length(yPos),numReps);
               count = 1;
               for jj=1:numStimuli
                       [minVal,~] = min(vepResponse(ii,jj,negParams(1):negParams(2)));
                       gaussResponse(count,1) = max(vepResponse(ii,jj,posParams(1):posParams(2)))...
                           -minVal;
                       %                        gaussResponse(1,count,2) = -min(tempResponse(ii,jj,kk,50:120));
                       
%                        x = centerVals(jj,1);y = centerVals(jj,2);
%                        xInd = find(xPos == x);yInd = find(yPos == y);
%                        VEPs(xInd,yInd,kk) = gaussResponse(count,1);
                       count = count+1;
               end
               
               % %                    figure();count = 1;edges = 0:50:1000;
               %                    meanVals = zeros(numStimuli,1);
               %                    count = 1;
               %                    for kk=length(yPos):-1:1
               %                        for jj=1:length(xPos)
               %                            temp = squeeze(VEPs(jj,kk,:));
               % %                            subplot(length(yPos),length(xPos),count);histogram(temp,edges);
               %                            stdev = std(temp);meanVal = mean(temp);
               % %                            title(sprintf('%3.0f  %3.0f',meanVal,stdev));
               %                            meanVals(count) = meanVal;
               %                            count = count+1;
               %                        end
               %                    end
               [littleAIC] = FitModel1(gaussResponse,xaxis,yaxis,centerVals);
               AIC(1) = AIC(1)+littleAIC;
               [littleAIC] = FitModel2(gaussResponse,xaxis,yaxis,centerVals);
               AIC(2) = AIC(2)+littleAIC;
               [littleAIC] = FitModel3(gaussResponse,xaxis,yaxis,centerVals);
               AIC(3) = AIC(3)+littleAIC;
               [littleAIC] = FitModel4(gaussResponse,xaxis,yaxis,centerVals);
               AIC(4) = AIC(4)+littleAIC;
               [littleAIC] = FitModel5(gaussResponse,xaxis,yaxis,centerVals);
               AIC(5) = AIC(5)+littleAIC;
               [littleAIC] = FitModel6(gaussResponse,xaxis,yaxis,centerVals);
               AIC(6) = AIC(6)+littleAIC;
               [littleAIC] = FitModel7(gaussResponse,xaxis,yaxis,centerVals);
               AIC(7) = AIC(7)+littleAIC;
               [littleAIC] = FitModel8(gaussResponse,xaxis,yaxis,centerVals);
               AIC(8) = AIC(8)+littleAIC;
               [littleAIC] = FitModel9(gaussResponse,xaxis,yaxis,centerVals);
               AIC(9) = AIC(9)+littleAIC;
           end
       end

end
numParameters
AIC
figure();plot(numParameters,AIC,'b','LineWidth',2);title('Retinotopic-Map Model Comparison: AIC');
xlabel('Number of Parameters');ylabel('AIC');
save('LFPRetinotopy_ModelComparison.mat','AIC','numParameters');
end

% finalIm = zeros(length(xaxis),length(yaxis));
% for ii=1:length(xaxis)
%     for jj=1:length(yaxis)
%         distance = DistFun(finalParameters(1,2:3),[ii,jj]);
%         finalIm(ii,jj) = hyperParameterFun([finalParameters(1,1),finalParameters(1,4),-120],distance);
%     end
% end
% figure();imagesc(finalIm');set(gca,'YDir','normal');

