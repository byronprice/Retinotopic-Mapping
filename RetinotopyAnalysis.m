% RetinotopyAnalysis.m
%  Use data from 2016 summer retinotopic mapping experiments with new
%   retinotopic region models (Gaussian distribution and binomial
%   distribution) ... try to figure which most accurately portrays the 
%   data
xaxis = 1:2500; yaxis = 1:1400;

files = dir('RetinoMap*.mat');
numFiles = length(files);

allGaussErrors = [];
allBinomErrors = [];

allGaussParams = [];
allBinomParams = [];
dataFile = [];
for zz=1:numFiles
   load(files(zz).name);
   
   try
       numChans = MapParams.numChans;
       centerVals = MapParams.centerVals;
       tempResponse = MapParams.Response;
       numStimuli = size(tempResponse,2);
       numReps = size(tempResponse,3);
       temp = MapParams.centerMass;
       
       centerMass = zeros(numChans,2);
       centerMass(:,1) = temp.x;
       centerMass(:,2) = temp.y;
       
       display(strcat('Running file-',files(zz).name));
       
       for ii=1:numChans
           if isnan(centerMass(ii,1)) == 0
               dataFile = [dataFile;files(zz).name];
               gaussResponse = zeros(1,numStimuli*numReps,2);
               binomResponse = zeros(1,numStimuli*numReps,2);
               count = 1;
               for jj=1:numStimuli
                   for kk=1:numReps
                       gaussResponse(1,count,1) = jj;
%                        gaussResponse(1,count,2) = max(tempResponse(ii,jj,kk,150:250))-min(tempResponse(ii,jj,kk,50:120));
                       gaussResponse(1,count,2) = min(tempResponse(ii,jj,kk,50:120));
                       binomResponse(1,count,1) = jj;
                       if gaussResponse(1,count,2) < -150
                           binomResponse(1,count,2) = 1;
                       end
                       count = count+1;
                   end
               end
               
               
               
               [gaussianParameters] = FitLFPGaussRetinoModel(gaussResponse,xaxis,yaxis,centerVals);
               allGaussParams = [allGaussParams;gaussianParameters];
               tempX = gaussianParameters(1,2);
               tempY = gaussianParameters(1,3);
               
               errorX = abs(tempX-centerMass(ii,1))/centerMass(ii,1);
               errorY = abs(tempY-centerMass(ii,2))/centerMass(ii,2);
               gaussError = (errorX+errorY)/2;display(gaussError);
               allGaussErrors = [allGaussErrors,gaussError];
               
%                [bernoulliParameters] = FitLFPretinoModel(binomResponse,xaxis,yaxis,0.2,centerVals);
%                allBinomParams = [allBinomParams;bernoulliParameters];
%                tempX = bernoulliParameters(1,2);
%                tempY = bernoulliParameters(1,3);
%                
%                errorX = abs(tempX-centerMass(ii,1))/centerMass(ii,1);
%                errorY = abs(tempY-centerMass(ii,2))/centerMass(ii,2);
%                binomError = (errorX+errorY)/2;display(binomError);
%                allBinomErrors = [allBinomErrors,binomError];
           end
       end
   catch
       continue;
   end
end

save('MaxMin_BiggestModelCheckResults.mat','allGaussParams','allBinomParams','allGaussErrors','allBinomErrors','dataFile');
% finalIm = zeros(length(xaxis),length(yaxis));
% for ii=1:length(xaxis)
%     for jj=1:length(yaxis)
%         distance = DistFun(finalParameters(1,2:3),[ii,jj]);
%         finalIm(ii,jj) = hyperParameterFun([finalParameters(1,1),finalParameters(1,4),-120],distance);
%     end
% end
% figure();imagesc(finalIm');set(gca,'YDir','normal');
