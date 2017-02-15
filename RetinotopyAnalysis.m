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

AIC = zeros(3,1);numParameters = [6,7,8];
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
%        filename = files(zz).name;filename = strcat(filename(1:end-4),'.fig');
%        openfig(filename);
       for ii=1:numChans
           if isnan(centerMass(ii,1)) == 0
%                dataFile = [dataFile;files(zz).name];
               gaussResponse = zeros(1,numStimuli*numReps,2);
               binomResponse = zeros(1,numStimuli*numReps,2);
               
               xPos = unique(centerVals(:,1));yPos = unique(centerVals(:,2));
               VEPs = zeros(length(xPos),length(yPos),numReps);
               count = 1;
               for jj=1:numStimuli
                   for kk=1:numReps
                       gaussResponse(1,count,1) = jj;
%                        gaussResponse(1,count,2) = max(tempResponse(ii,jj,kk,150:250))-min(tempResponse(ii,jj,kk,50:120));
                       gaussResponse(1,count,2) = min(tempResponse(ii,jj,kk,50:120));
                       
                       x = centerVals(jj,1);y = centerVals(jj,2);
                       xInd = find(xPos == x);yInd = find(yPos == y);
                       VEPs(xInd,yInd,kk) = gaussResponse(1,count,2);
                       binomResponse(1,count,1) = jj;
                       if gaussResponse(1,count,2) < -150
                           binomResponse(1,count,2) = 1;
                       end
                       count = count+1;
                   end
               end
               
%                figure();count = 1;edges = 0:0.5:10;
%                for kk=length(yPos):-1:1
%                    for jj=1:length(xPos)
%                       temp = log(abs(squeeze(VEPs(jj,kk,:))));
%                       subplot(length(yPos),length(xPos),count);histogram(temp,edges);
%                       stdev = std(temp);meanVal = mean(temp);
%                       title(sprintf('%3.0f  %3.0f',stdev,meanVal));
%                       count = count+1;
%                    end
%                end
                
%                [littleAIC] = FitModel1(gaussResponse,xaxis,yaxis,centerVals);
%                AIC(1) = AIC(1)+littleAIC;
%                
%                [littleAIC] = FitModel2(gaussResponse,xaxis,yaxis,centerVals);
%                AIC(2) = AIC(2)+littleAIC;

%                [littleAIC] = FitModel3(gaussResponse,xaxis,yaxis,centerVals);
%                AIC(3) = AIC(3)+littleAIC;
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
% figure();plot(numParameters,AIC,'b','LineWidth',2);title('Retinotopic-Map Model Comparison: AIC');
% xlabel('Number of Parameters');ylabel('AIC');
% save('LFPRetinotopy_ModelComparison.mat','AIC','numParameters');
end

% finalIm = zeros(length(xaxis),length(yaxis));
% for ii=1:length(xaxis)
%     for jj=1:length(yaxis)
%         distance = DistFun(finalParameters(1,2:3),[ii,jj]);
%         finalIm(ii,jj) = hyperParameterFun([finalParameters(1,1),finalParameters(1,4),-120],distance);
%     end
% end
% figure();imagesc(finalIm');set(gca,'YDir','normal');

function [AIC] = FitModel1(Response,xaxis,yaxis,centerVals)
%FitModel1

% hyperParameterFun = @(b,distX,distY) (b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4));

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 6;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
numRepeats = 1000;
maxITER = 500;
tolerance = 1e-3;
h = ones(numParameters,1)./100;

% parameters are 
%  1) b(1) - rise of map at center
%  2)+3) the position (x and y) of the center of mass
%  4) b(2) - standard deviation or spread of the map in x
%  5) b(3) - standard deviation or spread of map in y
%  6) rise for sigma, from the Gaussian likelihood, at map center
%  7) b(4) - peak negativity at edges of retinotopic region
%    b(1)+b(4) = peak negativity at retinotopic center of mass
%  8) rho - allows for elliptical contours
%      9) sigma at edges of retinotopic region
Bounds = [-1000,0;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,2000;1,1000;-1000,0];

% display('Steepest Ascent ...');
for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = squeeze(Response(zz,:,2));
    
    bigParameterVec = zeros(numRepeats,numParameters);
    bigLogLikelihoods = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions
    parfor repeats = 1:numRepeats
        gradientVec = zeros(numParameters,1);
        parameterVec = zeros(maxITER,numParameters);
        logLikelihood = zeros(maxITER,1);
        %parameterVec(1,:) = squeeze(bigParameterVec(repeats,:));
        for ii=1:numParameters
            parameterVec(1,ii) = Bounds(ii,1)+(Bounds(ii,2)-Bounds(ii,1)).*rand;
        end
        
        % for each starting position, do maxITER iterations
        for iter=2:maxITER
            % calculate likelihood at the current position in parameter
            %  space
            [logLikelihood(iter-1)] = GetLikelihood1(reps,parameterVec(iter-1,:),peakNegativity,flashPoints);
            
            if iter > 5
                check = sum(diff(logLikelihood(iter-2:iter-1)));
                if check < tolerance
                    break;
                end
            end
%             
%             % check if a random jump along one dimension improves the
%             %  likelihood
            temp = parameterVec(iter-1,:);
            randInd = random('Discrete Uniform',numParameters,1);
            temp(randInd) = Bounds(randInd,1)+(Bounds(randInd,2)-Bounds(randInd,1)).*rand;
            
            [likely] = GetLikelihood1(reps,temp,peakNegativity,flashPoints);
            if likely > logLikelihood(iter-1) %|| mod(iter,200) == 0
                parameterVec(iter-1,:) = temp';
                logLikelihood(iter-1) = likely;
            end
            
            % calculate the gradient by calculating the likelihood
            %  after moving over a small step h along each
            %  dimension
            for jj=1:numParameters
                tempParameterVec = parameterVec(iter-1,:);
                tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
                [gradLikelihoodplus] = GetLikelihood1(reps,tempParameterVec,peakNegativity,flashPoints);
                
                tempParameterVec = parameterVec(iter-1,:);
                tempParameterVec(jj) = tempParameterVec(jj)-h(jj);
                [gradLikelihoodminus] = GetLikelihood1(reps,tempParameterVec,peakNegativity,flashPoints);
                gradientVec(jj) = (gradLikelihoodplus-gradLikelihoodminus)./(2*h(jj));
            end
            
            % line search to get distance to move along gradient
            alpha = [0,1e-8,1e-6,1e-4,1e-2,1e0,1e1];
            lineSearchLikelihoods = zeros(length(alpha),1);
            lineSearchLikelihoods(1) = logLikelihood(iter-1);
            
            for ii=2:length(alpha)
                tempParameterVec = parameterVec(iter-1,:)'+gradientVec.*alpha(ii);
                [lineSearchLikelihoods(ii)] = GetLikelihood1(reps,tempParameterVec,peakNegativity,flashPoints);
            end
            
            newInds = find(imag(lineSearchLikelihoods)==0);
            newValues = zeros(length(newInds),1);
            for ii=1:length(newInds)
                newValues(ii) = lineSearchLikelihoods(newInds(ii));
            end
            [~,ind] = max(newValues);
            ind = newInds(ind);
            
            for jj=1:numParameters
                parameterVec(iter,jj) = max(Bounds(jj,1),min(parameterVec(iter-1,jj)+alpha(ind)*gradientVec(jj),Bounds(jj,2)));
            end
        end
        bigParameterVec(repeats,:) = parameterVec(iter-1,:);
        bigLogLikelihoods(repeats) = logLikelihood(iter-1);
    end
    [~,index] = max(bigLogLikelihoods);
    finalParameters(zz,:) = bigParameterVec(index,:);
    AIC = 2*numParameters-2*bigLogLikelihoods(index);
%     display(finalParameters(zz,:));
end

end

function [loglikelihood] = GetLikelihood1(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    distX = flashPoints(kk,1)-parameterVec(2);
    distY = flashPoints(kk,2)-parameterVec(3);
    b = [parameterVec(1),parameterVec(4),parameterVec(6)];
    mu = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(2)*b(2)))+b(3);
    stdev = parameterVec(5);%*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+parameterVec(9);
    loglikelihood = loglikelihood-(1/2)*log(2*pi*stdev*stdev)-(1/(2*stdev*stdev))*(peakNegativity(kk)-mu).^2;
%     summation = summation+(peakNegativity(kk)-mu).^2;
end

% loglikelihood = (-reps/2)*log(2*pi*parameterVec(6)*parameterVec(6))-...
%     (1/(2*parameterVec(6)*parameterVec(6)))*summation;
end

function [AIC] = FitModel2(Response,xaxis,yaxis,centerVals)
%FitModel2

% hyperParameterFun = @(b,distX,distY) (b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4));

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 7;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
numRepeats = 1000;
maxITER = 500;
tolerance = 1e-3;
h = ones(numParameters,1)./100;

% parameters are 
%  1) b(1) - rise of map at center
%  2)+3) the position (x and y) of the center of mass
%  4) b(2) - standard deviation or spread of the map in x
%  5) b(3) - standard deviation or spread of map in y
%  6) rise for sigma, from the Gaussian likelihood, at map center
%  7) b(4) - peak negativity at edges of retinotopic region
%    b(1)+b(4) = peak negativity at retinotopic center of mass
%  8) rho - allows for elliptical contours
%      9) sigma at edges of retinotopic region
Bounds = [-1000,0;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,2000;1,2000;1,1000;-1000,0];

% display('Steepest Ascent ...');
for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = squeeze(Response(zz,:,2));
    
    bigParameterVec = zeros(numRepeats,numParameters);
    bigLogLikelihoods = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions
    parfor repeats = 1:numRepeats
        gradientVec = zeros(numParameters,1);
        parameterVec = zeros(maxITER,numParameters);
        logLikelihood = zeros(maxITER,1);
        %parameterVec(1,:) = squeeze(bigParameterVec(repeats,:));
        for ii=1:numParameters
            parameterVec(1,ii) = Bounds(ii,1)+(Bounds(ii,2)-Bounds(ii,1)).*rand;
        end
        
        % for each starting position, do maxITER iterations
        for iter=2:maxITER
            % calculate likelihood at the current position in parameter
            %  space
            [logLikelihood(iter-1)] = GetLikelihood2(reps,parameterVec(iter-1,:),peakNegativity,flashPoints);
            
            if iter > 5
                check = sum(diff(logLikelihood(iter-2:iter-1)));
                if check < tolerance
                    break;
                end
            end
%             
%             % check if a random jump along one dimension improves the
%             %  likelihood
            temp = parameterVec(iter-1,:);
            randInd = random('Discrete Uniform',numParameters,1);
            temp(randInd) = Bounds(randInd,1)+(Bounds(randInd,2)-Bounds(randInd,1)).*rand;
            
            [likely] = GetLikelihood2(reps,temp,peakNegativity,flashPoints);
            if likely > logLikelihood(iter-1) %|| mod(iter,200) == 0
                parameterVec(iter-1,:) = temp';
                logLikelihood(iter-1) = likely;
            end
            
            % calculate the gradient by calculating the likelihood
            %  after moving over a small step h along each
            %  dimension
            for jj=1:numParameters
                tempParameterVec = parameterVec(iter-1,:);
                tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
                [gradLikelihoodplus] = GetLikelihood2(reps,tempParameterVec,peakNegativity,flashPoints);
                
                tempParameterVec = parameterVec(iter-1,:);
                tempParameterVec(jj) = tempParameterVec(jj)-h(jj);
                [gradLikelihoodminus] = GetLikelihood2(reps,tempParameterVec,peakNegativity,flashPoints);
                gradientVec(jj) = (gradLikelihoodplus-gradLikelihoodminus)./(2*h(jj));
            end
            
            % line search to get distance to move along gradient
            alpha = [0,1e-8,1e-6,1e-4,1e-2,1e0,1e1];
            lineSearchLikelihoods = zeros(length(alpha),1);
            lineSearchLikelihoods(1) = logLikelihood(iter-1);
            
            for ii=2:length(alpha)
                tempParameterVec = parameterVec(iter-1,:)'+gradientVec.*alpha(ii);
                [lineSearchLikelihoods(ii)] = GetLikelihood2(reps,tempParameterVec,peakNegativity,flashPoints);
            end
            
            newInds = find(imag(lineSearchLikelihoods)==0);
            newValues = zeros(length(newInds),1);
            for ii=1:length(newInds)
                newValues(ii) = lineSearchLikelihoods(newInds(ii));
            end
            [~,ind] = max(newValues);
            ind = newInds(ind);
            
            for jj=1:numParameters
                parameterVec(iter,jj) = max(Bounds(jj,1),min(parameterVec(iter-1,jj)+alpha(ind)*gradientVec(jj),Bounds(jj,2)));
            end
        end
        bigParameterVec(repeats,:) = parameterVec(iter-1,:);
        bigLogLikelihoods(repeats) = logLikelihood(iter-1);
    end
    [~,index] = max(bigLogLikelihoods);
    finalParameters(zz,:) = bigParameterVec(index,:);
    AIC = 2*numParameters-2*bigLogLikelihoods(index);
%     display(finalParameters(zz,:));
end

end

function [loglikelihood] = GetLikelihood2(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    distX = flashPoints(kk,1)-parameterVec(2);
    distY = flashPoints(kk,2)-parameterVec(3);
    b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(7)];
    mu = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4);
    stdev = parameterVec(6);%*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+parameterVec(9);
    loglikelihood = loglikelihood-(1/2)*log(2*pi*stdev*stdev)-(1/(2*stdev*stdev))*(peakNegativity(kk)-mu).^2;
%     summation = summation+(peakNegativity(kk)-mu).^2;
end

% loglikelihood = (-reps/2)*log(2*pi*parameterVec(6)*parameterVec(6))-...
%     (1/(2*parameterVec(6)*parameterVec(6)))*summation;
end

function [AIC] = FitModel3(Response,xaxis,yaxis,centerVals)
%FitModel2

% hyperParameterFun = @(b,distX,distY) (b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4));

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 8;
   
finalParameters = zeros(numChans,numParameters);
fisherInfo = zeros(numChans,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,numParameters);
numRepeats = 1000;
maxITER = 500;
tolerance = 1e-3;
h = ones(numParameters,1)./100;

% parameters are 
%  1) b(1) - rise of map at center
%  2)+3) the position (x and y) of the center of mass
%  4) b(2) - standard deviation or spread of the map in x
%  5) b(3) - standard deviation or spread of map in y
%  6) rise for sigma, from the Gaussian likelihood, at map center
%  7) b(4) - peak negativity at edges of retinotopic region
%    b(1)+b(4) = peak negativity at retinotopic center of mass
%  8) rho - allows for elliptical contours
%      9) sigma at edges of retinotopic region
Bounds = [-1000,0;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,2000;1,2000;1,1000;-1000,0;-11,11];

% display('Steepest Ascent ...');
for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    peakNegativity = squeeze(Response(zz,:,2));
    
    bigParameterVec = zeros(numRepeats,numParameters);
    bigLogLikelihoods = zeros(numRepeats,1);
    % repeat gradient ascent from a number of different starting
    %  positions
    parfor repeats = 1:numRepeats
        gradientVec = zeros(numParameters,1);
        parameterVec = zeros(maxITER,numParameters);
        logLikelihood = zeros(maxITER,1);
        %parameterVec(1,:) = squeeze(bigParameterVec(repeats,:));
        for ii=1:numParameters
            parameterVec(1,ii) = Bounds(ii,1)+(Bounds(ii,2)-Bounds(ii,1)).*rand;
        end
        
        % for each starting position, do maxITER iterations
        for iter=2:maxITER
            % calculate likelihood at the current position in parameter
            %  space
            [logLikelihood(iter-1)] = GetLikelihood3(reps,parameterVec(iter-1,:),peakNegativity,flashPoints);
            
            if iter > 5
                check = sum(diff(logLikelihood(iter-2:iter-1)));
                if check < tolerance
                    break;
                end
            end
%             
%             % check if a random jump along one dimension improves the
%             %  likelihood
            temp = parameterVec(iter-1,:);
            randInd = random('Discrete Uniform',numParameters,1);
            temp(randInd) = Bounds(randInd,1)+(Bounds(randInd,2)-Bounds(randInd,1)).*rand;
            
            [likely] = GetLikelihood3(reps,temp,peakNegativity,flashPoints);
            if likely > logLikelihood(iter-1) %|| mod(iter,200) == 0
                parameterVec(iter-1,:) = temp';
                logLikelihood(iter-1) = likely;
            end
            
            % calculate the gradient by calculating the likelihood
            %  after moving over a small step h along each
            %  dimension
            for jj=1:numParameters
                tempParameterVec = parameterVec(iter-1,:);
                tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
                [gradLikelihoodplus] = GetLikelihood3(reps,tempParameterVec,peakNegativity,flashPoints);
                
                tempParameterVec = parameterVec(iter-1,:);
                tempParameterVec(jj) = tempParameterVec(jj)-h(jj);
                [gradLikelihoodminus] = GetLikelihood3(reps,tempParameterVec,peakNegativity,flashPoints);
                gradientVec(jj) = (gradLikelihoodplus-gradLikelihoodminus)./(2*h(jj));
            end
            
            % line search to get distance to move along gradient
            alpha = [0,1e-8,1e-6,1e-4,1e-2,1e0,1e1];
            lineSearchLikelihoods = zeros(length(alpha),1);
            lineSearchLikelihoods(1) = logLikelihood(iter-1);
            
            for ii=2:length(alpha)
                tempParameterVec = parameterVec(iter-1,:)'+gradientVec.*alpha(ii);
                [lineSearchLikelihoods(ii)] = GetLikelihood3(reps,tempParameterVec,peakNegativity,flashPoints);
            end
            
            newInds = find(imag(lineSearchLikelihoods)==0);
            newValues = zeros(length(newInds),1);
            for ii=1:length(newInds)
                newValues(ii) = lineSearchLikelihoods(newInds(ii));
            end
            [~,ind] = max(newValues);
            ind = newInds(ind);
            
            for jj=1:numParameters
                parameterVec(iter,jj) = max(Bounds(jj,1),min(parameterVec(iter-1,jj)+alpha(ind)*gradientVec(jj),Bounds(jj,2)));
            end
        end
        bigParameterVec(repeats,:) = parameterVec(iter-1,:);
        bigLogLikelihoods(repeats) = logLikelihood(iter-1);
    end
    [~,index] = max(bigLogLikelihoods);
    finalParameters(zz,:) = bigParameterVec(index,:);
    AIC = 2*numParameters-2*bigLogLikelihoods(index);
%     display(finalParameters(zz,:));
end

end

function [loglikelihood] = GetLikelihood3(reps,parameterVec,peakNegativity,flashPoints)
loglikelihood = 0;
for kk=1:reps
    distX = flashPoints(kk,1)-parameterVec(2);
    distY = flashPoints(kk,2)-parameterVec(3);
    rho = parameterVec(8);rho = 2./(1+exp(-rho))-1;
    b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(7),rho];
    mu = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+b(4);
    stdev = parameterVec(6);%*exp(-(distX.^2)./(2*b(2)*b(2))-b(5)*distX*distY/(2*b(2)*b(3))-(distY.^2)./(2*b(3)*b(3)))+parameterVec(9);
    loglikelihood = loglikelihood-(1/2)*log(2*pi*stdev*stdev)-(1/(2*stdev*stdev))*(peakNegativity(kk)-mu).^2;
%     summation = summation+(peakNegativity(kk)-mu).^2;
end

% loglikelihood = (-reps/2)*log(2*pi*parameterVec(6)*parameterVec(6))-...
%     (1/(2*parameterVec(6)*parameterVec(6)))*summation;
end



