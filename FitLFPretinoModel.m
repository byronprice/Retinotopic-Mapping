function [finalParameters] = FitLFPretinoModel(Response,xaxis,yaxis,PrHitNoise,centerVals)
%FitLFPretinoModel.m
%   Use data from LFP retinotopic mapping experiment to fit a non-linear
%    model of that retinotopy (using hits/misses and a binomial
%    distribution)

%Created: 2016/12/12, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2017/01/17
% By: Byron Price

% DistFun = @(stimCenter,centerVals) (ceil(sqrt((stimCenter(1)-centerVals(:,1)).^2+(stimCenter(2)-centerVals(:,2)).^2))+1);
% hyperParameterFun = @(b,distToCenterMass) (b(1)*exp(-(distToCenterMass.^2)./(2*b(2)*b(2)))+b(3));
% hyperParameterFun = @(b,distX,distY) (b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4));

numChans = size(Response,1);
reps = size(Response,2);

numParameters = 5;
   
finalParameters = zeros(numChans,numParameters);
numRepeats = 10;
maxITER = 1000;
tolerance = 0.01;
h = [0.01,1,1,1,1];

% parameters are 
%  1) b(1) - probability of hit at center of mass
%  2)+3) the position (x and y) of the center of mass
%  4) b(2) - standard deviation or spread of the map in x
%  5) b(3) - standard deviation in y
Bounds = [0,1-PrHitNoise;min(xaxis),max(xaxis);min(yaxis),max(yaxis);1,2000;1,2000];

% display('Steepest Ascent ...');
for zz=1:numChans
%     display(sprintf('Running Data for Channel %d...',zz));
    flashPoints = centerVals(squeeze(Response(zz,:,1)),:);
    hitMiss = squeeze(Response(zz,:,2));
    
    gradientVec = zeros(numParameters,1);
    bigParameterVec = zeros(numRepeats,numParameters);
    % repeat gradient ascent from a number of different starting
    % positions
    for repeats = 1:numRepeats
        parameterVec = zeros(maxITER,numParameters);
        logLikelihood = zeros(maxITER,1);
        %parameterVec(1,:) = squeeze(bigParameterVec(repeats,:));
        for ii=1:numParameters
            parameterVec(1,ii) = Bounds(ii,1)+(Bounds(ii,2)-Bounds(ii,1)).*rand;
        end
        
        % for each starting position, do maxITER iterations
        for iter=2:maxITER
            % calcualte likelihood at the current position in parameter
            %  space
            
            [temp] = GetLikelihood(reps,parameterVec(iter-1,:),hitMiss,flashPoints,PrHitNoise);
            logLikelihood(iter-1) = temp;
            
            if iter > 250
                check = sum(diff(logLikelihood(iter-200:iter-1)));
                if check < tolerance
                    break;
                end
            end
            
            temp = parameterVec(iter-1,:);
            randInd = random('Discrete Uniform',numParameters,1);
            temp(randInd) = Bounds(randInd,1)+(Bounds(randInd,2)-Bounds(randInd,1)).*rand;
            [likely] = GetLikelihood(reps,temp,hitMiss,flashPoints,PrHitNoise);
            if likely > logLikelihood(iter-1)
                parameterVec(iter-1,:) = temp';
                logLikelihood(iter-1) = likely;
            end
            
            
            % calculate the gradient by calculating the likelihood
            %  after moving over a small step h along each
            %  dimension
            for jj=1:numParameters
                tempParameterVec = parameterVec(iter-1,:);
                tempParameterVec(jj) = tempParameterVec(jj)+h(jj);
                [gradLikelihood] = GetLikelihood(reps,tempParameterVec,hitMiss,flashPoints,PrHitNoise);
                gradientVec(jj) = (gradLikelihood-logLikelihood(iter-1))./h(jj);
            end
            
            % line search to get distance to move along gradient
            alpha = [0,1e-6,1e-4,1e-2,1e-1,1e0,1e1];
            lineSearchLikelihoods = zeros(length(alpha),1);
            lineSearchLikelihoods(1) = logLikelihood(iter-1);
            
            for ii=2:length(alpha)
                tempParameterVec = parameterVec(iter-1,:)'+gradientVec.*alpha(ii);
                [temp] = GetLikelihood(reps,tempParameterVec,hitMiss,flashPoints,PrHitNoise);
                lineSearchLikelihoods(ii) = temp;
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
    end
    finalParameters(zz,:) = median(bigParameterVec,1);
end

fisherInfo = zeros(numParameters,numParameters);
end

function [loglikelihood] = GetLikelihood(reps,parameterVec,hitMiss,flashPoints,PrHitNoise)
loglikelihood = 0;
for kk=1:reps
    distX = flashPoints(kk,1)-parameterVec(2);
    distY = flashPoints(kk,2)-parameterVec(3);
    b = [parameterVec(1),parameterVec(4),parameterVec(5),PrHitNoise];
    p = (b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4));
    loglikelihood = loglikelihood+log((p.^(hitMiss(kk))).*((1-p).^(1-hitMiss(kk))));
end
end

% FisherInfo = zeros(3,3);
% % [ sigma_xx sigma_xy sigma_xstd        ]
% % [ sigma_yx sigma_yy sigma_ystd        ] 
% % [ sigma_stdx sigma_stdy sigma_stdstd  ]
% 
% FisherInfo(1,1) = (reducedLikely(2,2+1,2)-2*reducedLikely(2,2,2)+...
%     reducedLikely(2,2-1,2))./(h_dist^2);
% FisherInfo(2,2) = (reducedLikely(2+1,2,2)-2*reducedLikely(2,2,2)+...
%     reducedLikely(2-1,2,2))./(h_dist^2);
% FisherInfo(3,3) = (reducedLikely(2,2,2+1)-2*reducedLikely(2,2,2)+...
%     reducedLikely(2,2,2-1))./(h_std^2);
% FisherInfo(1,2) = (reducedLikely(2+1,2+1,2)-reducedLikely(2+1,2-1,2)-...
%     reducedLikely(2-1,2+1,2)+reducedLikely(2-1,2-1,2))./(4*h_dist^2);
% FisherInfo(2,1) = FisherInfo(1,2);
% FisherInfo(1,3) = (reducedLikely(2,2+1,2+1)-reducedLikely(2,2+1,2-1)-...
%     reducedLikely(2,2-1,2+1)+reducedLikely(2,2-1,2-1))./(4*h_dist*h_std);
% FisherInfo(3,1) = FisherInfo(1,3);
% FisherInfo(2,3) = (reducedLikely(2+1,2,2+1)-reducedLikely(2+1,2,2-1)-...
%     reducedLikely(2-1,2,2+1)+reducedLikely(2-1,2,2-1))./(4*h_dist*h_std);
% FisherInfo(3,2) = FisherInfo(2,3);
% 
% covMat = inv(-FisherInfo);
% 
% maxX = xaxis(J);maxY = yaxis(I);maxSTD = STDS(K);
% errorX = sqrt(covMat(1,1));errorY = sqrt(covMat(2,2));
% errorSTD = sqrt(covMat(3,3));

