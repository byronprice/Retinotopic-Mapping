% RealTimeRetinoBayes_RecordComp.m

% To be used with a Plexon recording system and in conjunction with a 
%  stimulus computer that generates images on a screen. The two computers
%  will communicate via a tcp/ip server and also via a DAQ that sends
%  event times from the stimulus to the recording computer. On this side,
%  the data will pulled in real time from the acquisition software
%  (PlexControl). The VEPs recorded in response to individual stimulus
%  presentations will be analyzed in the following way: if the trace on an
%  individual trial crossing some threshold within a window from 60 to
%  120 milliseconds, then that presentation will get a value of 1,
%  otherwise 0. 

s = PL_InitClient(0);
if s == 0
   return
end

% setup server, send signal to stim computer that this computer is ready to
%  begin
display('Opening TCP/IP Server ...');
tcpipServer = tcpip('128.197.59.166',30000,'NetworkRole','server');
bufferSize = 50000;
set(tcpipServer,'OutputBufferSize',bufferSize);
fopen(tcpipServer);

p = PL_GetPars(s);
%sampleFreq = p(8);
sampleFreq = p(13);
stimLen = round(0.2*sampleFreq);
window = [round(0.06*sampleFreq),round(0.12*sampleFreq)];

startEXP = 254;
endEXP = 255;

startRUN = 252;
endRUN = 253;

endCHAN = 251;

tsChans = [6,8];  
adChans = [6,8]; % channels 22 and 24 have something on them, maybe continuous spiking activity
numChans = length(adChans);

% collect some data to get a baseline of the noise
%  wait for the signal from the stim computer that the experiment is ready
%  to begin

totalHeld = 30*sampleFreq;
D = zeros(totalHeld,numChans);
check = 0;
display('Obtaining estimate of noise ...');
index = 1;
tEvs = zeros(1,4);
while check == 0
    pause(1/1000);
    [~,tEvs] = PL_GetTS(s);
    [n,~,d] = PL_GetADV(s);

    D(index:index+n-1,:) = d(:,adChans);
    index = index+n;
    check = sum(tEvs(:,3) == startEXP);
    
    if index >= totalHeld-n
        index = 1;
    end
end

lastNonzero = find(D(:,1),1,'last');
D = D(1:lastNonzero,:).*1e6;
% approximate sigma = k*MAD 
k = 1.4826;
threshold = -(2*k).*mad(D,1,1);
display(sprintf('Noise threshold: %3.2f',threshold(1)));

%figure();plot(D(:,1)+threshold(1));hold on;
%plot(-threshold(1)*ones(lastNonzero,1));
% get probability of threshold crossings during noise events

N = 2000;
indeces = random('Discrete Uniform',lastNonzero-2*stimLen,[N,1]);

Pr = zeros(numChans,1);
for ii=1:numChans
    for jj=1:N
        VEP = D(indeces(jj):indeces(jj)+stimLen-1,ii);
        %saveVEP = VEP;
        %minInds = find(VEP<0);
       % maxInds = find(VEP>=0);
        %VEP(minInds) = -1;
        %VEP(maxInds) = 1;
        %tot = diff(VEP);
        [minVal,minInd] = min(VEP);
        if minInd > window(1) && minInd < window(2) && minVal < threshold(ii)
            Pr(ii) = Pr(ii)+1;
        end
    end
    if Pr(ii) == 0
        Pr(ii) = 1/N;
    end
end
Pr = Pr./N;

fwrite(tcpipServer,Pr,'double');

display(sprintf('Probability of VEP-like event: Channel 1- %3.2f, Channel 2- %3.2f',Pr(1),Pr(2)));
% this is the actual experiment
%  perform simple data analysis on the VEPs to determine whether or not a
%  significant VEP has occurred, relay the information to the stim computer

D = zeros(10*sampleFreq,2);
T = zeros(10*sampleFreq,4);
display('Beginning mapping experiment ...');
for ii=1:numChans
    display(sprintf('Mapping channel %d ...',ii));
   
    chanCheck = 0;
    tEvs = zeros(1,4);
    while chanCheck == 0
       [n, t, d] = PL_GetADV(s);
       [~, tEvs] = PL_GetTS(s);
       % tEvs contains the event timeStamps ... if tEvs(x,1) = 4 and tEvs(x,2)
       % = 257 , then tEvs(x,3) = strobed event word and tEvs(x,4) = timeStamp
       % this can be compared against t from GetADV and t(x,4) for all of the
       % events in which tEvs(x,2) = 6 (your desired channel)
       
       % continue while loop
       %  until startRUN event is strobed, perform analysis
       if sum(tEvs(:,3) == startRUN) > 0 
           runCheck = 0;
           dInd = 1;
           tInd = 1;
           tEvs = zeros(1,4);
           nEvs = 0;
           while runCheck == 0
              pause(1/1000);
              [n, t, d] = PL_GetADV(s);
              [nEvs,tEvs] = PL_GetTS(s);
              D(dInd:dInd+n-1,1) = d(:,adChans(ii));
              D(dInd,2) = t;
              dInd = dInd+n;
              
              if nEvs > 0 
                T(tInd:tInd+nEvs-1,:) = tEvs;
                tInd = tInd+nEvs;
              end
              runCheck = sum(tEvs(:,3) == endRUN);
           end
           t = D(find(D(:,2),1,'first'),2);
           timeStamps = t:1/sampleFreq:t+length(D)/sampleFreq-1;
           D(:,1) = D(:,1).*1e6;
           
           strobedEventNums = T(T(:,2)==257,3);
           strobedEventTimes = T(T(:,2) == 257,4);
           
           stimNum = strobedEventNums(strobedEventNums<250);
           stimOnset = strobedEventTimes(strobedEventNums<250); % stim onset time
           
           data = zeros(2,1);
           data(2) = stimNum; % event number
           [~,index] = min(abs(timeStamps-stimOnset));

           VEP = D(index:index+stimLen-1,1);

           [minVal,minInd] = min(VEP);
           if minInd > window(1) && minInd < window(2) && minVal < threshold(ii)
               data(1) = 1;
           end
           % send data over tcp/ip server
           fwrite(tcpipServer,data(:),'double');
           T(:) = 0;
           D(:) = 0;
           tEvs(:) = 0;
       end
       chanCheck = sum(tEvs(:,3) == endCHAN);
    end
end

fclose(tcpipServer);
display('Experiment over.');
% you need to call PL_Close(s) to close the connection
% with the Plexon server
PL_Close(s);
s = 0;

