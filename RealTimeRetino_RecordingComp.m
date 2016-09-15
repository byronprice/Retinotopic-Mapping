% RealTimeRetino_RecordingComp.m

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
window = [round(0.06*sampleFreq),round(0.12*sampleFreq)];

startEXP = 254;
endEXP = 255;

startRUN = 252;
endRUN = 253;

tsChans = [6,8];  
adChans = [6,8]; % channels 22 and 24 have something on them, maybe continuous spiking activity
numChans = length(adChans);

% collect some data to get a baseline of the noise
%  wait for the signal from the stim computer that the experiment is ready
%  to begin
threshold = zeros(numChans,1);

totalHeld = 30*sampleFreq;
D = zeros(totalHeld,numChans);
check = 0;
display('Obtaining estimate of noise ...');
index = 1;
while check == 0
    [~,tEvs] = PL_GetTS(s);
    [n,~,d] = PL_GetADV(s);
    
    D(index:index+n-1,:) = d(:,adChans);
    check = sum(tEvs(:,3) == startEXP);
    
    index = index+n;
    if index >= totalHeld-n
        index = 1;
    end
    pause(1);
end

lastNonzero = find(D(:,1),1,'last');
D = D(1:lastNonzero,:);
% approximate sigma = k*MAD 
k = 1.4826;
threshold = (k).*mad(D,1,1);

% get probability of threshold crossings during noise events
stimLen = window(2);
N = 2000;
indeces = random('Discrete Uniform',lastNonzero-2*stimLen,[N,1]);

Pr = zeros(numChans,1);
for ii=1:numChans
    for jj=1:N
        VEP = D(indeces(jj)+window(1):indeces(jj)+window(2)-1,ii);
        maxs = max(VEP,0);maxs = min(maxs,1);
        mins = min(VEP,0);mins = max(mins,-1);
        tot = maxs+mins;
        if sum(diff(tot)<=-2) >= 1
            Pr(ii) = Pr(ii)+1;
        end
    end
end
Pr = Pr./N;

% this is the actual experiment
%  perform simple data analysis on the VEPs to determine whether or not a
%  significant VEP has occurred, relay the information to the stim computer



display('Beginning mapping experiment ...');
for ii=1:numChans
    display(sprintf('Mapping channel %d ...',ii));
    check = 0;
    while check == 0
       pause(4);
       [~, tEvs] = PL_GetTS(s);
       % tEvs contains the event timeStamps ... if tEvs(x,1) = 4 and tEvs(x,2)
       % = 257 , then tEvs(x,3) = strobed event word and tEvs(x,4) = timeStamp
       % this can be compared against t from GetADV and t(x,4) for all of the
       % events in which tEvs(x,2) = 6 (your desired channel)
       
       [n, t, d] = PL_GetADV(s);
       
       % if no events were strobed, continue while loop
       % if events were strobed, perform analysis
       if sum(tEvs(:,2) == 257) > 0 
           timeStamps = t:1/sampleFreq:t+n/sampleFreq-1;
           d = d(:,adChans(ii))+threshold(ii);
           
           strobedEventInds = find(tEvs(:,2) == 257);
           numStrobes = length(strobedEventInds);
           
           display(numStrobes)
           data = zeros(numStrobes,2);
           stimOnset = 0;
           for jj=1:numStrobes
                   data(jj,2) = tEvs(strobedEventInds(jj),3); % event number
                   stimOnset = tEvs(strobedEventInds(jj),4); % strobe time / stimulus onset time
                   [~,index] = min(abs(timeStamps-stimOnset));
                   VEP = d(index+window(1):index+window(2)-1);
                   maxs = max(VEP,0);maxs = min(maxs,1);
                   mins = min(VEP,0);mins = max(mins,-1);
                   tot = maxs+mins;
                   if sum(diff(tot)<=-2) >= 1
                       data(jj,1) = 1;
                   end
           end
           % send data over tcp/ip server
           w = whos('data');
           dataSize = [size(data),w.bytes];
           fwrite(tcpipServer,dataSize,'double');
           fwrite(tcpipServer,data(:),'double');
       end
       check = sum(tEvs(:,3) == endRUN);
    end
end

fclose(tcpipServer);
display('Experiment over.');
% you need to call PL_Close(s) to close the connection
% with the Plexon server
PL_Close(s);
s = 0;
