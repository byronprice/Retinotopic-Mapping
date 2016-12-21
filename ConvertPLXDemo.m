% ConvertPLXDemo.m
%   Example Script to Get VEPs

EphysFileName = 'whatever'; % put your filename here, exclude .plx ending

% call .mat to .plx conversion, convert to volts, could filter after this
%  step
[ChanData,timeStamps,tsevs,svStrobed,sampleFreq] = ExtractSignal(EphysFileName);

% outputs are the LFP data for each channel, the timestamps for that LFP
% data, tsevs and svStrobed contain information regarding the strobed
% events and event names

[stimTriggeredResponse,VEPS,strobeTimes] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed,sampleFreq);

% stimTriggeredResponse contains the set of all LFP responses to the individual
%  stimulus presentations, as a cell array (500 milliseconds of data)
% VEPS is a matrix sized number_of_channels by number of stimuli by
%  stimulus length