# RetinotopicMapping
Code to display visual stimuli and analyze data to determine retinotopy of an LFP recording electrode. Use in conjunction with a Plexon electrophysiology system, Plexon's MATLAB Offline SDK, the Psychtoolbox, and Jeff Gavornik's StimulusSuite. The StimulusSuite is used to send event times through a DAQ to the Plexon system. In the future, I'd like to break free from that because it's too cumbersome to simply use for sending event times. Plexon outputs a .plx file.  In the MapRetinotopy.m file, I use the Offline SDK to convert a .plx to a .mat file.  This is done by the MyReadall.m function.  This function is nearly identical to Plexon's readall.m, but saves the .plx file to .mat with the same name (see below for file naming rules).
# 
The lab implants one microelectrode (3mm tip) into binocular primary visual cortex (one in each hemisphere for a total of two channels).  We then record the LFP using a Plexon system. Mice are positioned at 25 cm from a screen with a refresh rate of 80 Hz.  The specifics of your configuration can be adjusted in the Retinotopy.m file. The RetinotopyCallaway.m file and its supporting files, RetinotopyCallaway.frag.txt, RetinotopyCallaway.vert.txt, and MapRetinotopyCallaway.m are works in progress. We give our mice a unique identifier, e.g. 26881 corresponding to its cage ID and its personal ID within the cage.  So, we might have 26880, 26881, 26882, and 26883. 

#Steps:
1) Run Retinotopy(26881) to mouse #26881 while recording the LFP using Plexon system. Name the .plx file as RetinoDataDate_AnimalName.plx (e.g. July 5, 2016 is 20160705 so RetinoData20160705_26881.plx). The Retinotopy.m file will output a file with the stimulus parameters named RetinoStim20160705_26881.mat .
#
2) Run MapRetinotopy(26881,20160705) as long as the RetinoData file and the RetinoStim file are both on the MATLAB path.
#
3) This will output a figure and a file named RetinoMap20160705_26881.mat . 
#
4) Clear the workspace, open that file, view the figure, and select the channel with the best retinotopy. It should be a heat map with a clear center of greatest activity with a relatively circular decay around that center (like a 2D Gaussian). Write Channel = 1 if channel 1 is the best channel, or Channel = 2 if that is best (use the channel name at the top of the figure, ignoring the fact that some systems start at Channel 6 or what have you). Save all of the variables in the workspace as RetinoMap26881.mat . This is for use in subsequent analyses or stimulus generation, see Sequence-Learning repository for example.
#
5) If you would like to run the MapRetinotopy.m function for many animals, the best way is to use MapRetWrapper.m .  If you place all of the RetinoData* and RetinoStim files into a common folder, this can easily be done by simply calling MapRetWrapper('RetinoData*') .  If you want to only do so for a single animal, call MapRetWrapper('RetinoData*26881*') 

