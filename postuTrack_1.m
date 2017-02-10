% -----------------------------------------------------------------
%   Coded: Rafa Sebastian, University of Valencia 
%   Date: 06/03/2015
% -----------------------------------------------------------------
%   Code to post-process data detected by the uTrack software.
%   INPUT: .mat files produced by uTrack
%   OUTPUT: Intensity profiles of the detected events
% ------------------------------------------------------------------

% Set-up the detections params

clear all;
close all;

debug = 1;
% Detection parameters
ProfRadius = 2;       % Radius to calculate intensity profile
flagIntesity = 1;     % Sets =1(maximum gay level); =0 (averaged gray level) 
MaxLatDispl = 5;     % Maximum lateral displacement of a dot allowed
TotalDotsDetect = 0;  % Final number of dots (static and long enough) detected
savePlots = 0;        % Save plots after processing (=0 NO, =1 YES)

% CRITICAL USER PARAMETERS ------------------------------------------------
Ithres = 0.40;         % Events with Max Intensity < Ithres (%) of maximum intensity in the whole sequence, will be removed
MinDuration = 10;      % Minimum number of frames that lasts an event
MaxMissFrAllowed = 2; % Maximum number of frames where the spot "disappears"
CorrThres = 0.45;     % Minimum correlation between channels required
FitThres = 0.60;
% -------------------------------------------------------------------------


nEx = '1';             % Experiment number (for storing image seq.)
nC = '1';              % Cell number (for storing image seq.)
nCh = '1';             % Cell number (for storing image seq.)

filename = 'ReportDetection.txt';
% -------------------------------------------------------------------------
% Asking the number of channels
nChannels=input('Indicate the number of channels to be analysed: ');

% Reading channel with A3
[utrack_name1,utrack_folder] = uigetfile('*.mat','Select the .mat uTrack A3 file');
[image_name1,image_folder] = uigetfile({'*.tif';'*.tiff'},'Select the .tif A3 stack file');
inputfile_ch1 = [utrack_folder utrack_name1];

% Load utrack output data
fprintf('Reading .mat A3 file from uTrack...');
uDataCh1 = load(inputfile_ch1);
fprintf('[DONE]\n');
% Reading Channel with Dynamin
if (nChannels > 1) % If available:
    [utrack_name2,utrack_folder] = uigetfile('*.mat','Select the .mat uTrack Dynamin file');
    [image_name2,image_folder] = uigetfile({'*.tif';'*.tiff'},'Select the .tif Dynamin stack file');
    inputfile_ch2 = [utrack_folder utrack_name2];
    fprintf('Reading .mat Dynamin file from uTrack...');
    uDataCh2 = load(inputfile_ch2);
    fprintf('[DONE]\n');
end;
    
% Extracting A3 uTrack data
nFrames = length(uDataCh1.movieInfo);   % Number of frames in channel
nSpots = 0;
Framei = uDataCh1.movieInfo(1); 
for i=1:nFrames
   nSpots = length(Framei.xCoord) + nSpots;     
end
fprintf('\n Classifying %0.d spots from uTrack... \n ',nSpots);
fprintf(' Average # spots per frame = %2.d \n ',nSpots/nFrames);

% -------------------------------------------------------------------------
% REPORTING INFO READ
File_props = imfinfo([image_folder image_name1]);
fprintf('\n\n------------------------------------------------------------\n');
fprintf('File Information (A3 Channel):\n');
fprintf('   Stack length: %.0d \n',length(File_props));
fprintf('   Dimensions (px):%.0f x %.0f \n', File_props(1).Height,File_props(1).Width);
fprintf('   Resolution (px):%.3f x %.3f \n', File_props(1).XResolution,File_props(1).YResolution);
fprintf('   Bit Depth :%.0d bits \n', File_props(1).BitDepth);
% -------------------------------------------------------------------------


% Load original A3 Image stack
fprintf('Loading TIFF File [%s]\n',[image_folder image_name1]);
fprintf('  Removing Background...');
se = strel('disk',12);
for  i=1:length(File_props) 
    S(i).ima = imtophat(imread([image_folder image_name1],'Index',i),se);
    Io(i).ima = imread([image_folder image_name1],'Index',i);
end;
fprintf('[DONE]\n');

% If available: Load origial DYNAMIN Image stack
if (nChannels > 1) 
    fprintf('Loading TIFF File [%s]\n',[image_folder image_name2]);
    fprintf('  Removing Background...');
    for  i=1:length(File_props) 
        SD(i).ima = imtophat(imread([image_folder image_name2],'Index',i),se);
        IoD(i).ima = imread([image_folder image_name2],'Index',i);
    end;
    fprintf('[DONE]\n');
end;

[columnsInImage rowsInImage] = meshgrid(1:File_props(1).Width, 1:File_props(1).Height);

% -------------------------------------------------------------------------
% STEP 1:
% Add all tentative SPOTS seen in Frame 1:

ActiveEvents = zeros(10000,1); % Maximum spots allowed 10000
Framei = uDataCh1.movieInfo(1); 
nInitialSpots = length(Framei.xCoord);
s = nInitialSpots;
for (i = 1:s)
    Ispots(i).x = Framei.xCoord(i,1);
    Ispots(i).y = Framei.yCoord(i,1);
    Ispots(i).fr = 1;
    Ispots(i).Thres = 0;  % No Threshold passed yet
    % Calculate Mean intensity profile inside a circle of radius ProfRadius
    CirclePixels= (rowsInImage - round(Framei.yCoord(i,1))).^2 + (columnsInImage - round(Framei.xCoord(i,1))).^2 <= ProfRadius.^2;
    Ispots(i).I = mean(Io(1).ima(CirclePixels)); 
    Ispots(i).I2 = mean(S(1).ima(CirclePixels));
    ActiveEvents(i) = 1; % Activating spots detected in fr  1
    
    % Update info for Dymanin
    if (nChannels > 1)
        Ispots(i).DI = mean(IoD(1).ima(CirclePixels)); 
        Ispots(i).DI2 = mean(SD(1).ima(CirclePixels));
    end
end;

% -------------------------------------------------------------------------
% STEP 2: Checking for spatio-temporal coherence of spots (no lateral movement)
fprintf('\nChecking spatial coherence: \n ');

h = waitbar(0,'Checking spatial coherence'); 
for i = 2:nFrames
   waitbar(i/nFrames,h);
   Framei = uDataCh1.movieInfo(i); 
   
   % SPOTs detected in FRAME i:
   %  1 -> To be assigned to an EXISTING EVENTS
   %  2 -> To be created as NEW STARTING EVENTS
   nTentativeSpots = length(Framei.xCoord);
   
   nPrevSpots = length(Ispots); % Previous EVENTS detected so far (not closed yet)
   Updated = zeros(nPrevSpots,1); % For each new Frame_i "Reset Updated Spots"
   
   for j = 1:nTentativeSpots
       % Retrieve Coordinates [x y] of spot j at time = i
       x_j = Framei.xCoord(j,1);
       y_j = Framei.yCoord(j,1);
       
       % Check if [x_j y_j] matches a previous in Ispots()
       nPrevSpots = sum(ActiveEvents);
       Ind = find(ActiveEvents); % Indices of Events to check
       MinDist = MaxLatDispl + 1;
       MinInd = 0;
       
       %  1 -> Check distance to EXISTING EVENTS
       for k = 1:nPrevSpots     
           % if previous dot is close enough in time to this tentativespot 
           if ((i-Ispots(Ind(k)).fr(end)) <= MaxMissFrAllowed ) 
             % Retrieve Coordinates [x y] of spot k at time = i-1
              x_k = Ispots(Ind(k)).x(end);
              y_k = Ispots(Ind(k)).y(end);
              dist = sqrt((x_k-x_j)^2+(y_k-y_j)^2);
              if dist < MinDist
                MinInd = Ind(k); % Min distance from all previous Ispots to tentative [x_j y_j]
                MinDist = dist;
              end;
           else
               ActiveEvents(Ind(k)) = 0; % This event had no spots for > MaxMissFrAllowed
           end;
       end;
       
       if (MinInd > 0) % If the spot is close to a prev. one in space and time:
           
           if (Updated(MinInd) == 0) % if Event was not updated already in frame_i: 
               Ispots(MinInd).x(end+1) = x_j;
               Ispots(MinInd).y(end+1) = y_j;
               Ispots(MinInd).fr(end+1) = i;
               
               % Calculate Mean intensity profile inside a circle of radius ProfRadius
               CirclePixels = (rowsInImage - round(y_j)).^2 + (columnsInImage - round(x_j)).^2 <= ProfRadius.^2;
               Ispots(MinInd).I(end+1) = mean(Io(i).ima(CirclePixels));
               Ispots(MinInd).I2(end+1) = mean(S(i).ima(CirclePixels)); 
               
               Ispots(MinInd).Thres = 0;  % No Threshold passed yet
               Updated(MinInd) = j; % Spot MARKED as updated with j-th
               
               % Update info for Dymanin
               if (nChannels > 1)
                   Ispots(MinInd).DI(end+1) = mean(IoD(i).ima(CirclePixels)); 
                   Ispots(MinInd).DI2(end+1) = mean(SD(i).ima(CirclePixels));
               end
               
           else
               % If it was already updated already there is a problem -> 2 tentative
               % to same stored spot, THEN:
               % Check if [x_k y_k] is closer than previous STORED [x_k2 y_k2]
               k2 = Updated(MinInd);
               x_k2 = Framei.xCoord(k2,1);
               y_k2 = Framei.yCoord(k2,1);
               MinDist2 = sqrt((x_k2-x_j)^2+(y_k2-y_j)^2);
               if (MinDist < MinDist2) % We found a better match > Overwrite previous
                   Ispots(MinInd).x(end) = x_j;
                   Ispots(MinInd).y(end) = y_j;
                   Ispots(MinInd).fr(end) = i;
                   % Calculate Mean intensity profile inside a circle of radius ProfRadius
                   CirclePixels= (rowsInImage - round(y_j)).^2 + (columnsInImage - round(x_j)).^2 <= ProfRadius.^2;
                   Ispots(MinInd).I(end) = mean(Io(i).ima(CirclePixels)); 
                   Ispots(MinInd).I2(end) = mean(S(i).ima(CirclePixels));
                   Ispots(MinInd).Thres = 0;  % No Threshold passed yet
                   Updated(MinInd) = j; % Spot MARKED as updated with j-th
                   
                    % Update info for Dymanin
                   if (nChannels > 1)
                       Ispots(MinInd).DI(end) = mean(IoD(i).ima(CirclePixels)); 
                       Ispots(MinInd).DI2(end) = mean(SD(i).ima(CirclePixels));
                   end; %if_nchannels
               end; %if_MinDist
           end; %if_Updated
       else
           % 2 -> To be created as NEW STARTING EVENTS

           s = s + 1; % Increase number of EVENTS
           Updated(end+1) = 0; % Create a new index
           Ispots(s).x = x_j;
           Ispots(s).y = y_j;
           Ispots(s).fr = i;
           % Calculate Mean intensity profile inside a circle of radius ProfRadius
           CirclePixels= (rowsInImage - round(y_j)).^2 + (columnsInImage - round(x_j)).^2 <= ProfRadius.^2;
           Ispots(s).I = mean(Io(i).ima(CirclePixels));
           Ispots(s).I2 = mean(S(i).ima(CirclePixels));
           Ispots(s).Thres = 0;  % No Threshold passed yet
           ActiveEvents(s) = 1;
           
           % Update info for Dymanin
           if (nChannels > 1)
               Ispots(s).DI = mean(IoD(i).ima(CirclePixels)); 
               Ispots(s).DI2 = mean(SD(i).ima(CirclePixels));
           end           
           
       end;
           
   end;
   if (debug)
    fprintf('(%g/%g) Candidate spots = %g >  Active Events = %g \n',i,nFrames, nTentativeSpots,sum(ActiveEvents));
    fprintf('  Total potential events so far %g\n',s);
   end;
end
close(h);

% ----------------------------------------
% Checking Duration Thresholds:
fprintf('Checking minumum event duration: \n');
hh = waitbar(0,'Checking minumum event duration...'); 
k = 0; % Events Kept
d = 0; % Events deleted
for (i = 1:s)
    waitbar(i/s,hh);
    spot_i = Ispots(i);
    TimeLength = spot_i.fr(end)-spot_i.fr(1);
    
    % THRESHOLD 1: minimum temporal duration
    if (TimeLength >= MinDuration)
       % Detected in all Frames: keep
       %Save spot
       k = k + 1; % Events Kept
       Events(k) = spot_i;
       AllDuration(k) = TimeLength;
       MaxIntensity(k) = max(spot_i.I); % Without background
    else
       d = d + 1; % Events deleted
    end
end;


fprintf('    Total initial events: %.0d \n',s);
fprintf('    Events DELETED (duration < %.0d): %.0d \n',MinDuration,d);
close(hh);


% ----------------------------------------
% Checking ALL THRESHOLDS
fprintf('Checking ALL thresholds and assigning event type: \n');
hh = waitbar(0,'Checking ALL Thresholds');
MxMxIntensity = max(MaxIntensity);
d = 0;
d1 = 0;
d2 = 0;
s2 = length(Events); % Number of spots that PASSED (space+time) = s2
for (i = 1:s2)
    waitbar(i/s2,hh);
    
    % Check all thresholds in order
    if (MaxIntensity(i) >= Ithres*MxMxIntensity)
       Events(i).Thres = 1; % Pass Intensity Threshold
       
       % Check cross-correlation between channels:
       aux = corrcoef(Events(i).I2,Events(i).DI2);
       Mcor = aux(1,2); 
       if (Mcor >= CorrThres)
            Events(i).Thres = 2; % Pass correlation Threshold
           
           % Check Gaussian GOF Fitting of A3 channel:
            [curvefitA3,gofA3,outputA3] = fit(Events(i).fr',Events(i).I2','gauss2');
            if (gofA3.rsquare >= FitThres)
                Events(i).Thres = 3; % Pass Gaussian GOF Threshold
            else
                d2 = d2 + 1;
            end 
       else
           d1 = d1 + 1;
       end
    else
        d = d + 1;
    end
end;
close(hh);


fprintf('Total initial events: %g \n',s2);
fprintf('    Events with (GreyLevel < %.2d) = %g \n',Ithres*MxMxIntensity,d);
fprintf('    Events with (Correlation < %.2d) = %g \n',CorrThres,d1);
fprintf('    Events with (Gaussian Fitting < %.2d) = %g \n',FitThres,d2);
fprintf('    Events ACCEPTED : %g - (%g + %g + %g) = %g \n',s2,d,d1,d2,s2-d-d1-d2);

fid = fopen(filename,'w');
fprintf(fid,'File Information (A3 Channel):\n');
fprintf(fid,'   Stack length: %.0d \n',length(File_props));
fprintf(fid,'   Dimensions (px):%.0f x %.0f \n', File_props(1).Height,File_props(1).Width);
fprintf(fid,'   Resolution (px):%.3f x %.3f \n', File_props(1).XResolution,File_props(1).YResolution);
fprintf(fid,'   Bit Depth :%.0d bits \n', File_props(1).BitDepth);
fprintf(fid,'\n\n------------------------------------------------------------\n');
fprintf(fid,'Total initial events: %g \n',s2);
fprintf(fid,'    Events with (GreyLevel < %.2d) = %g \n',Ithres*MxMxIntensity,d);
fprintf(fid,'    Events with (Correlation < %.2d) = %g \n',CorrThres,d1);
fprintf(fid,'    Events with (Gaussian Fitting < %.2d) = %g \n',FitThres,d2);
fprintf(fid,'    Events ACCEPTED : %g - (%g + %g + %g) = %g \n',s2,d,d1,d2,s2-d-d1-d2);
fclose(fid);

Data = [];
Data = Events;
save(['EventsPASS_Ith(' Ithres ')MinD(' MinDuration ')_' image_name1(1:end-4) '.mat'],'Data');

