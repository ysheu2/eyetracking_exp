% Clear the workspace and the screen
close all;
clear all;
sca


% Seed the random number generator.
RandStream.create('mrg32k3a','seed',sum(100*clock));

% Set scanner true or false
scanner = true;

% Set eyetracker dummymode true or false
dummymode = 0; % set to 1 to initialize in dummymode

% Set eyetracker true or false
eyelink = true;

%% Unify and Surpress
KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests',1);
Screen('Preference','VisualDebugLevel', 0);
Screen('Preference', 'SuppressAllWarnings', 1);
warning('off','MATLAB:dispatcher:InexactMatch')

%% Important directories
expDir = pwd;
outDir = [pwd filesep 'Output' filesep]; %must create a folder name "Output" in your dir

%% Timing information in seconds

% mock time       
% T_syncwaste = 0;
% T_fix = [0.01 0.01];
% T_cue = 0.01;
% T_sample = 0.01;
% T_delay = [0.02 0.01]; 
% T_probe = 0.03;
% T_interval = 0.01;
% T_feed  = 0.2;
% T_end   = 2;

T_syncwaste = 0;
T_fix = [1 3]; % =ITI
T_cue = 1.5;
T_sample = 3;
T_delay = [4 6]; 
T_probe = 3;
T_interval = 0.01;
T_feed = 0.3;
T_end = 2;

T_nopDelay = [3 5 7 9];

totaltime1 = ((2+T_cue+T_sample+5+T_probe)*40)/60;
totaltime2 = ((2+T_cue+T_sample+6+T_probe)*8)/60;
totaltime = totaltime1+totaltime2;

%% Keyboard/Input Information

STS = 100; %standard text size
STM = 2/3;%standard text size multiplier
STS_cross = 64; %standard text size for fixation cross
ScaleDown = .9;
ScaleUp   = 1.2;

if scanner
    leftkey = 2; %253; %console 2 (match)
    rightkey = 3; %254; %console 1 (non-match), check that and make sure the right hand: use superlab 1.2
    synckey = 7; %223; %mri sync button
else
    leftStr = 'n';
    rightStr = 'm';
    killStr = 'q';
    nextStr = 'space';
    leftkey = KbName(leftStr);
    rightkey = KbName(rightStr);
    killkey = KbName(killStr);
    nextkey = KbName(nextStr);
end


%% Get the Subject ID, make sure you don't overwrite anything

subjID = input('Enter Subject ID: ', 's');
runNumber = input('Enter Run Number: ', 's');
%subjID = 123;
%numRuns = str2num(input('Enter number of trials: ', 's'));
outputFname = [outDir 'NIMH'];  %need to be short and no symbols
outputFname = [outputFname '_' subjID '_' runNumber '.mat'];
if exist(outputFname, 'file')
    disp('Output file already exists!');
    cd(expDir)
    return
end

%%%%%%%%%%
% STEP 1 %
%%%%%%%%%%
% Added a dialog box to set your own EDF file name before opening
% experiment graphics. Make sure the entered EDF file name is 1 to 8
% characters in length and only numbers or letters are allowed.
% Note: Octave does not support GUIs. replace lines below with
% %edfFile= 'DEMO.EDF'

if eyelink
    prompt = {'Enter tracker EDF file name (1 to 8 letters or numbers)'};
    dlg_title = 'Create EDF file';
    num_lines= 1;
    def     = {'DEMO'};
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    edfFile = answer{1};
    fprintf('EDFFile: %s\n', edfFile );
end
%% Set up screen options

%%%%%%%%%%
% STEP 2 %
%%%%%%%%%%
% This gives us a number for each of the screens attached to our computer.
screens = Screen('Screens');

% To draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen.
screenNumber = max(screens);

% Basis Colors
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);
gray = (white+black)/2;

% Opens window, offscreen, color gray
[window, theRect] = Screen(screenNumber, 'OpenWindow', gray, [], [], 2); %other script use PschImaging
%[window, theRect] = PsychImaging('OpenWindow', screenNumber, gray);
Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Get the center coordinate of the window
[xCenter, yCenter] = RectCenter(theRect);  

% Query the frame duration
FlipInt=Screen('GetFlipInterval', window);

% Set the text font and size
Screen('TextSize', window, STS);
Screen('TextFont', window, 'Helvetica');

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

baseRect = [0 0 130 130]; % we want all images to show up 200 x 200 pixels

% Here we set the size of the arms of our fixation cross
fixCrossDimPix = 20;
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 4;

% Center the left hand side squares on positions in the screen.
UpperLeftRect = CenterRectOnPointd(baseRect, screenXpixels * 0.40, yCenter + 110);
UpperMiddleRect = CenterRectOnPointd(baseRect,screenXpixels * 0.50, yCenter + 110);
UpperRightRect = CenterRectOnPointd(baseRect,screenXpixels * 0.60 , yCenter + 110);
 
% Do the same of the right hand side squares, but not concatonate these
% into a single matrix. This is bacause we will be drawing these both in a
% single line of code. For more details use Screen DrawRect?
LowerLeftRect = CenterRectOnPointd(baseRect,screenXpixels * 0.40, yCenter - 110);
LowerMiddleRect = CenterRectOnPointd(baseRect,screenXpixels * 0.50, yCenter - 110);
LowerRightRect = CenterRectOnPointd(baseRect,screenXpixels * 0.60, yCenter - 110); 


%% Set up eye tracker

if eyelink
    
    %%%%%%%%%%
    % STEP 3 %
    %%%%%%%%%%
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    % make necessary changes to calibration structure parameters and pass
    % it to EyelinkUpdateDefaults for changes to take affect
    
    el=EyelinkInitDefaults(window);
    
    % Disable key output to Matlab window:
    %ListenChar(2);
    
    % We are changing calibration to a black background with white targets,
    % no sound and smaller targets
    el.backgroundcolour = GrayIndex(el.window); %(BlackIndex(el.window);
    el.msgfontcolour  = WhiteIndex(el.window);
    el.imgtitlecolour = WhiteIndex(el.window);
    el.targetbeep = 0;
    el.calibrationtargetcolour= WhiteIndex(el.window);
    % for lower resolutions you might have to play around with these values
    % a little. If you would like to draw larger targets on lower res
    % settings please edit PsychEyelinkDispatchCallback.m and see comments
    % in the EyelinkDrawCalibrationTarget function
    el.calibrationtargetsize= 1;
    el.calibrationtargetwidth=0.5;
    
    EyelinkUpdateDefaults(el);
    
    %%%%%%%%%%
    % STEP 4 %
    %%%%%%%%%%
    
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end
    
    % open file to record data to
    res = Eyelink('Openfile', edfFile);
    if res~=0
        fprintf('Cannot create EDF file ''%s'' ', edffilename);
        cleanup;
        return;
    end
    
    % make sure we're still connected.
    if Eyelink('IsConnected')~=1 && ~dummymode
        cleanup;
        return;
    end
    
    %%%%%%%%%%
    % STEP 5 %
    %%%%%%%%%%
    
    % SET UP TRACKER CONFIGURATION
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    winWidth = screenXpixels; %2880 (mac 15")
    winHeight = screenYpixels; %1800 (mac 15)
    Eyelink('command', 'add_file_preamble_text ''Recorded by EyelinkToolbox''');
    % This command is crucial to map the gaze positions from the tracker to
    % screen pixel positions to determine fixation
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, winWidth-1, winHeight-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, winWidth-1, winHeight-1);
    % set calibration type.
    Eyelink('command', 'calibration_type = HV13');
    Eyelink('command', 'generate_default_targets = YES');
    % set parser (conservative saccade thresholds)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    % set EDF file contents
    % 5.1 retrieve tracker version and tracker software version
    [v,vs] = Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    vsn = regexp(vs,'\d','match');
    
    if v == 3 && str2double(vsn{1}) == 4 % if EL 1000 and tracker version 4.xx
        % remote mode possible add HTARGET ( head target)
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT,HTARGET');
    else
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end
    % allow to use the big button on the eyelink gamepad to accept the
    % calibration/drift correction target
    Eyelink('command', 'button_function 5 "accept_target_fixation"');
    
    %%%%%%%%%%
    % STEP 6 %
    %%%%%%%%%%
    
    % Hide the mouse cursor and Calibrate the eye tracker
    Screen('HideCursorHelper', window);
    
    % enter Eyetracker camera setup mode, calibration and validation
    EyelinkDoTrackerSetup(el);
end

%% Load the images

% Read images to matrix and make textures for faster display
imgPath = [pwd filesep 'fractals' filesep];
imgType = '*.png'; % change based on image type
images  = dir([imgPath imgType]);
for idx = 1:length(images)
    sequence{idx} = imread([imgPath images(idx).name]);
    texture{idx} = Screen('MakeTexture',window, sequence{idx});
end

imgPath = [pwd filesep 'fillers' filesep];
imgType = '*.png'; % change based on image type
images  = dir([imgPath imgType]);
for idx = 1:length(images)
    sequence_fillers{idx} = imread([imgPath images(idx).name]);
    texture_fillers{idx} = Screen('MakeTexture',window, sequence_fillers{idx});
end

%% Stimuli
combo = ['B' 'C' 'D' 'F' 'G' 'H' 'J' 'K' 'L' 'M' 'N' 'P' 'Q' 'R' 'S' 'T' 'V' 'W' 'X' 'Z'];
filling = ['%' '%' '%' '%'];
cueList = {'press yes/no', 'press yes/no', 'press yes/no', 'recite all letters', 'say yes/no', 'no probe' 'no probe'};
%% Create experiment arrays


condMatrixBase_main = [sort(repmat([1 1 2 2 3 3 4 4 5 5], 1, 4)); repmat([1 1 1 1 2 2 2 2], 1, 5); repmat([1 1 2 2 1 1 2 2], 1, 5);repmat([1 2 1 2 1 2 1 2], 1, 5)]; % 5 types of trials, 2 match/nonmatch, 2 delay periods, 2 ITI periods
%condMatrixBase = [sort(repmat([3 3 3 3 3], 1, 2)); repmat([1 2], 1, 5)]; % 5 types of trials, 2 match/nonmatch
%condMatrixBase = [sort(repmat([1 2 3], 1, 8)); repmat([1 1 1 1 2 2 2 2], 1, 3); repmat([1 1 2 2 1 1 2 2], 1, 3); repmat([1 2 1 2 1 2 1 2], 1, 3)]; % 5 types of trials, 2 match/nonmatch, 2 delay periods, 2 ITI periods

condMatrixBase_nop = [sort(repmat([6 6 6 6 6 6 6 7], 1, 1)); repmat([1 1 1 1 2 2 2 2], 1, 1); repmat([1 2 3 4 1 2 3 4], 1, 1);repmat([1 2 1 2 1 2 1 2], 1, 1)]; % 2 types of trials, 2 match/nonmatch, 4 delay periods, 2 ITI periods

% combine matrices
condMatrixBase = [condMatrixBase_main condMatrixBase_nop];

% Number of trials per condition.
trialsPerCondition = 1; % 1 trials * (6 * 2 * 2 * 2) = 48 trials in a run.  8 trials per condition

% Duplicate the condition matrix to get the full number of trials
condMatrix = repmat(condMatrixBase, 1, trialsPerCondition);

% Get the size of the matrix
[~, numTrials] = size(condMatrix);

% Randomise the conditions
condMatrixShuffled = Shuffle(condMatrix,1);

%% Find canonical responses

% This is a 2 row matrix the first row records the letters or images (I disable the recording of images, because it takes too much space) we presented,
% the second row records the displayed image or letter,

respMat = cell(2, numTrials);

%% Wait for the experiment to start (spacebar)

borderColor = gray;
wordColor = black;
fliptimectr = 0;
responsearr = 999*ones(numTrials,1); % record subject's response 
rtarr       = 999*ones(numTrials,1); % record RT (ptb's timestamp)
incorrectResp = 0;
correctResp = 0;
fliptimes   = zeros(10000,1);
ListenChar(2); %turn off the keyboard
HideCursor;

%start out Cedrus

if scanner
    try
        port = '/dev/cu.usbserial-00002014'; %the right upper USB on your macbook;
        handle = CedrusResponseBox('Open', port);
        CedrusResponseBox('SetConnectorMode', handle, 'ReflectiveSinglePulse');
    catch
        sca;
        rethrow(lasterror)
    end
    runText = 'Experiment will begin soon';
else
    runText = 'Press spacebar to begin';
end

try
    Screen('FillRect',window, gray);
    Screen('TextSize', window , STS);
    DrawFormattedText(window, runText, 'center', 'center', wordColor);
    [VBLTimestamp, startrt, fliptime]= Screen('Flip',window);
catch
    ListenChar(0);
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end

fliptimectr=fliptimectr+1;
fliptimes(fliptimectr)= GetSecs;

% Set up the trigger for cedrus box or press space for keyboard version

if scanner
    CedrusResponseBox('ClearQueues', handle);
    CedrusResponseBox('FlushEvents', handle);
    while 1
        evt = CedrusResponseBox('WaitButtons', handle);
        if ~isempty(evt) && evt.button == 7 %sync pulse happens
            break;
        end
        WaitSecs(0.1)
    end
else
    [KeyIsDown, endrt, KeyCode] = KbCheck;
    while 1
        [KeyIsDown, endrt, KeyCode]=KbCheck;
        if KeyCode(nextkey)
            break;
        end
        WaitSecs(0.1);
    end
end

% starttime=GetSecs; %record the start time, then later add each schedule time to it
% scheduletime=0;
% WaitSecs(T_syncwaste);
% scheduletime = scheduletime + T_syncwaste;

%% Experimental Loop

%%%%%%%%%%
% STEP 7 %
%%%%%%%%%%

% Now starts running individual trials;
% You can keep the rest of the code except for the implementation
% of graphics and event monitoring
% Each trial should have a pair of "StartRecording" and "StopRecording"
% calls as well integration messages to the data file (message to mark
% the time of critical events and the image/interest area/condition
% information for the trial)

regressor_fix = zeros(numTrials,1);
regressor_cue = zeros(numTrials,1);
regressor_sample = zeros(numTrials,1);
regressor_delay = zeros(numTrials,1);
regressor_probe = zeros(numTrials,1);

taskBegin = GetSecs;

for trial = 1:numTrials
    
    if eyelink
        % determine current trial type and send msg to edf
        if condMatrixShuffled(1,trial) == 1
            Eyelink('Message', 'highload');
        elseif condMatrixShuffled(1,trial) == 2;
            Eyelink('Message', 'lowload');
        elseif condMatrixShuffled(1,trial) == 3
            Eyelink('Message', 'pictures');
        elseif condMatrixShuffled(1,trial) == 4
            Eyelink('Message', 'recite');
        elseif condMatrixShuffled(1,trial) == 5
            Eyelink('Message', 'answer');
        elseif condMatrixShuffled(1,trial) == 6
            Eyelink('Message', 'noprobe');
        elseif condMatrixShuffled(1,trial) == 7
            Eyelink('Message', 'catch');
        end
        
        % STEP 7.1
        % Sending a 'TRIALID' message to mark the start of a trial in Data
        % Viewer.  This is different than the start of recording message
        % START that is logged when the trial recording begins. The viewer
        % will not parse any messages, events, or samples, that exist in
        % the data file prior to this message.
        Eyelink('Message', 'TRIALID %d', trial);
        
        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "TRIAL %d/%d"', trial, numTrials);
        
        % Before recording, we place reference graphics on the host display
        % Must be offline to draw to EyeLink screen
        Eyelink('Command', 'set_idle_mode');
        
        % clear tracker display and draw box at center
        Eyelink('Command', 'clear_screen 0');
        
        % STEP 7.2
        % Do a drift correction at the beginning of each trial
        % Performing drift correction (checking) is optional for
        % EyeLink 1000 eye trackers.
        % EyelinkDoDriftCorrection(el);
        
        % STEP 7.3
        % start recording eye position (preceded by a short pause so that
        % the tracker can finish the mode transition)
        % The paramerters for the 'StartRecording' call controls the
        % file_samples, file_events, link_samples, link_events availability
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.5);
        Eyelink('StartRecording');
        % record a few samples before we actually start displaying
        % otherwise you may lose a few msec of data
        WaitSecs(0.5);
    end
    % I am moving the starttime here because there are many waitsecs before
    % that.  The schedule time will be wrong if I put the following before
    % the trial loop.
    
    if trial == 1
        start = GetSecs;
        starttime=GetSecs; %record the start time, then later add each schedule time to it
        scheduletime=0;
        WaitSecs(T_syncwaste);
        scheduletime = scheduletime + T_syncwaste;
    end
    
    % STEP 7.4
    %% Fixation (aka ITI)

    % condition
    conditionNum = condMatrixShuffled(1, trial);
    matchNum = condMatrixShuffled(2, trial);
    delayNum = condMatrixShuffled(3, trial);
    ITINum = condMatrixShuffled(4, trial);

    try
        Screen('DrawLines', window, allCoords, lineWidthPix, black, [xCenter yCenter], 0 );
        [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
        if eyelink
            Eyelink('Message', 'fixation');
        end
    catch
        ListenChar(0);
        Screen('CloseAll');
        psychrethrow(psychlasterror);
    end
    
    % record onset of fixation
    regressor_fix(trial,1) = scheduletime;
    
    % Save the fliptime, schedule the next flip.
    fliptimectr = fliptimectr+1;
    fliptimes(fliptimectr) = fliptime;
    scheduletime = scheduletime + T_fix(ITINum);

    %% Cue display
    currCue = cueList{conditionNum};
    runText = currCue;
    
    try
        Screen('FillRect',window, borderColor);
        Screen('TextSize', window , STS);
        DrawFormattedText(window, runText, 'center', 'center', wordColor);
        [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
        if eyelink
            Eyelink('Message', 'cue');
        end
    catch 
        ListenChar(0);
        Screen('CloseAll');
        psychrethrow(psychlasterror);
    end
    
    % record onset of cue
    regressor_cue(trial,1) = scheduletime;
    
    % Save the fliptime, schedule the next flip.
    fliptimectr = fliptimectr+1;
    fliptimes(fliptimectr) = fliptime;
    scheduletime = scheduletime + T_cue;
    
    %% Sample display

    if conditionNum == 1 || conditionNum == 4 || conditionNum == 5 || conditionNum == 6 || conditionNum == 7 
        letter_set = combo(randperm(length(combo)));
        letterMatchSample = letter_set(1,1:6);  %matched response letters set
        letterNonMatchSample = letter_set(1,7:end);  %non-match response letters set
        respMat{1,trial} = num2cell(letter_set(1,1:6)); %save the sample to the first row of the matrix
        try
            Screen('DrawLines', window, allCoords, lineWidthPix, black, [xCenter yCenter], 0 );
            DrawFormattedText(window, letterMatchSample(1,1),screenXpixels * 0.40 -50, yCenter + 70,wordColor); %lower left
            DrawFormattedText(window, letterMatchSample(1,2),'center', yCenter + 70,wordColor); %lower center
            DrawFormattedText(window, letterMatchSample(1,3),screenXpixels * 0.60 -25, yCenter + 70,wordColor); %lower right
            DrawFormattedText(window, letterMatchSample(1,4),screenXpixels * 0.40 -50, yCenter - 140,wordColor); %upper left
            DrawFormattedText(window, letterMatchSample(1,5),'center', yCenter - 140,wordColor); %upper center
            DrawFormattedText(window, letterMatchSample(1,6),screenXpixels * 0.60 -25, yCenter - 140,wordColor); %upper right
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'sample');
            end
        catch
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
        
    elseif conditionNum == 2 % low load condition
        letter_set = combo(randperm(length(combo)));
        letterMatchSample = letter_set(1,1:2);  %matched response letters set
        letterMatchSampleFilling = [letterMatchSample filling];
        letterMatchSampleFilling = letterMatchSampleFilling(randperm(length(letterMatchSampleFilling)));
        letterNonMatchSample = letter_set(1,3:end);  %non-match response letters set
        respMat{1,trial} = num2cell(letter_set(1,1:2)); %save the sample to the first row of the matrix        
        try
            Screen('DrawLines', window, allCoords, lineWidthPix, black, [xCenter yCenter], 0 );
            DrawFormattedText(window, letterMatchSampleFilling(1,1),screenXpixels * 0.40 -50, yCenter + 70,wordColor);
            DrawFormattedText(window, letterMatchSampleFilling(1,2),'center', yCenter + 70,wordColor);
            DrawFormattedText(window, letterMatchSampleFilling(1,3),screenXpixels * 0.60 -25, yCenter + 70,wordColor);
            DrawFormattedText(window, letterMatchSampleFilling(1,4),screenXpixels * 0.40 -50, yCenter - 140,wordColor);
            DrawFormattedText(window, letterMatchSampleFilling(1,5),'center', yCenter - 140,wordColor);
            DrawFormattedText(window, letterMatchSampleFilling(1,6),screenXpixels * 0.60 -25, yCenter - 140,wordColor);
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'sample');
            end
        catch
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
    else
        fractal_set = texture(randperm(length(texture)));
        %respMat{1,trial}= fractal_set(1,1:3);
        imgTex = [fractal_set(1,1:3) texture_fillers(1,1:3)]; %only show 3 fractals
        imgTex = imgTex(randperm(length(imgTex)));   
        imagesMatchSample = fractal_set(1,1:3);
        imagesNonmatchSample = fractal_set(1,4:end);
        try
            Screen('FillRect',window,borderColor);
            Screen('DrawLines', window, allCoords, lineWidthPix, black, [xCenter yCenter], 0 );
            Screen('DrawTextures', window, imgTex{1}, [], UpperLeftRect, [], [], [], []);
            Screen('DrawTextures', window, imgTex{2}, [], UpperMiddleRect, [], [], [], []);
            Screen('DrawTextures', window, imgTex{3}, [], UpperRightRect, [], [], [], []);
            Screen('DrawTextures', window, imgTex{4}, [], LowerLeftRect, [], [], [], []);
            Screen('DrawTextures', window, imgTex{5}, [], LowerMiddleRect, [], [], [], []);
            Screen('DrawTextures', window, imgTex{6}, [], LowerRightRect, [], [], [], []);
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'sample');
            end
        catch
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
    end
    
    % record onset of sample
    regressor_sample(trial,1) = scheduletime;
    
    % Save the fliptime, schedule the next flip.
    fliptimectr = fliptimectr+1;
    fliptimes(fliptimectr) = fliptime;
    scheduletime = scheduletime + T_sample;
    
    %% delay
    try
        Screen('FillRect',window,borderColor);
        [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
        if eyelink
            Eyelink('Message', 'delay');
        end
    catch
        ListenChar(0);
        Screen('CloseAll');
        psychrethrow(psychlasterror);
    end
    
    % record onset of delay
    regressor_delay(trial,1) = scheduletime;
    
    if conditionNum == 6 || conditionNum == 7
        
        % Save the fliptime, schedule the next flip.
        fliptimectr = fliptimectr+1;
        fliptimes(fliptimectr) = fliptime;
        scheduletime = scheduletime + T_nopDelay(delayNum);
    else
        % Save the fliptime, schedule the next flip.
        fliptimectr = fliptimectr+1;
        fliptimes(fliptimectr) = fliptime;
        scheduletime = scheduletime + T_delay(delayNum);
    end
    
    %% Response phase (probe presented)
    
    % Regular, high load, finger press match or nonmatch
    if conditionNum == 1
        if matchNum == 1
            runText = letterMatchSample(randi(numel(letterMatchSample)));
        else
            runText = letterNonMatchSample(randi(numel(letterNonMatchSample)));
        end
        
        respMat{2,trial} = runText; %save the display letter
        
        try
            Screen('FillRect',window, borderColor);
            Screen('TextSize', window , STS);
            DrawFormattedText(window, runText, 'center', 'center', wordColor);
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'probe');
            end
        catch %#ok<CTCH>
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
        
    % Regular, low load, finger press match or nonmatch
        
    elseif conditionNum == 2
        if matchNum == 1
            runText = letterMatchSample(randi(numel(letterMatchSample)));
        else
            runText = letterNonMatchSample(randi(numel(letterNonMatchSample)));
        end
        
        respMat{2,trial} = runText; %save the display letter
        
        try
            Screen('FillRect',window, borderColor);
            Screen('TextSize', window , STS);
            DrawFormattedText(window, runText, 'center', 'center', wordColor);
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'probe');
            end
        catch 
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
    
    elseif conditionNum == 3
        if matchNum == 1
            imgTex = imagesMatchSample(randi(numel(imagesMatchSample)));
        else
            imgTex = imagesNonmatchSample(randi(numel(imagesNonmatchSample)));
        end
        
        try
        %respMat{2,trial} = imgTex; %save the display image
        Screen('FillRect',window, borderColor);
        dstRect = CenterRect(baseRect, theRect);
        Screen('DrawTexture', window, imgTex{1}, [], dstRect);
        [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
        if eyelink
            Eyelink('Message', 'probe');
        end
        catch
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end 
        
     % recite probe is a question mark
        
    elseif conditionNum == 4
        try
            respMat{2,trial} = '?'; %subject recite the letters
            Screen('FillRect',window, borderColor);
            Screen('TextSize', window , STS);
            DrawFormattedText(window, '?', 'center', 'center', wordColor);
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'probe');
            end
        catch %#ok<CTCH>
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
        
     % verbally say yes or no  
     
    elseif conditionNum == 5
        if matchNum == 1
            runText = letterMatchSample(randi(numel(letterMatchSample)));
        else
            runText = letterNonMatchSample(randi(numel(letterNonMatchSample)));
        end
        
        respMat{2,trial} = runText; %save the display letter
        
        try
            Screen('FillRect',window, borderColor);
            Screen('TextSize', window , STS);
            DrawFormattedText(window, runText, 'center', 'center', wordColor);
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'probe');
            end
        catch %#ok<CTCH>
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
        
    elseif conditionNum == 6
        try
            respMat{2,trial} = 'x'; %indicate no probe trial
            Screen('FillRect',window, borderColor);
            Screen('TextSize', window , STS);
            DrawFormattedText(window, 'stop', 'center', 'center', [1 0 0]);
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'probe');
            end
        catch %#ok<CTCH>
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
        
    elseif conditionNum == 7 %catch trial
        if matchNum == 1
            runText = letterMatchSample(randi(numel(letterMatchSample)));
        else
            runText = letterNonMatchSample(randi(numel(letterNonMatchSample)));
        end
        
        respMat{2,trial} = runText; %save the display letter
        
        try
            Screen('FillRect',window, borderColor);
            Screen('TextSize', window , STS);
            DrawFormattedText(window, runText, 'center', 'center', wordColor);
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'probe');
            end
        catch %#ok<CTCH>
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
        
    end
    
    % record onset of probe
    regressor_probe(trial,1) = scheduletime;

    % Save the fliptime, schedule the next flip.
    fliptimectr = fliptimectr+1;
    fliptimes(fliptimectr) = fliptime;
    scheduletime = scheduletime + T_probe; %comment this line if the resp is depends on the subject's speed
    
    %% Process input
    if scanner
        CedrusResponseBox('FlushEvents', handle);
    end
    
    while GetSecs < starttime+scheduletime-FlipInt
        if scanner
            evt = CedrusResponseBox('GetButtons', handle);
            [KeyIsDown, endrt, KeyCode]=KbCheck;
            if ~isempty(evt) && evt.button > 0 && evt.button < 5 && responsearr(trial, 1) == 999
                %resp_time2 = evt.ptbfetchtime; %Cedrus time, also test out evt.ptbtime
                if evt.button == leftkey %leftkey is 2
                    resp = 1; %match
                elseif evt.button == rightkey %rightkey is 3
                    resp = 2; %non-match
                else
                    resp = 5; %if any weird button press happens
                    %error('hint: response button error, check line 579')
                end
                responsearr(trial,1)=resp;
                rtarr(trial,1) = endrt-startrt; %record the response time
            end
        else
            [KeyIsDown, endrt, KeyCode] = KbCheck;
            if KeyCode(killkey)
                ShowCursor;
                Screen('CloseAll');
                %save(['QuitDump' date '.m']);
                return;
            end
            if responsearr(trial,1)==999 %don't let them change their response, only make one for the first time
                if KeyCode(leftkey)
                    responsearr(trial,1)=1;
                    rtarr(trial,1) =endrt-startrt;
                end
                if KeyCode(rightkey)
                    responsearr(trial,1)=2;
                    rtarr(trial,1)=endrt-startrt;
                end
            end
            WaitSecs(0.01);
        end
    end
    %% If test phase, accrue correct responses (also create feedback text)
    if conditionNum == 1 || conditionNum == 2 || conditionNum == 3
        if responsearr(trial,1) == condMatrixShuffled(2, trial);
            correctResp = correctResp + 1; %counts # of rewarded corr responses
            performance = 'Correct';
        else
            incorrectResp = incorrectResp + 1;%counts # of incorr responses
            if responsearr(trial,1) == 999;
                performance = 'Too slow';
            else
                performance = 'Incorrect';
            end
        end
    else
        performance = 'good';
    end
    tally = [correctResp; incorrectResp]'; %note that I transpose it
    
    %% Small interval
    if scanner == false
    try
        Screen('FillRect',window,borderColor);
        %Screen('DrawTexture', window, imgFixTex);
        [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
    catch %#ok<CTCH>
        ListenChar(0);
        Screen('CloseAll');
        psychrethrow(psychlasterror);
    end
    % Save the fliptime, schedule the next flip.
    fliptimectr = fliptimectr+1;
    fliptimes(fliptimectr) = fliptime;
    scheduletime = scheduletime + T_interval;
    %% Provide Feedback
        runText = performance;
        try
            Screen('FillRect',window, borderColor);
            Screen('TextSize', window , STS);
            DrawFormattedText(window, runText, 'center', 'center', wordColor);
            %CenterText2(window, runText, borderColor, 0, 0);
            [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
            if eyelink
                Eyelink('Message', 'feedback');
            end
        catch %#ok<CTCH>
            ListenChar(0);
            Screen('CloseAll');
            psychrethrow(psychlasterror);
        end
        % Save the fliptime, schedule the next flip.
        fliptimectr = fliptimectr+1;
        fliptimes(fliptimectr) = fliptime;
        scheduletime = scheduletime + T_feed;
    end
    %% STEP 7.6
    if eyelink
        % add 100 msec of data to catch final events and blank display
        WaitSecs(0.1);
        Eyelink('StopRecording');
        Screen('FillRect', window, el.backgroundcolour);
        Screen('Flip', window);
        
        % STEP 7.7
        % Send out necessary integration messages for data analysis
        % See "Protocol for EyeLink Data to Viewer Integration-> Interest
        % Area Commands" section of the EyeLink Data Viewer User Manual
        % IMPORTANT! Don't send too many messages in a very short period of
        % time or the EyeLink tracker may not be able to write them all
        % to the EDF file.
        % Consider adding a short delay every few messages.
        WaitSecs(0.001);
        
        % Send an integration message so that an image can be loaded as
        % overlay backgound when performing Data Viewer analysis.  This
        % message can be placed anywhere within the scope of a trial (i.e.,
        % after the 'TRIALID' message and before 'TRIAL_RESULT')
        % See "Protocol for EyeLink Data to Viewer Integration -> Image
        % Commands" section of the EyeLink Data Viewer User Manual.
        % Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', imgfile,  winWidth/2, winHeight/2);
        
        % interest areas
        %Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, rect1(1), rect1(2), rect1(3), rect1(4),'target');
        %Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 2, rect2(1), rect2(2), rect2(3), rect2(4),'distractor');
        % Send messages to report trial condition information
        % Each message may be a pair of trial condition variable and its
        % corresponding value follwing the '!V TRIAL_VAR' token message
        % See "Protocol for EyeLink Data to Viewer Integration-> Trial
        % Message Commands" section of the EyeLink Data Viewer User Manual
        % WaitSecs(0.001);
        
        %Eyelink('Message', '!V TRIAL_VAR index %d', trial);
        %Eyelink('Message', '!V TRIAL_VAR amplitude %d', perm(3)); % 5 or 10
        %Eyelink('Message', '!V TRIAL_VAR saccade %s', saccade);   % pro or anti
        %Eyelink('Message', '!V TRIAL_VAR feedback %s', feedbackpos); % left or right?
        
        % STEP 7.8
        %     Sending a 'TRIAL_RESULT' message to mark the end of a trial in
        %     Data Viewer. This is different than the end of recording message
        %     END that is logged when the trial recording ends. The viewer will
        %     not parse any messages, events, or samples that exist in the data
        %     file after this message.
        
        Eyelink('Message', 'TRIAL_RESULT 0');
    end
end

%% End screen

% End of Experiment; close the file first
% close graphics window, close data file and shut down tracker
if eyelink
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch %#ok<*CTCH>
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end
end

% this is for the scanner
if scanner
    CedrusResponseBox('CloseAll')
end
borderColor = white;

% try
%     Screen('FillRect',window,borderColor);
%     Screen('DrawTexture', window, imgEndTex);
%     [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
%     % add 100 msec of data to catch final events and blank display
% catch %#ok<CTCH>
%     ListenChar(0);
%     Screen('CloseAll');
%     psychrethrow(psychlasterror);
% end

try
    Screen('FillRect',window, borderColor);
    Screen('TextSize', window , STS);
    DrawFormattedText(window, 'end', 'center', 'center', wordColor);
    [vblTimestamp, startrt, fliptime] = Screen('Flip', window, starttime+scheduletime-FlipInt/2, 0, 0);
    % add 100 msec of data to catch final events and blank display
catch %#ok<CTCH>
    ListenChar(0);
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
% Save the fliptime, schedule the next flip.
fliptimectr = fliptimectr+1;
fliptimes(fliptimectr) = fliptime;
scheduletime = scheduletime + T_end;


if eyelink
    %%%%%%%%%%
    % STEP 9 %
    %%%%%%%%%%
    
    % run cleanup function (close the eye tracker and window).
    %cleanup;
    
    % Shutdown Eyelink:
    Eyelink('Shutdown');
    %
    % !./edf2asc 1.edf
    % %not really sure about saccade
    % !grep cue 1.asc > cue.txt
end

%% Save everything
save(outputFname)
ListenChar(0);
theEnd = GetSecs;
Screen('CloseAll');
%accuracy = correctResp/numTrials;