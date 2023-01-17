%function esc = MOTtask(subject)

esc = 0;

trl.exp_number = 1;
trl.sub_number = 1;
trl.blk_number = 4;
arduinoFlag = 0;

% sca
% ListenChar(0);
% keyboard

try

    %% set design
    rng('default'); % reset random number generator in case it is in legacy mode
    rng('shuffle'); % use current time as seed
    trl = setDesign(trl);
    
    %% subject data folder
    [esc dat] = setSubjectData(esc, trl);
    if esc, prepareExit; return; end

    % sync
    if arduinoFlag
        a = arduino();
        writeDigitalPin(a,'D13',1); %high
    end
    
    %% set ptb
    setPsychtoolbox;
    
    %% open window
    gfx.background = [128 128 128];
    screen = 0; % home: 1; 
    windowPtr = Screen('OpenWindow', screen, gfx.background); %, [0 0 1920 1080]   

    %% graphic
    [esc gfx] = setGFX(windowPtr, trl, dat, gfx);
    if esc, prepareExit; return; end
    
    %% mot 
    mot = setMot(windowPtr, trl, gfx);
    
    %% checker
    che = setChecker(windowPtr, gfx);
 
    %% instruction
    esc = showInstruction(esc, windowPtr, gfx, trl); 
    if esc, prepareExit; return; end
     
    %% loop thru all trials
    
    exp = GetSecs; time.E = GetSecs - exp;
    if arduinoFlag
        writeDigitalPin(a,'D13',0); %low
        WaitSecs(0.1);
        writeDigitalPin(a,'D13',1); %high
    end

    if trl.blk_number == 0 % block condition
        showWait = 1;
        if showWait %i1 == 1 | ~mod(i1-1,5)
            %% show wait
            timeWait = 180; % in s
            numFramesWait = ceil(timeWait/gfx.timeFrame);

            DrawFormattedText(windowPtr, ['Baseline'], 'center', gfx.cy, [0 255 255], [], [], [], 2.4);
            Screen('DrawingFinished', windowPtr);
            wait = Screen('Flip', windowPtr); wait1 = wait;

            for iw = (timeWait-1):-1:1
                DrawFormattedText(windowPtr, ['Baseline\n' num2str(iw) ' s'], 'center', gfx.cy, [0 255 255], [], [], [], 3);
                Screen('DrawingFinished', windowPtr);
                wait = Screen('Flip', windowPtr, [wait + 1 - gfx.timeFrame/2]);
            end
        end
    else % block condition
        for i1 = 1:size(mot.trajectory,4)

            tri = GetSecs; time.triE = tri - exp; time.triT = GetSecs - tri; % current time stamp; time relative to Experiment start; relative to Trials start; relative to preceding event

            trl.num = i1;

            showWait = 1;
            if showWait %i1 == 1 | ~mod(i1-1,5)
                %% show wait
                timeWait = 6% in s
                numFramesWait = ceil(timeWait/gfx.timeFrame);

                DrawFormattedText(windowPtr, ['Trial ' num2str(i1) '\n' num2str(timeWait) ' s'], 'center', gfx.cy, [0 255 255], [], [], [], 3);
                Screen('DrawingFinished', windowPtr);
                wait = Screen('Flip', windowPtr); wait1 = wait;

                for iw = (timeWait-1):-1:1
                    DrawFormattedText(windowPtr, ['Trial ' num2str(i1) '\n' num2str(iw) ' s'], 'center', gfx.cy, [0 255 255], [], [], [], 3);
                    Screen('DrawingFinished', windowPtr);
                    wait = Screen('Flip', windowPtr, [wait + 1 - gfx.timeFrame/2]);
                end

                time.waitE = wait1 - exp; time.waitT = wait - tri; %% QUES: `wait1` for `time.waitT`? 
            else
                time.waitE = 0; time.waitT = 0; wait1 = GetSecs;
            end

            %% show targets
            x1 = mot.trajectory(1,:,1,i1); % coordinates in frame 1 in trial i3
            y1 = mot.trajectory(2,:,1,i1);

            if trl.numTargets(i1)==0
                trl.targetIdx = nan;
                colorMat1 = mot.colorDots;
            else
                targetIdx = randperm(mot.numDots);
                trl.targetIdx = targetIdx( 1:trl.numTargets(i1) ); % 1 random indices if 1 target; 4 random indices if 4 targets
                colorMat1 = mot.colorDots;
                colorMat1(:,trl.targetIdx)= repmat( gfx.colorTargets,1,trl.numTargets(i1) ); % the default dot color at the indexed positions is replaced by the target color
            end
            timeTargets      = 2 % in s
            numFramesTargets = ceil(timeTargets/gfx.timeFrame);

            % 1st target frame
            Screen('DrawDots', windowPtr, [x1; y1], 2*mot.radiusDot, colorMat1, [], 3);
            %        Screen('FillRect', windowPtr, [0 255 255], gfx.fixXY);
            Screen('DrawingFinished', windowPtr);
            if showWait
                target = Screen('Flip', windowPtr, [wait + 1 - gfx.timeFrame/2]);
            else
                target = Screen('Flip', windowPtr);
            end

            % remaining target frames
            for it = 1:numFramesTargets-1
                Screen('DrawDots', windowPtr, [x1; y1], 2*mot.radiusDot, colorMat1, [], 3);
                %            Screen('FillRect', windowPtr, [0 255 255], gfx.fixXY);
                Screen('DrawingFinished', windowPtr);
                Screen('Flip', windowPtr);
            end

            time.targetE = target - exp; time.targetT = target - tri; time.waitP = target - wait1;

            %% loop thru trajecotry of a single trial
            % show dots and/or checker
            flick = 1;
            for i2 = 1:size(mot.trajectory,3)
                if trl.checkerPres(i1)==1
                    if ~mod(i2,che.numFramesPerFlick)
                        flick = 3 - flick;
                    end
                    if trl.checkerSide(i1) == 1
                        Screen('DrawTexture', windowPtr, che.tR(flick));
                    elseif trl.checkerSide(i1) == 2
                        Screen('DrawTexture', windowPtr, che.tL(flick));
                    end
                end

                x = mot.trajectory(1,:,i2,i1);
                y = mot.trajectory(2,:,i2,i1);
                % draw disks with 'DrawDots', requires Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                Screen('DrawDots', windowPtr, [x; y], 2*mot.radiusDot, [255 255 255], [], 3);

                %            Screen('FillRect', windowPtr, [0 255 255], gfx.fixXY);
                Screen('DrawingFinished', windowPtr);
                if i2 == 1
                    motion = Screen('Flip',windowPtr); % time of motion onset
                else
                    Screen('Flip',windowPtr);
                end
                %            imageArray=Screen('GetImage', windowPtr);

            end
            %save('imageArray', 'imageArray');
            time.motionE = motion - exp; time.motionT = motion - tri; time.targetP = motion - target;

            %% show probe
            xe = mot.trajectory(1,:,end,i1); % coordinates in frame end in trial i1
            ye = mot.trajectory(2,:,end,i1);

            colorMat2 = mot.colorDots;
            if trl.numTargets(i1)==0
                trl.probeIdx = nan;
            else
                if trl.probeMatch(i1) == 1
                    trl.probeIdx = targetIdx(1); % the first entry is always a target
                    colorMat2(:,targetIdx(1))= gfx.colorProbe; % the default dot color at the indexed positions is replaced by the target color
                elseif trl.probeMatch(i1) == 0
                    trl.probeIdx = targetIdx(end); % the last entry is always a distractor
                    colorMat2(:,targetIdx(end))= gfx.colorProbe; % the default dot color at the indexed positions is replaced by the target color
                end
            end

            timeProbe      = 2; % in s
            numFramesProbe = ceil(timeProbe/gfx.timeFrame);
            %for ip = 1:numFramesProbe
            Screen('DrawDots', windowPtr, [xe; ye], 2*mot.radiusDot, colorMat2, [], 3);
            %            Screen('FillRect', windowPtr, [0 255 255], gfx.fixXY);
            Screen('DrawingFinished', windowPtr);
            probe = Screen('Flip', windowPtr, [motion + size(mot.trajectory,3)*gfx.timeFrame - gfx.timeFrame/2]);
            %end
            % check assigned numbers
            % [trl.probeMatch(i1) trl.targetIdx trl.probeIdx]

            time.probeE = probe - exp; time.probeT = probe - tri; time.motionP = probe - motion;

            %% get response
            FlushEvents(['mouseUp'],['mouseDown'],['keyDown'],['autoKey'],['update']);
            [secs, dat.keyCode] = KbWait(-1, 2, probe+trl.RTmaxWait); % arg1: device1, arg2: 2=key press (not key release), arg2: wait only n secs % keyCode is a vector, its length is the number of keys on the keyboard
            dat.RT = secs - probe; % in s
            esc = checkEscape(dat.keyCode);
            if esc, prepareExit; return; end

            %% check response
            [dat] = analyseResp(trl, dat);

            %% show feedback
            colorMat3 = mot.colorDots;

            if trl.numTargets(i1)==0
                if dat.IDcor ~= 1 | dat.RTcor ~= 1
                    colorMat3 = repmat(gfx.colorFeedInc, 1, size(colorMat3,2));
                end
            else
                if dat.IDcor == 1 & dat.RTcor == 1
                    colorMat3(:,trl.probeIdx) = gfx.colorFeedCor;
                else
                    colorMat3(:,trl.probeIdx)= gfx.colorFeedInc;
                end
            end

            % show feedback
            Screen('DrawDots', windowPtr, [xe; ye], 2*mot.radiusDot, colorMat3, [], 3);
            %        Screen('FillRect', windowPtr, [0 255 255], gfx.fixXY);
            Screen('DrawingFinished', windowPtr);
            feedb = Screen('Flip', windowPtr, [probe + numFramesProbe*gfx.timeFrame - gfx.timeFrame/2] );

            time.feedbE = feedb - exp; time.feedbT = feedb - tri; time.probeP = feedb - probe;

            % sound feedback
            %         if dat.IDcor == 1 & dat.RTcor == 1
            %             Beeper(gfx.BeepFqH, 0, gfx.BeepDur);
            %         else
            %             Beeper(gfx.BeepFqL, gfx.BeepVol, gfx.BeepDur);
            %         end
            % [trl.probeMatch(i1) trl.targetIdx trl.probeIdx   dat.IDcor == 1 & dat.RTcor == 1]

            % last but 2nd feedback frame
            timeFeed      = 2; % in s
            numFramesFeed = ceil(timeFeed/gfx.timeFrame)-2; %-1 because this flip and the finish flip are two flips
            Screen('DrawDots', windowPtr, [xe; ye], 2*mot.radiusDot, colorMat3, [], 3);
            %        Screen('FillRect', windowPtr, [0 255 255], gfx.fixXY);
            Screen('DrawingFinished', windowPtr);
            Screen('Flip', windowPtr, [feedb + numFramesFeed*gfx.timeFrame - gfx.timeFrame/2] );

            % last feedback frame, to finish presentation
            Screen('DrawDots', windowPtr, [xe; ye], 2*mot.radiusDot, colorMat3, [], 3);
            %        Screen('FillRect', windowPtr, [0 255 255], gfx.fixXY);
            Screen('DrawingFinished', windowPtr);
            Screen('Flip', windowPtr);

            feedbEnd = GetSecs; time.feedbP = feedbEnd - feedb;

            %% save data
            dat = saveData(trl, dat, time);

            if mod(i1,trl.fbkT) == 0
                showFeedback(esc, windowPtr, gfx, trl, dat);
            end

        end % end of trial loop
    end % block condition

    expEnd = GetSecs - exp;
    if arduinoFlag
        writeDigitalPin(a,'D13',0); %low
        WaitSecs(0.1);
        writeDigitalPin(a,'D13',1); %high
    end
    save( dat.time_name , 'expEnd');

    
    %%
    DrawFormattedText(windowPtr, 'End', 'center', 'center', [0 255 255], [], [], [], 3);
    Screen('Flip', windowPtr);
    WaitSecs(1);
    % exit
    prepareExit;
    esc = 0;
    return;
catch
    % this "catch" section executes in case of an error in the "try" section above.
    prepareExit;
    psychrethrow(psychlasterror);
    esc = 0;
    return;
end % try catch end

%end % main function's end    

%% set design
function trl = setDesign(trl)

blockSize = 6;
bs        = ceil(blockSize/2);

probeMatch = [zeros(bs,1); ones(bs,1)];

numTargets = [0 2 3 4 5];
nt         = numTargets(randperm(length(numTargets)));

design = [];

for i = 1:length(numTargets)
    ntb = nt(i)*ones(blockSize,1);
    probeMatch = probeMatch(randperm(length(probeMatch)));
    design = [design; [ntb probeMatch]];
end
   
    trl.numTrials = length(design);

    trl.numTargets  = design(:,1);
    trl.probeMatch  = design(:,2);
    trl.checkerPres = zeros(trl.numTrials,1);
    trl.checkerSide = zeros(trl.numTrials,1);
    
    trl.RTmin           = 0.20; % in s
    trl.RTmax           = 2;
    trl.RTmaxWait       = 2;
    
    trl.fbkT  = 20;
end

%%
function setPsychtoolbox
    %HideCursor;
    KbName('UnifyKeyNames'); 
    Screen('Preference', 'Verbosity', 0);
    Screen('Preference', 'VisualDebugLevel',0);
    Screen('Preference', 'SkipSyncTests',1);
    Screen('Preference', 'TextRenderer', 1);
    Screen('Preference', 'SuppressAllWarnings', 1);
    %PsychDebugWindowConfiguration;
    KbName('UnifyKeyNames');
end

%% close windows and prepare for exit
function prepareExit
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    FlushEvents(['mouseUp'],['mouseDown'],['keyDown'],['autoKey'],['update']);
end

%% subject data folder
function [esc dat] = setSubjectData(esc, trl)
    data_folder = [pwd '\data'];
    if ~exist(data_folder, 'dir')
        mkdir(data_folder);
    end
    dat.sub_folder = [data_folder '\s' num2str(trl.sub_number)];
    if ~exist(dat.sub_folder, 'dir')
        mkdir(dat.sub_folder);
    else
        
    end 
    dat.log_name   = [dat.sub_folder '\s' num2str(trl.sub_number) 'b' num2str(trl.blk_number) '.mat'];    
    % fprintf('dat.file_name: %s\n', dat.log_name); % print file name
    if exist(dat.log_name, 'file') > 0 && trl.sub_number ~= 99 % if filename exists, exit, unless subject is 99        
        fprintf('file exists: %s\n', dat.log_name);
        dat.log_name = [];
        esc = 1;
    end

    dat.time_name   = [dat.sub_folder '\s' num2str(trl.sub_number) 'b' num2str(trl.blk_number) 'EndTime.mat'];
    if exist(dat.time_name, 'file') > 0 && trl.sub_number ~= 99 % if filename exists, exit, unless subject is 99
        fprintf('file exists: %s\n', dat.time_name);
        dat.time_name = [];
        esc = 1;
    end
    
    dat.blockdata =[];
    
end

%% save data
function dat = saveData(trl, dat, time)

    tri_line = [trl.exp_number trl.sub_number trl.num ...
                trl.numTargets(trl.num) trl.probeMatch(trl.num) trl.checkerPres(trl.num) trl.checkerSide(trl.num) ...
                dat.IDcor dat.RTcor dat.RT...
                time.triE time.waitE time.targetE time.motionE time.probeE ...
                time.triT time.waitT time.targetT time.motionT time.probeT ...
                time.waitP time.targetP time.motionP time.probeP time.feedbP
                ];
            
    blockdata = [dat.blockdata; tri_line];
    dat.blockdata = blockdata;
    save( [dat.log_name] , 'blockdata');

end

%% sets the screen's physical parameters  
function [esc gfx] = setGFX(windowPtr, trl, dat, gfx)
    Screen('TextFont', windowPtr, 'Arial');
    Screen('TextStyle', windowPtr , 0);
    gfx.timeFrame      = Screen('GetFlipInterval', windowPtr);

    % get the pixel width and height for a fullscreen output window
    [gfx.width, gfx.height] = Screen('WindowSize', 0);
    
    gfx.cx = gfx.width/2;
    gfx.cy = gfx.height/2;
    
    gfx.screenH = 27.4;
    gfx.screenV = 19.2;

    gfx.widthPixelPerCM  = gfx.width/gfx.screenH;
    gfx.heightPixelPerCM = gfx.height/gfx.screenV;
   
	gfx.subDistance     = 57; 
	gfx.deg2pix         = 1 / ((180/pi)*atan((gfx.screenH/gfx.height)/(gfx.subDistance)));
	gfx.mm2pix          = gfx.height / (10*gfx.screenH);
    scrH2               = gfx.height/2;
    scrV2               = gfx.width/2;
    
    gfx.colorDots = [255 255 255]';
    gfx.colorTargets = [0 0 255]';
    gfx.colorProbe = [0 255 255]';
    gfx.colorFeedCor = [0 255 0]';
    gfx.colorFeedInc = [255 0 0]';
    
    % sound feedback
    gfx.BeepDur = 0.1;% duration
    gfx.BeepVol = .6;  % volume
    gfx.BeepFqL = 130; % bottom frequency
    gfx.BeepFqH = 1050;% upper frequency

    % if abs(gfx.monitorFreq - 60) > 2, esc = 1; else esc = 0; end
    
    fixlength = 40; % pix
    if mod(fixlength,2) == 1
        fixlength = fixlength+1; % make sure that fixlength is even, otherwise the two lines will not cross at center
    end
    gfx.fixlength = fixlength;
    horiz = [gfx.cx-fixlength/2-1; gfx.cy; gfx.cx+fixlength/2; gfx.cy+1];
    vertic = [gfx.cx-1; gfx.cy-fixlength/2; gfx.cx; gfx.cy+fixlength/2+1];
    gfx.fixXY = [vertic, horiz];
    
    esc = 0;
end

%%
function mot = setMot(windowPtr, trl, gfx)

    %'DrawDots' requires: 
    Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % gfx.width gfx.height gfx.timeFrame
        
    mot.timeMotion = 10;
    mot.numFrames = round(mot.timeMotion/gfx.timeFrame);
    
    mot.radiusDot = round(0.7*gfx.widthPixelPerCM/2); %mvl 10, radius of drawn dot in pix
    mot.radiusDotA= round(1.0*gfx.widthPixelPerCM/2); r = mot.radiusDotA; % mvl 14, radius of area surrounding dot
                        v = round(4.0* gfx.widthPixelPerCM/ (1/gfx.timeFrame)); % mvl 2, initial speed (pixels per frame)
    mot.numDots=10;    n = mot.numDots;% number of circles to draw
    mot.colorDots = repmat(gfx.colorDots,1,n);
    
    fadeIn = 1;
    timeVChange = 1;
    numFramesCh = round(timeVChange/gfx.timeFrame);
    vChange     = linspace(0,v,numFramesCh);
    
    physics = 0; % set to 1 for realistic physics, 0 for constant speed
    dF = round(12.0*gfx.widthPixelPerCM/2); % mvl 100; dot Field is a square with d=2*dF in pixels

    xmax=gfx.cx+dF-r; % maximum x boundary
    xmin=gfx.cx-dF+r; % minimum x boundary
    ymax=gfx.cy+dF-r; % maximum y boundary
    ymin=gfx.cy-dF+r; % minimum y boundary
    
    trajectory = zeros(2, n, mot.numFrames, trl.numTrials);
    
    %%
    for ix1 = 1:trl.numTrials
        
        % find inital x,y that do not overlap
        o = 1; % overlap
        while nnz(o)>0
            % chose positions
            x=rand(1,n)*(2*dF-2*r) +xmin; % initial x positions
            y=rand(1,n)*(2*dF-2*r) +ymin; % initial y positions
            % check overlap
            idsq = (repmat(x',[1 n])-repmat(x,[n 1])).^2+ ...
                   (repmat(y',[1 n])-repmat(y,[n 1])).^2; % initial distance  
            o = idsq < 4*r^2 & triu(ones(n),1); 
        end 
    
        th=2*pi*rand(1,n); % initial headings
        dx=v*cos(th); % initial x velocity
        dy=v*sin(th); % initial y velocity
        
        for ix2 = 1:mot.numFrames
            % save x and y values
            % dim1 x,y; dim2 n dots; dim3 noFrames; dim4 noTrials
            % Goal: Screen('DrawDots', w, [x; y], 2*r, [255 255 255], [], 1); 
            trajectory(1,:,ix2,ix1)=x;
            trajectory(2,:,ix2,ix1)=y;
            
            % reflect off borders
            rx = x + dx > xmax | x + dx < xmin;
            dx(rx) = -dx(rx);
            ry = y + dy > ymax | y + dy < ymin;
            dy(ry) = -dy(ry);

            c=1;
            while any(c)
                dsq = (repmat(x',[1 n])-repmat(x,[n 1])).^2+ ...
                (repmat(y',[1 n])-repmat(y,[n 1])).^2; % old distance^2
                % fudge: temporarily disable collisions for overlapping objects
                dfudge = dsq < 4*r^2;
                dsq = (repmat((x+dx)',[1 n])-repmat(x+dx,[n 1])).^2+ ...
                (repmat((y+dy)',[1 n])-repmat(y+dy,[n 1])).^2; % new distance^2
                c = find(triu(ones(n)) & dsq < 4*r^2 & ~dfudge); % colliding pairs
                if any(c)
                    [i j]=ind2sub([n n],c);

                            thR=2*pi*rand(length(i),2);
                            dxR=v*cos(thR);
                            dyR=v*sin(thR);

                            dx([i j]) = dxR;
                            dy([i j]) = dyR;
                end
            end
 
            % linearly increase v during the inital frames
            if (ix2 < numFramesCh+1) & fadeIn == 1
                x=x+dx/v*vChange(ix2); % update positions
                y=y+dy/v*vChange(ix2);
            else
                x=x+dx; % update positions
                y=y+dy;
            end
            
            %check boundary crossing
            if1 = x<xmin;
            x(if1)=repmat(xmin,1,nnz(if1));
            if2 = x>xmax; 
            x(if2)=repmat(xmax,1,nnz(if2));
            if3 = y<ymin;
            y(if3)=repmat(ymin,1,nnz(if3));
            if4 = y>ymax; 
            y(if4)=repmat(ymax,1,nnz(if4));
            
        end
    end

    %%
    % save trajectory
    % do not use frames 1:100 because dots might overlap in these frames
    mot.trajectory = trajectory(:,:,1:end,:);
end

%% checker
function che = setChecker(windowPtr, gfx)
    rcycles = 8; % number of white/black circle pairs
    tcycles = 24; % number of white/black angular segment pairs (integer)
    flicker_freq = 4; % full cycle flicker frequency (Hz)
    flick_dur = 1/flicker_freq/2;
    che.numFramesPerFlick = round(flick_dur/gfx.timeFrame);
    
    hi_index=255;
    lo_index=0;
    bg_index =128;
    xysize = gfx.height;
    xylim = 2*pi*rcycles; % circumference of a r=8 unit circle
    [x,y] = meshgrid(-xylim:2*xylim/(xysize-1):xylim, ...
                     -xylim:2*xylim/(xysize-1):xylim);
    at = atan2(y,x);
    checks = ((1+sign(sin(at*tcycles)+eps) .* ...
                 sign(sin(sqrt(x.^2+y.^2))))/2) * (hi_index-lo_index) + lo_index;
             
    circle = x.^2 + y.^2 <= xylim^2; 
    checks = circle .* checks + bg_index * ~circle;
    % t(1) = Screen('MakeTexture', w, checks);
    % t(2) = Screen('MakeTexture', w, hi_index - checks); % reversed contrast
    
    % Change color of checker 
    indslow = find(checks==0); % 1: 80, 180 
    checks(indslow) = 50;      
    indshi = find(checks==255);
    checks(indshi) = 210;
    
    % ring checker
    ring = (x.^2 + y.^2 <= xylim^2) & (x.^2 + y.^2 >= xylim^2*(3.1/16) & ((x > x(1,round(xysize/2))+15) | (x < x(1,round(xysize/2))-15)) ); %3.1/16
    checksRing = ring .* checks + bg_index * ~ring;
    che.tR(1) = Screen('MakeTexture', windowPtr, checksRing);
    che.tR(2) = Screen('MakeTexture', windowPtr, hi_index - checksRing); % reversed contrast
    che.tL(1) = che.tR(1);
    che.tL(2) = che.tR(2);
    
%     % Right checker
%     moonR  = (x.^2 + y.^2 <= xylim^2) & (x.^2 + y.^2 >= xylim^2*(3.1/16) ) & (x > x(1,round(xysize/2)) );
%     checksMoonR = moonR .* checks + bg_index * ~moonR;
%     che.tR(1) = Screen('MakeTexture', windowPtr, checksMoonR);
%     che.tR(2) = Screen('MakeTexture', windowPtr, hi_index - checksMoonR); % reversed contrast
% 
%     % Right checker
%     moonL  = (x.^2 + y.^2 <= xylim^2) & (x.^2 + y.^2 >= xylim^2*(3.1/16) ) & (x < x(1,round(xysize/2)) );
%     checksMoonL = moonL .* checks + bg_index * ~moonL;
%     che.tL(1) = Screen('MakeTexture', windowPtr, checksMoonL);
%     che.tL(2) = Screen('MakeTexture', windowPtr, hi_index - checksMoonL); % reversed contrast

end

%% Instruction screen
function esc = showInstruction(esc, windowPtr, gfx, trl)
    Screen('TextSize', windowPtr , 20);
    
    % Text color as r g b - triplet
    if trl.sub_number == 99
        txtCol = [200 200 200];
    else
        txtCol = [0 255 255];
    end
    
    if trl.sub_number == 99
        strPractTrials = 'Practice\n';
    else
        strPractTrials = 'Experiment\n';
    end
    
    if 1
        strInfo = [
                    strPractTrials, ...
                    sprintf('Subject = %d\n\n', trl.sub_number), ...
                    'At the begnining of a trial 8 dots will be presented.\n', ...
                    '0, 2, 3, 4 or 5 dots will be blue. These are the targets (if any).\n', ...
                    'The remaining dots will be grey. These are distracters.\n\n', ...
                    'After a few seconds, the targets will turn grey and all dots will start to move.\n', ...
                    'Track the targets.\n', ...
                    'Keep tracking the assigned number of targets,', ...
                    'even when you are uncertain whether you are still tracking the orignial items. \n\n', ...
                    'At the end of a trial, all dots will stop to move and one dot will be higlighted in cyan.\n', ...
                    'If it was a target press ''Left Arrow'' on the number pad,\n', ...
                    'otherwise press ''Right Arrow''.\n\n', ...
                    'Press ''0'' to start.', ...
                ];
        DrawFormattedText(windowPtr, strInfo, 'center', 'center', txtCol, [], [], [], 3);
    end    
    Screen('Flip', windowPtr);
    
    % Wait for input
    WaitSecs(.5);
    fin = 0;
    zeroKey = KbName('0)'); F6 = KbName('F6'); space = KbName('SPACE');
    while ~fin & ~esc 
        FlushEvents(['mouseUp'],['mouseDown'],['keyDown'],['autoKey'],['update']);
        [secs, keyCode, deltaSecs] = KbWait; % KbWait waits for any key press, it waits infinite time % keyCode is a vector, its length is the number of keys on the keyboard
        esc = checkEscape(keyCode);

        if keyCode(zeroKey) | keyCode(F6) | keyCode(space)
            fin = 1;
        else
            fin = 0; 
        end
    end
end

function showFeedback(esc, windowPtr, gfx, trl, dat)
    Screen('TextSize', windowPtr , 20);
    idx = trl.num - trl.fbkT+1;
    corR = nnz(dat.blockdata(idx:end,8)==1);
    numR = length( dat.blockdata(idx:end,8) );
    perR = round((corR/numR)*100);    
    
    tF  = 6;
    numFramesF = ceil(tF/gfx.timeFrame);
    for it = 1:numFramesF
            DrawFormattedText(windowPtr, [num2str(perR) '% correct\n' '(' num2str(corR) '/' num2str(numR) ')'], 'center', gfx.cy, [0 255 255], [], [], [], 3);
            Screen('DrawingFinished', windowPtr);
            Screen('Flip', windowPtr);
    end
end


%% Check whether escape was pressed and exit if yes
function esc = checkEscape(keyCode)
    escapeKey = KbName('ESCAPE'); % returns key-code for ESCAPE
    if keyCode(escapeKey)         % check whether there was a 1 at the position 'escapeKey' in the vector keyCode 
        esc = 1;    
    else
        esc = 0;
    end
end

%% analyse response
function [dat] = analyseResp(trl, dat)
    numPad1 = KbName('LeftArrow');
    numPad2 = KbName('RightArrow');

    if dat.keyCode(numPad1)
        resp = 1; % match key
    elseif dat.keyCode(numPad2)
        resp = 2; % not match key
    elseif any(dat.keyCode)
        resp = 3;
    elseif ~any(dat.keyCode)
        resp = 4;
    end

    if trl.numTargets(trl.num)==0
        if resp == 4
            dat.IDcor = 1;
            dat.RTcor = 1;
        else
            dat.IDcor = 3;
            dat.RTcor = 3;
        end
    else
        % which key is the correct response in the current trial?
        if trl.probeMatch(trl.num) == 1 % match
            corResp = 1;
        elseif trl.probeMatch(trl.num) == 0 % not match
            corResp = 2;
        end

        % identiy press categories
        if resp == corResp % correct key
            dat.IDcor = 1;
        elseif resp ~= corResp & resp < 3 % opposite key
            dat.IDcor = 2;
        elseif resp == 3 % any other key
            dat.IDcor = 3;
        elseif resp == 4 % no key
            dat.IDcor = 4;
        end

        % time press cattegories
        if (dat.RT > trl.RTmin) & (dat.RT < trl.RTmax)
            dat.RTcor = 1;
        elseif (dat.RT <= trl.RTmin)
            dat.RTcor = 2;
        elseif (dat.RT >= trl.RTmax)
            dat.RTcor = 3;
        elseif (dat.RT > trl.RTmaxWait)
            dat.RTcor = 4;
        end
    end
 
end

