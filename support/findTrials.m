function [] = findTrials(filename)

%% check if file exists
% if exist([filename(1:end-6),'.mat'])
%     overwrite = questdlg([filename(1:end-6),'.mat already exists. Overwrite?'],'Existing file warning','Yes','No','No');
%     if ~strcmpi(overwrite,'Yes')
%         return
%     end
% end
%% define all my event codes

% these codes can be used as is. Presented in approx. the order they appear.
START = 1001;
PARADIGM = 32006;
PAUSECD = 1003; % I think?
SENDSTIM = 1060;

% these codes are added to the condition number 6xyy, where x are the
% following codes, yy is the stimulus number.
ARRAYON = 600;
FPOFF = 800;
SACINIT = 1000;
SACEND = 1400;

% back to standalone codes
REWON = 1030;

%% open file
try
    load(filename,'r_strobed','AD_heye','AD_veye')
catch
    try
        readall([filename(1:end-6),'.plx'])
        load(filename,'r_strobed','AD_heye','AD_veye')
    catch
        return
    end
end


%% dividing trial by events
AllEvts=r_strobed(r_strobed(:,2)==START,:);

j=1;
for i=1:length(AllEvts)-1, % for each trial
    trlcd_tmp=r_strobed(r_strobed(:,1)>AllEvts(i)-1 & r_strobed(:,1)<AllEvts(i+1),:); % relevant time stamps for current trial
    trlcd_tmp(trlcd_tmp(:,2)==PAUSECD,:)=[];
    if ~isempty(trlcd_tmp(trlcd_tmp(:,2)==REWON,1)), %was there a reward  
        % START
        try
            TrlCodes{j}=trlcd_tmp;
            TrlStart_ts(j,:)=trlcd_tmp(1,1); %trial start
            
            TARxy(j,:) = [trlcd_tmp(find(trlcd_tmp(:,2)==SENDSTIM)+1,2)-10000, trlcd_tmp(find(trlcd_tmp(:,2)==SENDSTIM)+2,2)-10000] / 10; % convert to degrees
            Cond(j,:) = trlcd_tmp(find(trlcd_tmp(:,2)==PARADIGM,1,'last')+1,2);
            if Cond(j,:) < 6000 || Cond(j,:) > 7000
                Cond(j,:) = [];
                continue
            else
                condCorr(j) = 6000;
            end
            
            % EVENTS
            Ton_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+ARRAYON), 1);
            FPoff_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+FPOFF), 1);
            Rewon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==REWON,1);
            SACon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+SACINIT),1);
            SACoff_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+SACEND),1);
            
            j=j+1;
        catch
        end
    end     
    clear trlcd_tmp; 
end

%% figure out the condition of each file
try
    numConds = 16;
    condition = rem(Cond-condCorr',numConds);
    condition(~condition)=numConds; % make condition 5/6 instead of "condition 0"
catch
    return
end

%% refining saccade times using eye velocity
SACCvel=[]; SACCon_ts=[];SACCoff_ts=[];

vel_crit = 50; % saccade if saccvel exceeds 50 degs / sec
max_vel = 1000;
eye_sample_rate = 1000; % this is based on sample rate of plexon for eye
nsamples = 3; % number of samples over which to estimate eye velocity
tm = nsamples/eye_sample_rate;

% nice function to get rid of HF noise
heye=smooth(AD_heye,50);
veye=smooth(AD_veye,50);

% compute saccade velocity across entire data
% NOTE for now blink removal not essential because trials with blink would not be rewarded
for i=2:length(AD_heye)-1
    displacement_deg = sqrt( ((heye(i-1) - heye(i+1)).^2) + ((veye(i-1) - veye(i+1)).^2) );
    SACCvel(i,:)=displacement_deg / tm;
end
SACCvel=[0; SACCvel];

for t=1:length(SACon_ts),     
    
    % finding saccade after reward
    y = 0; j = 0;
    while y<vel_crit && j<500
        y = SACCvel(SACon_ts(t)-j);
        j = j+1;
    end
    SACCon_ts = [SACCon_ts; SACon_ts(t)-j];
    
    % finding saccade off time
    y = vel_crit+1; j = 0;
    while y>vel_crit && j<200
        y = SACCvel(SACCon_ts(t)+j);
        j = j+1;
    end
    SACCoff_ts = [SACCoff_ts; SACCon_ts(t)+j];
end

cond1 = find(ismember(condition,[1,3,4,6,7,9,11,12,14,15]));
RF = repmat(TARxy(cond1(1),:),length(condition),1);

RFid = zeros(size(condition));
antiRFid = zeros(size(condition));

%% Hardcoding in the condition-target location relationship. Don't know how else to do this.
% Tar 1 in RF
condCheck = [1,3,7,9,11,15];
RFid(ismember(condition,condCheck)) = 1;

% Tar 1 opp RF
condCheck = [2,3,8,10,11,16];
antiRFid(ismember(condition,condCheck)) = 1;

% Tar 2 in RF
condCheck = [4,6,8,12,14,16];
RFid(ismember(condition,condCheck)) = 2;

% Tar 2 opp RF
condCheck = [5,6,7,13,14,15];
antiRFid(ismember(condition,condCheck)) = 2;

%% putting eye traces together
xpos = {};
ypos = {};
time = {};
dvect = {};
for n = 1:size(SACoff_ts,1)
    xpos = [xpos; AD_heye(TrlStart_ts(n):Rewon_ts(n))];
    ypos = [ypos; AD_veye(TrlStart_ts(n):Rewon_ts(n))];
    time = [time; TrlStart_ts(n):Rewon_ts(n)];
    vector = sqrt(AD_heye(TrlStart_ts(n):Rewon_ts(n)).^2+AD_veye(TrlStart_ts(n):Rewon_ts(n)).^2);
    dvect = [dvect; diff(vector)];
end

org_result = cell2struct([xpos,ypos,time,dvect,num2cell(condition),num2cell(RFid),num2cell(antiRFid),num2cell(FPoff_ts),num2cell(SACCon_ts),num2cell(SACCoff_ts),num2cell(Ton_ts), num2cell(RF)],...
    {'xpos','ypos','time','vectorvelocity','condition','RFid','antiRFid','FPoff','SacON','SacOFF', 'TarON','RFx','RFy'},2);

result = org_result;

trialind = 1;

save([filename(1:end-6),'.mat'],'org_result','result','trialind');


