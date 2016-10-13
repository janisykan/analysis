% in this version i'll try to see if the direction of the saccade after
% reward is in the general direction of RF.
% Also, this is the first version where i fixed the problem of low spks/s
% for individual blocks. It was 'cuz I included trials of another condition
% as "trials with no neural activity" in each condition.

function [convspks_delay_m, lfpOUT] = neural(file)

%% define all my event codes

% these codes can be used as is. Presented in approx. the order they appear.
START = 1001;
PARADIGM = 32006;
PAUSECD = 1003; % I think?
NSTIM = 196;
SENDSTIM = 1060;
CDEDIVIDER = 1062;
SENDSTIMEND = 1061;
SENDORI = 1070;
SENDORIEND = 1071;
SENDCOL = 1080;
SENDCOLEND = 1081;
WINDOPN = 800;
BCKLIT = 1027;

% these codes are added to the condition number 6xyy, where x are the
% following codes, yy is the stimulus number.
FPONSET = 200;
EYEINON = 400;
ARRAYON = 600;
FPOFF = 800;
SACEND = 1400;

% back to standalone codes
REWON = 1030;
REWOFF = 1040;

% add this to whatever basecode was there (ie. add -2000 to 6xyy to
% get 4xyy)
ABORT = -2000;

%% open file
channel = ['sig001',file.unit];
try
    load([file.path,file.name],'r_strobed',channel,'AD_heye','AD_veye','AD01')
catch
    try
        readall([file.path,file.name(1:end-6),'.plx'])
        load([file.path,file.name],'r_strobed',channel,'AD_heye','AD_veye','AD01')
    catch
        convspks_delay_m = 'corrupted';
        lfpOUT = 'corrupted';
        return
    end
end

AD_heye = -AD_heye;
AD_veye = -AD_veye;

%% dividing trial by events
AllEvts=r_strobed(r_strobed(:,2)==1001,:);% 1005

j=1;
for i=1:length(AllEvts)-1, % for each trial
    trlcd_tmp=r_strobed(r_strobed(:,1)>AllEvts(i)-1 & r_strobed(:,1)<AllEvts(i+1),:); % relevant time stamps for current trial
    trlcd_tmp(trlcd_tmp(:,2)==PAUSECD,:)=[];
    if ~isempty(trlcd_tmp(trlcd_tmp(:,2)==REWON,1)), %was there a reward  
        % START
        try
            TrlCodes{j}=trlcd_tmp;
            TrlStart_ts(j,:)=trlcd_tmp(1,1); %trial start
            
            RFxy(j,:) = [trlcd_tmp(find(trlcd_tmp(:,2)==SENDSTIM)+1,2)-10000, trlcd_tmp(find(trlcd_tmp(:,2)==SENDSTIM)+2,2)-10000] / 10; % convert to degrees
            Cond(j,:) = trlcd_tmp(find(trlcd_tmp(:,2)==PARADIGM,1,'last')+1,2);
            if Cond(j,:) < 6000 || Cond(j,:) > 7000
                Cond(j,:) = [];
                continue
            else
                condCorr(j) = 6000;
            end
            
            % EVENTS
            Ton_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==(Cond(j)+ARRAYON), 1);
            Rewon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==REWON,1);
            
            j=j+1;
        catch
        end
    end     
    clear trlcd_tmp; 
end

%% figure out the condition of each file
try
    numConds = length(file.condNames);
    condition = rem(Cond-condCorr',numConds);
    condition(~condition)=numConds; % make condition 5/6 instead of "condition 0"
catch
    convspks_delay_m = 'wrong format?';
    lfpOUT = 'wrong format?';
    return
end

for x = 1:numConds
    if sum(condition==x)<1
        convspks_delay_m = 'not enough';
        lfpOUT = 'not enough';
        keyboard
        return
    end
end

rewCond = condition;
if strcmpi(file.task,'prior04')
    if file.rew1 < file.rew2 % make rewCond relative to reward size, not color
        rewCond(condition==1) = 3;
        rewCond(condition==2) = 4;
        rewCond(condition==3) = 1;
        rewCond(condition==4) = 2;
        rewCond(condition==5) = 6;
        rewCond(condition==6) = 5;
        rewCond(condition==7) = 8;
        rewCond(condition==8) = 7;
    end
end

%% correct for vertical refresh (LCD 60 Hz, 16.67ms per refresh)
vertical_screen_deg = 50;
refresh_ms=16.67;
vertical_refresh_correction=((vertical_screen_deg - ((vertical_screen_deg / 2) - (RFxy(1,2) * .1))) / vertical_screen_deg) * refresh_ms; 
vertical_refresh_correction=ceil(vertical_refresh_correction);

%% plot align on target
xmin = -200;ymin=0;
xmax = 800; % time window to display spikes
sigma=5; nstd=sigma*5;
gauswin=normpdf(-nstd:nstd,0,sigma);
tt=xmin-nstd:xmax+nstd-1;
winsize = xmax - xmin;
spks_m=zeros(1,winsize);convspks_delay_m=[];
b=Ton_ts - vertical_refresh_correction;
for x=1:numConds
    eval(['c' num2str(x-1) ' =[];']);
    eval(['lfp' num2str(x-1) ' =[];']);
%     eval(['c' num2str(x-1) ' =zeros(1,winsize);']);
end
tick_st=0;

%% plotting rasters
figure
subplot(2,1,1)
hold on

% for finding saccade after reward
PUPvel=[]; SACCvel=[]; SACCon_ts=[];SACCoff_ts=[]; displacement_deg=[];

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

% filting lfp
d = designfilt('lowpassiir',...
    'PassbandFrequency',50,...
    'StopbandFrequency',55,...
    'DesignMethod','butter', ...         % Design method
    'MatchExactly','passband',...
    'SampleRate',1000);
% d = designfilt('bandpassfir',...
%         'FilterOrder',4,...
%         'StopbandFrequency1',1,...
%         'PassbandFrequency1',1.5,...
%         'PassbandFrequency2',50,...
%         'StopbandFrequency2',55,...
%        'DesignMethod','ls', ...         % Design method
%        'StopbandWeight1',1, ...         % Design method options
%        'PassbandWeight', 1, ...
%        'StopbandWeight2',1, ...
%        'SampleRate',1000);
lfpFilt = filtfilt(d,AD01);

for t=1:length(b),     
    
    % finding saccade after reward
    y = 0; j = 0;
    while y<vel_crit && j<500
        y = SACCvel(Rewon_ts(t)+j);
        j = j+1;
    end
    SACCon_ts = Rewon_ts(t)+j;
    
    % finding saccade off time
    y = vel_crit+1; j = 0;
    while y>vel_crit && j<200
        y = SACCvel(SACCon_ts+j);
        j = j+1;
    end
    SACCoff_ts = SACCon_ts+j;
    
    % find where saccade landed
    saccX = heye(SACCoff_ts);
    saccY = veye(SACCoff_ts);
    saccAngle = atan(saccY/saccX);
    fromRF = sqrt((saccX-RFxy(1,1))^2 + (saccY-RFxy(1,2))^2);
    
    % find angle of RF and window around it
    RFAngle = atan(RFxy(1,2)/RFxy(1,1));
    RFLength = sqrt(RFxy(1,1)^2+RFxy(1,2)^2);
    
    lowWinX = cos(RFAngle-deg2rad(45))*RFLength;
    lowWinY = sin(RFAngle-deg2rad(45))*RFLength;    
    highWinX = cos(RFAngle+deg2rad(45))*RFLength;
    highWinY = sin(RFAngle+deg2rad(45))*RFLength;
    
%     figure(1)
%     plot(-100:500,AD_heye(Rewon_ts(t)-100:Rewon_ts(t)+500));
%     hold on
%     plot(-100:500,AD_veye(Rewon_ts(t)-100:Rewon_ts(t)+500),'m');
%     plot([-100 500],[RFxy(1,1) RFxy(1,1)],'b')
%     plot([-100 500],[RFxy(1,2) RFxy(1,2)],'m')
%     plot(SACCoff_ts-Rewon_ts(t),saccX,'bo')
%     plot(SACCoff_ts-Rewon_ts(t),saccY,'mo')
%     title(['condition: ' num2str(rem(Cond(j)-condCorr(j),length(file.condNames)))])
%     hold off
%     
%     figure(2)
%     ezplot(@(x,y)myfun(x,y,RFxy(1,1),RFxy(1,2)),[RFxy(1,1)-8 RFxy(1,1)+8 RFxy(1,2)-8 RFxy(1,2)+8])
%     hold on
%     fill([0,highWinX,lowWinX,0],[0,highWinY,lowWinY,0],[1 0.6 0.6])
%     plot(AD_heye(Rewon_ts(t)-100:Rewon_ts(t)+500),AD_veye(Rewon_ts(t)-100:Rewon_ts(t)+500),'k')
%     plot(saccX,saccY,'rx');
%     hold off
%     axis([-30 30 -30 30])
%     axis square
    
%     if abs(saccAngle-RFAngle)<deg2rad(5)
%         continue
%     end
    
    eval(['spks = ' channel '(' channel '< b(t)+xmax & ' channel '> b(t)+xmin) - b(t)-xmin;']);
    tempspk =zeros(1,winsize);
    tempspk(1,spks)=1;
    eval(['c' num2str(rewCond(t)-1) '=[c' num2str(rewCond(t)-1) '; tempspk];']);
    
%     baseLFP = mean(lfpFilt(b(t)+(xmin-10):b(t)+(xmin+10)));
    baseLFP = mean(lfpFilt(b(t)-50:b(t)));
    lfp = lfpFilt(b(t)+xmin:b(t)+xmax-1)';
    lfp = lfp-baseLFP;
    
    eval(['lfp' num2str(rewCond(t)-1) '=[lfp' num2str(rewCond(t)-1) '; lfp];']);
%     d = designfilt('bandpassfir',...
%         'FilterOrder',2,...
%         'StopbandFrequency1',1,...
%         'PassbandFrequency1',1.5,...
%         'PassbandFrequency2',100,...
%         'StopbandFrequency2',100.5,...
%        'DesignMethod','ls', ...         % Design method
%        'StopbandWeight1',1, ...         % Design method options
%        'PassbandWeight', 1, ...
%        'StopbandWeight2',1, ...
%        'SampleRate',1000);
%     lfpFilt1 = filtfilt(d,lfp);
%     d = designfilt('bandstopfir',...
%         'FilterOrder',2,...
%         'StopbandFrequency1',55,...
%         'PassbandFrequency1',50,...
%         'PassbandFrequency2',65,...
%         'StopbandFrequency2',60,...
%        'DesignMethod','ls', ...         % Design method
%        'SampleRate',1000);
%     lfpFilt2 = filtfilt(d,lfpFilt1);
%     
%     figure
%     plot(lfp)
%     hold on
%     plot(lfpFilt)
%     hold off


end

for x=1:numConds
    eval(['curCond = c' num2str(x-1) ';']);
    if isempty(curCond)
        convspks_delay_m = 'not enough';
        lfpOUT = 'not enough';
        return
    end
end

for x=1:numConds
    eval(['curCond = c' num2str(x-1) ';']);
    for y=1:size(curCond,1)
        tick_st=tick_st+4;
        spks = find(curCond(y,:));
        if ~isempty(spks),
            for m=1:size(spks,2)
                h=plot([spks(m)+xmin spks(m)+xmin], [tick_st tick_st+2],file.colors{x});
                set(h,'linewidth',1.5)
            end
            clear h ;
        end
    end
end


%% plotting spike density function
subplot(2,1,2)
hold on
clear spks_m
for x=1:numConds
    eval(['spks_m(1,:)=mean(c' num2str(x-1) ',1);']);
    % convolve with gaussian
    convspks_delay_m(x,:)=conv(spks_m(1,:)*1000, gauswin);
    
    % spike density function
    h(x)=plot(tt,convspks_delay_m(x,:),file.colors{x});
    clear spks_m
    
    eval(['lfpOUT(x,:)=mean(lfp' num2str(x-1) ',1);']);
end


ymax1=ceil(max(max(convspks_delay_m,[],2)))+10;

if numConds == 5
    extraCond = nan(1,length(convspks_delay_m));
    convspks_delay_m = [extraCond;convspks_delay_m];
end

xlabel('time from target (ms) ')
ylabel('spikes/s ')
title([file.name(1:end-6) ' depth' num2str(file.depth) ', layer = ' file.layer ', x = ' num2str(RFxy(1,1)) ' , y = ' num2str(RFxy(1,2)) ', isolation: ' file.isolation]) 
legend(h,file.condNames)
set(gca,'XTick', xmin:100:xmax);
axis([xmin xmax 0 200])
% keyboard
close all


end

function z = myfun(x,y,RFx,RFy)
z = (x-RFx)^2/5^2+(y-RFy)^2/5^2-1;
end