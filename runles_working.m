% runles
% script for delay task only 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = runles_working()

clear all
close all

%% read
% fname = input(['filename: '],'s');    
[file,path]=uigetfile('C:\Users\Janis\Documents\Munoz Lab\Data\indy\delayed\','multiselect','on');
oldpath = cd;
cd(path)% the place I put the original .mat files after converted from plexon file
sheet = 1;

if ~iscell(file)
    [~,~,ext] = fileparts(file);
    if strcmp(ext,'.xlsx')
        [~,~,raw] = xlsread([path,file],sheet);
        relCell = cellfun(@findstr,repmat({'delayed'},[size(raw,1),1]),raw(:,2),'UniformOutput',false);
        relIndex = ~cellfun(@isempty,relCell);
%         relIndex = find(strcmp(raw(:,2),task));
        relevant = raw(relIndex,:);
        filenames = relevant(:,1);
        filenames = strcat(filenames,'-r.mat');
        unit = relevant(:,3);
        depth = relevant(:,5);
        tarx = relevant(:,6);
        tary = relevant(:,7);
        visual = nan(size(raw,1),1);
        motor = nan(size(raw,1),1);
        index = find(relIndex);
    else
        filenames = {file};
        unit = input('unit? ','s');
        unit = {unit};
    end
end

for n = 1:size(filenames,1) % for each file
    if depth{n} == -1000
        continue
    end
    fname = filenames{n};
    disp(['loading file ' fname])
    try
        load([path,fname]);
    catch
        try
            readall([path, fname(1:end-6),'.plx'])
            load([path, fname])
        catch
            disp([fname,' corrupted'])
            continue
        end
    end
    
    Comment % the Comment you wrote in the plexon file with important info
    disp('wait...')
    cd ..
    
    %% define
    %target distractor locations and codes
    % disp('RF locations in 10th of deg')
    % rfx = input('RF x?  ');
    % rfy = input('RF y?  ');
    %
    % eccentricity = (sqrt(abs(rfx^2) + abs(rfy^2)));
    
    %which unit we want to look at
    channel = ['sig001',unit{n}];
    
    % define basecodes for each task
    %task=input('task: (2)memory guided (4)gap (5)delay  ');
    task=5;
    
    if task==2
        Fon_basecode=[6222];
        Eyein_basecode = Fon_basecode+200;
        Ton_basecode = Eyein_basecode+200;
        Toff_basecode = Ton_basecode+200;
        Foff_basecode = Toff_basecode+200;
        Son_basecode = Foff_basecode+200;
        Eyein2_basecode = Son_basecode+200;
        Soff_basecode = Eyein2_basecode+200;
    end
    
    if task==4
        Fon_basecode=[6212];
        Eyein_basecode = Fon_basecode+200;
        Foff_basecode = Eyein_basecode+200;
        Ton_basecode = Foff_basecode+200;
        Son_basecode = Ton_basecode+200;
        Eyein2_basecode = Son_basecode+200;
        Soff_basecode = Eyein2_basecode+200;
    end
    
    if task==5 % delay
        Fon_basecode=[6252];
        Eyein_basecode = Fon_basecode+200;% 64xx
        Ton_basecode = Eyein_basecode+200;% 66xx
        Foff_basecode = Ton_basecode+200;% 68xx
        Son_basecode = Foff_basecode+200;% 70xx
        Eyein2_basecode = Son_basecode+200;% 72xx
        Soff_basecode = Eyein2_basecode+200;% 74xx
    end
    
    tnext_code=1060;
    cnext_code=1061;
    rewon=1030;
    rewend=1040;
    
    
    %% organize main events (set up for delay task only)
    
    AllEvts=r_strobed(r_strobed(:,2)==1001,:);% 1005
    j=1;
    clear RFxy Fixoff_ts Ton_ts Saccon_ts Rewon_ts TrlStart_ts
    for i=1:length(AllEvts)-1,
        trlcd_tmp=r_strobed(r_strobed(:,1)>AllEvts(i)-1 & r_strobed(:,1)<AllEvts(i+1),:); % relevant time stamps for current trial
        if ~isempty(trlcd_tmp(trlcd_tmp(:,2)==rewon,1)), %was there a reward
            
            basecode = trlcd_tmp(find(trlcd_tmp(:,2)==tnext_code)-1,2);
            if basecode~=Fon_basecode-200;
                continue
            end
            
            % START
            try
%                 TrlCodes{j}=trlcd_tmp;
                
                RFxy(j,:) = ([trlcd_tmp(find(trlcd_tmp(:,2)==tnext_code)+1,2), trlcd_tmp(find(trlcd_tmp(:,2)==tnext_code)+2,2)]-10000) / 10; % convert to degrees
%                 Tindex(j,:) = trlcd_tmp(find(trlcd_tmp(:,2)==tnext_code)+3,2);
                
                % EVENTS
%                 Fixon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Fon_basecode) & trlcd_tmp(:,2)<=max(Fon_basecode), 1);
%                 Fixon_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Fon_basecode) & trlcd_tmp(:,2)<=max(Fon_basecode), 2);
%                 Eyein_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Eyein_basecode) & trlcd_tmp(:,2)<=max(Eyein_basecode), 1);
%                 Eyein_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Eyein_basecode) & trlcd_tmp(:,2)<=max(Eyein_basecode), 2);
                Fixoff_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Foff_basecode) & trlcd_tmp(:,2)<=max(Foff_basecode), 1);
%                 Fixoff_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Foff_basecode) & trlcd_tmp(:,2)<=max(Foff_basecode), 2);
                Ton_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Ton_basecode) & trlcd_tmp(:,2)<=max(Ton_basecode), 1);
%                 Ton_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Ton_basecode) & trlcd_tmp(:,2)<=max(Ton_basecode), 2);
                Saccon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Son_basecode) & trlcd_tmp(:,2)<=max(Son_basecode), 1);
%                 Saccon_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Son_basecode) & trlcd_tmp(:,2)<=max(Son_basecode), 2);
%                 Saccoff_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Soff_basecode) & trlcd_tmp(:,2)<=max(Soff_basecode), 1);
%                 Saccoff_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)>=min(Soff_basecode) & trlcd_tmp(:,2)<=max(Soff_basecode), 2);
                Rewon_ts(j,:) = trlcd_tmp(trlcd_tmp(:,2)==rewon,1);
%                 Rewon_codes(j,:) = trlcd_tmp(trlcd_tmp(:,2)==rewon,2);
                TrlStart_ts(j,:)=trlcd_tmp(1,1); %trial start
                clear trlcd_tmp;
                j=j+1;
            catch
            end
        end
    end
    
    try
        Ntrl=length(TrlStart_ts);
    catch
        range = ['N' num2str(index(n))];
        msg = {'corrupted'};
        xlswrite([path file],msg,sheet,range)
        continue
    end
    
    % define some variables to come
    Targ_on=zeros(Ntrl,1);
    SACCon=zeros(Ntrl,1);
    SACCoff=zeros(Ntrl,1);
    SRT=zeros(Ntrl,1);
    SACCdur=zeros(Ntrl,1);
    SACCpvel=zeros(Ntrl,1);
    SACCamp=zeros(Ntrl,1);
    
    %% correct for vertical refresh (LCD 60 Hz, 16.67ms per refresh)
    vertical_screen_deg = 50;
    refresh_ms=16.67;
    vertical_refresh_correction=((vertical_screen_deg - ((vertical_screen_deg / 2) - (RFxy(1,2) * .1))) / vertical_screen_deg) * refresh_ms;
    vertical_refresh_correction=ceil(vertical_refresh_correction);
    
    % define variables needed for plotting neural
    sigma=5; nstd=sigma*5;
    gauswin=normpdf(-nstd:nstd,0,sigma);
    xminT = -100;
    xmaxT = 1000; % time window to display spikes
    ttT=xminT-nstd:xmaxT+nstd-1;
    winsizeT = xmaxT - xminT;
    
    xminS = -900;
    xmaxS = 200; % time window to display spikes
    ttS=xminS-nstd:xmaxS+nstd-1;
    winsizeS = xmaxS - xminS;
    
    spks_m=zeros(1,winsizeT);convspks_delayT_m=[];
    Nrast=0;
    b=Ton_ts - vertical_refresh_correction;
    cTon=zeros(1,winsizeT); % target oneset
    cSon=zeros(1,winsizeS); % saccade onset
    tick_st=0;
    
    %% auto mark
    PUPvel=[]; SACCvel=[]; SACCon_ts=[];SACCoff_ts=[]; displacement_deg=[];
    
    vel_crit = 50; % saccade if saccvel exceeds 50 degs / sec
    max_vel=1000;
    eye_sample_rate=1000;% this is based on sample rate of plexon for eye, even though EL desktop version at 500Hz
    nsamples=3;% number samples over which to estimate eye velocity is 3
    tm = nsamples / eye_sample_rate;
    
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
    
    automark=0;
    if automark, h = waitbar(0,'auto saccade marking...'); end
    
    for trl=1:Ntrl % for each trial
        if automark, waitbar(trl/Ntrl,h); end
        
        %%%%%% Edit saccade time stamps %%%%%%
        %%%% look from -30ms relative to current saccade onset to find more precise onset
        wd=30;
        j=1;y=0;
        while y < vel_crit
            y=SACCvel(Saccon_ts(trl)-wd+j);
            j=j+1;
        end
        SACCon_ts(trl)=Saccon_ts(trl)-wd+j;
        
        %%%% look from current saccade onset to find more precise saccade end
        j=1;y=vel_crit+1;
        while y > vel_crit
            y=SACCvel(SACCon_ts(trl)+j);
            j=j+1;
        end
        SACCoff_ts(trl)=SACCon_ts(trl)+j;
        
        SRT(trl,:)=SACCon_ts(trl)-Fixoff_ts(trl);
        SACCdur(trl,:)=SACCoff_ts(trl)-SACCon_ts(trl);
        SACCpvel(trl,:) = max(SACCvel(SACCon_ts(trl):SACCoff_ts(trl)));
        SACCamp(trl,:)=sqrt(  (heye(SACCon_ts(trl)) - heye(SACCoff_ts(trl))).^2   +  (veye(SACCon_ts(trl)) - veye(SACCoff_ts(trl))).^2  );
        
        %     TrlCodes{trl}
        
        %%% plot
%         if ~automark % needs a GUI to allow manual saccade marking
%             figure(1)
%             hold off
%             plot(heye(Ton_ts(trl):Rewon_ts(trl)));
%             hold on
%             plot(veye(Ton_ts(trl):Rewon_ts(trl)),'m');
%             plot([Fixoff_ts(trl)-Ton_ts(trl) Fixoff_ts(trl)-Ton_ts(trl)],[-20 +20],'g');
%             plot([SACCon_ts(trl)-Ton_ts(trl) SACCon_ts(trl)-Ton_ts(trl)],[-20 +20],'r');
%             plot([SACCoff_ts(trl)-Ton_ts(trl) SACCoff_ts(trl)-Ton_ts(trl)],[-20 +20],'r:');
%             plot([Rewon_ts(trl)-Ton_ts(trl) Rewon_ts(trl)-Ton_ts(trl)],[-20 +20],'y');
%             plot([0 Rewon_ts(trl)-Ton_ts(trl)], [tarx{n,1}/10 tarx{n,1}/10], 'k')
%             plot([0 Rewon_ts(trl)-Ton_ts(trl)], [tary{n,1}/10 tary{n,1}/10], 'k')
%             xlabel('Time from target onset (ms) ')
%             ylabel('Eye position (deg)')
%             h=legend('H eye', 'V eye', 'Fix off', 'Sacc on', 'Sacc end', 'Reward' ,'location', 'SouthWest');
%             %         input('');
%             
%             plot(heye(Ton_ts(trl):Rewon_ts(trl)),veye(Ton_ts(trl):Rewon_ts(trl)))
%             hold on
%             hO = ezplot(@(x,y)myfun(x,y,RFxy(1,1),RFxy(1,2)),[RFxy(1,1)-8 RFxy(1,1)+8 RFxy(1,2)-8 RFxy(1,2)+8]);
%             set(hO,'LineColor','r')
%             hO = ezplot(@(x,y)myfun(x,y,-RFxy(1,1),-RFxy(1,2)),[-RFxy(1,1)-8 -RFxy(1,1)+8 -RFxy(1,2)-8 -RFxy(1,2)+8]);
%             set(hO,'LineColor','g')
%         end
        
        %%%% neural stuff %%%%
        figure(2)
        subplot(2,1,1)
        hold on
        eval(['spks = ' channel '(' channel '< b(trl)+xmaxT & ' channel '> b(trl)+xminT) - b(trl)-xminT;']); % aligned on target onset
        cTon(trl,spks)=1;
        if ~isempty(spks),
            for m=1:size(spks,1)
                h=plot([spks(m)+xminT spks(m)+xminT], [tick_st tick_st+2],'b');
                set(h,'linewidth',1.5)
            end
            clear h ;
        end
        clear spks
        eval(['spks = ' channel '(' channel '< SACCon_ts(trl)+xmaxS & ' channel '> SACCon_ts(trl)+xminS) - SACCon_ts(trl)-xminS;']); % aligned on target onset
        cSon(trl,spks)=1;
        tick_st=tick_st+4;
        subplot(2,1,2)
        hold on
        if ~isempty(spks),
            for m=1:size(spks,1)
                h=plot([spks(m)+xminS spks(m)+xminS], [tick_st tick_st+2],'b');
                set(h,'linewidth',1.5)
            end
            clear h ;
        end
        
    end
    if automark, close(h); end
    
    figure(2)
    subplot(2,1,1)
    spks_m(1,:)=mean(cTon);
    convspks_delayT_m(1,:)=conv(spks_m(1,:)*1000, gauswin);
    h=plot(ttT,convspks_delayT_m(1,:));
    set(h,'color','k')
    clear spks_m h
    hold off
    
    subplot(2,1,2)
    spks_m(1,:)=mean(cSon);
    convspks_delayS_m(1,:)=conv(spks_m(1,:)*1000, gauswin);
    h=plot(ttS,convspks_delayS_m(1,:));
    set(h,'color','k')
    clear spks_m h
    hold off
    
    %% classifying neuron **
    baseWin = [-70 30];
    baseline = convspks_delayT_m(:,baseWin(1)-xminT+nstd:baseWin(2)-xminT+nstd);
    baseSD = std(baseline);
    
    visWin = [0 150];
    sigCount = 0;
    visual(index(n)) = 0;
    sigSpkV = nan;
    for compSpk = convspks_delayT_m(visWin(1)-xminT+nstd:visWin(2)-xminT+nstd)
        if compSpk >(mean(baseline)+3*baseSD)
            sigCount = sigCount+1;
            if sigCount >= 40
                visual(index(n)) = 1;
                break
            elseif sigCount == 1
                sigSpkV = compSpk;
            end
        else
            sigCount = 0;
            sigSpkV = nan;
        end
    end
    
    presacWin = [-150 -50];
    presac = convspks_delayS_m(:,presacWin(1)-xminS+nstd:presacWin(2)-xminS+nstd);
    presacSD = std(presac);
    
    sacWin = [-25 25];
    sigCount = 0;
    motor(index(n)) = 0;
    sigSpkM = nan;
    for compSpk = convspks_delayS_m(sacWin(1)-xminS+nstd:sacWin(2)-xminS+nstd)
        if compSpk >(mean(presac)+3*presacSD)
            sigCount = sigCount+1;
            if sigCount >= 40
                motor(index(n)) = 1;
                break
            elseif sigCount == 1
                sigSpkM = compSpk;
            end
        else
            sigCount = 0;
            sigSpkM = nan;
        end
    end
%     keyboard
    
end

visual = visual(2:end,1);
motor = motor(2:end,1);
xlswrite([path file],visual,sheet,'L2')
xlswrite([path file],motor,sheet,'M2')

%%

% %%%%% save %%%%%
% h = waitbar(1,['saving file ' fname]);  
% save([fname(1:end-2) '-m']) % appen 'm' for marked file
% close(h)

cd(oldpath)
end

function z = myfun(x,y,RFx,RFy)
z = (x-RFx)^2/5^2+(y-RFy)^2/5^2-1;
end



