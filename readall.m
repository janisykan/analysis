% this script reads all the spike timestamp and a/d info from a plx file into matlab
% variables.

function [] = readall(StartingFileName)

% Open a plx file
% this will bring up the file-open dialog
% StartingFileName = '';
%StartingFileName = 'c:\plexondata\NSSample.plx';
%StartingFileName = 'c:\plexondata\au101602a.plx';

if nargin<1
    StartingFileName = 'C:\Users\Janis\Documents\Munoz Lab\Data\indy\behavior\ic19c.plx';
end
[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(StartingFileName);

disp(['Opened File Name: ' OpenedFileName]);
disp(['Version: ' num2str(Version)]);
disp(['Frequency : ' num2str(Freq)]);
disp(['Comment : ' Comment]);
disp(['Date/Time : ' DateTime]);
disp(['Duration : ' num2str(Duration)]);
disp(['Num Pts Per Wave : ' num2str(NPW)]);
disp(['Num Pts Pre-Threshold : ' num2str(PreThresh)]);
% some of the information is only filled if the plx file version is >102
if ( Version > 102 )
    if ( Trodalness < 2 )
        disp('Data type : Single Electrode');
    elseif ( Trodalness == 2 )
        disp('Data type : Stereotrode');
    elseif ( Trodalness == 4 )
        disp('Data type : Tetrode');
    else
        disp('Data type : Unknown');
    end
        
    disp(['Spike Peak Voltage (mV) : ' num2str(SpikePeakV)]);
    disp(['Spike A/D Resolution (bits) : ' num2str(SpikeADResBits)]);
    disp(['Slow A/D Peak Voltage (mV) : ' num2str(SlowPeakV)]);
    disp(['Slow A/D Resolution (bits) : ' num2str(SlowADResBits)]);
end   

% get some counts
[tscounts, wfcounts, evcounts, slowcounts] = plx_info(OpenedFileName,1);

% tscounts, wfcounts are indexed by (unit+1,channel+1)
% tscounts(:,ch+1) is the per-unit counts for channel ch
% sum( tscounts(:,ch+1) ) is the total wfs for channel ch (all units)
% [nunits, nchannels] = size( tscounts )
% To get number of nonzero units/channels, use nnz() function

% gives actual number of units (including unsorted) and actual number of
% channels plus 1
[nunits1, nchannels1] = size( tscounts );   

% we will read in the timestamps of all units,channels into a two-dim cell
% array named allts, with each cell containing the timestamps for a unit,channel.
% Note that allts second dim is indexed by the 1-based channel number.
% preallocate for speed
allts = cell(nunits1, nchannels1);
for iunit = 0:nunits1-1   % starting with unit 0 (unsorted) 
    for ich = 1:nchannels1-1
        if ( tscounts( iunit+1 , ich+1 ) > 0 )
            % get the timestamps for this channel and unit 
            [nts, allts{iunit+1,ich}] = plx_ts(OpenedFileName, ich , iunit );
         end
    end
end

             
% get some other info about the spike channels
[nspk,spk_filters] = plx_chan_filters(OpenedFileName);
[nspk,spk_gains] = plx_chan_gains(OpenedFileName);
[nspk,spk_threshs] = plx_chan_thresholds(OpenedFileName);
[nspk,spk_names] = plx_chan_names(OpenedFileName);

%%
% added the following %
sigval = ['iabcdef'];
[rows,cols] = size(tscounts);
valid_units = [];
unsorted_units = [];
for x = 1:cols
    for y = 1:rows
        if tscounts(y,x) > 0
            eval([ spk_names((x-1),:) , sigval(y) '= round(allts{' num2str(y) ',' num2str(x-1), '} * 1000);']);      
            if sigval(y) ~= 'i'
                valid_units = [valid_units; [spk_names((x-1),:) , sigval(y)]];
            end
            if sigval(y) == 'i'
                unsorted_units=[unsorted_units; [spk_names((x-1),:) , sigval(y)]];
            end
        end
    end
end
valid_units
unsorted_units


%%
% get the a/d data into a cell array also.
% This is complicated by channel numbering.
% The number of samples for analog channel 0 is stored at slowcounts(1).
% Note that analog ch numbering starts at 0, not 1 in the data, but the
% 'allad' cell array is indexed by ich+1
numads = 0;
[u,nslowchannels] = size( slowcounts );   
if ( nslowchannels > 0 )
    % preallocate for speed
    allad = cell(1,nslowchannels);
    for ich = 0:nslowchannels-1
        if ( slowcounts(ich+1) > 0 )
			[adfreq, nad, tsad, fnad, allad{ich+1}] = plx_ad(OpenedFileName, ich);
            allad{ich+1}=[zeros(round(tsad*1000),1);allad{ich+1}];
			numads = numads + 1;
        end
    end
end

if ( numads > 0 )
    [nad,adfreqs] = plx_adchan_freqs(OpenedFileName);
    [nad,adgains] = plx_adchan_gains(OpenedFileName);
    [nad,adnames] = plx_adchan_names(OpenedFileName);

    % just for fun, plot the channels with a/d data
    iplot = 1;
    numPlots = 64;%min(nad, numads);
    for ich = 1:64;%nslowchannels
        [ nsamples, u ] = size(allad{ich});
       % if ( nsamples > 0 )
          %  subplot(numPlots/8,8,iplot); plot(allad{ich}); xlabel(adnames(ich,:));
            iplot = iplot + 1;
      %  end
        if iplot > numPlots
           break;
        end
    end
end

for i=1:nad,
    if (~isempty(allad{i}))
        if(~isempty(isspace(adnames(i,:))))
            adnames(i,isspace(adnames(i,:)))='_';
        end
        eval([adnames(i,1:length(deblank(adnames(i,:)))) '=allad{i};'])
    end
end


%%
% and finally the events
[u,nevchannels] = size( evcounts );  
if ( nevchannels > 0 ) 
    % need the event chanmap to make any sense of these
    [u,evchans] = plx_event_chanmap(OpenedFileName);
	for iev = 1:nevchannels
		if ( evcounts(iev) > 0 )
            evch = evchans(iev);
            if ( evch == 257 )
				[nevs{iev}, tsevs{iev}, svStrobed] = plx_event_ts(OpenedFileName, evch); 
			else
				[nevs{iev}, tsevs{iev}, svdummy] = plx_event_ts(OpenedFileName, evch);
            end
		end
	end
end
[nev,evnames] = plx_event_names(OpenedFileName);

for i=1:size(nevs,2),
    if (~isempty(tsevs{i}))
        if(~isempty(isspace(evnames(i,:))))
            evnames(i,isspace(evnames(i,:)))='_';
        end
        eval([evnames(i,1:length(deblank(evnames(i,:)))) '=round(tsevs{i} * 1000);'])
    end
end
r_strobed=[Strobed svStrobed];

p=11 %input('Photo diode channel: ');
x=13 %input('Horizontal eye channel: ');
y=14 %input('Vertical eye channel: ');

photo_diode_analog=allad{p};
conversion_factor = 0.025;% may need tweeking (check correspondence btw eye and cal points)
AD_heye=allad{x} * conversion_factor;
AD_veye=allad{y} * conversion_factor;

save([OpenedFileName(1:end-4) '-r.mat']);

% close
% clear

