function [eq, ffName] = binder(cPick, noffset, coffset, writePick)
%--------------------------------------------------------------------------
% binder.m
% Usage:
%   [eq, ffName] = binder(cPick, noffset, coffset, writePick);
% Function to take picks from cdataPlot and bind picks to
% certain events.  Writes UW-2 style pickfiles for input into spong. 
% Inputs:
%   cPick: Output from cdataPick.  A structure containing fields:
%           .pickTime
%           .station
%           .channel
%           .phase
%           .polarity
%           .uncertainty
%   noffset: number of seconds before or after the pick to consider part of
%            the same event (rule of thumb: radius of network/3 km/s)
%   coffset: number of seconds after the pick to consider a coda pick
%   <writePick>: (optional)boolean for writing pickfiles, 1 means write a UW-2
%              pickfile, 0 means do not write a pickfile (default is 1)
% Output:
%   eq: structure containing fields
%       type: type of pick
%       pick: time vector of pick time
%       unc: uncertainty of pick
%       pol: polarity of pick
%       sta: station of pick
%       chan: channel of pick
%       npick: number of picks
%   ffName: name of output file
%
% Written by WT on day 6 of waiting in Esso for a helicopter 8/12/07
% Modified by WT 10/10/07
% Modified by WT 10/10/08: Fixed case of one P or S pick, added option for
%                          not writing pickfiles
% Modified by WT 10/20/08: Fixed coda and detail cards to look for station
%                          instead of assuming first pick
%--------------------------------------------------------------------------

% Testing variables--------------------------------------------------------
% clear all; close all;
% !/bin/rm 07*s
test = 0;
% load 0710010855.pick.mat
% noffset = 3;
% coffset = 60;
%--------------------------------------------------------------------------
if nargin == 3
    writePick = 1;
end

if isfield(cPick, 'network')
    sncl = 1;
else
    sncl = 0;
end

cntr = 1;
%--------------------------------------------------------------------------
% Parse out cPick structure to get only P and S waves----------------------
IND = 1;
for i = 1 : length(cPick)
    if strcmp(cPick(i).phase, 'P') || strcmp(cPick(i).phase, 'S')
        type(IND) = cPick(i).phase;
        picks{IND} = cPick(i).pickTime;
        if isempty(cPick(i).uncertainty)
            unc(IND) = 0.05;
        else
            unc(IND) = cPick(i).uncertainty;
        end
        pol{IND} = cPick(i).polarity;
        sta{IND} = cPick(i).station;
        chan{IND} = cPick(i).channel;
        if sncl == 1
            net{IND} = cPick(i).network;
            loc{IND} = cPick(i).location;
        end
        IND = IND+1;
    end
end
     
numpicks = length(picks);
%--------------------------------------------------------------------------
% Bind P and S waves ------------------------------------------------------
while length(picks) > 1
    iind = 1;
    % Calculate the forward and backward times from a given pick
    for i = 1 : length(picks)
        if isempty(picks{i})
            eqDiff(i) = 9999999;
        elseif isempty(picks{1})
            eqDiff(i) = 999999;
        else
            eqDiff(i) = timediff(picks{1}, picks{i});
        end
    end
    % Find anything within noffset seconds
    rr = find (abs(eqDiff) < noffset);
 
    for i = 1 : length(rr)
        eq(cntr).type(i) = type(rr(i));
        eq(cntr).pick{i} = picks{rr(i)};
        eq(cntr).sta{i} = sta{rr(i)};
        eq(cntr).chan{i} = chan{rr(i)};
        if sncl == 1
            eq(cntr).net{i} = net{rr(i)};
            eq(cntr).loc{i} = loc{rr(i)};
        end
        eq(cntr).unc(i) = unc(rr(i));
        eq(cntr).pol(i) = pol(rr(i));
    end
    kk = find( abs(eqDiff) > noffset );
    kk2 = find( eqDiff == 999999 );
    ind1 = setxor(kk, kk2);
    kk3 = find( eqDiff == 9999999 );
    ind2 = setxor(ind1, kk3);
    
    kk = ind2;
    type = type(kk);
    picks = picks(kk);
    sta = sta(kk);
    chan = chan(kk);
    if sncl == 1
        net = net(kk);
        loc = loc(kk);
    end
    unc = unc(kk);
    pol = pol(kk);
    if test == 1
        D = Dpick;
        D = coralDemean(D);
        D = coralTaper(D);
        D = coralFilter(D, [1 20], 'bandpass', 2, 'minimum');
        figure (1); clf;
        %opt.y_offset = [1:length(eq(cntr).pick)];
        coralPlot(D);
        hold on;
        for i = 1 : length(eq(cntr).pick)
            for j = 1 : length(D)
                if strcmp(eq(cntr).sta(i), D(j).staCode)
                    ind = j;
                    break;
                end
            end
            pik = timediff(eq(cntr).pick{i}, D(i).recStartTime);
            plot(pik, -(ind-1), 'or');
        end
        pause;
    end     
    cntr = cntr + 1;
    iind = iind + 1;
    if iind > length(cPick)
        disp('Lots of iterations..., are you sure things are cool?');
        keyboard;
    end
    clear eqDiff D
end
if length(picks) == 1
    eq(cntr).type = type;
    eq(cntr).pick = picks;
    eq(cntr).sta = sta;
    eq(cntr).chan = chan;
    if sncl == 1
        eq(cntr).net = net;
        eq(cntr).loc = loc;
    end
    eq(cntr).unc = unc;
    eq(cntr).pol = pol;
    cntr = cntr + 1;
end
clear type picks sta chan unc pol
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Bind coda wave cards
% Parse out cPick structure
IND = 1;
for i = 1 : length(cPick)
    if strcmp(cPick(i).phase, 'C') 
        type(IND) = cPick(i).phase;
        picks{IND} = cPick(i).pickTime;
        unc(IND) = NaN;
        pol{IND} = cPick(i).polarity;
        sta{IND} = cPick(i).station;
        chan{IND} = cPick(i).channel;
        if sncl == 1
            net{IND} = cPick(i).network;
            loc{IND} = cPick(i).location;
        end
        IND = IND+1;
    end
end

% Calculate the forward and backward times from a given pick for a given eq
if exist('picks', 'var') && ~isempty(picks)
    for i = length(eq): -1 : 1
        for j = 1 : length(picks)
            % Try to find common station
            n = find(strcmp(sta{j},eq(i).sta) == 1);
            if length(n) == 0
                n = 1;
            elseif length(n) > 1
                n = n(1);
            end
            eqDiff(j) = timediff(picks{j}, eq(i).pick{n});
        end
        % Find anything within noffset seconds after initial pick
        rk = find (eqDiff >= 0);
        rr = find (eqDiff < coffset);
        cntr = length(eq(i).sta)+1;
        for ii = 1 : length(rr)
            if find(rr(ii) == rk)
                eq(i).type(cntr) = type(rr(ii));
                eq(i).pick{cntr} = picks{rr(ii)};
                eq(i).sta{cntr} = sta{rr(ii)};
                eq(i).chan{cntr} = chan{rr(ii)};
                if sncl == 1
                    eq(i).net{cntr} = net{rr(ii)};
                    eq(i).loc{cntr} = loc{rr(ii)};
                end
                eq(i).unc(cntr) = unc(rr(ii));
                eq(i).pol(cntr) = pol(rr(ii));
                cntr = cntr + 1;
            end
        end
        kk = find( eqDiff < 0 );
        type = type(kk);
        picks = picks(kk);
        sta = sta(kk);
        chan = chan(kk);
        if sncl == 1
            net = net(kk);
            loc = loc(kk);
        end
        unc = unc(kk);
        pol = pol(kk);
        if isempty(picks)
            break;
        end
        clear eqDiff;
    end
end
clear type picks sta chan unc pol
%--------------------------------------------------------------------------
% Bind comment cards
% Parse out cPick structure
IND = 1;
for i = 1 : length(cPick)
    if strcmp(cPick(i).phase, 'D')
        type(IND) = cPick(i).phase;
        picks{IND} = cPick(i).pickTime;
        unc(IND) = NaN;
        pol{IND} = cPick(i).polarity;
        sta{IND} = cPick(i).station;
        chan{IND} = cPick(i).channel;
        if sncl == 1
            net{IND} = cPick(i).network;
            loc{IND} = cPick(i).location;
        end
        IND = IND+1;
    end
end

% Calculate the forward and backward times from a given pick for a given eq
if exist('picks')
    for i = length(eq): -1 : 1
        for j = 1 : length(picks)
            % Try to find common station
            n = find(strcmp(sta{j},eq(i).sta) == 1);
            if length(n) == 0
                n = 1;
            elseif length(n) > 1
                n = n(1);
            end
            eqDiff(j) = timediff(picks{j}, eq(i).pick{n});
        end
        % Find anything within coffset seconds after initial pick
        rk = find (eqDiff >= 0);
        rr = find (eqDiff < coffset);
        cntr = length(eq(i).sta)+1;
        for ii = 1 : length(rr)
            if find(rr(ii) == rk)
                eq(i).type(cntr) = type(rr(ii));
                eq(i).pick{cntr} = picks{rr(ii)};
                eq(i).sta{cntr} = sta{rr(ii)};
                eq(i).chan{cntr} = chan{rr(ii)};
                if sncl == 1
                    eq(i).net{cntr} = net{rr(ii)};
                    eq(i).loc{cntr} = loc{rr(ii)};
                end
                eq(i).unc(cntr) = unc(rr(ii));
                eq(i).pol(cntr) = pol(rr(ii));
                cntr = cntr + 1;
            end
        end
        kk = find( eqDiff < 0 );
        type = type(kk);
        picks = picks(kk);
        sta = sta(kk);
        chan = chan(kk);
        if sncl == 1
            net = net(kk);
            loc = loc(kk);
        end
        unc = unc(kk);
        pol = pol(kk);
        if isempty(picks)
            break;
        end
        clear eqDiff;
    end
end
%-----------------------------------------------------------------------
% Write pickfiles----------------------------------------------------
% Construct filename
if exist('eq') == 1 && writePick == 1
for i = 1 : length(eq)
    % Get minimum time of picks for filename generation-----------------
    for j = 1 : length(eq(i).pick)
         timeC(j) = datenum(eq(i).pick{j}');
    end
    [m, ind] = min(timeC);
    Dseed = timeadd(eq(i).pick{ind}, -5);
    % Construct filename ------------------------------------------------
    yr = Dseed(1);
    sec = num2str(Dseed(6)); secF = str2num(sec(1));
    fName = sprintf('%04d%02d%02d%02d%02d%01ds', yr, Dseed(2), Dseed(3),...
        Dseed(4), Dseed(5), secF);
    ffName{i} = fName(3:end);
    sTime = Dseed;
    sTime(6) = 0;
    % Begin writing file ------------------------------------------------
    fid = fopen( fName(3:end), 'w' );
    acard = sprintf('%s', fName(1:end-2));
    fprintf( fid, 'A %s s\n', acard );
    for jj = 1 : length(eq(i).pick)
        tt(jj) = timediff(eq(i).pick{jj}, sTime);    % Subtract from filename time
        if strcmp(eq(i).type(jj), 'P') && length(eq(i).pol{jj}) == 0
            dotCard = sprintf('.%s.%s (P P _ %04.2f 0 %02.2f 0)',...
                    eq(i).sta{jj}, eq(i).chan{jj}, tt(jj), eq(i).unc(jj));
                fprintf( fid, '%s\n', dotCard );
        elseif strcmp(eq(i).type(jj), 'P') && length(eq(i).pol{jj}) ~= 0
            dotCard = sprintf('.%s.%s (P P %c %04.2f 0 %02.2f 0)',...
                    eq(i).sta{jj}, eq(i).chan{jj}, eq(i).pol{jj}, tt(jj), eq(i).unc(jj));
                fprintf( fid, '%s\n', dotCard );
        elseif strcmp(eq(i).type(jj), 'S') && length(eq(i).pol{jj}) == 0
            dotCard = sprintf('.%s.%s (P S _ %04.2f 0 %02.2f 0)',...
                    eq(i).sta{jj}, eq(i).chan{jj}, tt(jj), eq(i).unc(jj));
                fprintf( fid, '%s\n', dotCard );
        elseif strcmp(eq(i).type(jj), 'S') && length(eq(i).pol{jj}) ~= 0
            dotCard = sprintf('.%s.%s (P S %c %04.2f 0 %02.2f 0)',...
                    eq(i).sta{jj}, eq(i).chan{jj}, eq(i).pol{jj}, tt(jj), eq(i).unc(jj));
            fprintf( fid, '%s\n', dotCard );
        elseif strcmp(eq(i).type(jj), 'C')
            % Find station in earthquake structure
            for dd = 1 : length(eq(i).sta)
                if strcmp(eq(i).sta(dd), eq(i).sta(jj)) && strcmp(eq(i).type(dd), 'P')
                    hhh = dd;
                    break;
                end
            end
            % Subtract coda pick from P-phase pick
            if exist('hhh') && (hhh ~= jj)
                tc = timediff(eq(i).pick{jj}, eq(i).pick{hhh});
                dotCard = sprintf('.%s.%s (D %04.2f)',...
                    eq(i).sta{jj}, eq(i).chan{jj}, tc);
                fprintf( fid, '%s\n', dotCard );
		
            else
                disp(sprintf('Coda pick has no P-wave pick'));
            end
        elseif strcmp(eq(i).type(jj), 'D')
            dotCard = sprintf('C WT %s', eq(i).pol{jj});
            fprintf( fid, '%s\n', dotCard );
        end
    end
    fclose(fid);
    clear timeC tt
end
else
    ffName = {};
end

%-----------------------------------------------------------------------
% Stat reporting
disp(sprintf('Number of picks: %d', numpicks));
if exist('eq') == 1
	disp(sprintf('Number of events: %d', length(eq)));
	save binder.out.mat eq
else 
	disp(sprintf('Number of events: %d', 0));
	eq = struct([]);
	ffName = {};
end
