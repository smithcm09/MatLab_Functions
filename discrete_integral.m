%Make waveform object up here
stations = [get(w_clean, 'station')];
siteStruct = loadSiteTable('/raid/data/antelope/databases/PLUTONS/dbmerged'); %below is adding in station information
minEl = min(siteStruct.elev);
siteSta = siteStruct.sta;
staStruct = struct();
SECS2DAY = 60 * 60 * 24;
len = numel(w_clean);
for i=1:len
    for k = 1:numel(siteSta)
        if strcmp(stations{i}, siteSta{k})
            staStruct(i).sta = get(w_clean(i),'station');
            staStruct(i).lat = siteStruct.lat(k);
            staStruct(i).lon = siteStruct.lon(k);
            staStruct(i).elev = siteStruct.elev(k); %elevation of station in km
            %staStruct(i).dist = distance(eq(earthquake_number).lat, eq(earthquake_number).lon, siteStruct.lat(k), siteStruct.lon(k));
            [times, phasenames] = arrtimes(staStruct(i).dist, eq(earthquake_number).depth);
            staStruct(i);
            w_clean(i) = addfield(w_clean(i), 'ELEV', staStruct(i).elev);
            w_clean(i) = addfield(w_clean(i), 'LAT', staStruct(i).lat);
            w_clean(i) = addfield(w_clean(i), 'LON', staStruct(i).lon);
            w_clean(i) = addfield(w_clean(i), 'DISTANCE', distance(eq(earthquake_number).lat, eq(earthquake_number).lon, siteStruct.lat(k), siteStruct.lon(k)) );
        end
    end
end

%[Y,order] = sort([staStruct.dist]); 


%w_clean_sort = w_clean(order);

%%
SECSPERDAY = 60 * 60 * 24;
%index_values = load('index_values.mat');

index_valuez = [571 549 534 605 566 542 591 611 611 613 574 549 592 579 603 639 615 662 618 623 600]; %this is the index value of the peak for each station
n = 3; %this is because I'm doing 3-component seismic data - might still work for your infrasound
index = reshape(repmat(index_valuez(:).',n,1),1,[]); %reshape to make same length - make [571 571 571 549 549 549 ...]

%chan = get(w_clean_sort, 'channel')
%chan(1:3:numel(chan)) %run to check that this is all HHE
%%
%[index_valuez(1)-150:1:index_valuez(1)+60]
num_vals = 1:3:numel(w_clean);
for wavnum = 1:3:numel(w_clean)
        dataE=get(w_clean(wavnum),'data');
        dataN=get(w_clean(wavnum+1),'data');
        dataZ=get(w_clean(wavnum+2),'data');
        
        %to find min or max indices, use [M,I] = min(data(range of data in
        %which your min likely falls); where M is the minimum value, and I
        %is the index value, for use below. max() works the same way for
        %the maximum value.
        data_range = [index(wavnum)-150:1:index(wavnum)+60]; %specify range - you could determine by finding the peak of interest, then the minimum on each side (see above)
        N = numel(data_range);
        freqE = get(w_clean(wavnum), 'freq');
        dnumE(1) = datenum(get(w_clean(wavnum),'start'));
            for l = 2:numel(dataE)
                dnumE(l) = datenum((l/freqE)/SECSPERDAY+dnumE(1));
            end
            
        E_data = dataE(data_range); %narrows the data down to the range of interest
        N_data = dataN(data_range);
        Z_data = dataZ(data_range);
        
        int_E_data = sum(dataE); %calculating sums :)
        int_N_data = sum(dataN);
        int_Z_data = sum(dataZ);
end