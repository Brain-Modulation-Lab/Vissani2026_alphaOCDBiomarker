function TContr = get_MONOPATTERN(Tbip ,StringValue)

% Unique group identifiers
%groups = findgroups(Tbip.patient_id, Tbip.Hemisphere, Tbip.recording_id);
[GID, patientList, hemiList, sessList, recList] = findgroups(string(Tbip.patient_id), string(Tbip.Hemisphere), Tbip.session_id, Tbip.recording_id);


% Initialize output
TContr = table();

% Loop through seach group
for i = 1:max(GID)
    % Subset table for current group
    idx = GID == i;
    subT = Tbip(idx, :);
    ValueT = subT.(StringValue);

    % Initialize output for this group
    for chan = 0:3
        % Find all rows where chan is involved in the pair
        involvedIdx = (subT.Sensing_ChannelbipP == chan) | (subT.Sensing_ChannelbipN == chan);
        avgVal = mean(ValueT(involvedIdx));

        % For channels 1 and 2, check if SandwichhVal pair exists
        if chan == 1
            pairIdx = ( (subT.Sensing_ChannelbipP == 0 & subT.Sensing_ChannelbipN == 2) | ...
                        (subT.Sensing_ChannelbipP == 2 & subT.Sensing_ChannelbipN == 0) );
        elseif chan == 2
            pairIdx = ( (subT.Sensing_ChannelbipP == 1 & subT.Sensing_ChannelbipN == 3) | ...
                        (subT.Sensing_ChannelbipP == 3 & subT.Sensing_ChannelbipN == 1) );
        else
            pairIdx = false(height(subT), 1);
        end

        if any(pairIdx)
            SandwichVal = mean(ValueT(pairIdx));  % in case multiple rows match
            finalVal = max(avgVal, SandwichVal);

        else            
            SandwichVal = nan;
            finalVal = avgVal;
        end

        % Store result
        TContr = [TContr; 
            table(patientList(i), hemiList(i), sessList(i), recList(i), chan, SandwichVal, finalVal, ...
            'VariableNames', {'patient_id', 'Hemisphere', 'session_id','recording_id', 'Sensing_Channel',[StringValue,'_Sandwich'], [StringValue,'Contr']})];
    end
end