function [patterns, pattern_times, pattern_events] = makeOrderedLedPatterns(modalities, cfg)
    
    
    m_Ordered = orderLedModalities(modalities);
    
    times = findUniqueTimes(m_Ordered);
    events = findUniqueEvents(modalities);
    
    i_m=1; % index to modalities
    i_p=1;  % index to LED patterns
    nn=length(times)+length(events)-1; % this is not quite correct,
    % we may end up with too many patterns, and this migth lead to errors when
    % loading into the led controller (too many patterns)
    patterns=ledpattern(nn);
    pattern_times=[];
    pattern_events=[];
    
    for t=times
        %make a new LED pattern here; because we
        %handle changes to patterns, make it a copy of the previous one.
        if i_p > 1
            patterns(i_p).copyfrom(patterns(i_p-1));
        end
        
        m=m_Ordered(i_m);
        current_event = m.trigger_event;
        
        while m.trigger_delay==t
            if m.trigger_event ~= current_event
                % still the same trigger_delay, but because we are
                % after the next trigger, we make a new pattern.
                patterns(i_p+1).copyfrom(patterns(i_p));
                %%DEBUG
                %%patterns(i_p).dump;
                i_p = i_p+1;
                current_event=m.trigger_event;
            end
            
            % add this LED change to this pattern
            col=m.colour;
            onoff=m.trigger_action;
            patterns(i_p).set(m.Z,cfg.ledcolours{col},onoff);
            patterns(i_p).intensity(cfg.ledcolours{col},m.intensity); % hoop: range 0-255, sphere range 1-50
            pattern_times(i_p)=m.trigger_delay; %#ok<AGROW>
            pattern_events(i_p)=m.trigger_event; %#ok<AGROW>
            i_m = i_m+1;
            if i_m <= length(m_Ordered)
                m=m_Ordered(i_m);
            else
                break;
            end
        end
        %%DEBUG
        %%patterns(i_p).dump;
        i_p = i_p+1;
    end
end

function m_Ordered=orderLedModalities(modalities)
    
    m=findLedModalities(modalities);
    m_split = splitModalitiesOnOff(m);
    
    events = findUniqueEvents(m);
    
    m_Ordered=[];
    for ev=events
        m_byEv = findModalitiesByField(m_split,'trigger_event',ev);
        m_byAction = sortModalitiesByField(m_byEv,'trigger_action');
        m_byDelay = sortModalitiesByField(m_byAction,'trigger_delay');
        m_Ordered=[m_Ordered, m_byDelay]; %#ok<AGROW>
    end
    
end


function m = findLedModalities(modalities)
    m=findModalities(modalities, 'LED', 'SKY', 'LASER', 'IRLED');
end

function m = findModalities(modalities, varargin)
    found = false(size(modalities));
    for arg = varargin
        found = found | strcmpi({modalities.modality}, arg);
    end
    m = modalities(found);
end

function m = splitModalitiesOnOff(modalities)
    % split LED modalities into on and off events.
    % replace 'onevent', 'ondelay', 'offevent' and 'offdelay' by new fields
    % 'trigger_event', 'trigger_delay' and 'trigger_action'
    % 'trigger_action' is 0 for 'on' event and 1 for 'off' events
    % 'trigger_event' and 'trigger_delay' are copies of the original 'on*'
    % and 'off*' fields.
    
    % prepare output structures, replace
    replace_fields={'ondelay','offdelay','onevent','offevent'};
    m_on = rmfield(modalities,replace_fields);
    m_on(1).trigger_event=[];
    m_on(1).trigger_delay=[];
    m_on(1).trigger_action=[];
    
    m_off=m_on;
    
    for ii=1:length(modalities)
        m_on(ii).trigger_delay = modalities(ii).ondelay;
        m_on(ii).trigger_event = modalities(ii).onevent;
        m_on(ii).trigger_action = 1;
        
        m_off(ii).trigger_delay = modalities(ii).offdelay;
        m_off(ii).trigger_event = modalities(ii).offevent;
        m_off(ii).trigger_action = 0;
    end
    m = [m_on, m_off];
    
    cev={m.trigger_event};            % cell array with trigger_event
    nonempty=~cellfun('isempty',cev); % find the non-empty ones
    m=m(nonempty);                    % discard the empty modalities
end

function ev = findUniqueEvents(modalities)
    n=length(modalities);
    ev_on = cell(1,n);
    ev_off = cell(1,n);
    for ii=1:n
        ev_on{ii}=modalities(ii).onevent;
        ev_off{ii}=modalities(ii).offevent;
    end
    
    ev = unique([cell2mat(ev_on), cell2mat(ev_off)]);
    ev=sort(ev(~isnan(ev)));
end

function times = findUniqueTimes(modalities)
    times=unique(cell2mat({modalities.trigger_delay}));
end

function m = findModalitiesByField(modalities, fieldName, value)
    found=[];
    n=length(modalities);
    for ii=1:n
        if modalities(ii).(fieldName)==value
            found = [found, ii]; %#ok<AGROW>
        end
    end
    m=modalities(found);
end

function m = sortModalitiesByField(modalities, fieldName)
    keys=cell2mat({modalities.(fieldName)});
    [~,I]=sort(keys);
    m=modalities(I);
end
