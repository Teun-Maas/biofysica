classdef (Abstract) panaservo < handle

% 2018-01-02 GW: fixed maxint problem in convert_profile() function.
%            Now converting position values to int32 instead of int16.
% 

properties
    varmap;
    plc;

end

methods (Abstract, Access=protected)
    
    set_enable(this, value)
    
end


methods (Abstract)

    r = is_enabled(this)
    start(this)
    stop(this)
    print_status(this)
    result = write_profile(this, axis1, axis2)            
    [ax1, ax2] = read_profile_sv(this)
    [ax1, ax2] = read_profile_pv(this)

end

methods
    function this = panaservo(hostaddr, memorytable)
        this.plc = m2c_plc(hostaddr);
        this.varmap = m2c_memory_map(memorytable);
    end

    function delete(this)
        delete(this.varmap);
        delete(this.plc);
    end

    function enable(this)
        if ~this.is_enabled()
            warndlg(sprintf('Enable servo drives, then\npress OK to continue'),...
                'Warning','modal');
        end
    end
    
    function disable(this)
        this.set_enable(0);
    end
        
    function keepalive(this)
        this.plc.keepalive();
    end
    
    function idata = convert_profile(~, data, limit_r_vel, limit_r_pos)
        maxlen = 2000;

        data = reshape(data, [], 1); % make it a column vector
        len = length(data);
        if len > maxlen
           warning('profile too long, truncated to %d samples', maxlen);
           data = data(1:maxlen);
        end
        npad = maxlen-len;
        % convert to int32 data in 1/10 degrees, padded with zeros
        data = int32(round(10*data));
        r_pos = 10*range(data);  % r_pos in degrees
        if  nargin >=4
           if (r_pos(1) < limit_r_pos(1)) || (r_pos(2) > limit_r_pos(2))
              warning('profile position out of range');
           end
        end
% FIXME?
% ??? pad with one leading zero before d/dt? to keep buffer length == 2000
% padding is fine when running profile from only one buffer.

        idata = diff(data);
        idata = [0; idata];

        if nargin >= 3
            r_vel = max(abs(range(idata)));  % r_vel in degrees/second
            if (r_vel > limit_r_vel)
               warning('profile speed out of range');
            end
        end
        idata = int16([idata; zeros(npad, 1)]);
    end
    
    function result = clear_profile(this)
       result = this.write_profile([],[]);
    end


end  % methods

end  % classdef
