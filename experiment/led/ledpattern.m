classdef ledpattern
    % LEDPATTERN
    % Create a collection of stimulus patterns for the LEDCONTROLLER
    % and LEDCONTROLLER_PI classses.
    % s = ledpattern creates a single pattern
    % s = ledpattern(n) creates a vector of n patterns
    % s = ledpattern(m,n) creates a mxn matrix of patterns
    %
    % See also INTENSITY, SET, CLEAR, CLEAR_ALL, GET_LEDS, DELETE,
    % DUMP.
    
    % 20191025 GW Changed: ledpattern is not a handle class anymore to
    % enable by-value assigments
    properties (Access=protected)
        intensity_red = uint16(50);
        intensity_grn = uint16(50);
        leds_red = zeros(128, 1);
        leds_grn = zeros(128, 1);
    end

    methods
        function obj = ledpattern(m,n)
            % PA_LEDPATTERN constructor            
            if nargin ~= 0 % Allow nargin == 0 syntax
                if nargin < 2
                    n = 1;
                end
                obj(m, n) = ledpattern; % Preallocate object array
            end
        end
        
        function p=clone(this)
           p =  
        end

        function intensity(this, color, value)
            % INTENSITY(color, value)
            % Set the intensity of the red or green leds
            % if color starts with 'r' the intensity of the red leds is set
            % by manpulating the PWM output duty cycle. Otherwise the green
            % intensity is set.
            % The value should be in the range 1..50
            if color(1) == 'r'
                this.intensity_red = uint16(value);
            else
                this.intensity_grn  = uint16(value);
            end
        end

        function delete(obj) %#ok<INUSD>
            % DELETE - destructor
        end

        function set(this, lednr, color, value)
            % SET(lednr, color, value)
            % Turn the led with color on output lednr on or off, depending
            % on the value paramater. 0 is off, otherwise on. If not specified,
            % value is assumed to be 1
            %
            % Note that LEDNRs start at 0!
            if nargin < 4
                value = 1;
            end
            if color(1) == 'r'
                this.leds_red(lednr+1) = (value ~= 0);
            else
                this.leds_grn(lednr+1) = (value ~= 0);
            end
        end

        function clear(this, lednr, color)
            % CLEAR(lednr, color)
            % Shorthand for SET(lednr, color, 0);
            this.set(lednr, color, 0);
        end

        function clear_all(this)
            % CLEAR_ALL
            % Turn all leds off. 
						this.leds_red = zeros(128,1);
            this.leds_grn = zeros(128,1);
        end

        function [r, g, ir, ig] = get_leds(this)
            % [leds_red, leds_grn,intens_red, intens_grn] = GET_LEDS;
            % Returns the value of the internal variables holding led 
            % values and intensities.
            if nargout > 0
                r=this.leds_red;
            end
            if nargout > 1
                g=this.leds_grn;
            end
            if nargout > 2
                ir=this.intensity_red;
            end
            if nargout > 3
                ig=this.intensity_grn;
            end
        end
                
        function dump(this)
            % DUMP - show the content of the internal variables
            this.intensity_red
            this.intensity_grn
            this.leds_red' %#ok<NOPRT>
            this.leds_grn' %#ok<NOPRT>
        end
    end
end
