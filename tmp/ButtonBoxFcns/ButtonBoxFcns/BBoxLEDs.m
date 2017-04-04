function BBoxLEDs(Device, LEDVec)
% BBoxLEDs -- Light the LEDs on the ButtonBox. LEDVec is a 4-element
%   vector in which elements 1-4 correspond to LEDs 1-4. Device is the
%   System 3 device, normally RX6.

LEDVec= [0 0 LEDVec];   % shift right because bit 2 is LED1, etc.
LEDWord= 0;
LEDWord= sum(bitset(LEDWord,find(LEDVec)));
Device.SetTagVal('LEDWord',LEDWord);
