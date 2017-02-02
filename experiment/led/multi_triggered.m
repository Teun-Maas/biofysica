%
% example for using the led drivers with multiple triggers
% 
ledcontroller = pa_ledcontroller('fp-web2');
ledcontroller.trigger_enable(1);
ledcontroller.set(i,'red', 1);
ledcontroller.flush;

ledcontroller.set(i,'grn', 1);
ledcontroller.flush;

delete(ledcontroller);
