function r=have_biofpy
   % HAVE_BIOFPY - test if biofpy Python library is loaded
   % The library can be accessed through the GLOBAL variable biofpy
   %
   global biofpy
   r=~isempty(biofpy);
end
