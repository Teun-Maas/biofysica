function VS = spk_vectorstrength(T)
% SPK_VECTORSTRENGTH(T)
%
% determine vector strength
%
% e.g. JORIS, P.X. (2004) Neural Processing of Amplitude-Modulated Sounds. Physiol. Rev., 84, 541?577.

N	= numel(T);
VS	= sqrt( sum(cos(2*pi*T)).^2+sum(sin(2*pi*T)).^2 ) / N;

