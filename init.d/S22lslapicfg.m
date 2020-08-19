% print which lsl_api.cfg files we found
if ispc
  cfgs= { 'c:\etc\lsl_api\lsl_api.cfg', [ getenv('HOMEPATH') '\lsl_api\lsl_api.cfg'], getenv('LSLAPICFG') };
else
  cfgs= { '/etc/lsl_api/lsl_api.cfg', [ getenv('HOME') '/lsl_api/lsl_api.cfg'], getenv('LSLAPICFG') };
end

for ff=cfgs
  f=ff{1};
  if exist(f,'file')
      fprintf('LSL: found api configuration in %s\n',f);
  end
end
clear f ff cfgs
