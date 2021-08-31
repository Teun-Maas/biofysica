function Wn = spacing_log(Lofreq, Hifreq, NumFilter)

%Lofreq=350; Hifreq=5500; NumFilter=8;
filtrange1=[]; filtrange2=[]; filtcf=[];
for i=1:NumFilter
   filtcf=[filtcf Lofreq * (Hifreq/Lofreq)^((i-0.5)/NumFilter)];
   filtrange1=[filtrange1 Lofreq * (Hifreq/Lofreq)^((i-1)/NumFilter)];
   filtrange2=[filtrange2 Lofreq * (Hifreq/Lofreq)^((i)/NumFilter)];
end
Wn = [filtrange1 Hifreq];
return;
