vs=vs_servo;

clc;
tic
n=1;
for i=1:n
    home
    vs.print_status;
end
milliseconds_per_variable_lookup=1e3*toc/(n*37)
delete(vs);
