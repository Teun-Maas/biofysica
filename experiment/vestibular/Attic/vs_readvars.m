

map=m2c_memory_map('PLC Globale Variabelen.xlsx');

% try 
memory_map = {
    'X_SW4_S1_Home'  '40'  
    'X_SW4_S2_Home'  '41'  
};



address=map.lookup('X_SW4_S1_Home')

% catch ME
   % disp('something went wrong');
% 
% end

a1=map.X_SW4_S1_Home
a2=map.X_S2_At_Speed
a3=map.cmd_Run_Program_A
a4=map.Table_1A

delete(map);
