function [Bipo Aipo Z_f P_f] = Matching(R, Cst2)

Rsiz = size(R,2) - Cst2;
Cst    = R(1) + 1;

Bipo      = zeros(1, Rsiz/2);   
Aipo      = zeros(1, Rsiz/2 + Cst2);  
Aipo(1,1) = 1;
for i = 1:Rsiz/2           
% for i = 1:5             
    Bipo(1,i)   = R(i + 1);
%     Aipo(1,i + 1) = R((Rsiz/2) + i + 1);    %### For O < 5    
    if (i < Cst)                              %### For O >= 5
        Aipo(1,i + 1) = R((Rsiz/2) + i + 1);  %### For O >= 5
    end                                       %### For O >= 5        
end
Z_f = roots(Bipo);
P_f = roots(Aipo);

end




