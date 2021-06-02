function [Error] = Fitness(R, Hfilt, Wfilt, N, Ord, Logcoef, alpha, beta, gamma, Cst2)
T1     = 0;
T2     = 0;
Cst    = R(1) + 1;
Hfsiz  = size(Hfilt, 1);
Rsiz1  = size(R, 1);
Hfilt_N   = zeros(Hfsiz, 1);
q_N       = zeros(Hfsiz, 1);
F         = zeros(Hfsiz, 1);
Error     = zeros(Rsiz1, 1);

Rsiz2  = size(R,2) - Cst2;
Bipo   = zeros(Rsiz1, Rsiz2/2);        
Aipo   = zeros(Rsiz1, Rsiz2/2 + Cst2);      

Aipo(:,1) = 1;
for i = 1:Rsiz1
    for j = 1:Rsiz2/2              
        Bipo(i,j) = R(i,j + 1);
%         Aipo(i,1 + j) = R(i,(Rsiz2/2) + j + 1);       %### For O < 5
        if (j < Cst)                                    %### For O >= 5
            Aipo(i,1 + j) = R(i,(Rsiz2/2) + j + 1);     %### For O >= 5
        end                                             %### For O >= 5
    end
    
    [q,Q] = freqz(Bipo(i,:),Aipo(i,:),Hfsiz);


Hfilt_N = Hfilt .* N;     % output with input white noise
q_N     = q .* N;         % output with input white noise
    for j = 1:Hfsiz
     F(j,1) = (abs(Hfilt_N(j) - q_N(j)))^ 2;
    end
    P_f        = roots(Aipo);
    T1         = sum(abs(P_f(find(abs(P_f) >= 1))));
    T2         = sum(abs(P_f));
    Error(i,1) = (alpha*((1/Hfsiz)*sum(F))) + (beta*((R(1)/(Ord-1)) * (T1/T2)));
    Hfilt_N    = 0;
    q_N        = 0;
end
end

