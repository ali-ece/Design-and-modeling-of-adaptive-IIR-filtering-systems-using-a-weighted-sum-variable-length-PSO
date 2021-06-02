clc;
clear all;
close all
tic

%% IIR Filter Fitness

[O Hfilt Wfilt] = Main_IIR();

Logcoef = 1;
alpha   = 0.6; % 0.4 1 0.5 0.5
beta    = 0.4; % 0.6 0 0.0 0.1
gamma   = 0.0; % 0.0 1 0.5 0.4
Cst1    = 2; % number 1 for O<5 and 2 for >= for PSO function
Cst2    = 0; % number 1 for O<5,0 for >= for Fitness/IIRSOA function
nvar    = (2 * (O - 1)) + Cst1;
pnvar   = 0;
maxit   = 1000000;

xmin(1)           = 1;
xmax(1)           = (O - 1);
xmin(2:nvar)      = repmat(-2, 1, nvar - 1);
xmax(2:nvar)      = repmat(+2, 1, nvar - 1);
Xpos              = zeros(maxit,nvar);
Xfit              = zeros(maxit,1);
N                 = rand(size(Hfilt,1),1);
% gbest             = zeros(1,nvar);
gbest             = 0;
gbestfit          = inf;
meanfit           = zeros(maxit,1);
minfit            = zeros(maxit,1);

%% Main loop

for it = 1:maxit
    Xpos(it,:) = xmin + (xmax - xmin) .* rand(1,nvar); 
    Xpos(it,1) = round(Xpos(it,1));
    pnvar      = (2 * Xpos(it,1)) + Cst1;
    if pnvar < nvar
        Xpos(it,pnvar+1:nvar) = 0;
    end
    
    %---------------------------------------------------------------------
    Xfit(it)   = Fitness(Xpos(it,1:pnvar),Hfilt,Wfilt,N,O,Logcoef,alpha,beta,gamma,Cst2);
    %---------------------------------------------------------------------
    
    if Xfit(it) < gbestfit
        gbest    = Xpos(it,1:pnvar);
        gbestfit = Xfit(it);
    end
    meanfit(it) = sum(Xfit)/it;
    minfit(it)  = min(Xfit(1:it));
    D = it/maxit;
    if ( D == 1/5) | ( D == 2/5) | ( D == 3/5) | ( D == 4/5)
        it
    end
end

time = toc
gbest
gbestfit
Minimum_Fitness = minfit(it)
Mean_Fitness    = meanfit(it)
[Bsoa Asoa Z_f P_f] = Matching(gbest,Cst2)

% xx = [1,100,200,300,400,500];
% pdf(Xfit([1,100,200,300,400,500]),xx)

figure(1);
plot(1:maxit,meanfit,'r--','LineWidth',2)
ylim([-1 inf]);
% xlim([-1 inf]);
legend('MeanFitness')
xlabel('Iterations')
ylabel('Value')
title('Monte Carlo'); 
grid on

figure(2);
plot(1:maxit,minfit,'b--','LineWidth',2);
axis([0 inf -.01 inf])
grid on
legend('MinFitness')
xlabel('Iterations')
ylabel('Value')
title('Monte Carlo');  

figure(3);
zplane(Z_f,P_f); %%% Displays the poles and zeros of discrete-time systems.
legend('Zero','Pole');
xlabel('Real Part');
ylabel('Imaginary Part');
title('MonteCarlo: Pole-Zero Plot');




