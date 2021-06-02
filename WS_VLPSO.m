% for Rep=1:100
clc;
clear all;
close all
tic
%% IIR Filter Fitness

[O Hfilt Wfilt] = Main_IIR();

Logcoef = 1;
alpha   = 0.6;
beta    = 0.4;
gamma   = 0.0; 
Cst1    = 2; % number 1 for O<5 and 2 for upper for PSO function
Cst2    = 0; % number 1 for O<5,0 for upper for Fitness/IIRSOA function
nvar    = (2 * (O - 1)) + Cst1;
npop    = 50;
w       = 1;

maxit   = 200;
wdamp   = 0.99;
c1      = 2;
c2      = 2;
xmin(1) = 1;
xmax(1) = (O - 1);
xmin(2:nvar) = repmat(-2, 1, nvar - 1);
xmax(2:nvar) = repmat(+2, 1, nvar - 1);
dx             = xmax - xmin;
vmax           = 0.1 * dx;
empty_particle.position  = [];
empty_particle.velocity  = [];
empty_particle.cost      = [];
empty_particle.pbest     = [];
empty_particle.pbestcost = [];
empty_particle.nvar      = [];
empty_particle.pnvar     = [];

particle      = repmat(empty_particle,npop,1);
N             = rand(size(Hfilt,1),1);
gbest         = zeros(maxit,nvar);
gbestcost     = zeros(maxit,1);
gnvar         = zeros(maxit,1);
meanfit       = zeros(maxit,1);
Temp1particle = zeros(1,nvar); % The initial nvar is the maximom value of the first dimension
Temp1velocity = zeros(1,nvar);
Temp2particle = zeros(1,nvar);
Temp2velocity = zeros(1,nvar);
for it = 1:maxit
    if it == 1
        gbestcost(1) = inf;
        for i = 1:npop
            
            particle(i).position    = xmin + (xmax - xmin) .* rand(1,nvar); 
            particle(i).position(1) = round(particle(i).position(1));        
            particle(i).nvar        = (2 * (particle(i).position(1)))+Cst1;   
            if particle(i).nvar ~= nvar                                     
                dif1 = particle(i).nvar - nvar;                                
                if dif1 < 0                                                  
                    particle(i).position = particle(i).position(1:particle(i).nvar);
                else %dif1 > 0                                                  
                    
                    particle(i).position(nvar+1:nvar+dif1) = xmin(2:1+dif1)...
                        + (xmax(2:1+dif1) - xmin(2:1+dif1)) .* rand(1,dif1);       
                end                                                          
                dif1 = 0;
            end                                                              
            
            particle(i).velocity = zeros(1,particle(i).nvar);
            
            %%%%%%%%*********************************
 particle(i).cost = Fitness(particle(i).position,Hfilt,Wfilt,N,O,Logcoef,alpha,beta,gamma,Cst2);
            %%%%%%%%
            particle(i).pnvar     = particle(i).nvar;
            particle(i).pbest     = particle(i).position;%like line 81(gbest(it,1:gnvar)),nvar can be replace to term of line 75.
            particle(i).pbestcost = particle(i).cost;
            
            if particle(i).pbestcost < gbestcost(it)
                gnvar(it)                 = particle(i).nvar;
                gbest(it,1:gnvar(it))     = particle(i).pbest;
                gbest(it,gnvar(it)+1:end) = 0;
                gbestcost(it)             = particle(i).pbestcost;
            end
        end
    else
        gbest(it,:)   = gbest(it - 1,:);
        gbestcost(it) = gbestcost(it - 1);
        gnvar(it)     = gnvar(it - 1);
        
        for i = 1:npop
        
        dif3          = particle(i).pnvar - particle(i).nvar;
        Temp1particle = particle(i).position;
        Temp1velocity = particle(i).velocity;
 
        
            if dif3 ~= 0
                if dif3 < 0
                    Temp1particle = Temp1particle(1:particle(i).pnvar);
                    Temp1velocity = Temp1velocity(1:particle(i).pnvar);
                else %dif3 > 0
                    Temp1particle(particle(i).nvar+1:particle(i).nvar+dif3)...
                        = xmin(2:1+dif3) + (xmax(2:1+dif3) - xmin(2:1+dif3)) .* rand(1,dif3);
                    Temp1velocity(particle(i).nvar+1:particle(i).nvar+dif3)...
                        = min(max(rand(1,dif3),-vmax(2:1+dif3)),vmax(2:1+dif3));
                end
            end

            
        dif4          = gnvar(it) - particle(i).nvar;
        Temp2particle = particle(i).position;
        Temp2velocity = particle(i).velocity;
            
           if dif4 ~= 0
               if dif4 < 0
                   Temp2particle = Temp2particle(1:gnvar(it));
                   Temp2velocity = Temp2velocity(1:gnvar(it));
               else %dif4 > 0
                   Temp2particle(particle(i).nvar+1:particle(i).nvar+dif4)...
                       = xmin(2:1+dif4) + (xmax(2:1+dif4) - xmin(2:1+dif4)) .* rand(1,dif4);
                   Temp2velocity(particle(i).nvar+1:particle(i).nvar+dif4)...
                       = min(max(rand(1,dif4),-vmax(2:1+dif4)),vmax(2:1+dif4));

               end
           end
        Temp3pbest = particle(i).pbest - Temp1particle;
        Temp4gbest = gbest(it,1:gnvar(it)) - Temp2particle;
        
        dif5 = particle(i).pnvar - gnvar(it);

            if dif5 ~= 0
                if dif5 < 0
                    particle(i).position(particle(i).nvar+1:particle(i).nvar+abs(dif4))...
                        = xmin(2:1+abs(dif4)) + (xmax(2:1+abs(dif4))...
                        - xmin(2:1+abs(dif4))) .* rand(1,abs(dif4));
                    particle(i).position = particle(i).position(1:gnvar(it));
                    particle(i).nvar     = gnvar(it);
                    particle(i).velocity = Temp2velocity;
                    Temp3pbest(particle(i).pnvar + 1:gnvar(it)) = ...
                        min(max(rand(1,abs(dif5)),-vmax(2:1+abs(dif5))),vmax(2:1+abs(dif5)));
                else %dif5 > 0
                    particle(i).position(particle(i).nvar+1:particle(i).nvar+abs(dif3))...
                        = xmin(2:1+abs(dif3)) + (xmax(2:1+abs(dif3))...
                        - xmin(2:1+abs(dif3))) .* rand(1,abs(dif3));
                    particle(i).position = particle(i).position(1:particle(i).pnvar);
                    particle(i).nvar     = particle(i).pnvar;
                    particle(i).velocity = Temp1velocity;
                    Temp4gbest(gnvar(it) + 1:particle(i).pnvar) = ...
                        min(max(rand(1,dif5),-vmax(2:1+dif5)),vmax(2:1+dif5));
                end
            else
                particle(i).position = particle(i).position(1:gnvar(it));
                particle(i).nvar     = gnvar(it);
                particle(i).velocity = Temp2velocity;
            end

            
            particle(i).velocity = w * particle(i).velocity...
                + c1 * rand * Temp3pbest + c2 * rand * Temp4gbest;
            
            particle(i).velocity = min(max(particle(i).velocity,...
                -vmax(1:particle(i).nvar)),vmax(1:particle(i).nvar));
            
            particle(i).position = particle(i).position + particle(i).velocity;
            
            particle(i).position = min(max(particle(i).position,...
                xmin(1:particle(i).nvar)),xmax(1:particle(i).nvar));
            
            particle(i).position(1) = round(particle(i).position(1));
            Temp5nvar               = particle(i).nvar;
            particle(i).nvar        = (2 * (particle(i).position(1)))+Cst1;
            
            if particle(i).nvar ~= Temp5nvar
                dif2 = particle(i).nvar - Temp5nvar;
                if dif2 < 0

                    particle(i).position = particle(i).position(1:particle(i).nvar);
                    particle(i).velocity = particle(i).velocity(1:particle(i).nvar);
                else %dif2 > 0
                    particle(i).position(Temp5nvar+1:Temp5nvar+dif2)...
                        = xmin(2:1+dif2) + (xmax(2:1+dif2) - xmin(2:1+dif2)) .* rand(1,dif2);
                    particle(i).velocity(Temp5nvar+1:Temp5nvar+dif2)...
                        = min(max(rand(1,dif2),-vmax(2:1+dif2)),vmax(2:1+dif2));
                end
                Temp5nvar = 0;
                dif2      = 0;
            end
            
          %%%************************ 
particle(i).cost = Fitness(particle(i).position,Hfilt,Wfilt,N,O,Logcoef,alpha,beta,gamma,Cst2);
          %%%
            if particle(i).cost < particle(i).pbestcost
                particle(i).pnvar     = particle(i).nvar;
                particle(i).pbest     = particle(i).position;
                particle(i).pbestcost = particle(i).cost;

                if particle(i).pbestcost < gbestcost(it)
                    gnvar(it)                 = particle(i).nvar;
                    gbest(it,1:gnvar(it))     = particle(i).pbest;
                    gbest(it,gnvar(it)+1:end) = 0;
                    gbestcost(it)             = particle(i).pbestcost;
                end
            end
        end
    end
  
    disp(['Iter= ' num2str(it) ' // Best Cost = ' num2str(gbestcost(it))]);
    
    for k = 1:npop
        meanfit(it) = meanfit(it) + particle(k).cost;
    end
    meanfit(it) = meanfit(it) / npop;
    
    w = w * wdamp;
end

disp([ ' Best Solution = '  num2str(gbest(it,1:gnvar(it)))])
disp([ ' Best Fitness = '  num2str(gbestcost(it))])

% disp([ ' Time = '  num2str(toc)])
time=toc
% e = particle().pbest;
[Bsoa Asoa Z_f P_f] = Matching(gbest(it,1:gnvar(it)),Cst2)

figure(1);
plot(gbestcost,'r','LineWidth',2);
% hold on
% plot(meanfit,'.b','LineWidth',2);
legend('Bests')
xlabel('Iteration')
ylabel('Fitness')
title('WSPSO');

figure(2);
zplane(Z_f,P_f); %%% Displays the poles and zeros of discrete-time systems.
legend('Zero','Pole');
xlabel('Real Part');
ylabel('Imaginary Part');
title('Pole-Zero Plot');

% bmain = [0.1084 0.5419 1.0837 1.0837 0.5419 0.1084];                             %for case5                  
% amain = [1 0.9853 0.9738 0.3864 0.1112 0.0113];

% figure(1);
% impz(bmain,amain);
% hold on
% impz(Bsoa,Asoa);
% legend('Real Plant','WS-VLPSO');
% xlabel('Time (sec)');
% ylabel('Amplitude');
% title('Impulse Response');
% hold off

% figure(2);
% stepz(bmain,amain);
% hold on
% stepz(Bsoa,Asoa);
% legend('Real Plant','WS-VLPSO');
% xlabel('Time (sec)');
% ylabel('Amplitude');
% title('Step Response');
% hold off

% figure(3);
% H = tf([bmain],[amain],1,'variable','z^-1');
% % bode(H,{10^-4,10^1})
% bode(H)
% hold on
% H = tf([Bsoa],[Asoa],1,'variable','z^-1');
% bode(H)
% legend('Real Plant','WS-VLPSO');
% hold off
