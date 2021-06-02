
bmain = [0.1084 0.5419 1.0837 1.0837 0.5419 0.1084];    %for case5                  
amain = [1 0.9853 0.9738 0.3864 0.1112 0.0113];

figure(1);
plot(1:maxit,meanfit,'r--','LineWidth',2)
ylim([-1 inf]);
% xlim([-1 inf]);
legend('MeanFitness')
xlabel('Iterations')
ylabel('Value')
title('Monte-Carlo'); 
grid on

figure(2);
plot(1:maxit,minfit,'b--','LineWidth',2);
axis([0 inf -.01 inf])
grid on
legend('MinFitness')
xlabel('Iterations')
ylabel('Value')
title('Monte-Carlo');  

figure(3);
zplane(Z_f,P_f); %%% Displays the poles and zeros of discrete-time systems.
legend('Zero','Pole');
xlabel('Real Part');
ylabel('Imaginary Part');
title('Monte-Carlo: Pole-Zero Plot');

figure(4);
impz(bmain,amain);
hold on
impz(Bsoa,Asoa);
legend('Real Plant','Monte-Carlo');
xlabel('Time (sec)');
ylabel('Amplitude');
title('Impulse Response');
hold off

figure(5);
stepz(bmain,amain);
hold on
stepz(Bsoa,Asoa);
legend('Real Plant','Monte-Carlo');
xlabel('Time (sec)');
ylabel('Amplitude');
title('Step Response');
hold off

figure(6);
H = tf([bmain],[amain],1,'variable','z^-1');
% bode(H,{10^-4,10^1})
bode(H)
hold on
H = tf([Bsoa],[Asoa],1,'variable','z^-1');
bode(H)
legend('Real Plant','Monte-Carlo');
hold off
