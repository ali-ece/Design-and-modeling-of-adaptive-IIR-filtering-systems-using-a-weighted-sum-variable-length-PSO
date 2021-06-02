
bmain = [0.1084 0.5419 1.0837 1.0837 0.5419 0.1084];                             %for case5                  
amain = [1 0.9853 0.9738 0.3864 0.1112 0.0113];

figure(1);
impz(bmain,amain);
hold on
impz(Bsoa,Asoa);
legend('Real Plant','WS-VLPSO');
xlabel('Time (sec)');
ylabel('Amplitude');
title('Impulse Response');
hold off

figure(2);
stepz(bmain,amain);
hold on
stepz(Bsoa,Asoa);
legend('Real Plant','WS-VLPSO');
xlabel('Time (sec)');
ylabel('Amplitude');
title('Step Response');
hold off

figure(3);
H = tf([bmain],[amain],1,'variable','z^-1');
% bode(H,{10^-4,10^1})
bode(H)
hold on
H = tf([Bsoa],[Asoa],1,'variable','z^-1');
bode(H)
legend('Real Plant','WS-VLPSO');
hold off