function [O h w] = Main_IIR()

O = 5;
b = [0.1084 0.5419 1.0837 1.0837 0.5419 0.1084];      %for case4
a = [1 0.9853 0.9738 0.3864 0.1112 0.0113];           %for case4

[h,w]=freqz(b,a,200);

end

