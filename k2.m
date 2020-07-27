fs = 250000;
p1=13000; %kHz passband1
p2=23000; %kHz passband2
s1=15000; %kHz stopband1 (transition band=2kHz)
s2=21000; %kHz stopband2 (transition band=2kHz)
ps1=(p1+s1)/2;
ps2=(p2+s2)/2;
wc1=ps1*2*pi/fs;
wc2=ps2*2*pi/fs;

%Kaiser paramters
A = -20*log10(0.15);
freqT=2000;
wT=2*pi*freqT/fs; 
M1=ceil((A-8)/(2.285*wT)); %Window length for Kaiser Window
M=M1+25;

%Ideal bandpass impulse response of length "M"

alpha=(M-1)/2;
n=[0:1:(M-1)];
m=n-alpha+eps;
hd1=sin(wc1*m)./(pi*m);
hd2=sin(wc2*m)./(pi*m);
hd3=sin(pi*m)./(pi*m);

bs_ideal=hd3+hd1-hd2;

%Kaiser Window of length "M" with beta=0
kaiser_win = (kaiser(M,0))';

firbs = bs_ideal .* kaiser_win;

%magnitude response
[H,f] = freqz(firbs,1,4096, fs);
plot(f,abs(H))
grid