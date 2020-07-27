fs = 320000;
p1=21000; %Hz passband1
p2=31000; %Hz passband2
s1=19000; %Hz stopband1 (transition band=2kHz)
s2=33000; %Hz stopband2 (transition band=2kHz)

ps1=(p1+s1)/2;
ps2=(p2+s2)/2;
wc1=ps1*2*pi/fs;
wc2=ps2*2*pi/fs;

%Kaiser paramters
A = -20*log10(0.15);
freqT=2000;
wT=2*pi*freqT/fs; 
M1=ceil((A-8)/(2.285*wT)); %Window length for Kaiser Window
M=M1+36;

%Ideal bandpass impulse response of length "M"

alpha=(M-1)/2;
n=[0:1:(M-1)];
m=n-alpha+eps;
hd1=sin(wc1*m)./(pi*m);
hd2=sin(wc2*m)./(pi*m);
bp_ideal=hd2-hd1;

%Kaiser Window of length "M" with beta=0
kaiser_win = (kaiser(M,0))';

firbp = bp_ideal .* kaiser_win;

%magnitude response
[H,f] = freqz(firbp,1,4096, fs);
plot(f,abs(H))
grid