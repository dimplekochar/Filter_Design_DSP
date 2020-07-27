%Kaiser
fsamp=320000;
A = -20*log10(0.15);
freqT=2000;
wT=2*pi*freqT/fsamp;
pp1=(p1*1000)-(freqT/2);
pp2=(p2*1000)+(freqT/2);
M1=ceil((A-8)/(2.285*wT));
M=131;
B = fir1(M,[2*pp1/fsamp 2*pp2/fsamp],kaiser(M+1,0));
[H,f] = freqz(B,1,1024,fsamp);
%plot(f,abs(H))
%grid on;

%ideal BPF
for i=1:M+1
h(i)=(sin(w(4)*i)-sin(w(3)*i))/(pi*i);
end