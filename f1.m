p1=21; %kHz passband1
p2=31; %kHz passband2
s1=19; %kHz stopband1 (transition band=2kHz)
s2=33; %kHz stopband2 (transition band=2kHz)
fs=320; %sampling  frequency
w0=0;
ws1=2*pi*s1/fs;
wp1=2*pi*p1/fs;
wp2=2*pi*p2/fs;
ws2=2*pi*s2/fs;
wn=pi;
w(1)=w0;
w(2)=ws1;
w(3)=wp1;
w(4)=wp2;
w(5)=ws2;
w(6)=wn;
omega=tan(w./2); %bilinear transformation
omega0=sqrt(omega(3)*omega(4));
B=omega(4)-omega(3);
omegaL=((omega.^2)-(omega0^2))./(B.*omega);
del1=0.15;
del2=0.15;
D1=(1/((1-del1)^2))-1;
D2=(1/(del2^2))-1;
epsilon=sqrt(D1);
omegaLs=min(-omegaL(2), omegaL(5));
omegaLp=1;
Nmin=ceil(acosh(sqrt(D2/D1))/acosh(omegaLs/omegaLp));
N=Nmin;
d2=cosh((1/N)*asinh(1/epsilon));
d1=sinh((1/N)*asinh(1/epsilon));
j=1;
for k=1:2*N
    sig(k)=-sin((2*k-1)*pi/(2*N))*d1;
    ome(k)=cos((2*k-1)*pi/(2*N))*d2;
    if(sig(k)<0)
        p(j)=sig(k)+i*ome(k);
        j=j+1;
    end;
end
%scatter(sig,ome);
%axis equal;
%grid on;
%for k=1:numel(x)
%text(x(k),y(k),['(' num2str(x(k)) ',' num2str(y(k)) ')'], 'FontSize', 7);
%end
numm=p(1)*p(2)*p(3)*p(4)/sqrt(1+D1);
denn=poly([p(1), p(2), p(3), p(4)]);

%[bt,at] = lp2bp(numm,denn,omega0,B); %low pass to band pass transformation
%[num,den] = bilinear(bt,at,0.5); %bilinear transformation

syms s z;
lpf(s) = poly2sym(numm,s)/poly2sym(denn,s);    %analog lpf transfer function
bpf(s) = lpf((s*s +omega0*omega0)/(B*s));     %bandpass transformation
bpfz(z) = bpf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BPF
[ns, ds] = numden(bpf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete BPF
[nz, dz] = numden(bpfz(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;                                    %frequency response in dB

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 320e3);
plot(f,abs(H))
grid
%H=((2.45*10^-5)*s^4)/((s^4+0.1314*s^2+4.316*10^-3+0.0256*s^3+1.68*10^-3*s+0.0105*s^2)*(s^4+0.1314*s^2+4.316*10^-3+0.0619*s^3+4.066*10^-3*s+2.73*10^-3*s^2));
%num=[2.45*10^-5, 0, 0, 0, 0];
%den1=[1, 0.0256, 0.1314+0.0105, 0.00168, 0.004316];
%den2=[1, 0.0619, 0.1314+0.00273, 0.004066, 0.004316];
%den = conv(den1, den2);
%F=tf (num, den);
%[numd,dend] = bilinear(num,den,0.5);


