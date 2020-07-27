p1=13; %kHz passband1
p2=23; %kHz passband2
s1=15; %kHz stopband1 (transition band=2kHz)
s2=21; %kHz stopband2 (transition band=2kHz)
fs=250; %sampling  frequency
w0=0;
ws1=2*pi*s1/fs;
wp1=2*pi*p1/fs;
wp2=2*pi*p2/fs;
ws2=2*pi*s2/fs;
wn=pi;
w(1)=w0;
w(2)=wp1;
w(3)=ws1;
w(4)=ws2;
w(5)=wp2;
w(6)=wn;
omega=tan(w./2); %bilinear transformation
omega0=sqrt(omega(2)*omega(5));
B=omega(5)-omega(2);
omegaL=(B.*omega)./((omega0.^2)-(omega.^2));
del1=0.15;
del2=0.15;
D1=(1/((1-del1)^2))-1;
D2=(1/(del2^2))-1;
epsilon=sqrt(D1);
omegaLs=min(-omegaL(4), omegaL(3));
omegaLp=1;
Nmin=log(sqrt(D2/D1))/log(omegaLs/omegaLp);
N=ceil(Nmin);

lower=omegaLp*D1^(-1/(2*N));
higher=omegaLs*D2^(-1/(2*N));
omegaC=1.0850;
j=1;
for k=1:2*N
    pp(k)=omegaC*i*exp(-i*pi*(2*k-1)/(2*N));
    xx=real(pp(k));
    if xx<0
        p(j)=pp(k);
        j=j+1;
    end
end
%scatter(real(pp), imag(pp));
%x=real(pp);
%y=imag(pp);
%axis equal;
%grid on;
%for k=1:numel(x)
%text(x(k),y(k),['(' num2str(x(k)) ',' num2str(y(k)) ')'], 'FontSize', 7);
%end


numm=omegaC^N;
denn=poly([p(1), p(2), p(3), p(4), p(5), p(6)]);
%[bt,at] = lp2bs(numm,denn,omega0,B);%low pass to band stop transformation
%[num,den] = bilinear(bt,at,0.5); %bilinear transformation

syms s z;
lpf(s) = poly2sym(numm,s)/poly2sym(denn,s);        %analog LPF Transfer Function
bsf(s) = lpf((B*s)/(s*s + omega0*omega0));        %bandstop transformation
bsfz(z) = bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(bsfz(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;                                        %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 250e3);
plot(f,abs(H))
grid
