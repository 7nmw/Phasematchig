clc
clear all
close all
% lyamp=0.9; %микрон
% lyamp=[.6:.001:2];
lyams=[.70:.0001:.90]; %микрон
lyami=[1.3:.0001:1.7];
% lyamp=2./(1./lyami+1./lyams);

[lyami, lyams]=ndgrid(lyami,lyams);
% lyamp=2./(1./lyami+1./lyams);
lyamp=1./(1./lyami+1./lyams);

sigma=1;
delta=7.47; %переод модуляции кристалла
L=100; %длина кристалла
T=85; % температура кристалла
% T=[60:1:110];
% [T, T]=ndgrid(T,T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calc  no
a1=5.653;
a2=0.1185;
a3=0.2091;
a4=89.61;
a5=10.85;
a6=1.97*10^-2;
b1=7.941*10^-7;
b2=3.134*10^-8;
b3=-4.641*10^-9;
b4=-2.188*10^-6;

f=(T-24.5).*(T+570.82);

no_p=sqrt(a1+b1*f+(a2+b2*f)./(lyamp.^2-(a3+b3*f)^2)+(a4+b4*f)./(lyamp.^2-a5^2)-a6*lyamp.^2);
no_s=sqrt(a1+b1*f+(a2+b2*f)./(lyams.^2-(a3+b3*f).^2)+(a4+b4*f)./(lyams.^2-a5^2)-a6*lyams.^2);
no_i=sqrt(a1+b1*f+(a2+b2*f)./(lyami.^2-(a3+b3*f)^2)+(a4+b4*f)./(lyami.^2-a5^2)-a6*lyami.^2);

% plot(lyams,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calc ne
a1=5.756;
a2=0.0983;
a3=0.2020;
a4=189.3;
a5=12.52;
a6=1.32*10^-2;
b1=2.779*10^-6;
b2=5.763*10^-8;
b3=3.729*10^-8;
b4=1.415*10^-4;

ne_p=sqrt(a1+b1*f+(a2+b2*f)./(lyamp.^2-(a3+b3*f)^2)+(a4+b4*f)./(lyamp.^2-a5^2)-a6.*lyamp.^2);
ne_s=sqrt(a1+b1*f+(a2+b2*f)./(lyams.^2-(a3+b3*f)^2)+(a4+b4*f)./(lyams.^2-a5^2)-a6.*lyams.^2);
ne_i=sqrt(a1+b1*f+(a2+b2*f)./(lyami.^2-(a3+b3*f)^2)+(a4+b4*f)./(lyami.^2-a5^2)-a6.*lyami.^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ni=abs(ni);
% ns=abs(ns);
% figure
% plot(lyamp,ne_p)

% c=2.998e8 % миркометрах в секунду

% no_p=np;
phi=0;
n_eff_p=((cos(phi)^2./no_p.^2)+sin(phi)^2./ne_p.^2).^(-1/2);
n_eff_s=((cos(phi)^2./no_s.^2)+sin(phi)^2./ne_s.^2).^(-1/2);
n_eff_i=((cos(phi)^2./no_s.^2)+sin(phi)^2./ne_s.^2).^(-1/2);
% plot(lyamp,n_eff_p)


% %Calc dk
% kp=2*pi*real(np)./lyamp;
% ks=2*pi*real(ns)./lyams;
% ki=2*pi*real(ni)./lyami;
% ki=kp-ks;

kqpm=2*pi/delta;
% dk=2*kp-ks-ki;
% dk=2*pi*ni./lyami;
dk=1*2*pi*n_eff_p./lyamp-2*pi*no_s./lyams-2*pi*no_i./lyami-kqpm;

% plot(lyams,
% Z=sin(dk*L/2+eps)./(dk*L/2+eps);
% F=exp(-1i*dk.*L/2+eps).*sinc(dk.*L/2+eps);
F = (sinc(dk.*L/2+eps).*exp(1i*dk.*L/2+eps));
F=abs(F);
% F=sin(dk.*L/2)./(dk.*L/2);
% F=F.^2;
% wp=2*pi*c./lyamp;
% ws=2*pi*c./lyams;
% wi=2*pi*c./lyami;
% wp0=2*pi*c./lambdaP;
% 
sigma=.005; %um
lambdaP=.532; %um

% wp0=2*pi*c./lambdaP;
pump=exp((-(lyamp-lambdaP).^2)./(2*sigma.^2));

%function pump
figure
mesh(lyami,lyams,pump);
colormap('jet');
view(0,90);0
xlabel('\lambda_i'); ylabel('\lambda_s');

% pump_width = 0.05; %[um]
% p_pulsed =(2*sqrt(log(2))/sqrt(pi)/pump_width)^2*exp((-(2./((1./lyami+1./lyams))).^2)*4*log(2)/(pump_width^2));
% p=p_pulsed.*p_pulsed;
% F=F.*pump;
F=F.*conj(F);


%function phasemaching
figure
mesh(lyami,lyams,F);
set(gca,'view',[0 90]);%shading INTERP;
xlabel('\lambda_i'); ylabel('\lambda_s');title('spectral density')
