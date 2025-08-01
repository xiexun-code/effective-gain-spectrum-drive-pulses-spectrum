close all;
clear all;
belta2=-0.024;%ps^2/m
gama=0.005;%/Wm
gain0=0.524;%/m
sigema=0.01;%ps^2
E_sat=1e6;%pJ
L_GM=1;%m
D_DL=0.012;%ps^2
T_0=0.65;
delta_T=0.35;
P_sat=4;%W
k=0.7;
c = 299792.458;     

tol_round=200;
step_z=1e2; % length step
%nt =512;
nt=2^12; 
Tmax = 100;%30; %ps
deltah = 0.01; 
dtau = (2*Tmax)/nt; % time step

% mshape=9;

tau = (-nt/2:nt/2-1)*dtau; % temporal grid
omega = fftshift((pi/Tmax) * (-(nt/2):1:(nt/2-1))); 
delta_o=5;
%initial pulse
%uu = 10*sech(tau/5);

uu=sech(tau);
uu0=uu;
temp = fftshift(fft(uu)).*2*Tmax/sqrt(2*pi); % spectrum


fo=c/1600;
hhz = 1i*gama*deltah; % nonlinear phase factor
roundtrip=100000;
cen=zeros(1,roundtrip);
for i=1:roundtrip
D= exp(0.5*(1-sigema*(omega-delta_o).^2)*deltah);
f_temp=fft(uu).*D;
uu=ifft(f_temp);
[~,ind]=max(f_temp);
cen(i)=c/(omega(ind)+fo);
if i==1000
    uu1k=uu;
end
if i==2000
    uu2k=uu;
end
if i==30000
    uu30k=uu;
end
end
figure(1)
hold on;
uu_init=abs(fft(uu0).^2);
plot(c./(omega + fo),uu_init./max(uu_init),'-','LineWidth',1.5);
spect=abs(fft(uu1k).^2);
plot(c./(omega + fo),spect./max(spect),'--','LineWidth',1);
spect=abs(fft(uu2k).^2);
plot(c./(omega + fo),spect./max(spect),'--','LineWidth',1);
spect=abs(fft(uu30k).^2);
plot(c./(omega + fo),spect./max(spect),'--','LineWidth',1);
gain_s=1-sigema*(omega-delta_o).^2;
plot(c./(omega + fo),gain_s./max(gain_s),'--','LineWidth',3);
axis([1450,1650,0,1])
hold off;
ylabel('\fontname{Times New Roman} Normalized Int.');
xlabel('\fontname{Times New Roman} Wavelength (nm)');
legend('\fontname{Times New Roman} Initial pulse','\fontname{Times New Roman} Pulse after 1k RTs',...
    '\fontname{Times New Roman} Pulse after 2k RTs','\fontname{Times New Roman} Pulse after 30k RTs',...
    '\fontname{Times New Roman} Gain spectrum')
ax=gca;
ax.FontSize=14;

figure(2)
plot(1:roundtrip,cen);
hold on;
plot(1:roundtrip,ones(1,roundtrip)*c/(fo+delta_o),'--','LineWidth',3);
hold off;
ylabel('\fontname{Times New Roman}Central Wavelength (nm)');
xlabel('\fontname{Times New Roman}RTs')
legend('\fontname{Times New Roman} Spectral trace','\fontname{Times New Roman} Gain center')
ax=gca;
ax.FontSize=14;

% plot(tau,uu);