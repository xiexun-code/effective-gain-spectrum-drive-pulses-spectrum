clear all;
close all;
% clc;
        lambda_bw = 8; %gain bandwidth (1Thz~10nm)
        Psat =300;      %saturation level   pJ  
        GAIN=1200;     %small signal gain
        filt_band=10;   % bandwidth of the filter
        %% Simulation Parameters
        c = 299792.458;                      	% speed of light nm/ps
        lambda_pulse = 1820;                  	% pulse central lambda (nm)
        fo=c/lambda_pulse;                       % central pulse frequency (THz)

        %% Numerical Parameters 
        nt = 2^14;                              % number of spectral points
        time = 400;                            	% ps
        dt = time/nt;                           % ps
        t = -time/2:dt:(time/2-dt);             % ps
        df=1/(nt*dt);                           % frequencies separation (Thz)
        f=-(nt/2)*df:df:(nt/2-1)*df;            % frequencies vector (en THz)
%         lambda = c./(f + c/lambda_pulse);    	% lambdas vector (nm)
%         w = 2*pi*f;                             % angular frequencies vector (in THz)
        deltah = 0.00005;                        	% 5cm Initial longitudinal step (km)
        tol = 1e-5;                             %calcaulate tolerance
        roundtrip =10000;% set roundtrips

% two fiber
        expand_ind=0.5;	% estimated index of mode effective area
        TDF.Aeff = pi*2.2^2; %core area of TDF
        TDF.n2 = 2.0;                      	% Kerr coefficient (10^-20*m^2/W)
        TDF.gamma = 2*pi*TDF.n2/lambda_pulse/TDF.Aeff*1e4/expand_ind;	% W^-1 * km^-1
        TDF.alpha =1;                    	% atenuation coef. (km^-1)8 48
        TDF.betaw = [0 30 0 0];	% beta coefficients (ps^n/nm)
        amf1 = TDF;
        amf1.L = 0.005; 
        amf1.gss=GAIN;
        amf1.Psat = Psat;
        amf1.fbw = lambda_bw;
        amf1.fc = c/1850;

        SMF.Aeff = pi*4.1^2;%core area of SMF
        SMF.n2 = 2.0;                      	% Kerr coefficient (10^-20*m^2/W)
        SMF.gamma = 2*pi*SMF.n2/lambda_pulse/SMF.Aeff*1e4/expand_ind;	% W^-1 * km^-1
        SMF.betaw = [0 -50 0 0];
        SMF.alpha =1;                    	% atenuation coef. (km^-1)8 48
        smf1 = SMF;
        smf1.L = 0.01;   

% filter
   
        filter.lambda_c = 1750;
        filter.fc = c/filter.lambda_c;
        filter.f3dB = filt_band;
        filter.n = 1;
        %filter power index
        % coupler parameters                       
        ratio = 0.75;         % Output coupler 75:25 
        spec_z = zeros(roundtrip,nt); %space initial spectrum
        u_z = spec_z;

        %SA
        T1=1;%modulation depth
        SA_pow=15;%SA power

%store output
UU_prof_out=zeros(roundtrip,nt);
UU_spec_out=zeros(roundtrip,nt);
pos_cc=zeros(roundtrip);
u=15*sech(t*2);
eng=zeros(1,roundtrip);
cen=zeros(1,roundtrip);
for i=1:roundtrip
    u =  transfiber(u,dt,deltah,smf1,fo,tol);
    [u,gain_s] = gainfiber(u,dt,deltah,amf1,fo,tol);
    u=SA_modelock(u,T1,SA_pow);
    [u,filt_s] = gaussian(u,filter.f3dB,filter.fc,fo,df);
    [u,uout] = couple(u,0,ratio); 
    eng(i)=sum(abs(uout.^2))*dt;
    cen_THz=sum(abs(fft(uout)).^2.*fftshift(f+fo))./sum(abs(fft(uout)).^2);
    cen(i)=c/cen_THz;
    %CENTERALIZE
%             [~,pos]=max(abs(u));
%             pos_cc(i+1)=pos-nt/2+pos_cc(i);
%             up=[u,u,u];
%             u = up((pos+nt/2):(pos+3*nt/2-1));
            %stop the circle if output  is transient result
            UU_prof_out(i,:)=abs(uout).^2;
            UU_spec_out(i,:)=fftshift(abs(fft(uout)).^2);
            if sum(abs(uout).^2)<0.1
                break;
            end
            
              Z3_b2=abs(fftshift(fft(uout))).^2;
%               Z3_b2=log10(1000*Z3_b2+eps)*10;
               figure(1) 
               %time
               soliton_prof=UU_prof_out(i,:);
               subplot(1,2,1);plot(t,soliton_prof);%./max(soliton_prof),'r');
%                 axis([-inf,inf,0,200])
               %frequency
                subplot(1,2,2);plot(c./(f + fo),(Z3_b2-min(Z3_b2))./(max(Z3_b2)-min(Z3_b2)),'r');
                hold on;
                Gain_eff=gain_s.*filt_s;
                plot(c./(f + fo),Gain_eff/max(Gain_eff),'g');
                hold off;
%                 axis([1720,1800,0,1])
                figure(2)
%                 eng(i)=sum(abs(uout.^2))*dt;
%                 cen_THz=sum(abs(fft(uout)).^2.*fftshift(f+fo))./sum(abs(fft(uout)).^2);
%                 cen(i)=c/cen_THz;
                left_color = [0 0 1]; % blue
                right_color = [1 0 0]; % red
%                 set(fig,'defaultAxesColorOrder', [left_color; right_color]);
                yyaxis left;
                plot(1:i,eng(1:i),'b','linewidth',1);
                xlabel('RTs');
                ylabel('Energy (pJ)');
                yyaxis right;
                ylabel('Central Wave.(nm)');
                plot(1:i,cen(1:i),'r','linewidth',1);
%                 yticks([-20,-10,0,10,20])
%                 yticklabels([1760,1770,1780,1790,1800]);
                set(gca,'YDir','reverse')
%}

end
% '\fontname{Times New Roman}USA', 
figure,
set(figure,'defaultAxesColorOrder', [left_color; right_color]);
yyaxis left;
plot(1:i,eng(1:i),'b','linewidth',1);
xlabel('\fontname{Times New Roman}RTs');
ylabel('\fontname{Times New Roman}Energy (pJ)');   
yyaxis right;

plot(1:i,cen(1:i),'r','linewidth',1);                                                
ylabel('\fontname{Times New Roman}Central Wave.(nm)');             
ax=gca;                               
ax.FontSize=14;
set(gca,'YDir','reverse')
             




function [uu,nf]=transfiber(uu,dt,deltah,modfib,fo,~)
nt = length(uu);                        
w = fftshift(2*pi*(-nt/2:nt/2-1)/(dt*nt));    

% Fix parameter
betaw=modfib.betaw;
LIN= -modfib.alpha + 1i*betaw(2)*(w).^2/2-1i*betaw(3)*w.^3/6+1i*betaw(4)*w.^4/12;

% initialize simulated conditions
propagedlength=0;
Init_mark=1;
Stop_mark=0;
nf=zeros(1,length(nt))+1;
while propagedlength < modfib.L
     if (deltah + propagedlength) >= modfib.L
        deltah = modfib.L - propagedlength;
        Stop_mark=1;
     end

      if Init_mark
         Init_mark=0;
         D_f= exp((LIN)*deltah/2);
         f_temp = fft(uu).*D_f;
         uu=ifft(f_temp);
     else
        NON = uu.*exp(abs(uu).^2.*1i*modfib.gamma*deltah);
         if Stop_mark
            D_f=exp((LIN)*deltah/2);
            f_temp = fft(NON).*D_f;
            uu = ifft(f_temp);
            break;
        else
        DIS= exp((LIN)*deltah);
        f_temp = fft(NON).*DIS;
        uu = ifft(f_temp);
        propagedlength=deltah + propagedlength;
        end
     end
end

end

function [uu,nf]=gainfiber(uu,dt,deltah,modfib,fo,~)
nt = length(uu);                        
w = fftshift(2*pi*(-nt/2:nt/2-1)/(dt*nt));              
df = 1/(nt*dt);


% Fix parameter
betaw=modfib.betaw;
LIN= -modfib.alpha + 1i*betaw(2)*(w).^2/2-1i*betaw(3)*w.^3/6+1i*betaw(4)*w.^4/12;
gain_w = filter_lorentz_tf(uu,modfib.fbw,modfib.fc,fo,df);

% initialize simulated conditions
propagedlength=0;
Init_mark=1;
Stop_mark=0;
nf=zeros(1,length(nt))+1;
while propagedlength < modfib.L
     if (deltah + propagedlength) > modfib.L
        deltah = modfib.L - propagedlength;
        Stop_mark=1;
     end

      if Init_mark
         Init_mark=0;
            gain = gainsaturation(uu,modfib,dt).*gain_w;
            D_f = exp((LIN+fftshift(gain))*deltah/2); 
         f_temp = fft(uu).*D_f;
         uu=ifft(f_temp);
     else
         gain_sat = gainsaturation(uu,modfib,dt);
         gain=gain_sat.*gain_w;
         nf=nf.*exp(deltah*gain);
        NON = uu.*exp(abs(uu).^2.*1i*modfib.gamma*deltah);
        DIS= exp((LIN+fftshift(gain))*deltah); 
        if Stop_mark
            D_f=exp((LIN+fftshift(gain))*deltah/2);
            nf=nf.*exp(deltah/2*gain);
            f_temp = fft(NON).*D_f;
            uu = ifft(f_temp);
            break;
        else
        f_temp = fft(NON).*DIS;
        uu = ifft(f_temp);
        propagedlength=deltah + propagedlength;
        end
      end
end
end

function uu=SA_modelock(uu,T1,SA_pow)
T2=1-T1;% T1 modulation depth
uu=sqrt(T2+T1.*(1./(1+SA_pow./abs(uu).^2))).*uu;
end

function gain_sat = gainsaturation(uu,modfib,dt)
    gain_sat = modfib.gss/2./(1+dt*sum(abs(uu).^2)./modfib.Psat); 
end

function [uout1,uout2]=couple(uin1,uin2,ratio)
uout1 = sqrt(ratio)*uin1 + 1i*sqrt(1-ratio)*uin2;
uout2 = 1i*sqrt(1-ratio)*uin1 + sqrt(ratio)*uin2;

end


function tf = filter_lorentz_tf(ui,fbw,fc,fo,df)
N = size(ui,2);
% The Ui's corresponding frequency
f = (-(N/2)*df:df:(N/2-1)*df) + fo;      % frequencies vector (THz)
tf = (fbw)/2/pi./((f-fc).^2+(fbw/2)^2);
tf = tf/max(tf(:));
end

function [uout,tf2] = gaussian(uin,fbw,fc,fo,df)
    N = size(uin,2);
    % The Ui's corresponding frequency
    f = (-(N/2)*df:df:(N/2-1)*df) + fo;      % frequencies vector (THz)
    tf2 = exp(-log(sqrt(2))*(2/fbw*(f-fc)).^(2)); 
    tf2 = tf2/max(tf2(:));
    uout=ifft(fft(uin).*fftshift(tf2));
end

