clear all
close all
load("save2.mat");

[RT,nt]=size(UU_prof_out);
dt=time/nt;
t = -time/2:dt:(time/2-dt);
df=1/(nt*dt);
f=-(nt/2)*df:df:(nt/2-1)*df;
lambda_pulse= 1820;
c = 299792.458;
fo=c/lambda_pulse;
lambda=c./(f+fo);

UU_spec_out=UU_spec_out*time/sqrt(2*pi)*1e-12;

eng=sum(UU_prof_out.*dt,2)';

s_cen=zeros(RT,1);
for i=1:RT
    s_cen(i)=sum(UU_spec_out(i,:).*lambda)./sum(UU_spec_out(i,:));
end

delta_f=5;
delta_t=10;
delta_s=5;
initRT=2000;

figure(1),
FIG_s=log(sum(UU_spec_out,1));

cen_d=cen-mean(cen);
eng_d=eng-mean(eng);
filt_ind=5/100;
midf=fft(eng_d);
midf(floor(end*filt_ind):floor(end-end*filt_ind))=0;
eng_filt=real(ifft(midf));
midf=fft(cen_d);
midf(floor(end*filt_ind):floor(end-end*filt_ind))=0;
cen_filt=real(ifft(midf));
%              engM=lowpass(engM,0.01);
engnorm=eng_filt;%(engM-min(engM))/(  max(engM)-min(engM));
corr_ind1=conv(cen_d,eng_d(end:-1:1))./sqrt(sum(cen_d.^2)*sum(eng_d.^2));
corr_ind2=conv(cen_filt,eng_filt(end:-1:1))./sqrt(sum(cen_filt.^2)*sum(eng_filt.^2));
figure(1)
plot(corr_ind1);
hold on;%-0.77
plot(corr_ind2);
hold off;

UU_prof_out2=zeros(size(UU_prof_out));
for i=1:RT
    utemp=UU_prof_out(i,:);
    [~,pos]=max(abs(utemp));
    up=[utemp,utemp,utemp];
    UU_prof_out2(i,:)=up((pos+nt/2):(pos+3*nt/2-1));
end

figure(2)
                left_color = [0 0 1]; % blue
                right_color = [1 0 0]; % red
set(figure,'defaultAxesColorOrder', [left_color; right_color]);
yyaxis left;
plot(1:RT,eng(1:RT),'b','linewidth',1);
xlabel('\fontname{Times New Roman}RTs');
ylabel('\fontname{Times New Roman}Energy (pJ)');   
yyaxis right;

plot(1:RT,cen(1:RT),'r','linewidth',1);                                                
ylabel('\fontname{Times New Roman}Central Wave.(nm)');             
ax=gca;                               
ax.FontSize=14;
set(gca,'YDir','reverse')


figure(3)
S=max(UU_spec_out);
plot(c./(f + fo),log10(1000*S)*10,'k','linewidth',1)
xlabel('\fontname{Times New Roman} Wavelength (nm)');
ylabel('\fontname{Times New Roman} Intensity (dBm)');
box off
ax=gca;
ax.FontSize=14;
box off



figure(4)
Z2_b2=UU_prof_out2(initRT:delta_t:end,1:delta_s:end);
Z2=Z2_b2-min(min(Z2_b2));
Z2=Z2./max(max(Z2));
surf(Z2','YData',t(1:delta_s:end),'XData',initRT:delta_t:RT,'Edgecolor','none');
% axis([-inf,inf,-250,250]);
colormap(jet)
xlabel('\fontname{Times New Roman} RTs');ylabel('\fontname{Times New Roman} Time (ps)');
zticklabels([])
% zlabel('\fontname{Times New Roman} Norm. Intensity');
colorbar('Ticks',[0,50/200,100/200,150/200,1],'TickLabels',{'\fontname{Times New Roman} dBm','\fontname{Times New Roman} -60dBm',...
    '\fontname{Times New Roman} -20dBm','\fontname{Times New Roman} 13dBm'})
col=colorbar;
col.Title.String='\fontname{Times New Roman} Norm. Intensity';
box off
grid off;
ax=gca;
ax.FontSize=14;
view(89.9,75);

figure(5)
Z3_b2=10*log(1000*UU_spec_out(initRT:delta_t:end,1:delta_f:end));
maxz=max(max(Z3_b2));%-25
minz=min(min(Z3_b2));%0
Z3_b2(Z3_b2<-125)=-125;
Z3_b=Z3_b2-min(min(Z3_b2));
Z3_b=Z3_b./max(max(Z3_b));
surf(Z3_b','YData',c./(f(1:delta_f:end) + fo),'XData',initRT:delta_t:RT,'Edgecolor','none');
axis ([-inf,inf,1720,1900,-inf,inf])
set(gca,'yDir','normal')
xlabel('\fontname{Times New Roman} RTs');ylabel('\fontname{Times New Roman} Wavelength (nm)');
% zlabel('\fontname{Times New Roman} Intensity (dBm)')
axis tight;
bar=colormap(jet);
zticks([]);
ax=gca;
ax.FontSize=15;
box off
grid off;
view(89.9,75);
axis ([-inf,inf,1720,1900,-inf,inf])


col=colorbar('Ticks',[0,25/100,50/100,75/100,1],'TickLabels',{'\fontname{Times New Roman} -125 dBm','\fontname{Times New Roman} -100 dBm','\fontname{Times New Roman} -75 dBm',...
    '\fontname{Times New Roman} -50 dBm','\fontname{Times New Roman} -25 dBm'});
col.Title.String='\fontname{Times New Roman} Norm. Intensity';

figure(7),
% eng=sum(abs(UU_out(1:end,:)),2)*dt;
plot(1000:(RT+999),eng(:),'b','linewidth',1);
xlabel('RTs');
ylabel('Energy (pJ)');
ax=gca;
ax.FontSize=14;

figure,
                left_color = [0 0 1]; % blue
                right_color = [1 0 0]; % red
set(figure,'defaultAxesColorOrder', [left_color; right_color]);
yyaxis left;
plot(initRT:i,eng(initRT:i),'b','linewidth',1);
xlabel('\fontname{Times New Roman}RTs');
ylabel('\fontname{Times New Roman}Energy (pJ)');   
yyaxis right;
plot(initRT:i,cen(initRT:i),'r','linewidth',1);                                                
ylabel('\fontname{Times New Roman}Central Wave.(nm)');             
ax=gca;                               
ax.FontSize=14;
set(gca,'YDir','reverse')




%
init=4913;
count=30;
x=1:count;

x=ones(nt,1)*x;
figure,
for i=0:(count-1)
plot3(t,x(:,i+1),abs(UU_prof_out2(init+i,:)),'color',"#505050",'LineWidth',1.5);hold on;
end
hold off
axis([-50,50,-inf,inf,-inf,inf]);
ylabel('\fontname{Times New Roman} RTs');
xlabel('\fontname{Times New Roman} Time (ps)');
zlabel('\fontname{Times New Roman} Intensity (W)')
ax=gca;
ax.FontSize=14;

