function [E_hertz, E_sned, E_r1,E_std1,n1,z0,height]=E_sneddon_Chizhik_spm(ubic, xTrace, xRetrace, yTrace, yRetrace, xLabel, yLabel,cte,R,saturation,iter,thresh, saveq,folder_name)
% data = dlmread(file,'\t',1,0);
% [trace, retrace, scaleUnit, dataTypeDesc] = GetForceCurveData(this, ChannelNumber, UnitType)
% data= data(1:end,:);
nu=0.5;
beta=2;

% start_fc = 1;
% end_fc = length(data(:,1));
% data_z = dlmread('D:\BACKUP DISCO CATA\DATA AFM polymers\microgeles lbl 2\reversibility\analisis\40C_2a.003.pfc-5410_ForceCurveIndex_9826.spm - NanoScope Analysis.txt','\t',1,0);
% z = data_z(start_fc:end_fc,1); %extend

z=xTrace;
d=yTrace;

dr=yRetrace;
% dr = fliplr(data(start_fc:end_fc,4)')'; %retract
% dr = dr(d~=0);

[~,Inmax] = max(d);
ext = d;
ret = dr;
% ext = [dr(Inmax+1:end)' d(1:Inmax)']';
% ret = fliplr([d(Inmax+1:end)' dr(1:Inmax)'])';
% ext = [fliplr(dr(1:128-Inmax)') d(1:Inmax)']';
% ret = [dr(128-Inmax+1:end)' fliplr(d(Inmax+1:end)')]';
% ext = [fliplr(dr(1:55)') d(1:73)']';
% ret = [dr(56:end)' fliplr(d(74:end)')]';
% figure; plot(z,ext,'.k', z, ret, '.r')
d=ext(1:end);
dr=ret(1:end);
z=z(1:end);
% zf=z(end)

scale = length(d)/abs(z(end)-z(1)); %scale in units/nm

% problem_points = find(z==0 & d==0);
% 
% if length(problem_points)>=1
% indices =(1:end_fc);
% z=z(indices<problem_points(1,1));
% d=d(indices<problem_points(1,1));
% end

if saturation~=0
z = z(z<saturation); %extend
d = d(z<saturation); %extend
dr = dr(z<saturation); %retract
end

% if init~=0
% z = z(z>init); %extend
% d = d(z>init); %extend
% dr = dr(z>init); %extend
% end

% z = z(d~=0); %extend
% d = d(d~=0); %extend
% if isempty(dr(d~=0))
% figure; plot(z,d,'.k')
% else
% dr = dr(d~=0); %extend
if saveq==1
figure; plot(z,d,'.k', z, dr, '.g')
%     end
hold on
end

%%%
b0=1;
points = length(d);

for j=1:points-round(points/4)
    contact(j)= round(points/4)+j-1;
    p = polyfit(z(1:contact(j)),d(1:contact(j)),1);
    yresid1 = d(1:contact(j))- polyval(p,z(1:contact(j)));
    MSE1 = sum(yresid1.^2)/length(yresid1);
    [b,r,J,COVB,MSE2] = nlinfit(z(contact(j):end),d(contact(j):end)-d(contact(j)),@pseudohertz,b0);
    MSE(j)=50*MSE1+100*MSE2;
%     figure; plot(z,d,'.k'); hold on
%     plot(z(1:contact(i)), polyval(p,z(1:contact(i))),'b', z(contact(i):end), pseudohertz(b,z(contact(i):end))+d(contact(i)), 'r' )
%     waitforbuttonpress
%     close gcf
end    
% figure;plot(MSE)

% p1=polyfit(z(c_point:end), d(c_point:end),0);
% hold on
% plot(z(c_point:end),polyval(p1,z(c_point:end)),'-b')


c_point=contact(MSE==min(MSE));
z0 = z(c_point);
height = z(end)-z0;

if saveq==1
plot(z0,d(c_point),'*r')
title(['curve ' num2str(iter)])
ylabel('Deflection (nm)')
xlabel('z (nm)')
bb=gcf;
saveas(bb,[folder_name '\FDcurve' num2str(iter) '.tif']) 
close (gcf)
end




%%
d0 = mean(d(c_point-3:c_point));%calculo deflexión 0
d=d-d0; 
z=z-z(c_point);

h= z-d;

if saveq==1
figure; plot(cte*d(h>-20),h(h>-20),'.k')
ylabel('Indentation (nm)')
xlabel('Force (nN)')
hold on
plot(cte*d(c_point),h(c_point),'*r')
title(['curve ' num2str(iter)])
bb=gcf;
saveas(bb,[folder_name '\FIcurve' num2str(iter) '.tif']) 
% waitforbuttonpress
close (gcf)
end
%%
% if saveq==1
% figure; plot(cte*d(h>h(c_point+1)),h(h>h(c_point+1)).^(3/2),'.k')
% % lsline
% ylabel('Indentation^3^/^2 (nm)')
% xlabel('Force (nN)')
% hold on
% plot(cte*d(c_point),(h(c_point))^(3/2),'*r')
% hold on
% % lsline
% title(['curve ' num2str(iter)])
% bb=gcf;
% saveas(bb,[folder_name '\FI32curve' num2str(iter) '.tif']) 
% % waitforbuttonpress
% close (gcf)
% end

%%
% a=1;
% figure; loglog(cte*d(h>-20),h(h>-20),'.k');lsline
% figure; plot(log(h(10:c_point-a)), log(cte*d(10:c_point-a)),'.k')
% hold on
% pp=polyfit(log(h(20:c_point-a)), log(cte*d(20:c_point-a)),1);
% lsline
% plot(log(h(10:c_point-a)),polyval(pp,log(h(10:c_point-a))),'-r');

%%
p1=polyfit(h(10:c_point-5).^(3/2), cte*d(10:c_point-5),1);
E_hertz=10^(9)*3/4*(1-nu^2)*R^(-0.5)*p1(1);
p2=polyfit(h(10:c_point-5).^(2), cte*d(10:c_point-5),1);
E_sned=10^(9)*1/(2*beta)*(1-nu^2)*R^(-0.5)*p2(1);


if saveq==1
figure; 
subplot(1,2,1)
plot(cte*d(h>h(c_point+1)),h(h>h(c_point+1)).^(3/2),'.k')
hold on
plot(polyval(p1,h(10:c_point-5).^(3/2)),h(10:c_point-5).^(3/2),'-b')
% plot(p1(1)*h(10:c_point-5).^(3/2)+p1(2),h(10:c_point-5).^(3/2),'--g')
ylabel('Indentation^3^/^2 (nm)')
xlabel('Force (nN)')
hold on
plot(cte*d(c_point),(h(c_point))^(3/2),'*r')

subplot(1,2,2);
plot(cte*d(h>h(c_point+1)),h(h>h(c_point+1)).^(2),'.k')
hold on
plot(polyval(p2,h(10:c_point-5).^(2)),h(10:c_point-5).^(2),'-b')
ylabel('Indentation^2 (nm)')
xlabel('Force (nN)')
hold on
plot(cte*d(c_point),(h(c_point))^(2),'*r')
hold on
% lsline
title(['curve ' num2str(iter)])
bb=gcf;
saveas(bb,[folder_name '\FI2curve' num2str(iter) '.tif']) 
% waitforbuttonpress
close (gcf)
end
%%

% E = 10^(9)*(2*beta*R^0.5)^(-1)*(1-nu^2)*cte*(diff(flip(d(1:c_point)))./diff(flip(h(1:c_point))))./flip(h(2:c_point)).^0.5;
h_ord = flip(h(10:c_point));
d_ord = flip(d(10:c_point));
E= 10^(9)*(2*beta*R^0.5)^(-1)*(1-nu^2)*cte*(diff(d_ord)./diff(h_ord))./(h_ord(1:end-1)).^0.5;
% close gcf
% figure;
% title(['curve ' num2str(i)])
% semilogy(h(c_point:end-1), E,'ok')
% % plot(de, e-f0,'k', curve_x,4/3*E*(1-v^2)^(-1)*R^(1/2)*(curve_x-curve_x(1)).^(3/2), 'r')
% xlabel('Indentation (nm)')
% ylabel('Young Modulus (Pa)')
% title(['curve ' num2str(i)])
if saveq==1
figure; plot(h_ord(1:end-1), E,'o-k')
% plot(de, e-f0,'k', curve_x,4/3*E*(1-v^2)^(-1)*R^(1/2)*(curve_x-curve_x(1)).^(3/2), 'r')
xlabel('Indentation (nm)')
ylabel('Young Modulus (Pa)')
title(['curve ' num2str(iter)])
bb=gcf;
saveas(bb,[folder_name '\EIcurve' num2str(iter) '.tif']) 
close (gcf)
end
% 
% [~,I1] = min(abs(h-10));
% [~,I2] = min(abs(h-20));
% [~,I3] = min(abs(h-30));
% 
% E_r1 = mean(E(I1-c_point:I2-c_point));
% E_std1 = std(E(I1-c_point:I2-c_point));
% n1 = length(E(I1-c_point:I2-c_point));
% E_r2 = mean(E(I1-c_point:I3-c_point));
% E_std2 = std(E(I1-c_point:I3-c_point));
% n2 = length(E(I1-c_point:I3-c_point));

% [~,I0] = min(abs(h-range(1)));
% [~,If] = min(abs(h-range(2)));
[~,I0] = min(abs(h_ord-thresh));
E_r1 = mean(E(I0:end));
E_std1 = std(E(I0:end));
n1 = length(E(I0:end));

% E_r1 = mean(E(I0-c_point:end));
% E_std1 = std(E(I0-c_point:end));
% n1 = length(E(I0-c_point:end));

disp(['Eabs = ' num2str(E_r1/1000) '+/-' num2str(E_std1/1000) ' kPa, N = ' num2str(n1) 'Ehertz = ' num2str(E_hertz/1000) 'Esned = ' num2str(E_sned/1000)])
% disp(['Eabs2 = ' num2str(E_r2/1000) '+/-' num2str(E_std2/1000) ' kPa, N = ' num2str(n2)])

