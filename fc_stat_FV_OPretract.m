clear all
close all
clc

addpath (genpath('C:\Program Files\Bruker\AFM MATLAB Toolbox 1.90.15')) % add all necessary folders/subfolders to path
NSMU = NSMatlabUtilities();

ubic = 'D:\BACKUP DISCO CATA\DATA AFM polymers\Giussi_cop_nano\muestras BMA\agua\'; %C:\Users\catal\Desktop
pfc_file='bma350_0309.005';

NSMU.Open([ubic pfc_file]);
[imagePixel, forVolPixel] = NSMU.GetForceVolumeScanLinePixels();      
NumberOfCurves = NSMU.GetNumberOfForceCurves(); %Note: this should be equal to forVolPixel*forVolPixel
    
    %---------------------------------
    %Display FV height image
    [data, scaleUnit, dataTypeDesc] = NSMU.GetForceVolumeImageData(NSMU.METRIC);
    f1= figure;
    image(flipud(data),'CDataMapping','scaled');
    set(gca,'YDir','normal');
    axis('tight', 'square'); 
    colormap('Copper');
    hold on;
%     plot(1, 1, 's', 'MarkerSize', 10,'MarkerEdgeColor','b', 'MarkerFaceColor', 'b')   
%     plot(imagePixel, 1, 's', 'MarkerSize', 10,'MarkerEdgeColor','c', 'MarkerFaceColor','c')   
%     plot(imagePixel, imagePixel, 's', 'MarkerSize', 10,'MarkerEdgeColor','g', 'MarkerFaceColor','g')   
    %colorbar();
    
    title(strcat('Force Volume',{' '}, dataTypeDesc, ' Image')); 
    xLabel = NSMU.GetScanSizeLabel();
    xlabel(xLabel);
    
    
    % select region for analysis
    h = imrect;
    position = wait(h); % [xmin ymin width height]
    pos_fv=round(position);
    curves=[];
    for line=1:pos_fv(4)
        curve_1=forVolPixel*(pos_fv(2)-1)+forVolPixel*pos_fv(1)+forVolPixel*(line-1);
    curves=[curves curve_1:curve_1+pos_fv(3)-1];
    plot((pos_fv(1):pos_fv(1)+pos_fv(3)), (pos_fv(2)+line), 's', 'MarkerSize', 3,'MarkerEdgeColor','w', 'MarkerFaceColor', 'b')   
    end
    
    N_curves = length(curves);
    
%     curvas_sel=pos_fv(3)*pos_fv(4)
%     plot(1, 1, 's', 'MarkerSize', 3,'MarkerEdgeColor','b', 'MarkerFaceColor', 'b')   
%     plot(imagePixel/2, imagePixel/2, 's', 'MarkerSize', 3,'MarkerEdgeColor','c', 'MarkerFaceColor','c')   
%     plot(imagePixel, imagePixel, 's', 'MarkerSize', 3,'MarkerEdgeColor','g', 'MarkerFaceColor','g')   


%%

%cantilever parameters
cte = NSMU.GetSpringConstant(1);% 0.078941; %cantilever spring constant in N/m
angle = [19 2]; 
R=20; %tip radius in nm
nu=0.5;
% prompt={'Cantilever spring constant (N/m)'};
% name='K';
% numlines=1;
% defaultanswer={'0.06'};
% answer=inputdlg(prompt,name,numlines,defaultanswer);
% cte=str2num(answer{1});
rango=0.5; %rango en nN de los cortes para ajustar

%%
f2=figure;
% f2=figure;
b0=1;
xcurve = [];
ycurve = [];

%%
for i=1:N_curves
    curve_ind=curves(i);
    
  [xTrace, xRetrace, yTrace, yRetrace, xLabel, yLabel] = NSMU.CreateForceVolumeForceCurveZplot(curve_ind, NSMU.METRIC,0);
%   [xTraces, xRetraces, yTraces, yRetraces, xLabels, yLabels] = NSMU.CreateForceVolumeForceCurveZplot(curve_ind, NSMU.METRIC,1);
%     
% points = length(yTrace);
% 
xTrace=xTrace-yTrace;
xRetrace=xRetrace-yRetrace;
% figure; plot(xTracep-xTracep(end),xTraces,'.r')
% hold on; plot(xRetracep-xRetracep(1),xRetraces,'.k')
x_cont=xTrace;
y_cont=yTrace;
% 
% for i=1:points-round(points/4)
%     contact(i)= round(points/4)+i-1;
%     p = polyfit(x_cont(1:contact(i)),y_cont(1:contact(i)),1);
%     yresid1 = y_cont(1:contact(i))- polyval(p,x_cont(1:contact(i)));
%     MSE1 = sum(yresid1.^2)/length(yresid1);
%     [b,r,J,COVB,MSE2] = nlinfit(x_cont(contact(i):end),y_cont(contact(i):end)-y_cont(contact(i)),@pseudohertz,b0);
%     MSE(i)=100*MSE1+100*MSE2;
% %     plot(z(1:contact(i)), polyval(p,z(1:contact(i))),'b', z(contact(i):end), pseudohertz(b,z(contact(i):end))+d(contact(i)), 'r' )
% %     waitforbuttonpress
% %     close gcf
% end    
% 
% % figure;plot(MSE)
% c_point=contact(MSE==min(MSE));
[contactPointIndex] = FindContactPoint(x_cont, y_cont);
offset_def= mean(y_cont(1:10));

% if x_cont(contactPointIndex)<10 && y_cont(contactPointIndex)*cte<0.2
%   %%%%figures force vs separation 
    f3=figure;
    figure(f3)
    plot(xTrace-x_cont(contactPointIndex), cte*(yTrace-offset_def),'b-'); hold on
    plot(xRetrace-x_cont(contactPointIndex), cte*(yRetrace-offset_def),'m-'); hold on  
    xlabel('indentation (nm)');
    ylabel('Force (nN)');
hold on
% plot(x_cont(c_point), y_cont(c_point),'*r')
% hold on
plot(0, cte*(y_cont(contactPointIndex)-offset_def),'*g')
hold on

xTrace=xTrace-x_cont(contactPointIndex);
xRetrace=xRetrace-x_cont(contactPointIndex);

%%
clear I
clear f_step

for k=1:round((yRetrace(1)-offset_def)*cte/rango)-1
    f_step(k) = rango*(k);  
[~,I(k)] = min(abs((yRetrace(1:contactPointIndex)-offset_def)*cte-f_step(k)));
% I(k)=I(k)+contactPointIndex;
% if I(k)<length(xRetrace)
plot(xRetrace(I(k)),(yRetrace(I(k))-offset_def)*cte,'+g'); hold on;
% end
end


% ajustes y cálculo de E

tan_angle = tan(angle(1)*pi/180);
err_angle = angle(2)*pi/180;
E=zeros(1,k-1);
err_E=zeros(1,k-1);

for n=2:length(I)
    if I(n-1)~=I(n)
cf = fit(xRetrace(I(n):I(n-1)),yRetrace(I(n):I(n-1)),'poly1');
% cf2=fit(flip(xRetrace(I(n):I(n-1))),flip(yRetrace(I(n):I(n-1))),'poly1');
cf_coeff = coeffvalues(cf);
% cf_coeff2 = coeffvalues(cf)
% cf_confint = confint(cf);
plot(xRetrace(I(n):I(n-1)), (polyval(cf_coeff, xRetrace(I(n):I(n-1)))-offset_def)*cte,'-k');
% 
S(n-1) = cf_coeff(1);
% b = cf_coeff(2);
% S_uncert(n-1) = (cf_confint(2,1) - cf_confint(1,1))/2;
% plot(z(c_point:end),indent_all, 'og')
% hold on
% plot(indent_all+z(c_point),d(c_point:end), 'ob')


E(n-1)=0.5*(pi)^0.5*(1-0.5^2)*cte*S(n-1)/(2*tan_angle*(xRetrace(I(n))))*10^6;

% err_E(n-1) = E(n-1)*(0.1^2+(S_uncert(n-1)/S(n-1))^2+(scale*(z(end)-z(c_point)-d(end)+d(c_point)))^(-2)+(err_angle*(1+tan_angle^2)/tan_angle)^2)^0.5;
    end
end

bb=gcf;
saveas(bb,[ubic 'figures\FDcurve' num2str(i) '.tif']) 
close (gcf)

figure(f2)
plot (f_step(1:end-1),E,'*-k')
hold on


%%
E_mean= mean(E(end-4:end-1));
err_E=std(E(end-4:end-1));
disp(['Eabs = ' num2str(E_mean) '+/-' num2str(err_E) ' kPa'])
EM(i)=E_mean;
err_EM(i)=err_E;
end


 
%% Figures

figure; hist(EM)


%%

figure;
errorbar((1:length(EM)),EM,err_EM,'.-');
title('Young modulus r1 (E)')
xlabel('curve number')
ylabel('E (kPa)')
bb=gcf;
saveas(bb,[ubic 'Youngs1.fig']) 

%%
figure; hist(EM(EM>0),20)
title('Young modulus (E)')
xlabel('E(kPa)')
ylabel('counts')

disp(['Mean Young Modulus (kPa) = ' num2str(mean(EM(EM>0))) ' +/- ' num2str(std(EM(EM>0)))])
% disp(['Mean Young Modulus 10-30 nm(kPa) = ' num2str(mean_EM2) ' +/- ' num2str(err_EM2) ', std = ' num2str(sigma_EM2)])

%% save

save([ubic 'results_' pfc_file '_' date '_cte' num2str(cte) '.mat']) 

%%
figure
subplot(1,2,1)
    image(flipud(data),'CDataMapping','scaled');
    set(gca,'YDir','normal');
    axis('tight', 'square'); 
    colormap('Copper');
    hold on;
 A=zeros(imagePixel*imagePixel,1);
 for ii=1:i-1
     if EM(ii)>0
 A(curves(ii))=EM(ii);
     end
 end
 subplot(1,2,2)
 image(flipud(reshape(A,imagePixel,imagePixel)),'CDataMapping','scaled');
    set(gca,'YDir','normal');
    axis('tight', 'square'); 
    colormap('jet');

%%
figure; 
[counts,centers] = hist(EM(EM>0),40);
bar(centers, counts/N_curves*100)
xlabel('Young modulus(kPa)')
ylabel('Probability (%)')


%%
EM_nonano=EM(EM>0);
EM_nano = EM;

%%
xn=[0:0.05:2.2];
figure; 
[countsn,centersn] = hist(EM_nano/1000,xn);
bar(centersn, countsn/length(EM_nano)*100,'k')
xlabel('Young modulus(MPa)')
ylabel('Probability (%)')
xlim([0 2])
%%
% hold on
figure
x=[0:10:200];
[counts,centers] = hist(EM_nonano/1000,x);
bar(centers, counts/length(EM_nonano)*100,'g')
xlabel('Young modulus(MPa)')
ylabel('Probability (%)')
xlim([0 200])
