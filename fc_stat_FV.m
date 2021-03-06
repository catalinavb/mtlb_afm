% clear all
% close all
clc

addpath (genpath('C:\Program Files\Bruker\AFM MATLAB Toolbox 1.90.15')) % add all necessary folders/subfolders to path
NSMU = NSMatlabUtilities();

ubic = 'D:\BACKUP DISCO CATA\DATA AFM polymers\Giussi_cop_nano\18.05.18\'; %C:\Users\catal\Desktop
pfc_file='m2.004';

NSMU.Open([ubic pfc_file]);
[imagePixel, forVolPixel] = NSMU.GetForceVolumeScanLinePixels();      
NumberOfCurves = NSMU.GetNumberOfForceCurves(); %Note: this should be equal to forVolPixel*forVolPixel
    
    %---------------------------------
    %Display FV height image
    [data, scaleUnit, dataTypeDesc] = NSMU.GetForceVolumeImageData(NSMU.METRIC);
    
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
    pos_fv=round(position./2);
    curves=[];
    for line=1:pos_fv(4)
        curve_1=forVolPixel*(pos_fv(2)-1)+forVolPixel*pos_fv(1)+forVolPixel*(line-1);
    curves=[curves curve_1:curve_1+pos_fv(3)-1];
    plot((pos_fv(1)*2:pos_fv(1)*2+pos_fv(3)*2), 2*(pos_fv(2)+line), 's', 'MarkerSize', 3,'MarkerEdgeColor','w', 'MarkerFaceColor', 'b')   
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


%%
f1=figure;
% f2=figure;
b0=1;
xcurve = [];
ycurve = [];

for i=1:N_curves
    curve_ind=curves(i);
    
  [xTrace, xRetrace, yTrace, yRetrace, xLabel, yLabel] = NSMU.CreateForceVolumeForceCurveZplot(curve_ind, NSMU.METRIC,1); 
%     
% points = length(yTrace);
% 
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
    figure(f1);
    plot(xTrace, cte*(yTrace-offset_def),'b-'); hold on
    plot(xRetrace, cte*(yRetrace-offset_def),'m-'); hold on  
    xlabel(xLabel);
    ylabel('Force (nN)');
hold on
% plot(x_cont(c_point), y_cont(c_point),'*r')
% hold on
plot(x_cont(contactPointIndex), cte*(y_cont(contactPointIndex)-offset_def),'*g')
hold on


%   %%%%figures force vs separation^3/2  
%   figure(f2);
%   xTrace_corr=-xTrace+xTrace(contactPointIndex);
%   ind_cropped=find(xTrace_corr>0);
%   xx_plot=cte*yTrace(ind_cropped); %force
%   yy_plot=xTrace_corr(ind_cropped).^(3/2); %indentation^3/2
%   
%   plot(xx_plot, yy_plot,'k.-'); hold on
%     ylabel('indentation^3^/^2');
%     xlabel('Force(nN)');
%  hold on
%  
%  init=0.05;
%  endt=0.9;
%  
%     ps=polyfit(xx_plot(xx_plot<endt & xx_plot>init), yy_plot(xx_plot<endt & xx_plot>init),1);
%  
%     plot(xx_plot(xx_plot<endt & xx_plot>init),polyval(ps,xx_plot(xx_plot<endt & xx_plot>init)),'-r')
%     
% %     E_sned=10^(9)*2^0.5*(1-nu^2)/tan(angle(1))/ps(1)
%     E_sned=10^(9)*3/4*(1-nu^2)*R^(-0.5)/ps(1);
%     
%     EM(i).modulus=E_sned;
%     disp(['E = ' num2str( E_sned/1000) ' kPa'])



% height = xTrace(end)-z0;
% end
end


%
% figure;
% plot(xcurve/i,ycurve/i, '.-')
% % bb=gcf;
% saveas(gcf,[ubic '\E_vs_indentation.tif']) 
% save([ubic name], 'EM') 

%% Figures

% 
h=0;
for e=1:length(EM)
    if isempty(EM(e).modulus)==0
        h=h+1;
EM1(h) =EM(e).modulus/1000;
    end
end

figure; hist(EM1(EM1>=0),50)
% h=1;
% for e=1:length(EM)
%     if isempty(EM(e).modulus2)==0
%         h=h+1;
% EM2(h) =EM(e).modulus2;
% delta_EM2(h)=EM(e).err2;
%     end
% end

%%
% EM_kPa = EM0;%(EM0>0);
% delta_EM_kPa = delta_EM0;%(EM0>0);

EM1_kPa = EM1(EM1>0)/1000;
delta_EM1_kPa = delta_EM1(EM1>0)/1000;
% EM2_kPa = EM2(EM2>0)/1000;
% delta_EM2_kPa = delta_EM2(EM2>0)/1000;

% EM_kPa = EM_kPa(EM_kPa<100000);
% delta_EM_kPa = delta_EM_kPa(EM_kPa<100000);

figure;
errorbar((1:length(EM1_kPa)),EM1_kPa,delta_EM1_kPa,'.-');
title('Young modulus r1 (E)')
xlabel('curve number')
ylabel('E (kPa)')
bb=gcf;
saveas(bb,[ubic 'Youngs1.fig']) 

%%
figure;
hist(real(EM1_kPa(EM1_kPa<10000)),50)
title('Young modulus (E)')
xlabel('E (kPa)')
ylabel('counts')
h=gcf;
% saveas(h,[ubic 'histogram_250nm_20nm.fig']) 

% figure;
% hist(EM2_kPa,50)
% title('Young modulus (E)')
% xlabel('E (kPa)')
% ylabel('counts')
% h=gcf;
% saveas(h,[ubic 'histogram_Range2.fig']) 

mean_EM1=mean(real(EM1_kPa(EM1_kPa<10000)));
sigma_EM1=std(real(EM1_kPa(EM1_kPa<10000)));
err_EM1 = sigma_EM1/(length(EM1_kPa))^0.5;

mean_hertz =mean(Ef_hertz)
mean_sned = mean(Ef_sned)
% mean_EM2=mean(EM2_kPa);
% sigma_EM2=std(EM2_kPa);
% err_EM2 = sigma_EM2/(length(EM2_kPa))^0.5;
figure; hist(Ef_hertz,50)
title('Young modulus (E)')
xlabel('E hertz (kPa)')
ylabel('counts')
figure; hist(Ef_sned,50)
title('Young modulus (E)')
xlabel('E sned (kPa)')
ylabel('counts')
disp(['Mean Young Modulus (kPa) = ' num2str(mean_EM1) ' +/- ' num2str(err_EM1) ', std = ' num2str(sigma_EM1)])
% disp(['Mean Young Modulus 10-30 nm(kPa) = ' num2str(mean_EM2) ' +/- ' num2str(err_EM2) ', std = ' num2str(sigma_EM2)])

%% save

save([ubic 'results_' date '_cte' num2str(cte) 'nN_saturation' num2str(saturation) 'threshold' num2str(thresh) 'nm.mat']) 

  
% 25
% 186.0554 +/- 2.2085, std = 112.5448  
% 192.8448 +/- 2.0658, std = 105.2569  *
% 25 16/03
% 143.4302 +/- 1.0859, std = 57.5553  *
% 25 28/3
% Mean Young Modulus (kPa) = 516.4187 +/- 3.2668, std = 130.6313  *

% 40
% 445.6819 +/- 3.3169, std = 175.798  *
% 40 16/3
% 201.0392 +/- 0.44195, std = 23.2015   _010 *
% 40 28/3
% Mean Young Modulus (kPa) = 385.3745 +/- 2.1993, std = 87.9446  *

% 
% temp = [25 40 25 40 25 40];
% xaxis = (1:6);
% yng = [192.8448 445.6819 143.4302 201.0392 516.4187 385.3745];
% er_yng = [112.5 175.7 57.5 23.2 130.6 87.9];
% 
% figure; errorbar(xaxis, yng, er_yng)