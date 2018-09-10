% clear all
% close all
clc
conditions=3; %number of conditions to compare
addpath (genpath('C:\Program Files\Bruker\AFM MATLAB Toolbox 1.90.15')) % add all necessary folders/subfolders to path
NSMU = NSMatlabUtilities();

ubic = 'D:\BACKUP DISCO CATA\DATA AFM polymers\microgeles lbl 3\18.03.26\'; %C:\Users\catal\Desktop
pfc_file=['ugel-PAH9_40C_2nN.006P.pfc'];
f1=figure;
f2=figure;

for nn=1:conditions


NSMU.Open([ubic pfc_file(nn)]);
[imagePixel, forVolPixel] = NSMU.GetForceVolumeScanLinePixels();      
NumberOfCurves = NSMU.GetNumberOfForceCurves(); %Note: this should be equal to forVolPixel*forVolPixel
    
    %---------------------------------
    %Display FV height image
    f = figure();
    movegui(f,'northwest');
     
    [data, scaleUnit, dataTypeDesc] = NSMU.GetPeakForceCaptureImageData(NSMU.METRIC);
    image(flipud(data),'CDataMapping','scaled');
    set(gca,'YDir','normal');
    axis('tight', 'square'); 
    ax=gca;
%     ax.XTickLabel =[];
%     ax.YTickLabel =[];
    colormap('Copper');
    hold on;
    colorbar();    
    title(strcat('Peak Force',{' '}, dataTypeDesc, ' Image')); 
    xLabel = NSMU.GetScanSizeLabel();
    xlabel(xLabel)
    
    % select region for analysis
    h = imrect;
    position = wait(h); % [xmin ymin width height]
    pos_fv=round(position./2);
    curves=[];
    for line=1:pos_fv(4)
        curve_1=128*(pos_fv(2)-1)+128*pos_fv(1)+128*(line-1);
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

b0=1;
xcurve = [];
ycurve = [];

for i=20:30%N_curves
    curve_ind=curves(i);
    
  [xTrace, xRetrace, yTrace, yRetrace, xLabel, yLabel] = NSMU.CreatePeakForceForceCurveZplot(curve_ind, NSMU.METRIC,1); 
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


if x_cont(contactPointIndex)<10 && y_cont(contactPointIndex)*cte<0.2
%   %%%%figures def vs separation 
    figure(f1);
    plot(xTrace, cte*yTrace,'b-'); hold on
    plot(xRetrace, cte*yRetrace,'m-'); hold on  
    xlabel(xLabel);
    ylabel(yLabel);
hold on
% plot(x_cont(c_point), y_cont(c_point),'*r')
% hold on
plot(x_cont(contactPointIndex), cte*y_cont(contactPointIndex),'*g')
hold on


%   %%%%figures force vs separation^3/2  
  figure(f2);
  xTrace_corr=-xTrace+xTrace(contactPointIndex);
  ind_cropped=find(xTrace_corr>0);
  xx_plot=cte*yTrace(ind_cropped); %force
  yy_plot=xTrace_corr(ind_cropped).^(3/2); %indentation^3/2
  
  plot(xx_plot, yy_plot,'k.-'); hold on
    ylabel('indentation^3^/^2');
    xlabel('Force(nN)');
 hold on
 
 init=0.05;
 endt=0.9;
 
    ps=polyfit(xx_plot(xx_plot<endt & xx_plot>init), yy_plot(xx_plot<endt & xx_plot>init),1);
 
    plot(xx_plot(xx_plot<endt & xx_plot>init),polyval(ps,xx_plot(xx_plot<endt & xx_plot>init)),'-r')
    
%     E_sned=10^(9)*2^0.5*(1-nu^2)/tan(angle(1))/ps(1)
    E_sned=10^(9)*3/4*(1-nu^2)*R^(-0.5)/ps(1);
    
    EM(i).modulus=E_sned;
    disp(['E = ' num2str( E_sned/1000) ' kPa'])

% z0 = xTrace(c_point);
% 
% hTrace= xTrace-yTrace;
% hRetrace= xRetrace-yRetrace;

% ycurve = ycurve + yTrace(c_point-20,c_point+20);
% xcurve = xcurve + hTrace(c_point-20,c_point+20);

% plot((hTrace(15:end)-hTrace(c_point)),yTrace(15:end),'.b')
% plot((hTrace-hTrace(c_point)).^(3/2),yTrace,'b')
%   hold on
%   plot(xRetrace-xRetrace(c_point),yRetrace,'r')
%   hold on

%   plot(0,0,'*g')
% plot(xTrace(c_point),yTrace(c_point),'*g')

    
% height = xTrace(end)-z0;
end
end



end
