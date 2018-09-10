clear all
close all
clc

addpath (genpath('C:\Program Files\Bruker\AFM MATLAB Toolbox 1.90.15')) % add all necessary folders/subfolders to path
NSMU = NSMatlabUtilities();

ubic1 = 'D:\BACKUP DISCO CATA\DATA AFM polymers\Giussi_cop_nano\18.05.18\'; %C:\Users\catal\Desktop
pfc_file1='no nano.006';

ubic2 = ubic1;
pfc_file2='m2.004';

%%
NSMU.Open([ubic1 pfc_file1]);
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


%%

%cantilever parameters
cte = NSMU.GetSpringConstant(1);% 0.078941; %cantilever spring constant in N/m
angle = [19 2]; 
R=20; %tip radius in nm
nu=0.5;

%%
b0=1;
xcurve = [];
ycurve = [];

%%
f3=figure;

%%
for i=1:N_curves
    curve_ind=curves(i);
    
  [xTrace, xRetrace, yTrace, yRetrace, xLabel, yLabel] = NSMU.CreateForceVolumeForceCurveZplot(curve_ind, NSMU.METRIC,0);

xTrace=xTrace-yTrace;
xRetrace=xRetrace-yRetrace;

x_cont=xTrace;
y_cont=yTrace;
[contactPointIndex] = FindContactPoint(x_cont, y_cont);
offset_def= mean(y_cont(1:10));

%   %%%%figures force vs separation 
    figure(f3)
    plot(xTrace-x_cont(contactPointIndex), cte*(yTrace-offset_def),'b-'); hold on
%     plot(xRetrace-x_cont(contactPointIndex), cte*(yRetrace-offset_def),'m-'); hold on  
    xlabel('indentation (nm)');
    ylabel('Force (nN)');
hold on
% plot(x_cont(c_point), y_cont(c_point),'*r')
% hold on
% plot(0, cte*(y_cont(contactPointIndex)-offset_def),'*g')
hold on

xTrace=xTrace-x_cont(contactPointIndex);
xRetrace=xRetrace-x_cont(contactPointIndex);
end

%%
NSMU.Open([ubic2 pfc_file2]);
[imagePixel, forVolPixel] = NSMU.GetForceVolumeScanLinePixels();      
NumberOfCurves = NSMU.GetNumberOfForceCurves(); %Note: this should be equal to forVolPixel*forVolPixel
    
    %---------------------------------
    %Display FV height image
    [data2, scaleUnit, dataTypeDesc] = NSMU.GetForceVolumeImageData(NSMU.METRIC);
    f2= figure;
    image(flipud(data2),'CDataMapping','scaled');
    set(gca,'YDir','normal');
    axis('tight', 'square'); 
    colormap('Copper');
    hold on;
   
    title(strcat('Force Volume',{' '}, dataTypeDesc, ' Image')); 
    xLabel = NSMU.GetScanSizeLabel();
    xlabel(xLabel);
    
    
    % select region for analysis
    h = imrect;
    position = wait(h); % [xmin ymin width height]
    pos_fv=round(position);
    curves2=[];
    for line=1:pos_fv(4)
        curve_1=forVolPixel*(pos_fv(2)-1)+forVolPixel*pos_fv(1)+forVolPixel*(line-1);
    curves2=[curves2 curve_1:curve_1+pos_fv(3)-1];
    plot((pos_fv(1):pos_fv(1)+pos_fv(3)), (pos_fv(2)+line), 's', 'MarkerSize', 3,'MarkerEdgeColor','w', 'MarkerFaceColor', 'b')   
    end
    
    N_curves2 = length(curves2);


%%

%cantilever parameters
cte = NSMU.GetSpringConstant(1);% 0.078941; %cantilever spring constant in N/m
angle = [19 2]; 
R=20; %tip radius in nm
nu=0.5;

%%
b0=1;
xcurve = [];
ycurve = [];

%%
for i=1:N_curves2
    curve_ind=curves2(i);
    
  [xTrace, xRetrace, yTrace, yRetrace, xLabel, yLabel] = NSMU.CreateForceVolumeForceCurveZplot(curve_ind, NSMU.METRIC,0);

  xTrace=xTrace-yTrace;
xRetrace=xRetrace-yRetrace;

x_cont=xTrace;
y_cont=yTrace;
[contactPointIndex] = FindContactPoint(x_cont, y_cont);
offset_def= mean(y_cont(1:10));

%   %%%%figures force vs separation 
    figure(f3)
    plot(xTrace-x_cont(contactPointIndex), cte*(yTrace-offset_def),'k-'); hold on
%     plot(xRetrace-x_cont(contactPointIndex), cte*(yRetrace-offset_def),'r-'); hold on  
    xlabel('indentation (nm)');
    ylabel('Force (nN)');
hold on
% plot(0, cte*(y_cont(contactPointIndex)-offset_def),'*c')
hold on

xTrace=xTrace-x_cont(contactPointIndex);
xRetrace=xRetrace-x_cont(contactPointIndex);
end