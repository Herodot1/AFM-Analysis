clear all;

% Parameters:
% Number of points from the force distance curve used (starting the count 
% from maximal indentation)
NumPoints = 600;
% Order of the polynomial for offset correction:
Order = 1;
% Starting point for the offset correction (counting/direction as in 
% NumPoints). Needs to be smaller than NumPoints.
startOffset = 599;
% Diameter of indenter (in µm):
DiameterIndenter = 10^15;
% Young modulus of indenter in Pa:
YoungIndenter = 10^15;
% Fit range for the Hertz model: (in values of Force -> fitting of force x
% nN to y nN
FitRangeYoung = [0.3,1];
% Poisson ratio of sample:
my = 0.5;
% Fit range for the tension model: (in values of Force -> fitting of force x
% nN to y nN
FitRangeTension = [0.1,0.2];

% get current path and add it to the search path of matlab for additional
% .m -files
mPath = mfilename('fullpath');
Idx = max(strfind(mPath,filesep));
mPath = mPath(1:Idx);
addpath(mPath)

% Get data path:
Startpath=mPath;
FolderPath=uigetdir(Startpath, 'Chose the folder with the images');
cd(FolderPath);   
% Get all file names:
FileList = dir();
FileList([FileList.isdir]==1) = [];

% read in file with name "radius", containing the cells diameter (stupid me
% ...).
Diameter=dlmread('radius.txt');   
% load spring constant value:
SpringConstant = dlmread('k.txt');

% read in data files:
% as files are named as numbers of type  1.txt,2.txt, etc, determine n for
% reading. "-2" as k.txt and radius.txt are stored there as well.
n = size(FileList,1) -2;
Data = cell(n,1);
for k=1:n    
    Data{k}=dlmread(sprintf('%d.txt',k),'\t',1,0); %liest die textfiles 1 bis n aus und gibt ihnen die matrix data{k}    
end

% As the data might not have the same number of data points, crop it to
% a uniform size, removing some parts of the top part of the curve, above
% the sample. Afterwards, due a first offset correction.
CroppedDataY = NaN(NumPoints,n);
CroppedDataX = NaN(NumPoints,n); 
IdxOut = [];
count = 0;
% Curves to be kept:
KeptCurves = 1:size(Data,1);
for i = 1:size(Data,1)
    tmpX = Data{i}(:,1);
    tmpY = Data{i}(:,2);
    %Remove zeros:
    Removal = find(tmpY == 0);
    tmpX(Removal) = [];
    tmpY(Removal) = [];
    % corrected x:
    tmpX = abs(tmpX-max(tmpX));
    % adjust too high values (propably given in pm instead of nm:
    if  max(tmpX)-min(tmpX)  > 300
        tmpX = tmpX/1000;
    end
    
    
    if size(tmpY,1)>= size(CroppedDataY,1)
        % Get break points in curve:
        TF=ischange(tmpY(end-size(CroppedDataY,1)+1:end),'linear','MaxNumChanges',1); % Using ischange to find the "abrupt changes" in the graph
        % old:
        % TF=ischange(tmpY,'linear','MaxNumChanges',1); % Using ischange to find the "abrupt changes" in the graph
        Idx = find(TF == 1) -10;
        % offset correction:
        % s = 5;
        % [~,startOffset] = min(abs(tmpX - s))
        [tmpYCorrected] = OffsetCorrection(tmpX(end-size(CroppedDataY,1)+1:end),tmpY(end-size(CroppedDataY,1)+1:end),[size(CroppedDataY,1) - startOffset,Idx],Order);
        % old:
        % [tmpYCorrected] = OffsetCorrection(tmpX,tmpY,[length(tmpY) - startOffset,Idx],Order);
        % adjust too high values (propably given in pN instead of nN:
        if  max(tmpYCorrected)-min(tmpYCorrected)  > 100
            tmpYCorrected = tmpYCorrected/1000;
        end
        % Crop data:
        CroppedDataY(:,i) = tmpYCorrected(end-size(CroppedDataY,1)+1:end);
        CroppedDataX(:,i) = tmpX(end-size(CroppedDataY,1)+1:end);        
    else
        % If too few points, store index for later removal:
        count = count+1;
        IdxOut(count) = i;
    end
end
% Remove columns corresponding to curves with too low point number (lower
% than NumPoints)
KeptCurves(IdxOut) = [];
CroppedDataY(:,IdxOut) = [];
CroppedDataX(:,IdxOut) = [];
% Also delete associated diameters:
Diameter(IdxOut) = [];

% Denoise data by reconstructing a "noise free" curve using PCA:
[Reconstructed] =  ReconstructCurve(CroppedDataY');
Reconstructed = Reconstructed';
figure; plot(CroppedDataY)
figure; plot(Reconstructed')

% Calculate residuals of reconstructed and original curve:
Residuals = sum(abs((CroppedDataY-Reconstructed')),1);
% Find potential outliers, pointing to bad reconstructions, that might
% imply poor curves:
BadReconstructions =(Residuals> median(Residuals) + 3*mad(Residuals));
CroppedDataY(:,BadReconstructions) = [];
CroppedDataX(:,BadReconstructions) = [];
KeptCurves(BadReconstructions) = [];
% Repeat reconstruction step to improve reconstruction of remaining curves:
[Reconstructed] =  ReconstructCurve(CroppedDataY');
Reconstructed = Reconstructed';


% now use reconstructed curves to refine offset correction and get better
% fit values for offset, Young modulus and cortex tension determination:
% CortexThickness determination does currently not work particularily well,
% as I use for the Youngs modulus the one determined from the larger force
% spectrum. That is what I suppose at least. 
CorrectedReconstructed = NaN(size(Reconstructed));
% Allocate variables:
MaxIndentation = NaN(size(Reconstructed,1),1);
IndendationMaxFitForce = MaxIndentation;
Young = MaxIndentation;
RSquaredHertz = MaxIndentation;
Tension = MaxIndentation;
CortexThickness = MaxIndentation;
RSquaredTension = MaxIndentation;
RSquaredCortexThickness = MaxIndentation;

for i = 1:size(Reconstructed,1)
    % Current data set:
    tmpY = Reconstructed(i,:);
    % Range used for fitting:
    TF=ischange(tmpY,'variance','MaxNumChanges',1); % Using ischange to find the "abrupt changes" in the graph
    %endOffset = find(TF == 1) -150;
    %[~,startOffset] = min(abs(CroppedDataX - s));
    %startOffset = 1;
    %endOffset = startOffset+400;
    endOffset = find(TF == 1);
    [tmpYCorrected] = OffsetCorrection(CroppedDataX(:,i),tmpY,[length(tmpY) - startOffset,endOffset],Order);
    CorrectedReconstructed(i,:) = tmpYCorrected;     
    % Do Young modulus calculation:
    % first get last point where force is below zero -> this is probably
    % quite close to the contact point. Add 5 to be sure
    ContactPoint = 1 + find(tmpYCorrected<0,1,'last') + 5;
    % Maximal Indentation:
    MaxIndentation(i) = CroppedDataX(size(tmpYCorrected,2),i)-CroppedDataX(ContactPoint,i);   
    
    % Find fit range for given force interval:
    [~,startFitHertz] = min(abs(FitRangeYoung(1)-tmpYCorrected));
    [~,endFitHertz] = min(abs(FitRangeYoung(2)-tmpYCorrected)); 
    % Indentation for maximal fit force:
    IndendationMaxFitForce(i) = CroppedDataX(endFitHertz,i)-CroppedDataX(ContactPoint,i);
    % Contact point in µm:
    ContactPoint = CroppedDataX(ContactPoint,i);
    % Hertz fit:
    [tmpYoung,tmpRSquared] = HertzModel(CroppedDataX(:,i),tmpYCorrected,ContactPoint,[startFitHertz,endFitHertz],Diameter(i),DiameterIndenter,my,YoungIndenter);
    Young(i) = tmpYoung;
    RSquaredHertz(i) = tmpRSquared;      
    
    % Find fit range for tension model::
    [~,startFitTension] = min(abs(FitRangeTension(1)-tmpYCorrected));
    [~,endFitTension] = min(abs(FitRangeTension(2)-tmpYCorrected));
    FitRange = [startFitTension,endFitTension];
    % Calculate tension and cortex thickness. cortex thickness is probably
    % way too low:
    [tmpTension,tmpCortexThickness,tmpRSquaredTension,tmpRSquaredCortexThickness] = ...
        SurfaceTension(CroppedDataX(:,i),tmpYCorrected,ContactPoint,SpringConstant,Young(i),Diameter(i),FitRange);
    Tension(i) = tmpTension;
    CortexThickness(i) = tmpCortexThickness;
    RSquaredTension(i) = tmpRSquaredTension;
    RSquaredCortexThickness(i) = tmpRSquaredCortexThickness;
end
%figure; plot(CorrectedReconstructed')
% check for linearity of regime for hertz fit:
[Coefficient] = corrcoef(Young,abs(IndendationMaxFitForce).^(2/3))
nanmean(Young)
nanmean(Tension)


% save data:
folder = fullfile(strcat(FolderPath,'\Results')); 
% Checks the existence of the folder with the name "Results". If it
% does not exist it is created.
if (exist(folder) == 0)
    mkdir('Results');   
end
cd(folder)

%Save also as .mat file
save('Results.mat','Diameter','Young','Tension','RSquaredTension','RSquaredHertz','IndendationMaxFitForce','KeptCurves');
save('Parameters.mat','NumPoints','Order','startOffset','DiameterIndenter','YoungIndenter','FitRangeYoung','my','FitRangeTension')

a
% some plots:
% General plot:
[p,anovatab,stats] = anova1(youngAll');
[c,m,h,nms] = multcompare(stats,'display','on');
[p,anovatab,stats] = anova1(tensionAll');
[c,m,h,nms] = multcompare(stats,'display','on');
[p,anovatab,stats] = anova1(radiusAll');
[c,m,h,nms] = multcompare(stats,'display','on');

savename = 'Young U138 LN229 Glucose Deprivation';
% shape factor over density:
h=figure('units','normalized','outerposition',[0 0 1 1]);
set(h,'DefaultTextFontSize',24)
hold all
h = boxplot((youngAll'),{'U138 4.5g/l', 'U138 0.1g/l', 'LN229 2g/l', 'LN229 0.4g/l'},'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
%plot([3.81,3.81],[0.5,9.5],'--','LineWidth',2,'Color',[0.3,0.3,0.3])
set(h,{'linewidth'},{2})
xlabel('Young Modulus in Pa','FontSize',24)
set(gca,'fontsize', 24);
title('Young Modulus','FontSize',30);
saveas(gcf, sprintf('%s.png',savename));
savefig(sprintf('%s.fig',savename))
close all


savename = 'Tension U138 LN229 Glucose Deprivation';
% shape factor over density:
h=figure('units','normalized','outerposition',[0 0 1 1]);
set(h,'DefaultTextFontSize',24)
hold all
h = boxplot(10^6*tensionAll',{'U138 4.5g/l', 'U138 0.1g/l', 'LN229 2g/l', 'LN229 0.4g/l'},'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
%plot([3.81,3.81],[0.5,9.5],'--','LineWidth',2,'Color',[0.3,0.3,0.3])
set(h,{'linewidth'},{2})
xlabel('Tension in pN/µm','FontSize',24)
set(gca,'fontsize', 24);
title('Cortex Tension','FontSize',30);
saveas(gcf, sprintf('%s.png',savename));
savefig(sprintf('%s.fig',savename))
close all

savename = 'Diameter U138 LN229 Glucose Deprivation';
% shape factor over density:
h=figure('units','normalized','outerposition',[0 0 1 1]);
set(h,'DefaultTextFontSize',24)
hold all
h = boxplot(radiusAll',{'U138 4.5g/l', 'U138 0.1g/l', 'LN229 2g/l', 'LN229 0.4g/l'},'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
%plot([3.81,3.81],[0.5,9.5],'--','LineWidth',2,'Color',[0.3,0.3,0.3])
set(h,{'linewidth'},{2})
xlabel('Cell diameter in µm','FontSize',24)
set(gca,'fontsize', 24);
title('Cell Diameter','FontSize',30);
saveas(gcf, sprintf('%s.png',savename));
savefig(sprintf('%s.fig',savename))
close all








