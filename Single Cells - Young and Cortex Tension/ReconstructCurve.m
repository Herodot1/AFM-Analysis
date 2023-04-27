function [ReconstructedCurve] =  ReconstructCurve(InputDataY)
% Reconstructs curves - with less noise - using PCA.
% Input Y is a vector.
if size(InputDataY,1)>size(InputDataY,2)
    InputDataY = InputDataY';
end

% Do PCA: 
[coeff,score,latent,tsquared,explained] = pca(InputDataY','Centered','off');



% Evaluate the number of necessary 
% Take only these components explaining 99.9% of variance
CumExp = cumsum(explained);
IdxParam = find(CumExp >=99.99,1);

% Fo all possible combinations determine some key paramters, describing the
% deviation between all individual curves, their reconstructed counterparts
% and the remainder of the components.
% The idea behind is as follows: The chosen PCA components should allow for
% good reconstruction, while the remainders should ideally only contain
% noise.

ReconstructedCurve = zeros(size(InputDataY));
ReconstructedCurve = ReconstructedCurve';
Residuals = NaN(IdxParam,1);
MeanNoise = Residuals;
StdNoise = Residuals;
KurtosisNoise = Residuals;
SkewnessNoise = Residuals;
for CompNum = 1:IdxParam 
    ReconstructedCurve = ReconstructedCurve+score(:,CompNum).*coeff(:,CompNum)';       
    % Root mean squared error:
    %RMSE = sqrt((1./size(CroppedDataX,2)).*sum(abs(CroppedDataY-Reconstructed').^2,1));
    % Residuals -> pearson residuals:
    Residuals(CompNum) = sum(sum(   (InputDataY-ReconstructedCurve')));
    tmpNoise = zeros(size(InputDataY));
    tmpNoise = tmpNoise';
    for j = CompNum+1:size(InputDataY,1)
        tmpNoise = tmpNoise + score(:,j).*coeff(:,j)';
    end
    % Mean and std of higher order components:
    MeanNoise(CompNum) = mean(tmpNoise(:));
    StdNoise(CompNum) = std(tmpNoise(:));
    % Skewness and curtosis of higher order parameters:
    SkewnessNoise(CompNum) = skewness(tmpNoise(:),0);
    KurtosisNoise(CompNum) = kurtosis(tmpNoise(:),0);
end
% Normalize residuals to the number of measurements, as previously we have
% taken the sum:
Residuals = Residuals./size(InputDataY,2);
% Differential residuals:
DiffRes = -(diff(abs(Residuals)));
DiffMeanNoise = -(diff(abs(MeanNoise)));
DiffMeanStd = -(diff(abs(StdNoise)));
DiffKurtosis = -(diff(abs(KurtosisNoise)));
DiffSkewness = -(diff(abs(SkewnessNoise)));
% Find the adequate number of PCA components.
% Get break points in curve:
Idx = NaN(5,1);
met = 'linear';
TF=ischange(DiffRes,met,'MaxNumChanges',1); % Using ischange to find the "abrupt changes" in the graph
if sum(TF) > 0
    Idx = 1 + find(TF == 1,1); % +1 due to the determination via "diff" -> gives an offset of 1
end
TF=ischange(DiffMeanNoise,met,'MaxNumChanges',1);
if sum(TF) > 0
    Idx(2) = 1 + find(TF == 1,1);
end
TF=ischange(DiffMeanStd,met,'MaxNumChanges',1);
if sum(TF) > 0
    Idx(3) = 1 + find(TF == 1,1);
end
TF=ischange(DiffKurtosis,met,'MaxNumChanges',1);
if sum(TF) > 0
    Idx(4) = 1 + find(TF == 1,1);
end
TF=ischange(DiffSkewness,met,'MaxNumChanges',1);
if sum(TF) > 0
    Idx(5) = 1 + find(TF == 1,1);
end
Idx(Idx == 0) = [];

% Take the mean and round to take the final number of components:
Idx = ceil(nanmean(Idx));

% if only one parameter is needed anyway, all the above does not make
% really sense, thus just take the first component:
if IdxParam == 1 && sum(isnan(Idx)) == 5
   Idx = 1; 
end

% Reconstruct all curves using the said number of pca components:
% Reconstructed/denoised curves:

ReconstructedCurve = zeros(size(InputDataY));
ReconstructedCurve = ReconstructedCurve';
for CompNum = 1:Idx
    ReconstructedCurve = ReconstructedCurve+score(:,CompNum).*coeff(:,CompNum)';       
end