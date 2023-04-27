function [Young,RSquared] = HertzModel(InputDataX,InputDataY,ContactPoint,FitRange,DiameterSample,DiameterIndenter,my,YoungIndenter)
% Input:
% InputDataX/Y = vector of input data
% ContactPoint = contact point of cantilever with cell (actually not needed for calculation)
% FitRange = Range to fit
% DiameterIndenter/Sample = Diameter of indenter/sample ^^
% my = poisson ratio of indented material
% YoungIndenter = Young modulus of indenter

% Reduce input data to region of interest:
InputDataY = InputDataY(FitRange(1):FitRange(2));
InputDataX = InputDataX(FitRange(1):FitRange(2));
% Transform y-data so a linear fit can be used:
InputDataY = InputDataY.^(2/3);
% Set contact point to zero height:
InputDataX = -InputDataX+ContactPoint;
% Make them the same dimension:
InputDataX = InputDataX(:);
InputDataY = InputDataY(:);
% Fit Hertzmodel to data:
[Hertzmodell,G]=fit(InputDataX,InputDataY,'poly1'); 
% Get Young modulus and goodness of fit:
RSquared = G.rsquare;
% Reduced radius (denote: I assume the input is given as a diameter in µm):
RadiusReduced = 0.5*DiameterSample*0.5*DiameterIndenter ./(0.5*DiameterSample + 0.5*DiameterIndenter);
% Calculate Reduced Young`s modulus:
YoungReduced = (0.75*Hertzmodell.p1.^(1.5)) ./ sqrt(RadiusReduced*10^(-6));
% Calculate actual Young´s modulus:
Young = 1./ (1./(YoungReduced*(1-my^2)) - 1/YoungIndenter);