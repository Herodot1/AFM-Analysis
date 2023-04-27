% Berechnung der Oberflächenspannung und des Young Moduls nach:
% Actomyosin Cortical Mechanical Properties in Nonadherent Cells Determined by Atomic Force Microscopy
% Alexander X.Cartagena-Rivera1Jeremy S.Logue12Clare M.Waterman2Richard S.Chadwick1

function [Tension,CortexThickness,RSquaredTension,RSquaredCortexThickness] = SurfaceTension(InputDataX,InputDataY,ContactPoint,SpringConstant,YoungSample,DiameterSample,FitRange)


% Reduce input data to region of interest:
InputDataY = InputDataY(FitRange(1):FitRange(2));
InputDataX = InputDataX(FitRange(1):FitRange(2));
% Transform y-Data into true deflection in µm, given that force is in nN
% and k in N/m:
InputDataY =InputDataY./(SpringConstant*10^3);
% Set contact point to zero height:
InputDataX = -InputDataX+ContactPoint;
% Make them the same dimension:
InputDataX = InputDataX(:);
InputDataY = InputDataY(:);

% Fit tension model:
[TensionModel,G]=fit(InputDataY ,InputDataX,'poly1');  
RSquaredTension = G.rsquare;
% in N/m (unit of spring constant k):
Tension = SpringConstant/(pi*(TensionModel.p1-1));


% Estimate Cortex Size from input young modulus, etc. in µm.
% Both inputs given in µm;
TensionPlot = (SpringConstant/pi)./((InputDataX./InputDataY)-1);
[CortexThicknessModel,G]=fit(TensionPlot ,sqrt(InputDataY),'poly1'); 
RSquaredCortexThickness = G.rsquare;
% Cortex thickness estimate in mm:
CortexThickness = pi*0.5*DiameterSample./(2*SpringConstant*YoungSample*(CortexThicknessModel.p1^2)); 
% transform to µm;
CortexThickness = CortexThickness * 10^3;
% seems unreasonably low ...


