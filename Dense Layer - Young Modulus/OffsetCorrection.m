function [OutputDataY] = OffsetCorrection(InputDataX,InputDataY,FitRange,Order)

% Do offset correction using a polynom of order "Order"
% InputData = Vector used for offset correction
% FitRange = region used for fitting

if size(InputDataY,2) == 1
   InputDataY = InputDataY'; 
end
if size(InputDataX,2) == 1
   InputDataX = InputDataX'; 
end

offset = polyfit(InputDataX([FitRange(1):FitRange(2)]),InputDataY([FitRange(1):FitRange(2)]),Order);
% Get reconstructed offset - polynomial parts:
tmp = offset(1:end-1)' .*InputDataX; 
tmp = sum(tmp,1);
% Correct data, including the constant offset:
OutputDataY = InputDataY - offset(end) - tmp;
OutputDataY = sum(OutputDataY,1);

