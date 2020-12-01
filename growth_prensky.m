function maxGR = growth_prensky(Xdata, Ydata)
 
%set window size to appropriate interval
window = size(Xdata, 2) / max(Xdata);
newYdata = log(Ydata./min(Ydata));
smoothYdata = CalculateSlidingAvg(newYdata,window);

% now we should work with the smoothedYdata to calculate the max GR (max
% derivative)
%approximate derivative based on average increments of x (approx. .25hr)
timeDer = diff(Xdata);
h = mean(timeDer);

%first calculate the 1st derivative of the new Y data
firstDerivative = diff(smoothYdata) / h;

%Smooth this vector
firstDerivativeSmoothed = smoothdata(firstDerivative, 'movmean', h);

%calculate the max and locate it
maxDerivative = max(firstDerivativeSmoothed);
location = find(firstDerivativeSmoothed== maxDerivative);
%create tangent line to area on original curve where max occurs, extract
%slope
if location <= 2
    slope = polyfit(Xdata(location: location + 2), newYdata(location: location + 2), 1);
elseif location == length(Xdata) | location == (length(Xdata) - 1)
    slope = polyfit(Xdata(location - 2: location), newYdata(location - 2: location), 1);
else 
    slope = polyfit(Xdata(location - 2: location + 2), newYdata(location - 2: location + 2), 1);
end 
maxGR = slope;
if maxGR < 1E-10;
    maxGR = 0;
end 

end

