%       This script calculates growth rate, lag time, and geometric lag time
%       using the four independent methods described in the manuscript
%       "Conjugation dynamics depend on both plasmid acquisition cost and
%       fitness cost"
%          - AUTHOR: H PRENSKY
%          - DATE: 2020-12-01
%          - MEHTODS:
%               1. modified logistic method
%               2. modified gompertz method
%               3. baranyi method
%               4. prensky method


clear all, close all, clc

% load demo data for analysis
load sample_data.mat 
% define number of biological replicates
reps = size(OD,1); 
% open a new figure to plot logistic, gompertz, and baranyi fits
figure; hold on

% define collection matrixes, columns correspond to method number
growth_rates = zeros(reps, 4);
lag_times = zeros(reps, 4);
geo_lag_times = zeros(reps, 3);

% loop through each data replicate
for q = 1:reps
    % define y data, which are OD measurements
    ydata = OD(q,:); 
    
    %% Method 1 -- Modified Logistic
    method = 1;
    % plot OD as a function of time
    subplot(3,1,method), hold on, title('Logistic')
    xlabel('Time'), ylabel('OD_6_0_0')
    ydata1 = log(ydata./min(ydata)); %log-transforming the demo data
    plot(time,ydata1,'c','linewidth',4.0)
    
    % calculate growth rate and lag time using modified logistic equation
    %%%% x = [max, growthrate, lagtime] = [0.5 .3 15]
    x0 = [0.5 .3 15];
    % least squares fitting using the modified logistic equation
    x = lsqcurvefit(@growth_modifiedLogistic,x0,time,ydata1);
    ydatafit1 = growth_modifiedLogistic(x,time);
    plot(time,ydatafit1,':k','linewidth',4.0)
    growth_rates(q,method) = x(2); %Calculate growth rate and save to matrix
    lag_times(q,method) = x(3); %Calculate lag time and save to matrix
    % calculate lag time using logistic growth & geometric lag function
    output = calculate_geometric_lagTime(time,ydatafit1);
    geo_lag_times(q, method) = output; %save data to matrix
    
    
    
    %% Method 2 -- Modified Gompertz
    method = 2;
    % plot OD as a function of time
    subplot (3, 1, method), hold on, title('Gompertz')
    xlabel('Time'), ylabel('OD_6_0_0')
    ydata2 = log(ydata./min(ydata)); %log transform the demo data
    plot(time, ydata2, 'c', 'linewidth', 4.0)
    % calculate growth rate and lag time using modified logistic equation
    %%%% x = [max, growthrate, lagtime] = [0.5 .3 15]
    x0 = [0.5, 0.3, 15];
    % least squares fitting using the modified gompertz equation
    x = lsqcurvefit(@growth_modifiedGompertz, x0, time, ydata2);
    ydatafit2 = growth_modifiedGompertz(x, time);
    plot(time, ydatafit2, ':k', 'linewidth', 4.0)
    growth_rates(q,method) = x(2); %Calculate growth rate and save to matrix
    lag_times(q,method) = x(3); %Calculate lag time and save to matrix
    % calculate lag time using logistic growth & geometric lag function
    output = calculate_geometric_lagTime(time,ydatafit2);
    geo_lag_times(q, method) = output; %save data to matrix
    
    
    
    %% Method 3 -- Baranyi
    method = 3;
    % plot OD as a function of time
    subplot (3, 1, method), hold on, title('Baranyi')
    xlabel('Time'), ylabel('OD_6_0_0')
    ydata3 = log(ydata./min(ydata)); %log transform the demo data
    plot(time, ydata3, 'c', 'linewidth', 4.0)
    % calculate growth rate and lag time using modified logistic equation
    %%%% x = [growth rate, lag time, min, max] = [0.5 15 0 1]
    x0 = [0.5, 15, 0, 1];
    x = lsqcurvefit(@growth_baranyi, x0, time, ydata3);
    ydatafit3 = growth_baranyi(x, time);
    plot(time, ydatafit3, ':k', 'linewidth', 4.0)
    growth_rates(q,method) = x(1); %Calculate growth rate and save to matrix
    lag_times(q,method) = x(2); %Calculate lag time and save to matrix
    % calculate lag time using logistic growth & geometric lag function
    output = calculate_geometric_lagTime(time,ydatafit3);
    geo_lag_times(q, method) = output; %save data to matrix
    
    
    %% Method 4 -- Prensky
    method = 4;
    % This method is not a curve-fitting equation, rather, it calculates
    % growth rate and lag time using a tangent line to the maximum slope. This
    % section will calculate growth rate and lag time and save the data to
    % column 4 of the collection matrices
    
    tangent_line = growth_prensky(time, ydata); 
    growth_rates(q, method) = tangent_line(1);
    lag_times(q, method) = (tangent_line(2) / tangent_line(1)) * -1;
end


