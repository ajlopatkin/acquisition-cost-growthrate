function output = calculate_geometric_lagTime(time,ODfit)

ind_maxU = find(diff(ODfit)./diff(time) == max(diff(ODfit)./diff(time)));
if length(ind_maxU) > 1
    ind_maxU = ind_maxU(1);
end

ind_maxU2 = [ind_maxU-1 ind_maxU+1];
if ind_maxU2(1) == 0
    ind_maxU2 = ind_maxU2 + 1;
end

time_line = time(ind_maxU2);
OD_line = ODfit(ind_maxU2);
slope = (OD_line(2) - OD_line(1)) / (time_line(2) - time_line(1));

y = slope .* (time - time(ind_maxU)) +  ODfit(ind_maxU);
hold on
plot(time,y)

ind_LT = find(y>=0,1,'first');
output = time(ind_LT);
scatter(output,y(ind_LT),100,'filled','markerfacecolor',[.6 .6 .6], 'markeredgecolor','k')

end