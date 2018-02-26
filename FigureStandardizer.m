function FigureStandardizer(xlhand, ylhand, axishand)

%use this in code to set axes size/font and remove box
% FigureStandardizer(get(gca,'xlabel'), get(gca,'ylabel'), gca);

AxisLabelSize = 22;
AxisTickSize = 20;


set(xlhand,'fontsize',AxisLabelSize,'FontName','Arial');

set(ylhand,'fontsize',AxisLabelSize,'FontName','Arial');

set(axishand, 'fontsize',AxisTickSize,'FontName','Arial');

box off

end