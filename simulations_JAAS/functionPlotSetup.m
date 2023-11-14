function functionPlotSetup(squareLength,RP_positions,AP_positions,TP_POSITIONS)
[X,Y]=meshgrid(0:100:squareLength);
figure; hold on; 
plot(X,Y,":k");
plot(Y,X,":k");
RP_positions_X = real(RP_positions);
RP_positions_Y = imag(RP_positions);
AP_positions_X = real(AP_positions);
AP_positions_Y = imag(AP_positions);
TP_POSITIONS_X = real(TP_POSITIONS);
TP_POSITIONS_Y = imag(TP_POSITIONS);
plot(RP_positions_X,RP_positions_Y,"xb",'MarkerSize',10, 'LineWidth',1.5);
plot(TP_POSITIONS_X,TP_POSITIONS_Y,"xm",'MarkerSize',10, 'LineWidth',1);
plot(AP_positions_X,AP_positions_Y,"^r",'MarkerFaceColor',"r",'MarkerSize',9, 'LineWidth',2);
title('Positioning Simulation Setup');
xlabel('Horizontal Distance (m)')
ylabel('Vertical Distance (m)')