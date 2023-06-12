set (0 ,'defaultTextFontSize', 18) % Default Font Size
set (0 ,'defaultAxesFontSize', 18) % Default Font Size
set (0 ,'defaultAxesFontName','Times') % Default Font Type
set (0 ,'defaultTextFontName','Times') % Default Font Type
set (0 ,'defaultFigurePaperPositionMode','auto') % Default Plot position
set (0 ,'DefaultFigurePaperType','< custom >') % Default Paper Type
set (0 ,'DefaultFigurePaperSize',[14.5 7.7]) % Default Paper Size


% Create an environment with free-space path loss
s = qd_simulation_parameters ;
s.center_frequency = 2.4e9;

l = qd_layout(s);
l.no_tx = 9; 
tx_x_offset = 5;
tx_y_offset = 40;
tx_per_row = sqrt(l.no_tx);
l.tx_position(3,:,1,1) = 50; %50 m height of bs
for x = 1:tx_per_row
     for y = 1:tx_per_row
        l.tx_position(1,(x-1)*tx_per_row+y,1,1) = tx_x_offset + (x-1)*10 ;
        l.tx_position(2,(x-1)*tx_per_row+y,1,1) = tx_y_offset + (y-1)*10 ;
     end
end
tx_ant_2600_Mhz = qd_arrayant( '3gpp-3d', 1 , 1 , s.center_frequency(1) , 3) ;
l.tx_array(1,:) = tx_ant_2600_Mhz;


l.no_rx = 25; 
rx_x_offset = 40;
rx_y_offset = 30;
rx_per_row = sqrt(l.no_rx);
for x = 1:rx_per_row
     for y = 1:rx_per_row
        l.rx_position(1,(x-1)*rx_per_row+y,1,1) = rx_x_offset + (x-1)*10 ;
        l.rx_position(2,(x-1)*rx_per_row+y,1,1) = rx_y_offset + (y-1)*10 ;
     end
end
l.rx_array(1,:) = qd_arrayant('omni') ; 

l.rx_position(3,:,1,1) = 1.5; %50 m height of bs
indoor_rx = l.set_scenario ('3GPP_38.901_UMa' ,[] ,[] ,0) ; % Set the scenario
l.rx_position (3, ~indoor_rx ) = 1.5 ;  % Set outdoor - users to 1.5 m height


sample_distance = 5;
x_min = -20; % Area to be samples in [ m ]
x_max = 100;
y_min = -20;
y_max = 100;
rx_height = 1.5; % Mobile terminal height in [ m ]
tx_power = 0; % Tx - power in [ dBm ] per antenna element
i_freq = 1; 

[ map , x_coords , y_coords ] = l.power_map ( '3GPP_38.901_UMa_LOS', 'quick' , sample_distance , x_min , x_max , y_min , y_max , rx_height , tx_power , i_freq ) ;
Pow = [];
P_db = [];
total_pow = zeros(size(sum(map{i},4)));

for i = 1:l.no_tx
    Pow(:,:,i) = sum ( map {i} , 4 ) ;
    P_db(:,:,i) = 10*log10(Pow(:,:,i));
    total_pow = total_pow + Pow(:,:,i);
end


total_P_db = 10*log10(total_pow);
x_pos = unique(l.rx_position(1,:));
y_pos = unique(l.rx_position(2,:));

x_indices = [];
y_indices = [];
for x = 1:length(x_pos)
    x_indices(x) = find(x_coords == x_pos(x));
end

for y = 1:length(y_pos)
    y_indices(y) = find(y_coords == y_pos(y));
end

RSS = [];

for x_pos = 1:length(x_indices)
        for y_pos = 1:length(y_indices)
            x_index = x_indices(x_pos);
            y_index = y_indices(y_pos);
            for tx_pos = 1:l.no_tx
                RSS(x_pos,y_pos,tx_pos) = P_db(x_index,y_index,tx_pos);
            end
        end
end

% Plot the results
l.visualize([],[],0);                                   % Show BS and MT positions on the map
hold on; imagesc( x_coords, y_coords, total_P_db ); hold off  % Plot the antenna footprint
axis([x_min,x_max,y_min,y_max]);
caxis( max(total_P_db(:)) + [-20 0] );                        % Color range
colmap = colormap;
colormap( colmap*0.5 + 0.5 );                           % Adjust colors to be "lighter"
set(gca,'layer','top')                                  % Show grid on top of the map
colorbar('south')
title('Received power [dBm] for 2.4 GHz band')



