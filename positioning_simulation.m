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

[ map , map_x_coords , map_y_coords ] = l.power_map ( '3GPP_38.901_UMa_LOS', 'quick' , sample_distance , x_min , x_max , y_min , y_max , rx_height , tx_power , i_freq ) ;
Pow = [];
P_db = [];
total_pow = zeros(size(sum(map{1},4)));

for i = 1:l.no_tx
    Pow(:,:,i) = sum(map{i},4) ;
    P_db(:,:,i) = 10*log10(Pow(:,:,i));
    total_pow = total_pow + Pow(:,:,i);
end


total_P_db = 10*log10(total_pow);
RP_x_coords = unique(l.rx_position(1,:));
RP_y_coords = unique(l.rx_position(2,:));

x_indices = [];
y_indices = [];
for x = 1:length(RP_x_coords)
    x_indices(x) = find(map_x_coords == RP_x_coords(x));
end

for y = 1:length(RP_y_coords)
    y_indices(y) = find(map_y_coords == RP_y_coords(y));
end

[RSS_pos_grid_X, RSS_pos_grid_Y] = meshgrid(RP_x_coords,RP_y_coords);
RSS_pos_grid_X = reshape(RSS_pos_grid_X',l.no_rx,1);   %RSS_pos_grid_X --> 25x1
RSS_pos_grid_Y = reshape(RSS_pos_grid_Y',l.no_rx,1);   %RSS_pos_grid_Y --> 25x1
RSS = []; %RSS --> 5x5x9 vector
total_RSS = []; %total_RSS --> 5x5 vector

for x_pos = 1:length(x_indices)
        for y_pos = 1:length(y_indices)
            x_index = x_indices(x_pos);
            y_index = y_indices(y_pos);
            for tx_pos = 1:l.no_tx
                RSS(x_pos,y_pos,tx_pos) = P_db(x_index,y_index,tx_pos); 
            end
            total_RSS(x_pos,y_pos) = total_P_db(x_index,y_index); 
        end
end

%{
% Plot the results
l.visualize([],[],0);                                   % Show BS and MT positions on the map
hold on; imagesc( map_x_coords, map_y_coords, total_P_db ); hold off  % Plot the antenna footprint
axis([x_min,x_max,y_min,y_max]);
caxis( max(total_P_db(:)) + [-20 0] );                        % Color range
colmap = colormap;
colormap( colmap*0.5 + 0.5 );                           % Adjust colors to be "lighter"
set(gca,'layer','top')                                  % Show grid on top of the map
colorbar('south')
title('Received power [dBm] for 2.4 GHz band')

%}


%test point at (55,45)

x_indx_tp = find(map_x_coords == 75);
y_indx_tp = find(map_y_coords == 65);

RSS_tp = []; %RSS_tp --> 1x9 vector
for tx_pos = 1:l.no_tx
    RSS_tp(tx_pos) = P_db(x_indx_tp,y_indx_tp,tx_pos);
end

%algo
%offline
    no_subrgns = (rx_per_row-1)^2; %a 5x5 matrix has 16 SRs
    RPs_in_a_subrgn = []; % RPs_in_a_subrgn --> 16x4 vector; 4 vertices of 16 subregions
    no_of_RPs_in_a_subrgn = 4;
    for subrgn_idx = 1:no_subrgns
         vertex1 = subrgn_idx + floor((subrgn_idx-1)/4);  %formula -> SR1 has vertices 1,2,6,7; SR5 has vertices 6,7,11,12
         vertex2 = vertex1 + 1;
         vertex3 = vertex1 + rx_per_row;
         vertex4 = vertex2 + rx_per_row;
         RPs_in_a_subrgn(subrgn_idx,1) = vertex1;
         RPs_in_a_subrgn(subrgn_idx,2) = vertex2;
         RPs_in_a_subrgn(subrgn_idx,3) = vertex3;
         RPs_in_a_subrgn(subrgn_idx,4) = vertex4;
    end
    
%online
%APS ==> AP Subset
%RPS ==> RP Subset
RSS_fngprnt = reshape(RSS,l.no_rx,l.no_tx); %RSS_fngprnt --> 25x9 vector

%positioning with seed_ap
[RSS_tp_sorted, RSS_tp_sorted_idx] = sort(RSS_tp,'descend');  %RSS_tp --> 9x1 vector

seed_ap_idx = RSS_tp_sorted_idx(1); %seed_ap_idx --> 1x1 value
seed_ap_tp_rss = RSS_tp_sorted(1); %seed_ap_fngprnt --> 1X1 value
seed_ap_fngprnt = RSS_fngprnt(:,seed_ap_idx);   %seed_ap_fngprnt --> 1X25 vector
rss_diff_with_seed_ap =  abs(seed_ap_fngprnt-seed_ap_tp_rss); %rss_diff_with_seed_ap --> 25x1 vector

best_subrgn = find_best_subrgn(rss_diff_with_seed_ap,RPs_in_a_subrgn,no_subrgns,no_of_RPs_in_a_subrgn);  %index of best subregion
[RPS_RSS,RPS_x,RPS_y] = find_RPS(best_subrgn,seed_ap_fngprnt,RSS_pos_grid_X,RSS_pos_grid_Y,RPs_in_a_subrgn);  %RPS_RSS, RPS_x, RPS_y --> (APSx4) vector; for seed ap --> (1x4) vector
estimated_pos = wknn_positioning(RPS_RSS,seed_ap_tp_rss,1,RPS_x',RPS_y',no_of_RPs_in_a_subrgn);
estimated_pos

%positioning with JAAS
no_AP_in_APS = 3;
prev_itrn_subrgn = 0;
itrn_count = 1;
while (prev_itrn_subrgn ~= best_subrgn)
    itrn_count
    prev_itrn_subrgn = best_subrgn;
    APS = find_APS(best_subrgn,no_AP_in_APS,RPs_in_a_subrgn,RSS_tp,RSS_fngprnt); %1x3 vector
    APS_fngprnt = [];
    for i = 1:length(APS)
        APS_fngprnt(:,i) = RSS_fngprnt(:,APS(i)); %APS_fngprnt --> 25xno_AP_in_APS i.e. 25x3 vector
    end
    APS_RSS_diff_with_tp = abs(APS_fngprnt - RSS_tp(APS)); %RSS_tp --> 1x9 vector; APS_RSS_diff_with_tp--> 25x3 vector
    avg_APS_RSS_diff_with_tp = mean(APS_RSS_diff_with_tp,2); %avg_APS_RSS_diff_with_tp --> 9x1 vector; 2 ==> mean along the rows
    
    best_subrgn = find_best_subrgn(avg_APS_RSS_diff_with_tp,RPs_in_a_subrgn,no_subrgns,no_of_RPs_in_a_subrgn);  %index of best subregion
    [RPS_RSS,RPS_x,RPS_y] = find_RPS(best_subrgn,APS_fngprnt,RSS_pos_grid_X,RSS_pos_grid_Y,RPs_in_a_subrgn);  %RSS_fngprnt, RPS_x, RPS_y --> (APSx4) vector; for seed ap --> (1x4) vector
    estimated_pos = wknn_positioning(RPS_RSS,RSS_tp(APS),no_AP_in_APS,RPS_x',RPS_y',no_of_RPs_in_a_subrgn);
    itrn_count = itrn_count+1;
    if(itrn_count > 5) %if even after 5 iterations convergence of subregion is not reached, then break
        break;
    end
    estimated_pos
    best_subrgn
end

function best_subrgn = find_best_subrgn(rss_diff,RPs_in_a_subrgn,no_subrgns,no_of_RPs_in_a_subrgn) %rss_diff --> 1x25 vector
    %global RPs_in_a_subrgn
    %global no_subrgns
    %global no_of_RPs_in_a_subrgn
    subrgn_rss_sum = zeros(1,no_subrgns);
    for subrgn_idx = 1:no_subrgns
        for vertex = 1:no_of_RPs_in_a_subrgn
            RP_idx_in_subrgn = RPs_in_a_subrgn(subrgn_idx,vertex);
            subrgn_rss_sum(subrgn_idx) = subrgn_rss_sum(subrgn_idx) + rss_diff(RP_idx_in_subrgn);
        end
    end 
    [subrgn_sum_sorted, subrgn_sum_sorted_idx] = sort(subrgn_rss_sum,'ascend');
    best_subrgn = subrgn_sum_sorted_idx(1); 
end

function [RPS_RSS,RPS_x,RPS_y] = find_RPS(best_subrgn,aps_fngprnt,RSS_pos_grid_X,RSS_pos_grid_Y,RPs_in_a_subrgn)
    %global RPs_in_a_subrgn
    RPS = RPs_in_a_subrgn(best_subrgn,:); %RPS --> 4x1 vector of 4 vertices of best subregion
    RPS_RSS = aps_fngprnt(RPS,:);  %RPS_RSS --> 4xAPS vector
    RPS_x = RSS_pos_grid_X(RPS);   % 4x1 vector
    RPS_y = RSS_pos_grid_Y(RPS);   % 4x1 vector 
    % ex: RPS_x = [50;60;50;60] and RPS_y = [60;60;70;70]
    % coordinates of 4 vertices of best subregion are : (50,60),(60,60),(50,70) & (60,70) 
end

function APS = find_APS(best_subrgn,no_AP_in_APS,RPs_in_a_subrgn,RSS_tp,RSS_fngprnt)
    %global RPs_in_a_subrgn % RPs_in_a_subrgn --> 16x4 vector
    %global RSS_tp % RSS_tp --> 1x9 vector of RSS of 9APs at TP
    %global RSS_fngprnt %RSS_fngprnt --> 25x9 vector; RSS values of 9 APs at 25RPs 
    RPS = RPs_in_a_subrgn(best_subrgn,:); %RPS --> 4x1 vector of 4 vertices of best subregion
    ap_fngprnt_for_rps = RSS_fngprnt(RPS,:);  %ap_fngprnt_for_rps --> 4x9 vector; 4=RPS vertices;9=no of APs
    diff_rss = abs(ap_fngprnt_for_rps-RSS_tp);  % diff_rss --> 4x9 vector; 
    avg_diff_rss = mean(diff_rss,1); %avg_diff_rss --> 1x9 vector; mean of all APs at a RP; 1 => mean along the columns
    [avg_diff_rss_sorted, avg_diff_rss_sorted_idx] = sort(avg_diff_rss,'ascend');
    APS = avg_diff_rss_sorted_idx(1,1:no_AP_in_APS);
end

function estimated_pos = wknn_positioning(APS_RSS_of_subrgn_RPs, RSS_tp, no_APs_in_APS,RPS_x,RPS_y,no_of_RPs_in_a_subrgn)
    %global no_of_RPs_in_a_subrgn
    RSS_k = APS_RSS_of_subrgn_RPs; %APS_RSS_of_subrgn_RPs --> 4xAPS; 4x1 for seed AP
    inv_diff_rss = [];
    for m = 1:no_APs_in_APS 
        diff_rss = abs(RSS_k(:,m) - RSS_tp(m)); %diff_rss --> 4x1 vector; diff of 4 RSS values of the RPS with TP RSS values for an AP
        inv_diff_rss(m,:) = 1./diff_rss; %1x4 vector; 
    end
    %inv_diff_rss --> APSx4 vector
    w = [];
    for k = 1:no_of_RPs_in_a_subrgn
        w(k) = mean(inv_diff_rss(:,k)); % 1x4 vector
    end
    sigma_wk = sum(w); %1x1 sum of 4 avg values 
    estimated_pos_x = sum(w.*RPS_x)/sigma_wk;
    estimated_pos_y = sum(w.*RPS_y)/sigma_wk;
    estimated_pos = [estimated_pos_x,estimated_pos_y];
end

