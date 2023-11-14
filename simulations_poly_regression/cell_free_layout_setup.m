function [RP_positions,TP_positions,AP_positions] = cell_free_layout_setup(RP_positions_per_row,squareLength,L,minDistanceUE2AP,minDistanceAP2AP)
inter_RP_dist = squareLength/(RP_positions_per_row-1); %1000/20 => each RP spaced apart by 50
RP_positions_x = 1:RP_positions_per_row;
RP_positions_y = 1:RP_positions_per_row;
[RP_positions_X,RP_positions_Y] = meshgrid(RP_positions_x,RP_positions_y);
RP_positions_X = (RP_positions_X-1)*inter_RP_dist;
RP_positions_Y = (RP_positions_Y-1)*inter_RP_dist ;
RP_positions = RP_positions_X + 1i*RP_positions_Y;
RP_positions = RP_positions.'; % A.' will transpose a complex matrix. A' will provide complex conjugate. A' is used in functionFindRSS to find combining vector
TP_positions = RP_positions + (inter_RP_dist/2)*(1+1i); %TP placed at centre of the subregion
TP_positions = TP_positions(1:(RP_positions_per_row-1),1:(RP_positions_per_row-1));

%initialize infinite distance to ensure second condition inside for loop is satisfied at start
AP_positions = 2*squareLength*(ones(1,L) + 1i*ones(1,L)); 
%Random AP locations with uniform distribution


for AP_idx = 1:L
    AP_pos_rand = (rand + 1i*rand) * squareLength;
    %if the dist is less than 10m between any two APs or RP positions and any AP then reposition the AP
    while( (min(abs(RP_positions - AP_pos_rand), [], "all") < minDistanceUE2AP) || ...
           (min(abs(AP_positions - AP_pos_rand), [], "all") < minDistanceAP2AP) || ...
           (min(abs(TP_positions - AP_pos_rand), [], "all") < minDistanceUE2AP) ) 
        AP_pos_rand = (rand + 1i*rand) * squareLength;
    end
    AP_positions(AP_idx) =  AP_pos_rand;
end


