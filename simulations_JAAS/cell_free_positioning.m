rng(0);
diary myDiaryFile25RP_9AP_APS_5_RP_type1_shadowing_noise;
Tempr = 290; %simulation at 290K
noiseFigure = 9; %UT noise figure = 9 dB
k =  1.38065e-20; %k = boltzmann constant in mJ/K
B = 20e6; %communication channel bandwidth = 20MHz
noiseVariancedBm = 10*log10(k*Tempr*B) + noiseFigure;
%center frequency of 2GHz
fc = 2e9;
% N = Number of antennas per AP
N = 4;
sigma_sf = 6;

%Height of BS (in meters)
h_BS = 10;
%Height of UT (in meters)
h_UT = 1.5;
%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

squareLength = 100; %each side of a square simulation area 
%minimum distance between BSs and UEs
minDistanceUE2AP = 5.5; % (d2D > 5.5) => (d3D > 10m) for h_UT = 1m and h_BS = 10m; d3D must be > 10 for the UE to be in Fraunhofer region of the BS 
minDistanceAP2AP = 10;

L = 9 %total number of APs in simulation area
RP_SIM_TYPE = 1 %whether there are 4 RP in RPS or 6 RP
ITRN_COUNT_MAX = 10
nbrOfSetups = 1; %this controls the number of AP-RP layouts i.e. number of setups with random UE and AP locations

RP_positions_per_row = 5
[RP_positions,TP_positions,AP_positions] = cell_free_layout_setup(RP_positions_per_row,squareLength,L,minDistanceUE2AP,minDistanceAP2AP);

functionPlotSetup(squareLength,RP_positions,AP_positions,TP_positions);

K = numel(RP_positions); % K = total number of RPs in the simulation area
%Total uplink transmit power per UE (mW)
p = 100; %power is in mw
beta_fngprnt = zeros(K,L); %25x9
%comments correspond to a 25RP(5x5), 9AP setup

% Go through all RP positions
% For fingerprint database building, only one UE is placed at the RP location and RSS is measured.  

for n = 1:nbrOfSetups
    for AP_idx = 1:L
        for RP_idx = 1:numel(RP_positions)
            %disp(['Running RP' num2str(RP_idx) ' and AP ' num2str(AP_idx)]);
            d_2D = abs(RP_positions(RP_idx) - AP_positions(AP_idx));
            d_3D = sqrt((h_BS-h_UT)^2 + d_2D^2);     
            PL = 35.3*log10(d_3D) + 22.4 + 21.3*log10(fc/1e9);
            beta_fngprnt(RP_idx,AP_idx) = -PL ;%no shadowing in Offline phase owing to spatial averaging in practical considerations
        end %for AP_idx = 1:L 
    end %for RP_idx = 1:numel(RP_positions)

RSS_fngprnt_mW = N*p*db2pow(beta_fngprnt) + N*db2pow(noiseVariancedBm);
RSS_fngprnt_dB = 10*log10(RSS_fngprnt_mW/100); %25x9

num_tp_points = numel(TP_positions);
beta_fngprnt_tp = zeros(num_tp_points,L); %16x9 matrix
PL_fngprnt_tp = zeros(num_tp_points,L); %16x9 matrix
shadowing = zeros(num_tp_points,L);
for AP_idx = 1:L
    for TP_idx = 1:num_tp_points
        shadowing(TP_idx,AP_idx) = sigma_sf*randn;
        d_2D = abs(TP_positions(TP_idx)-AP_positions(AP_idx));
        d_3D = sqrt((h_BS-h_UT)^2 + d_2D^2);     
        PL   = 35.3*log10(d_3D) + 22.4 + 21.3*log10(fc/1e9);
        PL_fngprnt_tp(TP_idx,AP_idx) = -PL;
        beta_fngprnt_tp(TP_idx,AP_idx) = -PL + shadowing(TP_idx,AP_idx); %shadowing noise present in online phase as it can't be eliminated due to unavailability of spatial info of the UE
    end %for AP_idx = 1:L
end %TP_idx = 1:num_tp_points

RSS_tp_mW = N*p*db2pow(beta_fngprnt_tp) + N*db2pow(noiseVariancedBm); %RSS_tp_mW --> 16x9 vector
RSS_tp_dB = 10*log10(RSS_tp_mW/100); %16x9, with TPs in the centre of the 16 subregions
%RSS_tp_dB = pow2db(RSS_tp_mW);
[RPs_in_a_subrgn, no_of_RPs_in_a_subrgn, no_subrgns] = OfflineAlgo.functionOfflineSetup(RP_positions_per_row, RP_SIM_TYPE);
cluster_size = 3 % different sizes of clusters were considered for offline k-means; 
AP_cluster_idx = OfflineAlgo.cluster_AP_kmeans(AP_positions,cluster_size);
% RPs_in_a_subrgn : 16x4, 16 subregions, 4 RP numbers in each subregion; 24x6 , 24 subregions, 6 RPs in each subregion for RP_SIM_TYPE = 2
% no_of_RPs_in_a_subrgn = 4
% no_subrgns = 16
itrn_check = zeros(1,num_tp_points); %variable to monitorhow many iterations it takes for algorithm termination for each TP 
best_subrgn_seed = zeros(1,num_tp_points); %variable to monitor how the seed subrgn fares wrt the final selected subregion
for TP_idx = 1:num_tp_points % num_tp_points=16
    per_tp_RSS = RSS_tp_dB(TP_idx,:);
    %positioning with seed_ap
    [RSS_tp_sorted, RSS_tp_sorted_idx] = sort(per_tp_RSS,'descend');  %per_tp_RSS --> 9x1 vector of the TP; RSS_tp_sorted --> 1x9 vector; RSS_tp_sorted_idx --> 1x9 vector
    seed_ap_idx = RSS_tp_sorted_idx(1); %seed_ap_idx --> 1x1 value
    seed_ap_rss_at_tp = RSS_tp_sorted(1); %seed_ap_rss_at_tp --> 1X1 value
    seed_ap_fngprnt = RSS_fngprnt_dB(:,seed_ap_idx);   %seed_ap_fngprnt --> 25x1 vector
    rss_diff_with_seed_ap =  abs(seed_ap_fngprnt - seed_ap_rss_at_tp); %rss_diff_with_seed_ap --> 25x1 vector
    
    best_subrgn_seed(TP_idx) = OnlineAlgo.find_best_subrgn(rss_diff_with_seed_ap,RPs_in_a_subrgn,no_subrgns,no_of_RPs_in_a_subrgn);  %index of best subregion
    [RPS_RSS,RPS_x,RPS_y] = OnlineAlgo.find_RPS_RSS(best_subrgn_seed(TP_idx),seed_ap_fngprnt,RPs_in_a_subrgn,RP_positions);  %RPS_RSS --> 4xAPS vector, RPS_x --> 4x1 vector, RPS_y --> 4x1 vector; 
    %RPS_RSS --> (RPs_in_a_subrgn x APS) vector ex 4x3 for3APs in APs; for seed ap --> (RPs_in_a_subrgnx1) vector ex 4x1 vector
    estimated_pos_seed = OnlineAlgo.wknn_positioning(RPS_RSS,seed_ap_rss_at_tp,1,RPS_x,RPS_y,no_of_RPs_in_a_subrgn);
    TP_COORDS_seed = [real(TP_positions(TP_idx)) imag(TP_positions(TP_idx))];
    dist_coords = [estimated_pos_seed;TP_COORDS_seed];
    estimated_err_seed(TP_idx) = pdist(dist_coords,'euclidean');
    
    best_subrgn = best_subrgn_seed(TP_idx);
    itrn_count = 1;
    prev_itrn_subrgns = zeros(1,ITRN_COUNT_MAX);
    while (~ismember(best_subrgn,prev_itrn_subrgns)) %do JAAS until you get the same subregion
        prev_itrn_subrgns(itrn_count) = best_subrgn;
        APS = OnlineAlgo.find_APS(best_subrgn,RPs_in_a_subrgn,per_tp_RSS,RSS_fngprnt_dB,AP_cluster_idx,cluster_size); %1x3 vector
        APS_fngprnt = [];
        for i = 1:length(APS)
            APS_fngprnt(:,i) = RSS_fngprnt_dB(:,APS(i)); %APS_fngprnt --> 25xno_AP_in_APS i.e. 25x3 vector for APS=3
        end
        APS_RSS_diff_with_tp = abs(APS_fngprnt - per_tp_RSS(APS)); %per_tp_RSS --> 1x9 vector; per_tp_RSS(APS) --> 1x3 vector; APS_RSS_diff_with_tp--> 25x3 vector
        best_subrgn = OnlineAlgo.find_best_subrgn(APS_RSS_diff_with_tp,RPs_in_a_subrgn,no_subrgns,no_of_RPs_in_a_subrgn);  %index of best subregion
        itrn_count = itrn_count+1;
        if(itrn_count > ITRN_COUNT_MAX) %if even after 5 iterations convergence of subregion is not reached, then break
            break;
        end        
    end %while (~ismember(best_subrgn,prev_itrn_subrgns))
    itrn_check(TP_idx) = itrn_count;
    [RPS_RSS,RPS_x,RPS_y] = OnlineAlgo.find_RPS_RSS(best_subrgn,APS_fngprnt,RPs_in_a_subrgn,RP_positions);  %RPS_RSS, RPS_x, RPS_y --> (APSx4) vector; for seed ap --> (1x4) vector
    estimated_pos = OnlineAlgo.wknn_positioning(RPS_RSS,per_tp_RSS(APS),length(APS),RPS_x,RPS_y,no_of_RPs_in_a_subrgn);
    TP_COORDS = [real(TP_positions(TP_idx)) imag(TP_positions(TP_idx))];
    dist_coords = [estimated_pos;TP_COORDS];
    estimated_err(TP_idx) = pdist(dist_coords,'euclidean');
       
end %for TP_idx = 1:num_tp_points
fprintf("err_avg = %d\n",mean(estimated_err));
fprintf("err_avg_seed = %d\n",mean(estimated_err_seed));
fprintf("itrn_avg = %d\n",mean(itrn_check));

end %for n = 1:nbrOfSetups
diary off

