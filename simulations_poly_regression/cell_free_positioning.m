rng(3);
diary myDiaryFile_poly_regr_high_noise_TP9_n2;
Tempr = 290; %simulation at 290K
noiseFigure = 9; %UT noise figure = 9 dB
k =  1.38065e-20; %k = boltzmann constant in mJ/K
B = 20e6; %communication channel bandwidth = 20MHz
noiseVariancedBm = 10*log10(k*Tempr*B) + noiseFigure;
%center frequency of 2GHz
fc = 2e9;
% N = Number of antennas per AP
N = 4;
sigma_sf = 5;

%Height of BS (in meters)
h_BS = 10;
%Height of UT (in meters)
h_UT = 1.5;
%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

squareLength = 100; %each side of a square simulation area 
%minimum distance between BSs and UEs
minDistanceUE2AP = 5.5; % (d2D > 5.5) => (d3D > 10m) for h_UT = 1m and h_BS = 10m
minDistanceAP2AP = 10;

L = 9 %total number of APs in simulation area
RP_SIM_TYPE = 1 %whether there are 4 RP in RPS or 6 RP
ITRN_COUNT_MAX = 10
nbrOfSetups = 1; %this controls the number of AP-RP layouts i.e. number of setups with random UE and AP locations

RP_positions_per_row = 5
[RP_positions,TP_positions,AP_positions] = cell_free_layout_setup(RP_positions_per_row,squareLength,L,minDistanceUE2AP,minDistanceAP2AP);

functionPlotSetup(squareLength,RP_positions,AP_positions,TP_positions);

K = numel(RP_positions);
%Total uplink transmit power per UE (mW)
p = 100;
beta_fngprnt = zeros(K,L); %25x9
d_2D = zeros(K,L); %25x9
d_3D = zeros(K,L); %25x9
%comments correspond to a 25 (5x5)RP, 9AP setup

% Go through all RP positions
% For fingerprint database building, only one UE is placed at the RP location and RSS is measured.  

for n = 1:nbrOfSetups
    for AP_idx = 1:L
        for RP_idx = 1:numel(RP_positions)
            %disp(['Running RP' num2str(RP_idx) ' and AP ' num2str(AP_idx)]);
            d_2D(RP_idx,AP_idx) = abs(RP_positions(RP_idx) - AP_positions(AP_idx));
            d_3D(RP_idx,AP_idx) = sqrt((h_BS-h_UT)^2 + d_2D(RP_idx,AP_idx)^2);     
            PL = 35.3*log10(d_3D(RP_idx,AP_idx)) + 22.4 + 21.3*log10(fc/1e9);
            beta_fngprnt(RP_idx,AP_idx) = -PL ;%no shadowing in Offline phase owing to spatial averaging in practical considerations
        end %for AP_idx = 1:L 
    end %for RP_idx = 1:numel(RP_positions)

RSS_fngprnt_mW = N*p*db2pow(beta_fngprnt) + N*db2pow(noiseVariancedBm);
RSS_fngprnt_dB = 10*log10(RSS_fngprnt_mW/100); %25x9

fngprnt_poly = polynomial_regression(RSS_fngprnt_dB,d_3D);

num_tp_points = numel(TP_positions);
beta_fngprnt_tp = zeros(num_tp_points,L); %16x9 matrix
PL_fngprnt_tp = zeros(num_tp_points,L); %16x9 matrix
shadowing = zeros(num_tp_points,L);
d_2D_TP = zeros(1,num_tp_points);
for AP_idx = 1:L
    for TP_idx = 1:num_tp_points
        shadowing(TP_idx,AP_idx) = sigma_sf*randn;
        d_2D_TP(TP_idx) = abs(TP_positions(TP_idx)- AP_positions(AP_idx));
        d_3D_TP = sqrt((h_BS-h_UT)^2 + d_2D_TP(TP_idx)^2);     
        PL   = 35.3*log10(d_3D_TP) + 22.4 + 21.3*log10(fc/1e9);
        PL_fngprnt_tp(TP_idx,AP_idx) = -PL;
        beta_fngprnt_tp(TP_idx,AP_idx) = -PL + shadowing(TP_idx,AP_idx);
    end %for AP_idx = 1:L
end %TP_idx = 1:num_tp_points

RSS_tp_mW = N*p*db2pow(beta_fngprnt_tp) + N*db2pow(noiseVariancedBm); %RSS_tp_mW --> 16x9 vector
RSS_tp_dB = 10*log10(RSS_tp_mW/100); %16x9, with TPs in the centre of the 16 subregions

[RPs_in_a_subrgn, no_of_RPs_in_a_subrgn, no_subrgns] = OfflineAlgo.functionOfflineSetup(RP_positions_per_row, RP_SIM_TYPE);
cluster_size = 3; %floor(0.80*length(AP_positions))
AP_cluster_idx = OfflineAlgo.cluster_AP_kmeans(AP_positions,cluster_size);
% RPs_in_a_subrgn : 16x4, 16 subregions, 4 RP numbers in each subregion; 24x6 , 24 subregions, 6 RPs in each subregion for RP_SIM_TYPE = 2
% no_of_RPs_in_a_subrgn = 4
% no_subrgns = 16
itrn_check = zeros(1,num_tp_points);
best_subrgn_seed = zeros(1,num_tp_points);
APS_size =  3
APS_array = nchoosek((1:L),APS_size);
num_APS_combinations = nchoosek(L,APS_size);
best_subrgn = zeros(1,num_APS_combinations);
dist_diff_avg = zeros(1,num_APS_combinations);
for TP_idx = 1:10%1:num_tp_points % num_tp_points=16
    per_tp_RSS = RSS_tp_dB(TP_idx,:); %per_tp_RSS --> 9x1 vector of the TP
    per_tp_est_dists = polyval(fngprnt_poly,per_tp_RSS);
    per_tp_est_dists_2D = sqrt(per_tp_est_dists.^2 - (h_BS-h_UT)^2);
    for APS_idx = 1:num_APS_combinations
        APS_idx;
        APS = APS_array(APS_idx,:);
        APS_dists = [];
        APS_dists = d_2D(:,APS); %APS_dists --> 25xno_AP_in_APS i.e. 25x3 vector for APS=3
        [best_subrgn(APS_idx),dist_diff_avg(APS_idx)] = OnlineAlgo.find_best_subrgn(per_tp_est_dists_2D(APS),APS_dists,RPs_in_a_subrgn,no_subrgns,no_of_RPs_in_a_subrgn);  %index of best subregion
    end
    
    best_subrgn_final = mode(best_subrgn);
    best_subrgn_indices = find(best_subrgn == mode(best_subrgn));
    least_diff_of_best_subrgn = sort(dist_diff_avg(best_subrgn_indices));
    APS_final = APS_array(find(dist_diff_avg == least_diff_of_best_subrgn(1),1),:);
    APS_dists_final = d_2D(:,APS_final);
    [RPS_dists,RPS_x,RPS_y] = OnlineAlgo.find_RPS_dists(best_subrgn_final,APS_dists_final,RPs_in_a_subrgn,RP_positions);  %RPS_RSS, RPS_x, RPS_y --> (APSx4) vector
    estimated_pos = OnlineAlgo.wknn_positioning(RPS_dists,per_tp_est_dists_2D(APS),length(APS),RPS_x,RPS_y,no_of_RPs_in_a_subrgn);
    TP_COORDS = [real(TP_positions(TP_idx)) imag(TP_positions(TP_idx))];
    dist_coords = [estimated_pos;TP_COORDS];
    estimated_err(TP_idx) = pdist(dist_coords,'euclidean');
    
    x = unique(best_subrgn);
    num_x = [];
    for i = x
        num_x = [num_x, numel(find(best_subrgn==i))];
    end
    TP_idx
    ["subrgn" x;"rpt_cnt" num_x]
    selected_subrgn = best_subrgn_final
    fprintf("==================================");
end %for TP_idx = 1:num_tp_points

%{
for i = 1:num_tp_points
    TP_COORDS = [real(TP_positions(i)) imag(TP_positions(i))];
    fprintf("(%d,%d) %d\n",TP_COORDS(1),TP_COORDS(2),estimated_err(i))
end
%}
best_subrgn_final;
APS_final;
fprintf("err_avg = %d\n",mean(estimated_err));
end %for n = 1:nbrOfSetups
diary off

