classdef OnlineAlgo
    methods ( Static = true )
        function [best_subrgn,dist_diff_avg] = find_best_subrgn(APS_est_dist,APS_distances,RPs_in_a_subrgn,no_subrgns,no_of_RPs_in_a_subrgn) %APS_est_dist --> 1XAPS vector; 1x1 value for seed AP
            %fprintf("inside function find_best_subrgn\n");
            all_subrgns_dist_diffs = zeros(no_subrgns,no_of_RPs_in_a_subrgn); %16x4 vector
            for i = 1:size(APS_est_dist,2)
                each_AP_est_dist = APS_est_dist(i);
                each_AP_dist_diff = abs(APS_distances(:,i)' - each_AP_est_dist);
                for subrgn_idx = 1:no_subrgns
                        RP_indices_in_subrgn = RPs_in_a_subrgn(subrgn_idx,:);
                        subrgn_dist_diff = each_AP_dist_diff(RP_indices_in_subrgn);
                        %subrgn_rss_diff_without_outliers = OnlineAlgo.rss_diff_remove_outliers(subrgn_rss_diff);
                        all_subrgns_dist_diffs(subrgn_idx,:) = all_subrgns_dist_diffs(subrgn_idx,:) + subrgn_dist_diff;% + subrgn_rss_diff_without_outliers';
                end 
            end
            all_subrgns_dist_diffs;
            all_subrgns_dist_diff_avg = mean(all_subrgns_dist_diffs,2);
            [subrgn_dist_diff_avg_sorted, subrgn_dist_diff_sorted_idx] = sort(all_subrgns_dist_diff_avg,'ascend');
            best_subrgn = subrgn_dist_diff_sorted_idx(1);
            dist_diff_avg = subrgn_dist_diff_avg_sorted(1);
            %fprintf("out of function find_best_subrgn\n");
        end

        function [RPS_dists,RPS_x,RPS_y] = find_RPS_dists(best_subrgn,aps_dists,RPs_in_a_subrgn,RP_positions)
            %fprintf("inside function find_RPS_RSS\n");
            % RPs_in_a_subrgn --> 16x4, 16 subrgns, 4 RPs of the subrgn
            % RP_positions --> 5x5
            % aps_dists --> 25xAPS
            RPS = RPs_in_a_subrgn(best_subrgn,:); %RPS --> 1x4 vector of 4 vertices of best subregion
            RPS_dists = aps_dists(RPS,:);  %RPS_RSS --> 4xAPS vector
            RPS_xy = RP_positions(RPS);   % 4x1 vector
            RPS_x = real(RPS_xy);   % 4x1 vector 
            RPS_y = imag(RPS_xy);   % 4x1 vector
            % ex: RPS_x = [50;60;50;60] and RPS_y = [60;60;70;70]
            % coordinates of 4 vertices of best subregion are : (50,60),(60,60),(50,70) & (60,70) 
            %fprintf("out of function find_RPS_RSS\n");
        end


         function APS = find_APS(best_subrgn,RPs_in_a_subrgn,per_tp_est_dists,d_2D,AP_cluster_idx,cluster_size)
            %fprintf("inside function find_APS\n");
            % RPs_in_a_subrgn --> 16x4 vector
            % RSS_tp --> 1x9 vector of RSS of 9APs at TP
            % RSS_fngprnt --> 25x9 vector; All RSS values of 9 APs at 25RPs 
            RPS = RPs_in_a_subrgn(best_subrgn,:); %RPS --> (1x4) vector of 4 vertices of best subregion
            ap_dist_for_rps = d_2D(RPS,:);  %ap_fngprnt_for_rps --> 4x9 vector; 4=RPS vertices;9=no of APs
            diff_dist = abs(ap_dist_for_rps-per_tp_est_dists);  % diff_rss --> 4x9 vector; 
            avg_diff_dists = mean(diff_dist,1); %avg_diff_dists --> 1x9 vector; mean of all APs at a RP; 1 => mean along the columns
            [avg_diff_dists_sorted, avg_diff_dists_sorted_idx] = sort(avg_diff_dists,'ascend');
            APS = avg_diff_dists_sorted_idx(1:9);
         end

    %in WKNN positioning, instead of RSS, the dist diff is from AP to TP is %taken as weights
    function estimated_pos = wknn_positioning(APS_RP_dists, APS_tp_dist, no_APs_in_APS,RPS_x,RPS_y,no_of_RPs_in_a_subrgn)
            %fprintf("inside function wknn_positioning\n");
            dists_k = APS_RP_dists; %APS_RP_dists --> 4xAPS; 4x1 for seed AP
            inv_diff_dists = []; %APSx4 vector
            M = no_APs_in_APS;
            for m = 1:M 
                diff_dists = abs(dists_k(:,m) - APS_tp_dist(m)); %diff_dists --> 4x1 vector; diff of 4 RSS values of the RPS with TP RSS values for an AP
                inv_diff_dists(m,:) = 1./diff_dists; %1x4 vector; 
            end
            inv_diff_dists;
            w = []; %4x1 vector
            for k = 1:no_of_RPs_in_a_subrgn
                w(k) = mean(inv_diff_dists(:,k)); % 1x4 vector
            end
            w;
            sigma_wk = sum(w); %1x1 sum of 4 avg values 
            estimated_pos_x = sum(w.*RPS_x)/sigma_wk;
            estimated_pos_y = sum(w.*RPS_y)/sigma_wk;
            estimated_pos = [estimated_pos_x,estimated_pos_y];
            %fprintf("out of function wknn_positioning\n");
        end
    end %methods ( Static = true )
end %classdef OnlineAlgo