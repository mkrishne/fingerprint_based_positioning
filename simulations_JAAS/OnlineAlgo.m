classdef OnlineAlgo
    methods ( Static = true )
        function best_subrgn = find_best_subrgn(rss_diff,RPs_in_a_subrgn,no_subrgns,no_of_RPs_in_a_subrgn) %rss_diff --> 25XAPS vector; 25x1 vector for seed AP
            %fprintf("inside function find_best_subrgn\n");
            all_subrgns_rss = zeros(no_subrgns,no_of_RPs_in_a_subrgn); %16x4 vector
            for i = 1:size(rss_diff,2)
                each_AP_rss_diff = rss_diff(:,i);
                for subrgn_idx = 1:no_subrgns
                        RP_indices_in_subrgn = RPs_in_a_subrgn(subrgn_idx,:);
                        subrgn_rss_diff = each_AP_rss_diff(RP_indices_in_subrgn); 
                        all_subrgns_rss(subrgn_idx,:) = all_subrgns_rss(subrgn_idx,:) + subrgn_rss_diff'; %Recalculation Using APS
                end 
            end
            all_subrgns_rss_mean = sum(all_subrgns_rss,2);
            [subrgn_mean_sorted, subrgn_mean_sorted_idx] = sort(all_subrgns_rss_mean,'ascend');
            best_subrgn = subrgn_mean_sorted_idx(1); %Best Subregion Identification 
            %fprintf("out of function find_best_subrgn\n");
        end
        
       function [RPS_RSS,RPS_x,RPS_y] = find_RPS_RSS(best_subrgn,aps_fngprnt,RPs_in_a_subrgn,RP_positions)
            %fprintf("inside function find_RPS_RSS\n");
            % RPs_in_a_subrgn --> 16x4, 16 subrgns, 4 RPs of the subrgn
            % RP_positions --> 5x5
            % aps_fngprnt --> 25xAPS
            RPS = RPs_in_a_subrgn(best_subrgn,:); %RPS --> 1x4 vector of 4 vertices of best subregion
            RPS_RSS = aps_fngprnt(RPS,:);  %RPS_RSS --> 4xAPS vector
            RPS_xy = RP_positions(RPS);   % 4x1 vector
            RPS_x = real(RPS_xy);   % 4x1 vector 
            RPS_y = imag(RPS_xy);   % 4x1 vector
            % ex: RPS_x = [50;60;50;60] and RPS_y = [60;60;70;70]
            % coordinates of 4 vertices of best subregion are : (50,60),(60,60),(50,70) & (60,70) 
            %fprintf("out of function find_RPS_RSS\n");
        end

      function APS = find_APS(best_subrgn,RPs_in_a_subrgn,RSS_tp,RSS_fngprnt,AP_cluster_idx,cluster_size)
            % RPs_in_a_subrgn --> 16x4 vector
            % RSS_tp --> 1x9 vector of RSS of 9APs at TP
            % RSS_fngprnt --> 25x9 vector; All RSS values of 9 APs at 25RPs 
            RPS = RPs_in_a_subrgn(best_subrgn,:); %RPS --> (1x4) vector of 4 vertices of best subregion
            ap_fngprnt_for_rps = RSS_fngprnt(RPS,:);  %ap_fngprnt_for_rps --> 4x9 vector; 4=RPS vertices;9=no of APs
            diff_rss = abs(ap_fngprnt_for_rps-RSS_tp);  % diff_rss --> 4x9 vector; 
            avg_diff_rss = mean(diff_rss,1); %avg_diff_rss --> 1x9 vector; mean of all APs at a RP; 1 => mean along the columns
            [avg_diff_rss_sorted, avg_diff_rss_sorted_idx] = sort(avg_diff_rss,'ascend');
            avg_diff_rss_zcores = zscore(avg_diff_rss_sorted);
            APS = avg_diff_rss_sorted_idx(1:5); %5 APs with lowest diff; %always 1:L is giving best results
            
            %[avg_diff_rss_sorted, avg_diff_rss_sorted_idx] = sort(avg_diff_rss,'descend');
            %APS = avg_diff_rss_sorted_idx(1:12);
                        
            %{
            APS = zeros(1,cluster_size);
            for i = 1:cluster_size
                idx = find(AP_cluster_idx == i);
                [~,APS_idx] = sort(avg_diff_rss(idx),'ascend'); %take the lowest value from each cluster
                APS(i) = idx(APS_idx(1));
            end
            %}
        end

        %In this function the diff values were clustered using kmeans in
        %online phase
        function APS = find_APS_kmeans(best_subrgn,no_AP_in_APS,RPs_in_a_subrgn,RSS_tp,RSS_fngprnt)
            %fprintf("inside function find_APS_kmeans\n");
            RPS = RPs_in_a_subrgn(best_subrgn,:); %RPS --> 1X4 vector of 4 vertices of best subregion
            ap_fngprnt_for_rps = RSS_fngprnt(RPS,:);  %ap_fngprnt_for_rps --> 4x9 vector; 4=RPS vertices;9=no of APs
            diff_rss = abs(ap_fngprnt_for_rps-RSS_tp);  % diff_rss --> 4x9 vector; 
            total_no_APs = size(diff_rss,2);
            for each_AP = (1:total_no_APs)
                diff_rss(:,each_AP) = OnlineAlgo.rss_diff_remove_outliers(diff_rss(:,each_AP));
            end
            diff_rss;
            kmeans_ap_samples = diff_rss'; %kmeans_ap_samples --> 9x4 vector
            [idx,C,sumD,D] = kmeans(kmeans_ap_samples,no_AP_in_APS);
            % C = no_AP_in_APS cluster centroid locations in the no_AP_in_APS-by-4 matrix. 
            % e.g. if no_AP_in_APS=3, then C is a 3x4 matrix, where each row gives the coordinates of  each of the 3 cluster centroid location
            % sumD =  within-cluster sums of point-to-centroid distances in the no_AP_in_APS-by-1 vector sumd
            % e.g sumD is a 3x1 matrix providing the sum of all points in each cluster to its centroid
            % D = 9xno_AP_in_APS matrix providing distances from each point to every centroid 
            % e.g. D is a 9x3 matrix, where D(1,1) gives distance of point1 (i.e. AP1 diff) to centroid 1,
            % D(1,2) gives distance of point1 to centroid 2, 
            % D(2,3) gives distance of point2 to centroid 3
            
            APS = [];
            for i = (1:no_AP_in_APS)
                AP_cluster_indices = find(idx == i);
                dist_points_to_centroid = D(AP_cluster_indices,i);
                [dist, dist_idx] = sort(dist_points_to_centroid,'ascend');
                APS(i) = AP_cluster_indices(dist_idx(1)); %select one AP from every cluster; a AP with closest to the centroid
            end
            APS;
            %fprintf("out of function find_APS_kmeans\n");
        end

    function estimated_pos = wknn_positioning(APS_RSS_of_subrgn_RPs, RSS_tp_of_APS, no_APs_in_APS,RPS_x,RPS_y,no_of_RPs_in_a_subrgn)
            %fprintf("inside function wknn_positioning\n");
            RSS_k = APS_RSS_of_subrgn_RPs; %APS_RSS_of_subrgn_RPs --> 4xAPS; 4x1 for seed AP
            %RSS_tp_of_APS = 1xAPS vector of the RSS from TP at each AP
            inv_diff_rss = []; %APSx4 vector
            M = no_APs_in_APS;
            for m = 1:M 
                diff_rss = abs(RSS_k(:,m) - RSS_tp_of_APS(m)); %diff_rss --> 4x1 vector; diff of 4 RSS values of the RPS with TP RSS values for an AP
                inv_diff_rss(m,:) = 1./diff_rss; %1x4 vector; 
            end
            w = []; %4x1 vector
            for k = 1:no_of_RPs_in_a_subrgn
                w(k) = mean(inv_diff_rss(:,k)); % 1x4 vector
            end
            sigma_wk = sum(w); %1x1 sum of 4 avg values 
            estimated_pos_x = sum(w.*RPS_x)/sigma_wk;
            estimated_pos_y = sum(w.*RPS_y)/sigma_wk;
            estimated_pos = [estimated_pos_x,estimated_pos_y];
            %fprintf("out of function wknn_positioning\n");
        end
    end %methods ( Static = true )
end %classdef OnlineAlgo