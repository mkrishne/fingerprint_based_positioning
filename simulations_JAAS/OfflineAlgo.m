classdef OfflineAlgo
    methods ( Static = true )
        function [RPs_in_a_subrgn, no_of_RPs_in_a_subrgn, no_subrgns] = functionOfflineSetup(RP_positions_per_row, RP_SIM_TYPE)
            if(RP_SIM_TYPE == 1)
                no_subrgns = (RP_positions_per_row-1)^2; %a 5x5 matrix has 16 SRs
                no_subrgns_per_row = sqrt(no_subrgns);
                no_of_RPs_in_a_subrgn = 4;
                RPs_in_a_subrgn = zeros(no_subrgns,no_of_RPs_in_a_subrgn); % RPs_in_a_subrgn --> 16x4 vector; 4 vertices of 16 subregions
                for subrgn_idx = 1:no_subrgns
                     vertex1 = subrgn_idx + floor((subrgn_idx-1)/no_subrgns_per_row);  %formula -> SR1 has vertices 1,2,6,7; SR5 has vertices 6,7,11,12
                     vertex2 = vertex1 + 1;
                     vertex3 = vertex1 + RP_positions_per_row;
                     vertex4 = vertex2 + RP_positions_per_row;
                     RPs_in_a_subrgn(subrgn_idx,1) = vertex1;
                     RPs_in_a_subrgn(subrgn_idx,2) = vertex2;
                     RPs_in_a_subrgn(subrgn_idx,3) = vertex3;
                     RPs_in_a_subrgn(subrgn_idx,4) = vertex4;
                end
            
             elseif(RP_SIM_TYPE == 2)  
                no_minor_subrgns = (RP_positions_per_row-1)^2; %a 5x5 matrix has 16 minor SRs
                no_minor_subrgns_per_row = sqrt(no_minor_subrgns);
                no_major_subrgns = 2*(no_minor_subrgns_per_row-1)*no_minor_subrgns_per_row; % 5x5 matrix has 24 major SRs i.e. 2 smaller subregion areas
                no_of_RPs_in_mjr_subrgn = 6;
                RPs_in_mjr_subrgn = zeros(no_major_subrgns,no_of_RPs_in_mjr_subrgn); % RPs_in_mjr_subrgn --> 16x4 vector; 4 vertices of 16 subregions
                no_major_subrgns_type1 = no_major_subrgns/2;
                for subrgn_idx_type1 = 1:no_major_subrgns_type1
                     vertex1 = subrgn_idx_type1 + floor((subrgn_idx_type1-1)/(no_minor_subrgns_per_row-1))*2;  %formula -> major SR type1 has vertices 1,2,3,6,7,8; SR5 has vertices 7,8,9,12,13,14
                     vertex2 = vertex1 + 1;
                     vertex3 = vertex2 + 1;
                     vertex4 = vertex1 + RP_positions_per_row;
                     vertex5 = vertex2 + RP_positions_per_row;
                     vertex6 = vertex3 + RP_positions_per_row;
                     RPs_in_mjr_subrgn(subrgn_idx_type1,1) = vertex1;
                     RPs_in_mjr_subrgn(subrgn_idx_type1,2) = vertex2;
                     RPs_in_mjr_subrgn(subrgn_idx_type1,3) = vertex3;
                     RPs_in_mjr_subrgn(subrgn_idx_type1,4) = vertex4;
                     RPs_in_mjr_subrgn(subrgn_idx_type1,5) = vertex5;
                     RPs_in_mjr_subrgn(subrgn_idx_type1,6) = vertex6;
                end
                no_major_subrgns_type2 = no_major_subrgns/2;
                for subrgn_idx_type2 = 1:no_major_subrgns_type2
                     vertex1 = subrgn_idx_type2 + floor((subrgn_idx_type2-1)/no_minor_subrgns_per_row);  %formula -> major SR type2 has vertices 1,2,6,7,11,12; SR5 has vertices 7,8,12,13,17,18
                     vertex2 = vertex1 + 1;
                     vertex3 = vertex1 + RP_positions_per_row;
                     vertex4 = vertex2 + RP_positions_per_row;
                     vertex5 = vertex3 + RP_positions_per_row;
                     vertex6 = vertex4 + RP_positions_per_row;
                     RPs_in_mjr_subrgn(no_major_subrgns_type1+subrgn_idx_type2,1) = vertex1;
                     RPs_in_mjr_subrgn(no_major_subrgns_type1+subrgn_idx_type2,2) = vertex2;
                     RPs_in_mjr_subrgn(no_major_subrgns_type1+subrgn_idx_type2,3) = vertex3;
                     RPs_in_mjr_subrgn(no_major_subrgns_type1+subrgn_idx_type2,4) = vertex4;
                     RPs_in_mjr_subrgn(no_major_subrgns_type1+subrgn_idx_type2,5) = vertex5;
                     RPs_in_mjr_subrgn(no_major_subrgns_type1+subrgn_idx_type2,6) = vertex6;
                end
                RPs_in_a_subrgn = RPs_in_mjr_subrgn;
                no_of_RPs_in_a_subrgn = no_of_RPs_in_mjr_subrgn;
                no_subrgns = no_major_subrgns;
          elseif(RP_SIM_TYPE == 3)
              no_subrgns = (RP_positions_per_row-2)^2; %a 5x5 matrix has 9 minor SRs
              no_subrgns_per_row = sqrt(no_subrgns);
              no_of_RPs_in_a_subrgn = 9;
              RPs_in_a_subrgn = zeros(no_subrgns,no_of_RPs_in_a_subrgn); % RPs_in_a_subrgn --> 16x4 vector; 4 vertices of 16 subregionss
              for subrgn_idx = 1:no_subrgns
                     vertex1 = subrgn_idx + 2*floor((subrgn_idx-1)/no_subrgns_per_row);  %formula -> SR1 has vertices 1,2,3,6,7,8,11,12,13; SR4 has vertices 6,7,8,11,12,13,16,17,18
                     vertex2 = vertex1 + 1;
                     vertex3 = vertex2 + 1;
                     vertex4 = vertex1 + RP_positions_per_row;
                     vertex5 = vertex2 + RP_positions_per_row;
                     vertex6 = vertex3 + RP_positions_per_row;
                     vertex7 = vertex4 + RP_positions_per_row;
                     vertex8 = vertex5 + RP_positions_per_row;
                     vertex9 = vertex6 + RP_positions_per_row;
                     RPs_in_a_subrgn(subrgn_idx,1) = vertex1;
                     RPs_in_a_subrgn(subrgn_idx,2) = vertex2;
                     RPs_in_a_subrgn(subrgn_idx,3) = vertex3;
                     RPs_in_a_subrgn(subrgn_idx,4) = vertex4;
                     RPs_in_a_subrgn(subrgn_idx,5) = vertex5;
                     RPs_in_a_subrgn(subrgn_idx,6) = vertex6;
                     RPs_in_a_subrgn(subrgn_idx,7) = vertex7;
                     RPs_in_a_subrgn(subrgn_idx,8) = vertex8;
                     RPs_in_a_subrgn(subrgn_idx,9) = vertex9;
                end
          end %if(RP_SIM_TYPE == 1)
      end %functionOfflineSetup
      function AP_cluster_idx = cluster_AP_kmeans(AP_positions, cluster_size)
            AP_cluster_idx = kmeans([real(AP_positions); imag(AP_positions)]',cluster_size);
      end
   end %methods ( Static = true )
end %classdef OnlineAlgo