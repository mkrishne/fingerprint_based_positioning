function [fngprnt_poly] = polynomial_regression(RSS_fngprnt_dB,d_3D)
    poly_x = [];
    poly_y = [];
    test_poly_x = [];
    test_poly_y = [];
    L = size(RSS_fngprnt_dB,2);
    num_RPs = size(RSS_fngprnt_dB,1);
    for AP_idx = 1:(L-1)
        for RP_idx = 1:num_RPs
            poly_x = [poly_x,RSS_fngprnt_dB(RP_idx,AP_idx)];
            poly_y = [poly_y,d_3D(RP_idx,AP_idx)];
        end
    end

    for AP_idx = L
        for RP_idx = 1:num_RPs
            test_poly_x = [test_poly_x,RSS_fngprnt_dB(RP_idx,AP_idx)];
            test_poly_y = [test_poly_y,d_3D(RP_idx,AP_idx)];
        end
    end

    fngprnt_poly = polyfit(poly_x,poly_y,5);
    f1 = polyval(fngprnt_poly,test_poly_x);
    polyfit_err = zeros(1,numel(test_poly_y));
    for i = 1:numel(test_poly_y)
        polyfit_err(i) = pdist([f1(i);test_poly_y(i)],'euclidean');
    end
    fprintf("polyfit err_avg = %d\n",mean(polyfit_err));




