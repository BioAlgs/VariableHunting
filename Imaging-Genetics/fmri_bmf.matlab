sig_dir = dir('/home/stat/xinxing/extdata/cm/ADNI/final_region/*.csv')
sig_files = strcat({'/home/stat/xinxing/extdata/cm/ADNI/final_region/'},{sig_dir.name});

[~, nregion]= size(sig_files)

for region=1:nregion
    ymat_stuc= importdata(sig_files{region});

    [pathstr,sig_name,ext] = fileparts(sig_files{region});

    ymat_raw = ymat_stuc.data';
   

    [ns , nt] = size(ymat_raw);
    % bmat_out = zeros(ns, r);
    % zero_ind = find(sum(ymat_raw,2)==0)

    ymat_raw  = ymat_raw - sum(ymat_raw,2)/nt * ones(1,nt);
    yvar = var(ymat_raw, 0, 1);
    ymat = ymat_raw*diag(1./sqrt(yvar));


    r = 20
    bmat = (dec2bin(0:2^(r)-1)-48)';

    ymat(find(sum(ymat,2)==0),:) = [];
    [n,t] = size(ymat);

    %Fisrt find r row and column
    [U,S,V] = svd(ymat);
    umat = U(:,1:r);

    umat_r = [];
    umat_rind = [1];
    for i=1:r
        if i==1
            umat_r = umat(1,:);
        else
            proj_cosin = [];
            for j=i:n
                proj = umat_r' * inv(umat_r*umat_r') * umat_r * umat(j,:)';
                proj_cosin = [proj_cosin, norm(proj)/ norm(umat(j,:))];
            end
            [M,I] = min(proj_cosin);
            umat_rind = [umat_rind, I];
            umat_r = [umat_r ; umat(I+i-1,:)];
        end
    end

    bmat = 2*((dec2bin(0:2^(r)-1)-48)' - 0.5);

    tmat = umat * inv(umat_r) * bmat;
    tmat_0 = 2*((tmat > 0) - 0.5);
    score = sqrt(sum((tmat-tmat_0).^2, 1));
    [~, I] = sort(score);

    rt = 2;
    bmat_1 = tmat_0(:,I(1:(rt*r)));

    snp_file_name = strcat('/extdata/xx/project/fmri-snp/snp_std_1/', sig_name, '-snp.txt')

    dlmwrite(snp_file_name, bmat_1, ' ');
 

end
