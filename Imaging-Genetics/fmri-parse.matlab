
%%Read the fmri data for AD
file_name_ad =  dir('/home/stat/xinxing/extdata/xx/fmri-snp/ADNI_Common_Data/AD/*_*_*');
files = strcat({'/home/stat/xinxing/extdata/xx/fmri-snp/ADNI_Common_Data/AD/'},{file_name_ad.name});
file_name_lmci = dir('/home/stat/xinxing/extdata/xx/fmri-snp/ADNI_Common_Data/LMCI/*_*_*');
files = [files strcat({'/home/stat/xinxing/extdata/xx/fmri-snp/ADNI_Common_Data/LMCI/'},{file_name_lmci.name})];
file_name_emci = dir('/home/stat/xinxing/extdata/xx/fmri-snp/ADNI_Common_Data/EMCI/*_*_*');
files = [files strcat({'/home/stat/xinxing/extdata/xx/fmri-snp/ADNI_Common_Data/EMCI/'},{file_name_emci.name})];
file_name_normal = dir('/home/stat/xinxing/extdata/xx/fmri-snp/ADNI_Common_Data/Normal/*_*_*');
files = [files strcat({'/home/stat/xinxing/extdata/xx/fmri-snp/ADNI_Common_Data/Normal/'},{file_name_normal.name})];

[~, sample_name, ~]  = fileparts(files{1});

sig = importdata('/home/stat/xinxing/extdata/cm/ADNI/final_result/002_S_0413.csv')

atlas = importdata('/home/stat/xinxing/extdata/xx/fmri-snp/ADNI_Common_Data/atlasStd2mm_whole_b_signals.txt');

[nregion,~] = size(unique(atlas));
altas_id = unique(atlas);
[~, nsample] = size(files);

for i= 1:nregion
    ind = find(atlas == altas_id(i));
    sig_alt_i = [];
    for j= 1:nsample
        sig_file_name=strcat(files{j},'/filtered_func_data2std.detrend_whole_b_signals.txt');
        dmat=importdata(sig_file_name);
        sig_alt_i = [sig_alt_i; mean(dmat(ind,:),1)];
    end
    out_file_name = strcat('/home/stat/xinxing/extdata/xx/fmri-snp/fmri-atlas/', num2str(altas_id(i)), '.txt')
    dlmwrite(out_file_name, sig_alt_i, ' ');
end

for region=2:100
    out_file_name = strcat('/home/stat/xinxing/extdata/xx/fmri-snp/fmri-atlas/', num2str(altas_id(region)), '.txt');
    ymat= importdata(out_file_name);
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

    snp_file_name = strcat('/home/stat/xinxing/extdata/xx/fmri-snp/fmri-atlas/', num2str(altas_id(region)), '-snp.txt')

    dlmwrite(snp_file_name, bmat_1, ' ');
end

sample_names = []
for i=1:76
    [~, sample_name, ~]  = fileparts(files{i});
    % if i==67
    %     continue;
    % end
    sample_names = [sample_names;sample_name];
end
dlmwrite('/home/stat/xinxing/extdata/xx/fmri-snp/fmri-atlas/samplename.txt', sample_names, '')


plot(xmat(5:6,:)')
saveas(gcf,'1.png')


sig_dir = dir('/home/stat/xinxing/extdata/cm/ADNI/final_region/*.csv')
sig_files = strcat({'/home/stat/xinxing/extdata/cm/ADNI/final_region/'},{sig_dir.name});

[~, nregion]= size(sig_files)

sig = importdata('/home/stat/xinxing/extdata/cm/ADNI/final_region/1.csv')
for i= 1:nregion
    ind = find(atlas == altas_id(i));
    sig_alt_i = [];
    for j= 1:nsample
        sig_file_name=strcat(files{j},'/filtered_func_data2std.detrend_whole_b_signals.txt');
        dmat=importdata(sig_file_name);
        sig_alt_i = [sig_alt_i; mean(dmat(ind,:),1)];
    end
    out_file_name = strcat('/home/stat/xinxing/extdata/xx/fmri-snp/fmri-atlas/', num2str(altas_id(i)), '.txt')
    dlmwrite(out_file_name, sig_alt_i, ' ');
end
plot(sig.data(:,1:10))
saveas(gcf, '1.png')

phen = readtable('/home/stat/xinxing/extdata/cm/ADNI/finalreport.csv')

sig = importdata('/home/stat/xinxing/extdata/cm/ADNI/final_region/5.csv')
[nt , ns] = size(sig.data)

dmat_raw = sig.data';
dmat  = dmat_raw - sum(dmat_raw,2)/nt * ones(1,nt);



figure
hold on
for k=1:20
    p(k)=plot(dmat(k,:)); 
    if  strcmp(phen.DX_bl{k}, 'AD')
        set(p(k),'Color', 'red')
    elseif strcmp(phen.DX_bl{k}, 'LMCI')
        set(p(k),'Color', 'red')
    elseif strcmp(phen.DX_bl{k},'EMCI')
        set(p(k),'Color', 'black')
    elseif strcmp(phen.DX_bl{k},'CN')
        set(p(k),'Color', 'black')
    end

end

saveas(gcf, '1.png')


sig = importdata('/home/stat/xinxing/extdata/cm/ADNI/final_region/5.csv')
dmat_raw = sig.data';
dmat  = dmat_raw - sum(dmat_raw,2)/nt * ones(1,nt);
ph_id = {'AD', 'LMCI', 'EMCI', 'CN'}
figure
hold on
for k=1:4
    ind = find(strcmpi(phen.DX_bl, ph_id{k}))
    
    p(k)=plot(sum(dmat(ind,:))); 
    if  strcmp(ph_id{k}, 'AD')
        set(p(k),'Color', 'red')
    elseif strcmp(ph_id{k}, 'LMCI')
        set(p(k),'Color', 'red')
    elseif strcmp(ph_id{k},'EMCI')
        set(p(k),'Color', 'black')
    elseif strcmp(ph_id{k},'CN')
        set(p(k),'Color', 'black')
    end

end

saveas(gcf, '1.png')
