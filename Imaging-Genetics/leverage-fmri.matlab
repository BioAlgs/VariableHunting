run('/data/xinxing/software/spams-matlab/spams-matlab.m')

cd ~/extdata/xx/project/fmri-leverage
signame='~/extdata/xx/project/fmri-leverage/1.EMOTION.sig.txt';
asig = importdata(signame);
pvec = [0.01, 0.05, 0.1, 0.2];
res = [];
for i=1:4
    param.mode = 0;
    param.K = 50;
    param.lambda = 0.1;
    param.numThreads = -1;
    param.iter=1000;
    %param.batchsize = 400;
    D0=mexTrainDL(asig,param);
    alpha=mexLasso(asig,D0,param);

    %%Run dictionary learning with leverage sampling
    amat = alpha';
    acet = inv(amat' * amat);
    [n p] = size(amat);
    levec = zeros(n,1);
    for k = 1:n
        levec(k) = amat(k,:) * acet *amat(k,:)';
    end 
    levec = levec./sum(levec);
    rid_lev =  datasample(1:n,floor(n*pvec(i)),'Replace', false,'Weight', levec)
    asig_lev = asig(:,rid_lev);
    param.mode = 0;
    param.K = 50;
    param.lambda = 0.1;
    param.numThreads = -1;
    param.iter=1000;
    %param.batchsize = 400;
    D0_lev=mexTrainDL(asig_lev,param);
    [r1,r2_lev,r3_lev] = comparetwomatrix(D0,D0_lev)


    %%Run dictionary learning with uniform sampling
    [ntim, nvex] = size(asig);
    rid = randi(nvex,floor(nvex*pvec(i)) ,1);
    asig_uni = asig(:,rid);
    param.mode = 0;
    param.K = 50;
    param.lambda = 0.1;
    param.numThreads = -1;
    param.iter=1000;
    %param.batchsize = 400;
    D0_uni=mexTrainDL(asig_uni,param);

    [r1,r2_uni,r3_uni] = comparetwomatrix(D0,D0_uni)

    res = [res; [r2_lev, r2_uni, r3_lev,r3_uni]];
end



