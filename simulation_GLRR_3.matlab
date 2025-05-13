clear; %clc;

path(path,'/data/xx/project/variablehunting/GLRR/MlabFunctions'); % Bios Server

%% set some constants and hyperparameters
r = 2;  % rank of low-rank matrix
MCMCpara.nBurnin = 1e3; % number of burn-in samples
MCMCpara.nCollect = 1e3; % number of samples after burn-in
MCMCpara.a0 = 1e-6; % page 980, 1st gamma hyperparameter of tau_delta
MCMCpara.b0 = 1e-6; % page 980, 2nd gamma hyperparameter of tau_delta
MCMCpara.c0 = 1e-6; % page 980, 1st gamma hyperparameter of tau_u
MCMCpara.d0 = 1e-6; % page 980, 2nd gamma hyperparameter of tau_u
MCMCpara.e0 = 1e-6; % page 980, 1st gamma hyperparameter of tau_v,d
MCMCpara.f0 = 1e-6; % page 980, 2nd gamma hyperparameter of tau_v,d
MCMCpara.g0 = 1e-6; % 1st gamma hyperparameter, for the precision of the prior of B_p in equation (5)
MCMCpara.h0 = 1e-6; % 2nd gamma hyperparameter, for the precision of the prior of B_p in equation (5)
OutputFileName = 'Result';  %filename of output mat and csv file

%% input files

n = 200
m = 50
q = 20
p = 500

load('siglecellsim.RData')
[n, ~] = size(X);
nMissVec = [0 0.2 0.4 0.6 0.8].*20;
errVec = [1 2 5]./sqrt(m);

FPMat = zeros(length(nMissVec), length(errVec));
FNMat = zeros(length(nMissVec), length(errVec));

parfor ii=1:length(nMissVec)
	FPVec = zeros(1, length(errVec));
	FNVec = zeros(1, length(errVec));
	for jj = 1:length(errVec)

		% X = binornd(1,0.2, n, p);
		beta = normrnd(0,1,q,m);

		XT = X(:,1:q);
		sigma = errVec(jj);
		Y = XT * beta + normrnd(0, sigma ,n,m);


		nMiss = nMissVec(ii)
		idxOT = randsample(1:q, q-nMiss);
		XO = [XT(:,idxOT)  X(:,(q+1):p)];
		betaO = [beta(idxOT,:); zeros(p-q, m)];

		W = ones(n,1);


		%% RUN LGLRR %%
		TypeCov='Factor'; % Wish, Diag, Factor, default is iid
		%t0 = cputime; % track time

		tic
		Output = GLRRadj(XO,W,Y,90,5,MCMCpara,TypeCov); %Wish, Diag, Factor, iid(default)
		toc






		Bh = Output.PostB;
		Gh = Output.PostG;
 

		T0Vec = min(abs(Bh(:))) : 0.05 : max(abs(Bh(:)));
		score = zeros(length(T0Vec),5);
		nfold = 5;
		%%%Cross Validation Part
		for i = 1:length(T0Vec)
			Indices = crossvalind('Kfold', n, nfold);
			for j = 1:nfold
				XOj =XO(find(Indices~=j),:);
				Yj =Y(find(Indices~=j));
				Wj = ones(length(find(Indices~=j)),1);
				Outputj = GLRRadj(XOj,Wj,Yj,80,5,MCMCpara);
				Bhj = Outputj.PostB;
				Bhj(find(abs(Bhj)<T0Vec(i))) = 0;
				Ghj = Output.PostG;

				XOjtest = XO(find(Indices==j),:);
				YOjtest = Y(find(Indices==j));
				Wjtest = ones(length(find(Indices==j)),1);
				Ytest = XOjtest*Bhj + Wjtest*Ghj;% + W*B0;% + Z*b;
				score(i,j) = norm(Ytest-YOjtest);
				disp([i j score(i,j)])
			end
		end
			idx = find(mean(score,2)==min(mean(score,2)));
			T0Opt = T0Vec(1,idx);

			FPVec(1,jj)=length(find(abs(Bh)>T0Opt & betaO==0))
			FNVec(1,jj)=length(find(abs(Bh)<T0Opt & betaO~=0))
	end
	FPMat(ii,:) = FPVec;
	FNMat(ii,:) = FNVec;
end











Perf = PerformGLRR(Y, Yh, cov(XO), r, Bh);   
%time = (cputime - t0)/60;
%disp(['Time = ', num2str(time/60), ' hours']);

%% %%%%%%% Saving %%%%%%%%%
save ([OutputFileName,'.mat'],'-struct','Output');
csvwrite([OutputFileName,'.csv'], Perf);