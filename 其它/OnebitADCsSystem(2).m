addpath('../main/');

IterationTimes = 1000;

NT=16;      %发送天线数
NR=128;      %接受天线数
Eb=-10:5:10;     %信噪比

ntest = 10;         % Number of Monte-Carlo tests
numiter = 1000;      % Number of iterations
nx = 16;           % Number of input components
nz = 128;  % Number of output components
newdat = 1; 	    % Set =1 to generate new data each test

sparseRat = 1; 
QPSK_IN = 1;
GAUSS_OUT = 2;
SER_PERF = 3;
GAUSS_MAT = 4;
inDist = QPSK_IN; 
outDist = GAUSS_OUT; 
perfMetric = SER_PERF;
matrixType = GAUSS_MAT;

opt = GampOpt();        % default parameters
opt.nit = numiter;      % number of iterations
%opt.pvarMin = 1e-10;
%opt.xvarMin = 1e-10;
opt.adaptStep = 1;
%opt.stepWindow = 0;								
opt.uniformVariance = 0;	
opt.tol = -1;          % stopping tolerance 
opt.stepIncr = 1.1;	% multiplicative step increase
opt.stepMax = 0.3;	% max step 			
opt.step = 0.01;	% initial step 
%opt.stepMin = 0.01;	% min step
%opt.stepTol = -1;       % stopping tolerance based on step
opt.pvarStep = 1;	% apply step to pvar? 					
opt.varNorm = 0;	% normalize the variances? 
% opt.xvar0auto = true

xmax = 1;       % Max value
x0 = [-xmax xmax]';
px0 = [1 1]'/sqrt(2);
inputEst0 = DisScaEstim(x0, px0);
[xmean0,xvar0] = inputEst0.dist.meanVar();

if (sparseRat < 1)
    inputEst = SparseScaEstim( inputEst0, sparseRat );
else
    inputEst = inputEst0;
end

if (perfMetric == SER_PERF)
    x0 = inputEst.getPoints();
end

% Get mean and variance after sparsification
[xmean0s, xvar0s] = inputEst.estimInit;

% if (outDist == GAUSS_OUT) || (outDist == CGAUSS_OUT)
%     snr = 10;                       % SNR 
%     wmean = 0;                      % Noise mean
%     wvar = 10.^(-0.1*snr)*sparseRat*(abs(xmean0)^2+xvar0);    % Noise variance
%     if (matrixType == DFT_MAT_BLK)||(matrixType == DFT_MAT_RAND)
%       wvar = wvar*(nx/nz);	% DFT_MAT normalized differently
%     end;
%     if (outDist == GAUSS_OUT) && ((matrixType == DFT_MAT_BLK) ...
%                 ||(DFT_MAT_RAND)||(matrixType == CGAUSS_MAT)),
%         outDist = CGAUSS_OUT;
%     warning('Forcing outDist=CGAUSS_OUT because matrix is complex-valued')
%     end;
% elseif (outDist == POISSON_OUT)
%     snr = 20;  % Output SNR
%     poisson_scale = 10^(0.1*snr);
% elseif (outDist == LOGIT_OUT)
%     logitScale = 10;
% end

metricGAMP = nan(numiter, ntest);

for jj=1:IterationTimes
    
    ip = randn(2*NT,1)>0; %等概率产生0和1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
        
        Y = sqrt(1/NT)*H*X+(10^(-Eb(ii)/20))*N;%利用上面的信号x，信道H，noise n计算出输出信号Y
%       Y = H*X+(10^(-Eb(ii)/20))*N;
        Y_hat = [real(Y);imag(Y)];
        H_hat = [real(H),-imag(H);imag(H),real(H)];
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        wvar = 10^(-Eb(ii)/10);
        w = N_hat;
        
            % IterationTimes = 1000; % 发送的符号数目
      

    for itest = 1:ntest

        if (newdat)
            % Generate random input vector
            x = [real(X);imag(X)];
       if (matrixType == GAUSS_MAT)
            a0 = 0;
            A = H_hat;
            Aop = MatrixLinTrans(A);
       end
       z = Aop.mult(x);

      if (outDist == GAUSS_OUT)
            y = z + N_hat;
            outputEst = AwgnEstimOut(y, wvar);
      end
      snr_actual = 20*log10(norm(z)/norm(w));
      end%newdat
    [xhat, xvar, rhat, rvar, shatFinal, svarFinal,zhatFinal,zvarFinal, estHist] = ...
                    gampEst(inputEst, outputEst, Aop, opt);
            xhatGAMP = xhat;
            xhatTot = estHist.xhat;

        if (perfMetric == SER_PERF)

            metricMeth(itest,imeth) = measSER(xhat, x, x0);
        %             measSER(xhat, x, x0)
            if (meth == GAMP_METH)
        %                  measSER(xhatTot, x, x0)
        %                 metricGAMP(:,itest) = measSER(xhatTot, x, x0);
                SER = measSER(xhatTot, x, x0);
                metricGAMP(1:size(SER,1),itest) = measSER(xhatTot, x, x0);
                metricGAMP(size(SER,1)+1:end,itest) = metricGAMP(size(SER,1),itest);
            end
            % Display results
            fprintf(1,'it=%d %s ser=%f\n', itest, methStr{meth}, metricMeth(itest,imeth));

        end
end
        
%         GR_T = [];
%         GR = [real(H.'),imag(H.');-imag(H.'),real(H.')]';
%         GR_NEW = mat2cell(GR,ones(1,2*NR),2*NT);
%         GR_OLD = GR_NEW;
%         YR = [sign(real(Y)); sign(imag(Y))];
%         for kk = 1:2*NR
%             GR_NEW{kk,1}= GR_NEW{kk,1}*YR(kk);
%             GR_T = cat(2,GR_T,GR_NEW{kk,1}');
%         end
%         G_R = GR_NEW;
%         p = 10^(Eb(ii)/10);

  end
    

 
    
end

% 
SER_0 = (1/NT)*errorNumber/(2*NT*IterationTimes);
SER_5 = (1/NT)*errorNumber_5/(2*NT*IterationTimes);
SER_10 = (1/NT)*errorNumber_10/(2*NT*IterationTimes);
SER_f5 = (1/NT)*errorNumber_f5/(2*NT*IterationTimes);
SER_f10 = (1/NT)*errorNumber_f10/(2*NT*IterationTimes);
% 
SER = [SER_f10,SER_f5,SER_0,SER_5,SER_10];
x = [-10,-5,0,5,10];
semilogy(x,SER);
