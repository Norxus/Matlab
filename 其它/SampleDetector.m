clc
clear
IterationTimes = 10000; % 发送的符号数目
NT=4;      %发送天线数
NR=32;      %接受天线数
Eb=-10:5:10;     %信噪比
E_f10 = [];
E_f5 = [];
E_0 = [];
E_5 = []; 
E_10 = [];
errorNumber =0;
errorNumber_5 =0;
errorNumber_10 =0;
errorNumber_f5 =0;
errorNumber_f10 =0;
 l=0;
 l1=0;
 l2=0;
 l3=0;
 l4=0;
 M=4;
 X_S = [1+1i,1-1i,-1+1i,-1-1i];

for jj=1:IterationTimes
    ip = randn(2*NT,1)>0; %等概率产生0和1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    
    X_transmit = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)];
    X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
        p = 10^(Eb(ii)/10);
        Y = sqrt(p)*H*X+N;%利用上面的信号x，信道H，noise n计算出输出信号Y
%       Y = sqrt(2p)*H*X_transmit+(10^(-Eb(ii)/20))*N;
        Y_hat = [real(Y);imag(Y)];
        H_hat = [real(H),-imag(H);imag(H),real(H)];
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        
        GR_T = [];
        GR = [real(H.'),imag(H.');-imag(H.'),real(H.')]';
        GR_NEW = mat2cell(GR,ones(1,2*NR),2*NT); %Gr,n
        GR_OLD = GR_NEW;
        YR = [sign(real(Y)); sign(imag(Y))];
        r_zf = [YR(1:NR,1)+1i*YR(NR+1:2*NR,1)];
        for kk = 1:2*NR
            GR_NEW{kk,1}= GR_NEW{kk,1}*YR(kk);
            GR_T = cat(2,GR_T,GR_NEW{kk,1}');
        end
        G_R = GR_NEW;
        S_zf = (H.')*pinv(H*(H.'))*r_zf;
        S_zf_normlize = sqrt(NT)*(S_zf/norm(S_zf));
        S_zf_normlize_Second = [];
        
     for i_mid = 1:NT
%           S_zf_normlize(i_mid) = S_zf_normlize(i_mid)+1i*S_zf_normlize(i_mid+NT);
          mid = [];
          for i_inmid = 1:M
              midNumber = abs(S_zf_normlize(i_mid)-X_S(i_inmid));
              mid=[mid,midNumber];
          end
          S_zf_normlize(i_mid) = X_S(find(mid == min(mid)));
          mid(find(mid == min(mid))) = inf;
          S_zf_normlize_Second = [S_zf_normlize_Second;X_S(find(mid == min(mid)))];
     end
     S_zf_normlizeNeareast = [S_zf_normlize,S_zf_normlize_Second]; 
     
     X_candidateMatrix = [];
     for i_x_candidateIndex = 1:2
        X_candidate = [];
        X_candidateStartElement = S_zf_normlizeNeareast(1,i_x_candidateIndex);
        X_candidate = [X_candidate,X_candidateStartElement];
        for i_x_candidateSecondIndex = 1 : 2
%            IndexCopy(i_x_candidateIndex) = 0;
            if(length(X_candidate) == NT)
             X_candidate = X_candidate(1);
            end

             X_candidate = [X_candidate,S_zf_normlizeNeareast(2,i_x_candidateSecondIndex)];
            
          for i_x_candidateThirdIndex = 1 : 2
            if(length(X_candidate) == NT)
             X_candidate = X_candidate(1:2);
            end
             X_candidate = [X_candidate,S_zf_normlizeNeareast(3,i_x_candidateThirdIndex)];
            for i_x_candidateFourthIndex = 1 : 2

            if(length(X_candidate) == NT)
             X_candidate = X_candidate(1:3);
            end
             X_candidate = [X_candidate,S_zf_normlizeNeareast(4,i_x_candidateFourthIndex)];
             X_candidateMatrix = [X_candidateMatrix,X_candidate.'];

            end
         end
       end
     end
     
    X_candidateMatrix_R = [real(X_candidateMatrix);imag(X_candidateMatrix)];
    X_estimate = [];
    for k = 1:size(X_candidateMatrix_R,2)
        X_Two_nML = 1;
     for i = 1:2*NR
            X_Two_nML = X_Two_nML* normcdf(sqrt(2*p)*G_R{i,1}*X_candidateMatrix_R(:,k),0,1);
     end
      X_estimate = [X_estimate,X_Two_nML];
    end 
    X_TwoStage_FinalOutput = X_candidateMatrix(:,find(X_estimate == max(X_estimate)));
   if(ii == 1)
        E_f10 = [E_f10,X_TwoStage_FinalOutput];  
   end
      
   if(ii == 2)
        E_f5 = [E_f5,X_TwoStage_FinalOutput];  
   end
      
   if(ii == 3)
        E_0 = [E_0,X_TwoStage_FinalOutput];  
   end
      
   if(ii == 4)
        E_5 = [E_5,X_TwoStage_FinalOutput];  
   end
      
   if(ii == 5)
        E_10 = [E_10,X_TwoStage_FinalOutput];  
   end
 end
    
            
        l =NT-length(find(E_0(:,jj)==X_transmit));
        errorNumber = errorNumber + l;


       
        l1 =NT-length(find(E_5(:,jj)==X_transmit));
        errorNumber_5 = errorNumber_5 + l1;
            
        
                                                                                                                                                                                                                                                                                 
        l2 = NT-length(find(E_10(:,jj)==X_transmit));
        errorNumber_10 = errorNumber_10 + l2;

        
        
        l3 = NT-length(find(E_f5(:,jj)==X_transmit));
        errorNumber_f5 = errorNumber_f5 + l3;

        
        
        l4 = NT-length(find(E_f10(:,jj)==X_transmit));
        errorNumber_f10 = errorNumber_f10 + l4;
    
    end



SER_0 = (1/NT)*errorNumber/(2*NT*IterationTimes);
SER_5 = (1/NT)*errorNumber_5/(2*NT*IterationTimes);
SER_10 = (1/NT)*errorNumber_10/(2*NT*IterationTimes);
SER_f5 = (1/NT)*errorNumber_f5/(2*NT*IterationTimes);
SER_f10 = (1/NT)*errorNumber_f10/(2*NT*IterationTimes);

SER_0 = errorNumber/(NT*IterationTimes);
SER_5 = errorNumber_5/(NT*IterationTimes);
SER_10 = errorNumber_10/(NT*IterationTimes);
SER_f5 = errorNumber_f5/(NT*IterationTimes);
SER_f10 = errorNumber_f10/(NT*IterationTimes);

SER = [SER_f10,SER_f5,SER_0,SER_5,SER_10];
x = [-10,-5,0,5,10];
semilogy(x,SER,'-o');
grid on
