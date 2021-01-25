function [SER] = ZF(M,NT,NR,IterationTimes,Eb)
% M为调制方式（进制），NT为发射端个数，NR为接收端个数，IterationTimes为迭代次数，
% Eb为设定的SNR范围（目前只支持5个不同的信噪比进行比较）
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
if M == 2
      X_S = [1,-1];
elseif M == 4
     X_S = [1+1i,1-1i,-1+1i,-1-1i]/sqrt(2);
elseif M ==16
    X_S = [3+1i,3-1i,3-3i,3+3i,-3-3i,-3+3i,-3-i,-3+i,-1+i,-1-i,-1+3i,-1-3i,1+3i,1-3i,1+i,1-i]/sqrt(10);
else
    error('不支持该调制方式');
end

for jj=1:IterationTimes
    if M == 2
        ip = randn(NT,1)>0; %等概率产生0和1
        X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
        X = X_hat;
    elseif M == 4
         ip = randn(2*NT,1)>0; %等概率产生0和1
         X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
         X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
    elseif M == 16
        ip = randn(4*NT,1)>0; %等概率产生0和1
        X_bpsk = 2*ip-1; % 0 -> -1; 1 -> 1
        P = zeros(2*NT,4*NT);
        P_1 = zeros(1,4*NT);
        P_1(1) = 2;
        P_1(2) = 1;
        P(1,:) = P_1;
        for cc = 2:2*NT
            P_temp = circshift(P_1,2*(cc-1));
            P(cc,:) = P_temp;
        end
        X_hat = P*X_bpsk;
        X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(10);
    end
   
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
        p = 10^(Eb(ii)/10);
        Y = sqrt(p)*H*X+N;%利用上面的信号x，信道H，noise n计算出输出信号Y
%       Y = sqrt(2p)*H*X_transmit+(10^(-Eb(ii)/20))*N;
        Y_hat = [sign(real(Y));sign(imag(Y))];
        Y = [Y_hat(1:NR)+1i*Y_hat(NR+1:2*NR)];
        
        if M == 2
            H_hat = [real(H);imag(H)];
        else
            H_hat = [real(H),-imag(H);imag(H),real(H)];
        end
        
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        


      Y_mid_OneBit = pinv(H)*Y;  
      Y_mid_OneBit = sqrt(NT)*Y_mid_OneBit/norm(Y_mid_OneBit);
      for i_mid = 1:NT
          mid = [];
          for i_inmid = 1:M
              midNumber = abs(Y_mid_OneBit(i_mid)-X_S(i_inmid));
              mid=[mid,midNumber];
          end
          Y_mid_OneBit(i_mid) = X_S(find(mid == min(mid)));
      end
      
       if(ii == 1)
        E_f10 = [E_f10,Y_mid_OneBit(1:NT)];  
      end
      
      if(ii == 2)
        E_f5 = [E_f5,Y_mid_OneBit(1:NT)];  
      end
      
      if(ii == 3)
        E_0 = [E_0,Y_mid_OneBit(1:NT)];  
      end
      
      if(ii == 4)
        E_5 = [E_5,Y_mid_OneBit(1:NT)];  
      end
      
      if(ii == 5)
        E_10 = [E_10,Y_mid_OneBit(1:NT)];  
       end
  end
    
        
        
        l =NT-length(find(E_0(:,jj)==X));
        errorNumber = errorNumber + l;


       
        l1 =NT-length(find(E_5(:,jj)==X));
        errorNumber_5 = errorNumber_5 + l1;
            
        
                                                                                                                                                                                                                                                                                 
        l2 = NT-length(find(E_10(:,jj)==X));
        errorNumber_10 = errorNumber_10 + l2;

        
        
        l3 = NT-length(find(E_f5(:,jj)==X));
        errorNumber_f5 = errorNumber_f5 + l3;

        
        
        l4 = NT-length(find(E_f10(:,jj)==X));
        errorNumber_f10 = errorNumber_f10 + l4;
    
    end


SER_0 = errorNumber/(NT*IterationTimes);
SER_5 = errorNumber_5/(NT*IterationTimes);
SER_10 = errorNumber_10/(NT*IterationTimes);
SER_f5 = errorNumber_f5/(NT*IterationTimes);
SER_f10 = errorNumber_f10/(NT*IterationTimes);
% 
SER = [SER_f10,SER_f5,SER_0,SER_5,SER_10];
end

