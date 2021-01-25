clc
close all
clear all

syms u

IterationTimes = 30; % 发送的符号数目
NT=2;      %发送天线数
NR=6;      %接受天线数
Eb=4;     %信噪比
E_4 = [];
E_6 = [];
E_8 = [];
E_10 = []; 
E_12 = [];
E_14 = [];
E_16 = [];

errorNumber_4 =0;
errorNumber_6 =0;
errorNumber_8 =0;
errorNumber_10 =0;
errorNumber_12 =0;
errorNumber_14 =0;
errorNumber_16 =0;
 l4=0;
 l6=0;
 l8=0;
 l10=0;
 l12=0;
 l14=0;
 l16=0;
 M=4;
 X_S = [1+1i,1-1i,-1+1i,-1-1i];
%  X_S = [1,-1];
A = [0,1,2,3];
% A = [0,1];
 [x2,x1] = ndgrid(A,A);
 B = [x2(:),x1(:)]';
 X_candidate = zeros(NT,(M)^NT);
 X_transmit = zeros(NT,1);
 for l_i = 1:16
       l(l_i) = B(1,l_i)*1+B(2,l_i)*4;
%        l(l_i) = B(1,l_i)*1+B(2,l_i)*2;
       X_candidate_index = B(:,l_i);
       X_candidate_index_1 = X_candidate_index(1)+1;
       X_candidate_index_2 = X_candidate_index(2)+1;
       X_candidate(:,l_i) = [X_S(X_candidate_index_1),X_S(X_candidate_index_2)];
 end
 
for jj=1:IterationTimes
      jj
%     ip = randn(2*NT,1)>0; %等概率产生0和1
%     X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    X_index = round(rand(1,1)*15)+1;
    

     X_transmit = X_candidate(:,X_index);

    
    X = X_transmit/sqrt(2);
%     X = X_transmit;
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
        p = 10^(Eb(ii)/10);
        Y = sqrt(p)*H*X+N;%利用上面的信号x，信道H，noise n计算出输出信号Y
        Y_hat = [real(Y);imag(Y)];
        H_hat = [real(H),-imag(H);imag(H),real(H)];
%         H_hat = [real(H);imag(H)];
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        X_candidate_hat = [real(X_candidate);imag(X_candidate)];
%         X_candidate_hat = X_candidate;

        YR = [sign(real(Y)); sign(imag(Y))]; %r
        
        C = zeros(2*NR,(M)^NT);
        Probility = zeros(2*NR,(M)^NT);
        for i_c = 1 : (M)^NT
            for j_c = 1 : 2*NR
                C(j_c,i_c) = sign(H_hat(j_c,:)*X_candidate_hat(:,i_c));
            end
            for j_p = 1 : 2*NR
                 Lbound = double(2*abs(H_hat(j_p,:)*X_candidate_hat(:,i_c)));
                 errProbility = double((1/(2*pi))*int(exp(-(u^2)/2),Lbound,inf));
                if(C(j_p,i_c) == YR(j_p))
                  Probility(j_p,i_c)  =  1 - errProbility;
                else
                  Probility(j_p,i_c)  =  errProbility;  
                end
            end
        end
      [MLD_value,MLD_Index] = max(prod(Probility));  
      

 
        Y_mid_OneBit = X_candidate(:,MLD_Index);

   
      
       if(ii == 1)
        E_4 = [E_4,Y_mid_OneBit];  
      end
      
      if(ii == 2)
        E_6 = [E_6,Y_mid_OneBit];  
      end
      
      if(ii == 3)
        E_8 = [E_8,Y_mid_OneBit];  
      end
      
      if(ii == 4)
        E_10 = [E_10,Y_mid_OneBit];  
      end
      
      if(ii == 5)
        E_12 = [E_12,Y_mid_OneBit];  
      end
       
      if(ii == 6)
        E_14 = [E_14,Y_mid_OneBit];  
      end
      
      if(ii == 7)
        E_16 = [E_16,Y_mid_OneBit];  
      end
  end
    
        
        l4 =NT-length(find(E_4(:,jj)==X_transmit));
        errorNumber_4 = errorNumber_4 + l4;


       
        l6 =NT-length(find(E_6(:,jj)==X_transmit));
        errorNumber_6 = errorNumber_6 + l6;
            
        
                                                                                                                                                                                                                                                                                 
        l8 = NT-length(find(E_8(:,jj)==X_transmit));
        errorNumber_8 = errorNumber_8 + l8;

        
        
        l10 = NT-length(find(E_10(:,jj)==X_transmit));
        errorNumber_10 = errorNumber_10 + l10;

        
        
        l12 = NT-length(find(E_12(:,jj)==X_transmit));
        errorNumber_12 = errorNumber_12 + l12;
        
        l14 = NT-length(find(E_14(:,jj)==X_transmit));
        errorNumber_14 = errorNumber_14 + l14;
        
        l16 = NT-length(find(E_16(:,jj)==X_transmit));
        errorNumber_16 = errorNumber_16 + l16;
    
    end

SER_4 = errorNumber_4/(NT*IterationTimes);
SER_6 = errorNumber_6/(NT*IterationTimes);
SER_8 = errorNumber_8/(NT*IterationTimes);
SER_10 = errorNumber_10/(NT*IterationTimes);
SER_12 = errorNumber_12/(NT*IterationTimes);
SER_14 = errorNumber_14/(NT*IterationTimes);
SER_16 = errorNumber_16/(NT*IterationTimes);
% 
SER = [SER_4,SER_6,SER_8,SER_10,SER_12,SER_14,SER_16];
x = [4,6,8,10,12,14,16];
semilogy(x,SER,'-o');
grid on
