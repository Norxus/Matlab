clc
clear
IterationTimes = 10; % 发送的符号数目
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
A = [1/sqrt(2);-1/sqrt(2)];
[x8,x7,x6,x5,x4,x3,x2,x1] = ndgrid(A,A,A,A,A,A,A,A);
X_candidate = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:)]'; %所有排列组合情况
x_yalmip = (2*binvar(2*NT,1)-1)/sqrt(2);
Constraints = [];
options = [];
% options = sdpsettings('solver','cplex','verbose',2);
% options.cplex= cplexoptimset('cplex');
% options.cplex.mip.tolerances.mipgap = 0.01;

for jj=1:IterationTimes
    
    ip = randn(2*NT,1)>0; %等概率产生0和1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    
    X_transmit = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)];
    X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
        p = 10^(Eb(ii)/10);
        Y = H*X+(10^(-Eb(ii)/20))*N;%利用上面的信号x，信道H，noise n计算出输出信号Y
%       Y = sqrt(2p)*H*X_transmit+(10^(-Eb(ii)/20))*N;
        Y_hat = [real(Y);imag(Y)];
        H_hat = [real(H),-imag(H);imag(H),real(H)];
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        
        GR_T = [];
        GR = zeros(2*NR,2*NT);
        YR = [sign(real(Y)); sign(imag(Y))];
        for kk = 1:2*NR
            GR(kk,:)= H_hat(kk,:)*YR(kk);
        end
        Objective = sum(abs(GR*x_yalmip) - GR*x_yalmip);
   
        sol = optimize(Constraints, Objective, options);
   
        solution = value(x_yalmip);
        Y_mid_OneBit = [solution(1:NT) + 1i*solution(NT+1:2*NT)];

   
      
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
x = [-10,-5,0,5,10];
semilogy(x,SER,'-o');
grid on
