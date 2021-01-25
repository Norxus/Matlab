clc
close all
clear all

IterationTimes = 100000; % 发送的符号数目
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
x_yalmip = (2*binvar(2*NT,1)-1)/sqrt(2);
y = sdpvar(2*NR,1)
options= sdpsettings;
options.solver = 'intlinprog';
% options = optimoptions(@intlinprog,'OutputFcn',@savemilpsolutions,'PlotFcn',@optimplotmilp);


for jj=1:IterationTimes
    jj
    ip = randn(2*NT,1)>0; %等概率产生0和1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    
    X_transmit = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)];
    X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道

    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
        p = 10^(Eb(ii)/10);
        Y = H*X+(10^(-Eb(ii)/20))*N;%利用上面的信号x，信道H，noise n计算出输出信号Y
        H_hat = [real(H),-imag(H);imag(H),real(H)];
        
        GR = zeros(2*NR,2*NT);
        YR = [sign(real(Y)); sign(imag(Y))];
        for kk = 1:2*NR
            GR(kk,:)= H_hat(kk,:)*YR(kk);
        end
%         Objective = sum(abs(GR*x_yalmip) - GR*x_yalmip);
   
        Constraints = [y<=GR*x_yalmip,y<=0];
%         Constraints = [y<=GR*x_yalmip];
        Objective = -sum(y);
%         Objective = sum(y.^2);
        sol = optimize(Constraints, Objective, options);
        
   
        solution = round(value(x_yalmip));
 
        Y_mid_OneBit = [solution(1:NT) + 1i*solution(NT+1:2*NT)];%SER
%         Y_mid_OneBit = solution;%BER

   
      
       if(ii == 1)
        E_f10 = [E_f10,Y_mid_OneBit];  
      end
      
      if(ii == 2)
        E_f5 = [E_f5,Y_mid_OneBit];  
      end
      
      if(ii == 3)
        E_0 = [E_0,Y_mid_OneBit];  
      end
      
      if(ii == 4)
        E_5 = [E_5,Y_mid_OneBit];  
      end
      
      if(ii == 5)
        E_10 = [E_10,Y_mid_OneBit];  
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


SER_0 = errorNumber/(NT*jj);
SER_5 = errorNumber_5/(NT*jj);
SER_10 = errorNumber_10/(NT*jj);
SER_f5 = errorNumber_f5/(NT*jj);
SER_f10 = errorNumber_f10/(NT*jj);
% 
SER = [SER_f10,SER_f5,SER_0,SER_5,SER_10];
x = [-10,-5,0,5,10];
semilogy(x,SER,'-o');
grid on
