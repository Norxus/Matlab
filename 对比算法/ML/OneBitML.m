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
[x4,x3,x2,x1] = ndgrid(X_S,X_S,X_S,X_S);
 X_S_temp = [x4(:),x3(:),x2(:),x1(:)]'/sqrt(2);
 X_ML = [real(X_S_temp);imag(X_S_temp)];

for jj=1:IterationTimes
    jj
    
    ip = randn(2*NT,1)>0; %等概率产生0和1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    
    X_transmit = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)];
    X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    E = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 错误矩阵
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道

    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
        p = 10^(Eb(ii)/10);
        Y = sqrt(p)*H*X+N;%利用上面的信号x，信道H，noise n计算出输出信号Y
%       Y = sqrt(2p)*H*X_transmit+(10^(-Eb(ii)/20))*N;
        Y_hat = [sign(real(Y));sign(imag(Y))];
        H_hat = [real(H),-imag(H);imag(H),real(H)];
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        G = zeros(2*NT,2*NR);
        x_ML_detection = zeros(1,size(X_ML,2));
        
        for i = 1 : 2*NR
           G(:,i) =  H_hat(i,:)*Y_hat(i);
        end
        
        for i = 1 : size(X_ML,2)
           x_ML_detection(i) = prod(normcdf(sqrt(2*p)*(G')*X_ML(:,i),0,1));
        end
        
        Y_mid_OneBit  = X_ML(:,find(x_ML_detection == max(x_ML_detection)));
        Y_mid_OneBit = [ Y_mid_OneBit(1:NT)+1i* Y_mid_OneBit(NT+1:2*NT)];
       
      
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
len = length(find(E_0(1,:)~=0));

SER_0 = errorNumber/(NT*len);
SER_5 = errorNumber_5/(NT*len);
SER_10 = errorNumber_10/(NT*len);
SER_f5 = errorNumber_f5/(NT*len);
SER_f10 = errorNumber_f10/(NT*len);
% 
SER = [SER_f10,SER_f5,SER_0,SER_5,SER_10];
x = [-10,-5,0,5,10];
semilogy(x,SER,'-o');
grid on
