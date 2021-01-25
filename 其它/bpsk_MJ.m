clc
clear
IterationTimes = 1000; % 发送的符号数目
NT=2;      %发送天线数
NR=16;      %接受天线数
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
A = [1;-1];
[x2,x1] = ndgrid(A,A);
X_candidate = [x1(:) x2(:)]'; %所有排列组合情况
% [x8,x7,x6,x5,x4,x3,x2,x1] = ndgrid(A,A,A,A,A,A,A,A);
% X_candidate = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:)]';
% [x22,x21,x20,x19,x18,x17,x16,x15,x14,x13,x12,x11,x10,x9,x8,x7,x6,x5,x4,x3,x2,x1] = ndgrid(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A);
% X_candidate = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:) x9(:) x10(:) x11(:) x12(:) x13(:) x14(:) x15(:) x16(:) x17(:) x18(:) x19(:) x20(:) x21(:) x22(:)]';

for jj=1:IterationTimes
    jj
    
    ip = randn(NT,1)>0; %等概率产生0和1
    X = 2*ip-1; % 0 -> -1; 1 -> 1
    
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
        Y = H*X+(10^(-Eb(ii)/20))*N;%利用上面的信号x，信道H，noise n计算出输出信号Y
        Y_hat = [real(Y);imag(Y)];
        H_hat = [real(H);imag(H)];
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        
        GR = zeros(2*NR,NT);
        YR = [sign(real(Y)); sign(imag(Y))];
        for kk = 1:2*NR
            GR(kk,:)= H_hat(kk,:)*YR(kk);
        end
        
        
       for i = 1:size(X_candidate,2)
           Z(:,i) = GR*X_candidate(:,i);
       end
        dist = sum(Z-abs(Z),1);
        [maxvalue indest] = max(dist);
        X_Output = X_candidate(:,indest);
   
      Y_mid_OneBit = X_Output(1:NT);
      

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
