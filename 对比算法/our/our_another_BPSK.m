clc
clear
IterationTimes = 100000; % 发送的符号数目
NT=2;      %发送天线数
NR=2;      %接受天线数
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
X_candidate = [x1(:) x2(:)]';
x_yalmip = (2*binvar(2*NT,1)-1)/sqrt(2);
Constraints = [];
options = [];


for jj=1:IterationTimes
    jj
    ip = randn(NT,1)>0; %等概率产生0和1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    
    X_transmit = X_hat;
    X = X_hat;
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
        p = 10^(Eb(ii)/10);
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
        
        dist = sum(Z<0);;
        [maxvalue indest] = min(dist);
        X_Output = X_candidate(:,indest);
   
      Y_mid_OneBit = X_Output;
      

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
    
% 
%         l =2*NT-length(find(E_0(:,jj)==X_hat));
%         errorNumber = errorNumber + l;
% 
% 
%        
%         l1 =2*NT-length(find(E_5(:,jj)==X_hat));
%         errorNumber_5 = errorNumber_5 + l1;
%             
%         
%                                                                                                                                                                                                                                                                                  
%         l2 = 2*NT-length(find(E_10(:,jj)==X_hat));
%         errorNumber_10 = errorNumber_10 + l2;
% 
%         
%         
%         l3 =2*NT-length(find(E_f5(:,jj)==X_hat));
%         errorNumber_f5 = errorNumber_f5 + l3;
% 
%         
%         
%         l4 =2*NT-length(find(E_f10(:,jj)==X_hat));
%         errorNumber_f10 = errorNumber_f10 + l4;
        
        
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


% 
% SER_0 = (1/NT)*errorNumber/(2*NT*IterationTimes);
% SER_5 = (1/NT)*errorNumber_5/(2*NT*IterationTimes);
% SER_10 = (1/NT)*errorNumber_10/(2*NT*IterationTimes);
% SER_f5 = (1/NT)*errorNumber_f5/(2*NT*IterationTimes);
% SER_f10 = (1/NT)*errorNumber_f10/(2*NT*IterationTimes);

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
