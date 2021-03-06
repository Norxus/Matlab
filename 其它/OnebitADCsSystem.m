clc
clear
IterationTimes = 100; % 发送的符号数目
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
A = [1;-1];
[x8,x7,x6,x5,x4,x3,x2,x1] = ndgrid(A,A,A,A,A,A,A,A);
X_candidate = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:)]'; %所有排列组合情况
x_yalmip = 2*binvar(2*NT,1)-1;
Constraints = [];
options = [];
% Objective = norm(y-H*x)^2;

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
        GR_NEW = mat2cell(GR,ones(1,2*NR),2*NT);
        GR_OLD = GR_NEW;
        YR = [sign(real(Y)); sign(imag(Y))];
        for kk = 1:2*NR
            GR_NEW{kk,1}= GR_NEW{kk,1}*YR(kk);
            GR_T = cat(2,GR_T,GR_NEW{kk,1}');
        end
        G_R = GR_NEW;
        
%         X_0 =  sqrt(NT)*(GR_T*ones(2*NR,1)/norm(GR_T*ones(2*NR,1)));
%         f = [];
%         k_final = 0.01;
%         c_final = 10^(-3);
%         for TT = 1:2*NR
%              temp = (1/sqrt(2*pi))*(exp((-p)*abs(G_R{TT,1}*X_0)^2))/normcdf(sqrt(2*p)*G_R{TT,1}*X_0,0,1);
%              f = [f,temp];
%         end
%         X_1 = X_0 + k_final*GR_T*f';
%         if(norm(X_1)^2 > NT)
%                X_1 = sqrt(NT)*(X_1/norm(X_1));
%         end
%         R(:,1) = X_0;
%         R(:,2) = X_1;
%         CN = 2;
%         while(norm(R(:,CN) - R(:,CN-1))>= c_final*norm(R(:,CN-1)))
%           f_while = [];
%           for T = 1:2*NR
%              temp = (1/sqrt(2*pi))*(exp((-p)*abs(G_R {T,1}*R(:,CN-1))^2))/normcdf(sqrt(2*p)*G_R{T,1}*R(:,CN-1),0,1);
%              f_while = [f_while,temp];
%           end
%             CN=CN+1;
%             R(:,CN) = R(:,CN-1) +k_final*GR_T*f_while';
%             if(norm(R(:,CN))^2 > NT)
%                 R(:,CN) = sqrt(NT)*(R(:,CN)/norm(R(:,CN)));
%             end
%         end

        X_ML = 1;
         X_estimate = [];
    for k = 1:size(X_candidate,2)
        X_ML = 1;
         for i = 1:2*NR
                X_ML = X_ML* normcdf(sqrt(2*p)*G_R{i,1}*X_candidate(:,k),0,1);
         end
      X_estimate = [X_estimate,X_ML];
    end
    find(X_estimate == max(X_estimate))
      X_Output = X_candidate(:,find(X_estimate == max(X_estimate)));
      Y_mid_OneBit = [X_Output(1:NT) + 1i*X_Output(NT+1:end)]
      
%       Y_mid_OneBit = sqrt(NT)*(1/norm(R(:,CN)))*(R(:,CN));
% 
%       for i_mid = 1:NT
%           Y_mid_OneBit(i_mid) = Y_mid_OneBit(i_mid)+1i*Y_mid_OneBit(i_mid+NT);
%           mid = [];
%           for i_inmid = 1:M
%               midNumber = abs(Y_mid_OneBit(i_mid)-X_S(i_inmid));
%               mid=[mid,midNumber];
%           end
%           Y_mid_OneBit(i_mid) = X_S(find(mid == min(mid)));
%       end
%       FinalOutput = [real(Y_mid_OneBit(1:NT));imag(Y_mid_OneBit(1:NT))];
      
%       if(ii == 1)
%         E_f10 = [E_f10,FinalOutput];  
%       end
%       
%       if(ii == 2)
%         E_f5 = [E_f5,FinalOutput];  
%       end
%       
%       if(ii == 3)
%         E_0 = [E_0,FinalOutput];  
%       end
%       
%       if(ii == 4)
%         E_5 = [E_5,FinalOutput];  
%       end
%       
%       if(ii == 5)
%         E_10 = [E_10,FinalOutput];  
%       end
%        if(ii == 1)
%         E_f10 = [E_f10,Y_mid_OneBit(1:NT)];  
%       end
%       
%       if(ii == 2)
%         E_f5 = [E_f5,Y_mid_OneBit(1:NT)];  
%       end
%       
%       if(ii == 3)
%         E_0 = [E_0,Y_mid_OneBit(1:NT)];  
%       end
%       
%       if(ii == 4)
%         E_5 = [E_5,Y_mid_OneBit(1:NT)];  
%       end
%       
%       if(ii == 5)
%         E_10 = [E_10,Y_mid_OneBit(1:NT)];  
%        end
%   end

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


% 
% SER_0 = (1/NT)*errorNumber/(2*NT*IterationTimes);
% SER_5 = (1/NT)*errorNumber_5/(2*NT*IterationTimes);
% SER_10 = (1/NT)*errorNumber_10/(2*NT*IterationTimes);
% SER_f5 = (1/NT)*errorNumber_f5/(2*NT*IterationTimes);
% SER_f10 = (1/NT)*errorNumber_f10/(2*NT*IterationTimes);

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
