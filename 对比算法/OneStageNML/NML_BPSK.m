clc
clear
IterationTimes = 100000; % 发送的符号数目
NT=16;      %发送天线数
NR=64;      %接受天线数
Eb=-10:5:10;     %信噪比
E_f10 = [];
E_f5 = [];
E_0 = [];
E_5 = []; 
E_10 = [];
E_15 = [];
errorNumber =0;
errorNumber_5 =0;
errorNumber_10 =0;
errorNumber_f5 =0;
errorNumber_f10 =0;
errorNumber_f15 =0;
 l=0;
 l1=0;
 l2=0;
 l3=0;
 l4=0;
 l5=0;
 M=2;
 X_S = [1,-1];

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
        Y = sqrt(p)*H*X+N;%利用上面的信号x，信道H，noise n计算出输出信号Y
%       Y = sqrt(2p)*H*X_transmit+(10^(-Eb(ii)/20))*N;
        Y_hat = [real(Y);imag(Y)];
        H_hat = [real(H);imag(H)];
        N_hat = [real(10^(-Eb(ii)/20)*N);imag(10^(-Eb(ii)/20)*N)];
        
        GR_T = [];
        GR = [real(H);imag(H)];
        GR_NEW = mat2cell(GR,ones(1,2*NR),NT);
        GR_OLD = GR_NEW;
        YR = [sign(real(Y)); sign(imag(Y))];
        for kk = 1:2*NR
            GR_NEW{kk,1}= GR_NEW{kk,1}*YR(kk);
            GR_T = cat(2,GR_T,GR_NEW{kk,1}');
        end
        G_R = GR_NEW;

        X_0 =  sqrt(NT)*(GR_T*ones(2*NR,1)/norm(GR_T*ones(2*NR,1)));
        f = [];
        k_final = 0.01;
        c_final = 10^(-3);
        for TT = 1:2*NR
             temp = (1/sqrt(2*pi))*(exp((-p)*abs(G_R{TT,1}*X_0)^2))/normcdf(sqrt(2*p)*G_R{TT,1}*X_0,0,1);
             f = [f,temp];
        end
        X_1 = X_0 + k_final*GR_T*f';
        if(norm(X_1)^2 > NT)
               X_1 = sqrt(NT)*(X_1/norm(X_1));
        end
        R(:,1) = X_0;
        R(:,2) = X_1;
        CN = 2;
        while(norm(R(:,CN) - R(:,CN-1))>= c_final*norm(R(:,CN-1)))
%             norm(R(:,CN) - R(:,CN-1))
          f_while = [];
          for T = 1:2*NR
             temp = (1/sqrt(2*pi))*(exp((-p)*abs(G_R {T,1}*R(:,CN))^2))/normcdf(sqrt(2*p)*G_R{T,1}*R(:,CN),0,1);
             f_while = [f_while,temp];
          end
            CN=CN+1;
            R(:,CN) = R(:,CN-1) +k_final*GR_T*f_while';
%             if(CN>100||prod(isnan(R(:,CN))))
%                 break;
%             end
            if(norm(R(:,CN))^2 > NT)
                R(:,CN) = sqrt(NT)*(R(:,CN)/norm(R(:,CN)));
            end
        end
% 
%             if(CN>100||prod(isnan(R(:,CN))))
%                 break;
%             end 
      Y_mid_OneBit = sqrt(NT)*(1/norm(R(:,CN)))*(R(:,CN));

      for i_mid = 1:NT
%           Y_mid_OneBit(i_mid) = Y_mid_OneBit(i_mid)+1i*Y_mid_OneBit(i_mid+NT);
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
       
%        if(ii == 6)
%         E_15 = [E_15,Y_mid_OneBit(1:NT)];  
%        end
    end
    
%             if(CN>100||prod(isnan(R(:,CN))))
%                 E_f10(:,jj) = zeros(NT,1);  
%                 E_f5(:,jj) = zeros(NT,1); 
%                 E_0(:,jj) = zeros(NT,1);
%                 E_5(:,jj) = zeros(NT,1); 
%                 E_10(:,jj) = zeros(NT,1); 
% %                 E_15(:,jj) = zeros(NT,1);
%                 continue;
%             end
    
        
        
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
        
%         l5 = NT-length(find(E_15(:,jj)==X_transmit));
%         errorNumber_f15 = errorNumber_f15 + l5;
    
end
len = length(find(E_0(1,:)~=0))

SER_0 = errorNumber/(NT*len);
SER_5 = errorNumber_5/(NT*len);
SER_10 = errorNumber_10/(NT*len);
SER_f5 = errorNumber_f5/(NT*len);
SER_f10 = errorNumber_f10/(NT*len);
% 
SER = [SER_f10,SER_f5,SER_0,SER_5,SER_10];
x = -10:5:10;
semilogy(x,SER,'-o');
grid on
