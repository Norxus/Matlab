clc
clear
IterationTimes = 3000; % 发送的符号数目
NT=8;      %发送天线数
NR=192;      %接受天线数
Eb=-10:5:0;     %信噪比
E_f10 = [];
E_f5 = [];
E_0 = [];
E_5 = []; 
E_10 = [];
CNIteration = [];
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
 M=16;
 X_S = [3+1i,3-1i,3-3i,3+3i,-3-3i,-3+3i,-3-i,-3+i,-1+i,-1-i,-1+3i,-1-3i,1+3i,1-3i,1+i,1-i]/sqrt(10);

for jj=1:IterationTimes
     jj   
    Eb_error=[];
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
    
    X_transmit = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)];
    X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(10);
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
        
        R = zeros(2*NT,120);
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
          f_while = [];
          for T = 1:2*NR
             temp = (1/sqrt(2*pi))*(exp((-p)*abs(G_R {T,1}*R(:,CN-1))^2))/normcdf(sqrt(2*p)*G_R{T,1}*R(:,CN-1),0,1);
             f_while = [f_while,temp];
          end
            CN=CN+1;
            R(:,CN) = R(:,CN-1) +k_final*GR_T*f_while';
            if(CN>40)
                break;
            end
            if(norm(R(:,CN))^2 > NT)
                R(:,CN) = sqrt(NT)*(R(:,CN)/norm(R(:,CN)));
            end
        end
        
        
%          CNIteration = [CNIteration,CN];  
%         if(CN>50||prod(isnan(R(:,CN))))
%                 Eb_error = [Eb_error,ii];
% %                 CNIteration((jj-1)*5+ii)=0; 
%                 continue;
%             end
       
      Y_mid_OneBit = sqrt(NT)*(1/norm(R(:,CN)))*(R(:,CN));

      for i_mid = 1:NT
          Y_mid_OneBit(i_mid) = Y_mid_OneBit(i_mid)+1i*Y_mid_OneBit(i_mid+NT);
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
%       
%       if(ii == 4)
%         E_5 = [E_5,Y_mid_OneBit(1:NT)];  
%       end
%       
%       if(ii == 5)
%         E_10 = [E_10,Y_mid_OneBit(1:NT)];  
%        end
    end

%             if(length(Eb_error)>0||prod(isnan(R(:,CN))))
% 
%                 for i_error = 1 : length(Eb_error)
%                     if(Eb_error(i_error) == 1)
%                         E_f10(:,jj) = zeros(NT,1);
%                     elseif(Eb_error(i_error) == 2)
%                         E_f5(:,jj) = zeros(NT,1);
%                     elseif(Eb_error(i_error) == 3)
%                         E_0(:,jj) = zeros(NT,1);
%                      elseif(Eb_error(i_error) == 4)
%                         E_5(:,jj) = zeros(NT,1);
%                      elseif(Eb_error(i_error) == 5)
%                         E_10(:,jj) = zeros(NT,1);
%                     end
%                 end
%             end
            
        if(prod(E_0(:,jj) == zeros(NT,1)))
           errorNumber = errorNumber + 0; 
        else
           l =NT-length(find(E_0(:,jj)==X_transmit/sqrt(10)));
           errorNumber = errorNumber + l;
        end


%         if(prod(E_5(:,jj) == zeros(NT,1)))
%            errorNumber_5 = errorNumber_5 + 0; 
%         else
%            l1 =NT-length(find(E_5(:,jj)==X_transmit/sqrt(10)));
%            errorNumber_5 = errorNumber_5 + l1;
%         end
%             
%         
%         if(prod(E_10(:,jj) == zeros(NT,1)))
%            errorNumber_10 = errorNumber_10 + 0; 
%         else                                                                                                                                                                                                                                                                         
%            l2 = NT-length(find(E_10(:,jj)==X_transmit/sqrt(10)));
%            errorNumber_10 = errorNumber_10 + l2;
%          end

        
        if(prod(E_f5(:,jj) == zeros(NT,1)))
           errorNumber_f5 = errorNumber_f5 + 0; 
        else
           l3 = NT-length(find(E_f5(:,jj)==X_transmit/sqrt(10)));
           errorNumber_f5 = errorNumber_f5 + l3;
        end

        
        if(prod(E_f10(:,jj) == zeros(NT,1)))
           errorNumber_f10 = errorNumber_f10 + 0; 
        else
           l4 = NT-length(find(E_f10(:,jj)==X_transmit/sqrt(10)));
           errorNumber_f10 = errorNumber_f10 + l4;
        end
    
    end
len1 = length(find(E_f10(1,:)~=0));
len2 = length(find(E_f5(1,:)~=0));
len3 = length(find(E_0(1,:)~=0));
% len4 = length(find(E_5(1,:)~=0));
% len5 = length(find(E_10(1,:)~=0));

SER_0 = errorNumber/(NT*len3);
% SER_5 = errorNumber_5/(NT*len4);
% SER_10 = errorNumber_10/(NT*len5);
SER_f5 = errorNumber_f5/(NT*len2);
SER_f10 = errorNumber_f10/(NT*len1);
% 
SER = [SER_f10,SER_f5,SER_0];
x = [-10,-5,0];
semilogy(x,SER,'-o');
grid on
