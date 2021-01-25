clc
close all
clear all

IterationTimes = 10000; % 发送的符号数目
NT=2;      %发送天线数
NR=6;      %接受天线数
Eb=0:5:30;     %信噪比


%总错误码字个数
errorNumber_4 =0;
errorNumber_6 =0;
errorNumber_8 =0;
errorNumber_10 =0;
errorNumber_12 =0;
errorNumber_14 =0;
errorNumber_16 =0;

%每次迭代的错误码字数量
 l4=0;
 l6=0;
 l8=0;
 l10=0;
 l12=0;
 l14=0;
 l16=0;
 M=4; %调制方式M-QAM
 X_S = [1+1i,1-1i,-1+1i,-1-1i]; %可以发送的所有码元
 sigma = 30;
 T = 8;  %训练数据个数
 Td = 512;%检测数据个数
 Tt = T*(M)^NT;%训练数据个数
 Tu = 10*Tt;


 
for jj=1:IterationTimes
      jj
    ip = randn(2*NT,Td)>0; %等概率产生0和1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    X = [X_hat(1:NT,:)+1i*X_hat(NT+1:2*NT,:)]/sqrt(2);%传输向量
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    H_hat = [real(H),-imag(H);imag(H),real(H)];


    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y
%%
        Y = H*X + (10^(-Eb(ii)/20))*N;
        YR = [sign(real(Y));sign(imag(sign(Y)))];%待检测向量
        
        %数据训练阶段
        [x2,x1] = ndgrid(X_S,X_S);
        x_train_all = [x2(:),x1(:)]'/sqrt(2);%所有可能的传输向量
        x_D_all = zeros(2*NR,(M)^NT*T);%所有可能传输向量计算得到的最终接受向量
        x_D_prosibility = zeros(2*NR,(M)^NT);%错误概率参数
        x_D_exception = zeros(2*NR,(M)^NT);%用于对比的c码字
        x_DL_temp = zeros(2*NR,T);
        x_D_prosibility_calculate = zeros(2*NR,(M)^NT*T);
       
       %标记信号集D
            for j = 1 : (M)^(NT)%发送向量所有可能的情况
                for l = 1 : T %每种可能的情况都会重复发送T次
                    Z_train = (10^(-Eb(ii)/20))*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);%每一种情况对应不同的噪声向量.
                    y_temp = H*x_train_all(:,j) + Z_train;
                    y_1_kjl = [sign(real(y_temp));sign(imag(y_temp))];
                    x_D_all(:,l+(j-1)*T) = y_1_kjl;
                end %这个for循环结束后，就相当于训练完一种情况了
                x_DL_temp = x_D_all(:,(j-1)*T+1:j*T);
                x_D_exception(:,j) = sign(sum( x_DL_temp,2));%计算用于对比的c码字
                temp_exception = x_D_exception(:,j);
                temp_exception(find(temp_exception == 0)) = 1;
                x_D_exception(:,j) = temp_exception;%下面这三行主要是将sum得到的0换成1，因为sign函数无法进行此项替换
                
                for n = 1 : 2*NR %然后用这个for循环计算错误概率
                    errnumber_train = 0;
                    for l = 1 : T
                        if(temp_exception(n) ~= x_DL_temp(n,l))
                            errnumber_train = errnumber_train + 1;
                        end
                    end
                x_D_prosibility(n,j) = (errnumber_train+1)/(T+2);
                end
            end
            

     
      X_finalDetection = zeros(2*NT,Td);
      theta_exception_update = x_D_exception;
%       
      for k = 1 : Td %得到两个参数之后，进行检测
          detection_different = theta_exception_update - YR(:,k);%相减，等于0的数据就就以为这对应位置相等，将每一个待检测向量与所有可能的情况进行对比
          detection_posibility = zeros(2*NR,(M)^NT);
          theta_posibility_update = x_D_prosibility;
          for i = 1 : (M)^NT
              for j = 1 : 2*NR
                  if(detection_different(j,i) ~= 0)
                      detection_posibility(j,i) = theta_posibility_update(j,i);
                  else
                      detection_posibility(j,i) = 1 - theta_posibility_update(j,i);
                  end
              end
          end

         [value,detectionIndex] = max(prod(detection_posibility));%竖乘并进行找到相应的索引
          X_finalDetection_mid = [real(x_train_all(:,detectionIndex));imag(x_train_all(:,detectionIndex))];
          %按照顺序进行训练，那么对比码字的顺序也是相同的，
          X_finalDetection(:,k) = X_finalDetection_mid;
%         X_finalDetection(:,k) = x_train_all(:,detectionIndex);
      end
      
       if(ii == 1)
        E_4 = X_finalDetection;  
      end
      
      if(ii == 2)
        E_6 = X_finalDetection;  
      end
      
      if(ii == 3)
        E_8 = X_finalDetection;  
      end
      
      if(ii == 4)
        E_10 = X_finalDetection;  
      end
      
      if(ii == 5)
        E_12 = X_finalDetection;  
      end
       
      if(ii == 6)
        E_14 = X_finalDetection;  
      end
      
      if(ii == 7)
        E_16 = X_finalDetection;  
      end
 end
    
        
        l4 =length(find(E_4~=X_hat/sqrt(2)));
        errorNumber_4 = errorNumber_4 + l4;


       
        l6 =length(find(E_6~=X_hat/sqrt(2)));
        errorNumber_6 = errorNumber_6 + l6;
            
        
                                                                                                                                                                                                                                                                                 
        l8 =length(find(E_8~=X_hat/sqrt(2)));
        errorNumber_8 = errorNumber_8 + l8;

        
        
        l10 =length(find(E_10~=X_hat/sqrt(2)));
        errorNumber_10 = errorNumber_10 + l10;

        
        
        l12 =length(find(E_12~=X_hat/sqrt(2)));
        errorNumber_12 = errorNumber_12 + l12;
        
        l14 =length(find(E_14~=X_hat/sqrt(2)));
        errorNumber_14 = errorNumber_14 + l14;
        
        l16 =length(find(E_16~=X_hat/sqrt(2)));
        errorNumber_16 = errorNumber_16 + l16;
    
    end

SER_4 = errorNumber_4/(2*NT*IterationTimes*Td);
SER_6 = errorNumber_6/(2*NT*IterationTimes*Td);
SER_8 = errorNumber_8/(2*NT*IterationTimes*Td);
SER_10 = errorNumber_10/(2*NT*IterationTimes*Td);
SER_12 = errorNumber_12/(2*NT*IterationTimes*Td);
SER_14 = errorNumber_14/(2*NT*IterationTimes*Td);
SER_16 = errorNumber_16/(2*NT*IterationTimes*Td);
% 
SER = [SER_4,SER_6,SER_8,SER_10,SER_12,SER_14,SER_16];
x = 0:5:30;
semilogy(x,SER,'-o');
grid on
