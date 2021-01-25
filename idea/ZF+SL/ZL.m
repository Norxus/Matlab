clc
close all
clear all

IterationTimes = 100; % 发送的符号数目
NT=2;      %发送天线数
NR=8;      %接受天线数
Eb=5:5:25;     %信噪比
E_f10 = [];
E_f5 = [];
E_0 = [];
E_5 = []; 
E_10 = [];


%总错误码字个数
errorNumber_f10 =0;
errorNumber_f5 =0;
errorNumber_0 =0;
errorNumber_5 =0;
errorNumber_10 =0;


%每次迭代的错误码字数量
 l_f10=0;
 l_f5=0;
 l_0=0;
 l_5=0;
 l_10=0;
 M=4; %调制方式M-QAM
 L = 3;
 T = 500;
 Td = T - (M)^NT*L


 X_S = [1+1i,1-1i,-1+1i,-1-1i]; %可以发送的所有码元

 
for jj=1:IterationTimes
      jj
      
    ip = randn(2*NT,Td)>0; %等概率产生0和1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    X_transmit = [X_hat(1:NT,:)+1i*X_hat(NT+1:2*NT,:)];
    X = X_transmit/sqrt(2);
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0均值高斯白噪声
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % 瑞利衰落信道
    y_training = zeros(NR,(M)^NT*L);
    
       
    for ii = 1:length(Eb) %计算多个信噪比情况下的输出Y 
        %Tt训练阶段
            y_tranining_complex = zeros(NR,(M)^NT*L);
            [x2,x1] = ndgrid(X_S,X_S,X_S,X_S);
            x_training = [x2(:),x1(:)]'/sqrt(2);
            for i = 1 : (M)^NT
                for j = 1 : L
                    Z_train = sqrt(NT)*(10^(-(Eb(ii)+3)/20))*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);
                    y_training_complex(:,(i-1)*L+j) = H*x_training(:,i)+Z_train;
                end               
            end

            y_training = [sign(real(y_training_complex));sign(imag(y_training_complex))];
            factor = (1:2*NR)';
            y_training_sum = sum(y_training.*factor);
            y_unique = zeros(1,(M)^NT*L);
            y_prosibility =  zeros(1,(M)^NT*L);
            y_ID = zeros(1,(M)^NT*L);
            E_y = zeros(2*NR,(M)^NT);

             for i = 1 : (M)^NT
                Linterval = (i - 1)*L + 1;
                Uinterval =  i*L;
                temp = y_training(:,Linterval:Uinterval);
                unique_temp = unique(y_training_sum(:,Linterval:Uinterval));
                y_unique(Linterval:Linterval+length(unique_temp)-1) =  unique_temp; 
                y_prosibility(Linterval:Linterval+length(unique_temp)-1) = histc(y_training_sum(Linterval:Uinterval),unique_temp)/L;
                [y_isornot,y_Id]=ismember(unique_temp,y_training_sum(Linterval:Uinterval));
                y_ID(Linterval:Linterval+length(unique_temp)-1) = y_Id;
                E_y(:,i)=sum(temp(:,y_Id).*y_prosibility(Linterval:Linterval+length(unique_temp)-1),2);
             end
        % Td数据检测阶段
        for kk = 1:Td
           
            Y = H*X(:,kk)+ sqrt(NT)*(10^(-Eb(ii)/20))*N;     
            
            YR = [sign(real(Y)); sign(imag(Y))]; 
            vectorDetection = zeros(1,(M)^NT);
            for i = 1 : (M)^NT
                detectionData = norm(YR - E_y(:,i),2);
                vectorDetection(i) = detectionData;
            end
            [minValue,minIndex] = min(vectorDetection);
            X_finalDetection = x_training(:,minIndex);
       
       if(ii == 1)
        E_f10 = [E_f10,X_finalDetection];  
        l_f10 =NT-length(find(X_finalDetection==X(:,kk)));
        errorNumber_f10 = errorNumber_f10 + l_f10;
      end
      
      if(ii == 2)
        E_f5 = [E_f5,X_finalDetection]; 
        l_f5 =NT-length(find(X_finalDetection==X(:,kk)));
        errorNumber_f5 = errorNumber_f5 + l_f5;
      end
      
      if(ii == 3)
        E_0 = [E_0,X_finalDetection];  
        l_0 = NT-length(find(X_finalDetection==X(:,kk)));
        errorNumber_0 = errorNumber_0 + l_0;
      end
      
      if(ii == 4)
        E_5 = [E_5,X_finalDetection];
        l_5 = NT-length(find(X_finalDetection==X(:,kk)));
        errorNumber_5 = errorNumber_5 + l_5;
      end
      
      if(ii == 5)
        E_10 = [E_10,X_finalDetection];  
        l_10 = NT-length(find(X_finalDetection==X(:,kk)));
        errorNumber_10 = errorNumber_10 + l_10;
      end
   end
      
 end
             
end

SER_4 = errorNumber_f10/(NT*Td*IterationTimes);
SER_6 = errorNumber_f5/(NT*Td*IterationTimes);
SER_8 = errorNumber_0/(NT*Td*IterationTimes);
SER_10 = errorNumber_5/(NT*Td*IterationTimes);
SER_12 = errorNumber_10/(NT*Td*IterationTimes);

% 
SER = [SER_4,SER_6,SER_8,SER_10,SER_12];
x = 5:5:25;
semilogy(x,SER,'-o');
grid on

