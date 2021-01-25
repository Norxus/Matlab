clc
close all
clear all
syms u

IterationTimes = 3000; % ���͵ķ�����Ŀ
NT=6;      %����������
NT_2 = 1;
NT_1 = NT - NT_2;
NR=16;      %����������
Eb=-10:5:10;     %�����
E_4 = [];
E_6 = [];
E_8 = [];
E_10 = []; 
E_12 = [];
E_14 = [];
E_16 = [];

%�ܴ������ָ���
errorNumber_4 =0;
errorNumber_6 =0;
errorNumber_8 =0;
errorNumber_10 =0;
errorNumber_12 =0;
errorNumber_14 =0;
errorNumber_16 =0;

%ÿ�ε����Ĵ�����������
 l4=0;
 l6=0;
 l8=0;
 l10=0;
 l12=0;
 l14=0;
 l16=0;
 M=4; %���Ʒ�ʽM-QAM
 X_S = [1+1i,1-1i,-1+1i,-1-1i]; %���Է��͵�������Ԫ

deltaX = zeros(1,length(Eb));
b = zeros(5,length(Eb));
q = zeros(4,length(Eb));

for i = 1 : length(Eb)
    deltaX(i) = ((NT/(10^(Eb(i)/10)))+NT)/2;
    b(1,i) = -inf;
    b(2,i) = norminv(normcdf(-inf)+0.25)*sqrt(deltaX(i));
    b(3,i) = norminv(normcdf(-inf)+0.25+0.25)*sqrt(deltaX(i));
    b(4,i) = norminv(normcdf(-inf)+0.25+0.25+0.25)*sqrt(deltaX(i));
    b(5,i) = inf;
    f = (u/sqrt(2*pi*deltaX(i)))*exp(-u^2/(2*deltaX(i)));
    q(1,i) = double(int(f,u,b(1,i),b(2,i))/(normcdf(b(2,i),0,sqrt(deltaX(i)))-normcdf(b(1,i),0,sqrt(deltaX(i)))));
    q(2,i) = double(int(f,u,b(2,i),b(3,i))/(normcdf(b(3,i),0,sqrt(deltaX(i)))-normcdf(b(2,i),0,sqrt(deltaX(i)))));
    q(3,i) = double(int(f,u,b(3,i),b(4,i))/(normcdf(b(4,i),0,sqrt(deltaX(i)))-normcdf(b(3,i),0,sqrt(deltaX(i)))));
    q(4,i) = double(int(f,u,b(4,i),b(5,i))/(normcdf(b(5,i),0,sqrt(deltaX(i)))-normcdf(b(4,i),0,sqrt(deltaX(i)))));
end

 
for jj=1:IterationTimes
      jj
    ip = randn(2*NT,1)>0; %�ȸ��ʲ���0��1
    X_hat = 2*ip-1; % 0 -> -1; 1 -> 1
    X_transmit = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)];
    X = [X_hat(1:NT)+1i*X_hat(NT+1:2*NT)]/sqrt(2);
    N = [randn(NR,1)+1i*randn(NR,1)]/sqrt(2); % 0��ֵ��˹������
    H = [randn(NR,NT)+1i*randn(NR,NT)]/sqrt(2); % ����˥���ŵ�
    
  %�����ָ�׶�
    I1 = [];
    I2 = (1:NT)';
    for i_divide = 1 : NT_1
        D_chordal = [];
        for i_argmax = 1 : length(I2)
            I1 = [I1,I2(i_argmax)];
            I2_temp = zeros(length(I2),1);
            I2_i_argmax = I2(i_argmax);
            I2(i_argmax) = [];
            H1 = H(:,I1);
            H2 = H(:,I2);
            U1 = orth(H1);
            U2 = orth(H2);
            d_chordal = norm(U1*(U1')-U2*(U2'),'fro')/sqrt(2);
            D_chordal = [D_chordal,d_chordal];
            I1(end) = [];
            I2_temp(1:i_argmax-1) = I2(1:i_argmax-1);
            I2_temp(i_argmax) = I2_i_argmax;
            I2_temp(i_argmax+1:end) = I2(i_argmax:end);
            I2  = I2_temp;
        end
        k = find(D_chordal == max(D_chordal));
        I1 = [I1,I2(k)];
        I2(k) = [];
    end
    H1 = H(:,I1);
    H2 = H(:,I2);
    H_2R = [real(H2),-imag(H2);imag(H2),real(H2)];
    W = null(H_2R')';
    
    X1 = X(I1);
    X2 = X(I2);

    for ii = 1:length(Eb) %���������������µ����Y
%%
        Y = (H1*X1+H2*X2)+ sqrt(NT)*(10^(-Eb(ii)/20))*N;
        Y_temp = [real(Y); imag(Y)];
        Y_temp(find(Y_temp>b(4,ii))) = q(4,ii);
        Y_temp(find(Y_temp<b(2,ii))) = q(1,ii);
        Y_temp(find(Y_temp>b(2,ii)&Y_temp<=b(3,ii))) = q(2,ii);
        Y_temp(find(Y_temp<=b(4,ii)&Y_temp>b(3,ii))) = q(3,ii);
        YR = Y_temp;
        YR_1 = W*YR;
        
        %����ѵ���׶�
        L = 20;  %ѵ�����ݸ���
        X_factor = (randperm(40,2*(NR-NT_2)))';
        [x1] = ndgrid(X_S);
        x_2 = [x1(:)]'/sqrt(2);
        [x5,x4,x3,x2,x1] = ndgrid(X_S,X_S,X_S,X_S,X_S);
        x_1 = [x5(:),x4(:),x3(:),x2(:),x1(:)]'/sqrt(2);
        xk_temp = [];
        xk_temp_false = [];
        xk_2 = [];
        x_all = zeros(2*(NR-NT_2),(M)^NT*L);
        x_all_false = zeros(2*(NR-NT_2),(M)^NT*L);
        x_j_factor = 1;
        
%         %�������п��ܵĽ����ź�
%         for l = 1 : L
%             x_k_factor = 1;%��������
%             Z_train = sqrt(NT)*(10^(-Eb(ii)/20))*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);
%             for k_x1 = 1 : size(x_1,2)
%                 for j_x2 = 1 : size(x_2,2)
%                     y_temp = H1*x_1(:,k_x1)+H2*x_2(:,j_x2) + Z_train;
%                     y_1_kjl = [sign(real(y_temp)); sign(imag(y_temp))];
%                     y_1_kjl = w*y_1_kjl;
%                     y_1_kjl_false = X_factor.*y_1_kjl;
%                     if(k_x1 == x_k_factor)
%                        xk_temp = [xk_temp,y_1_kjl] ;
%                        xk_temp_false = [xk_temp_false,y_1_kjl_false] ;
%                        if(j_x2 == size(x_2,2))
%                          Lindex = 1 + (k_x1-1)*(M^(NT_2)*L)+(l-1)*size(x_2,2);
%                          Uindex = Lindex + M^(NT_2)-1;
%                          x_all(:,Lindex:Uindex) = xk_temp;
%                          x_all_false(:,Lindex:Uindex) = xk_temp_false;
%                          x_k_factor = x_k_factor + 1;
%                          xk_temp = [];
%                          xk_temp_false = [];
%                        end                    
%                     end                 
%                 end
%             end              
%         end       
       %���ϴ���Ӧ��ȫ����ȷ
       
       %�������п��ܵĽ����ź�
%        x_all = zeros(2*(NR-NT_2),(M)^NT*L);
       x_all_sign = zeros(2*NR,(M)^NT*L);
       x_k_factor = 1;%��������
        for k_x1 = 1 : size(x_1,2)
            for j_x2 = 1 : size(x_2,2)
                for l = 1 : L
                Z_train = sqrt(NT)*(10^(-Eb(ii)/20))*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);
                y_temp = (H1*x_1(:,k_x1)+H2*x_2(:,j_x2)) + Z_train;
                y_1_kjl_temp = [real(y_temp); imag(y_temp)];
                y_1_kjl_temp(find(y_1_kjl_temp>b(4,ii))) = q(4,ii);
                y_1_kjl_temp(find(y_1_kjl_temp<b(2,ii))) = q(1,ii);
                y_1_kjl_temp(find(y_1_kjl_temp>b(2,ii)&y_1_kjl_temp<=b(3,ii))) = q(2,ii);
                y_1_kjl_temp(find(y_1_kjl_temp<=b(4,ii)&y_1_kjl_temp>b(3,ii))) = q(3,ii);
                y_1_kjl = y_1_kjl_temp;
                x_all_sign(:,x_k_factor) = y_1_kjl;
%                 y_1_kjl = W*y_1_kjl;
%                 x_all(:,x_k_factor) = y_1_kjl;
                x_k_factor = x_k_factor + 1;
                end
            end
        end
        
%         x_all_sum = sum(x_all_false,1);
%         interval = M^(NT_2)*L;
% %         �������Ч�����ĸ��ʺ�����
%         x_unique = zeros(1,(M)^NT*L);
%         x_prosibility =  zeros(1,(M)^NT*L);
%         x_ID = zeros(1,(M)^NT*L);
%         for i_calculate = 1 : size(x_1,2)
%             Linterval = (i_calculate - 1)*interval + 1;
%             Uinterval =  i_calculate*interval;
%             unique_temp = unique(x_all_sum(:,Linterval:Uinterval));
%             x_unique(Linterval:Linterval+length(unique_temp)-1) =  unique_temp; 
%             x_prosibility(Linterval:Linterval+length(unique_temp)-1) = histc(x_all_sum(Linterval:Uinterval),unique_temp)/interval;
%             [x_isornot,x_Id]=ismember(unique_temp,x_all_sum(Linterval:Uinterval));
%             x_ID(Linterval:Linterval+length(unique_temp)-1) = x_Id;
%         end
%         y_effective = [ x_unique;x_prosibility;x_ID];
        interval = M^(NT_2)*L;
        y_effective = zeros(2*NR+1,(M)^NT*L);
        for i_calculate = 1 : size(x_1,2)
            Linterval = (i_calculate - 1)*interval + 1;
            Uinterval =  i_calculate*interval;
          [y_effective_temp] = SearchSame(x_all_sign(:,Linterval:Uinterval));
               
            y_effective(:,Linterval:Linterval+size(y_effective_temp,2)-1) = y_effective_temp;
        end
        
        y_receivetraining = x_all;
        
        %��������y_tk
        y_exception = zeros(2*(NR-NT_2),(M)^NT_1);
        for i_y_exception = 1 : (M)^NT_1
            yLinterval = (i_y_exception - 1)*interval + 1;
            yUinterval =  i_y_exception*interval;
%             y_procibility = y_effective(2,yLinterval:yUinterval);
             y_procibility = y_effective(2*NR+1,yLinterval:yUinterval)/interval;

             y_procibility(find(y_procibility==0))=[];
%             y_index = y_effective(3,yLinterval:yUinterval);
%             y_index(find(y_index==0))=[];
%             if(norm(sum(y_procibility)-1)>10^(-14))
%                 1
%             end
%             y_receivetraining_temp = y_receivetraining(1:2*(NR-NT_2),yLinterval:yUinterval);
              y_receivetraining_temp = y_effective(1:2*NR,yLinterval:yUinterval);
              y_receivetraining_temp(:,find(y_receivetraining_temp(1,:)==0)) = [];
%             y_effective_temp = y_receivetraining_temp(:,y_index);
              y_receivetraining_temp = W*y_receivetraining_temp;
            for i_y_product = 1 : length(y_procibility)
                y_receivetraining_temp(:,i_y_product) =y_procibility(i_y_product).* y_receivetraining_temp(:,i_y_product); 
            end
            y_exception(:,i_y_exception) = sum(y_receivetraining_temp,2);
        end
        
        %�Ե�һ���������ļ��
        
        vectorDetection = zeros(1,(M)^NT_1);
        for i_vectorDetector = 1 : (M)^NT_1
            detectionData = norm(YR_1 - y_exception(:,i_vectorDetector),2);
            vectorDetection(i_vectorDetector) = detectionData;
        end
        [minValue,minIndex] = min(vectorDetection);
         
         X_detectionOne = x_1(:,minIndex); 
%          X_two_tranining_false = zeros(2*NR, size(x_2,2)*L);
         X_two_tranining = zeros(2*NR, size(x_2,2)*L);
         for j_two = 1 : size(x_2,2)
             for l_two = 1 : L
                Z_train = sqrt(NT)*10^(-Eb(ii)/20)*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);
                    y_temp = H1*X_detectionOne+H2*x_2(:,j_two) + Z_train;
%                     y_2_jl_temp = [sign(real(y_temp)); sign(imag(y_temp))];
                    y_2_jl_temp = [real(y_temp); imag(y_temp)];
                    y_2_jl_temp(find(y_2_jl_temp>b(4,ii))) = q(4,ii);
                    y_2_jl_temp(find(y_2_jl_temp<b(2,ii))) = q(1,ii);
                    y_2_jl_temp(find(y_2_jl_temp>b(2,ii)&y_2_jl_temp<=b(3,ii))) = q(2,ii);
                    y_2_jl_temp(find(y_2_jl_temp<=b(4,ii)&y_2_jl_temp>b(3,ii))) = q(3,ii);
                    y_2_jl = y_2_jl_temp;
                    Lindex_two = l_two + (j_two-1)*L ;
                    X_two_tranining(:,Lindex_two) = y_2_jl;            
            end
         end
         
%          X_two_tranining_sum =  sum(X_two_tranining_false,1);   
%          x_twounique = zeros(1,size(x_2,2)*L);
%          x_twoprosibility =  zeros(1,size(x_2,2)*L);
%          x_twoID = zeros(1,size(x_2,2)*L);
%          
%          for j_two = 1 : size(x_2,2)
%             Linterval = (j_two - 1)*L + 1;
%             Uinterval =  j_two*L;
%             unique_temp = unique(X_two_tranining_sum(:,Linterval:Uinterval));
%             x_twounique(:,Linterval:Linterval+length(unique_temp)-1) =  unique_temp; 
%             x_twoprosibility(:,Linterval:Linterval+length(unique_temp)-1) = histc(X_two_tranining_sum(:,Linterval:Uinterval),unique_temp)/L;
%             [x_twoisornot,x_twoId]=ismember(unique_temp,X_two_tranining_sum(:,Linterval:Uinterval));
%             x_twoID(:,Linterval:Linterval+length(unique_temp)-1) = x_twoId;
%          end
%          y_twoeffective = [x_twounique;x_twoprosibility;x_twoID];
%          y_tworeceivetraining =  X_two_tranining;
        y_twoeffective = zeros(2*NR+1,size(x_2,2)*L);
        for i_calculate = 1 : size(x_2,2)
            Linterval = (i_calculate - 1)*L + 1;
            Uinterval =  i_calculate*L;
          [y_twoeffective_temp] = SearchSame( X_two_tranining(:,Linterval:Uinterval));
               
            y_twoeffective(:,Linterval:Linterval+size(y_twoeffective_temp,2)-1) = y_twoeffective_temp;
        end
        
         %��������2������
        y_twoexception = zeros(2*NR,size(x_2,2));
        for i_y_twoexception = 1 : (M)^NT_2
            yLinterval = (i_y_twoexception - 1)*L + 1;
            yUinterval =  i_y_twoexception*L;
%             y_procibility = y_effective(2,yLinterval:yUinterval);
             y_procibility = y_twoeffective(2*NR+1,yLinterval:yUinterval)/L;

             y_procibility(find(y_procibility==0))=[];
              y_tworeceivetraining_temp = y_twoeffective(1:2*NR,yLinterval:yUinterval);
              y_tworeceivetraining_temp(:,find(y_tworeceivetraining_temp(1,:)==0)) = [];
            for i_y_product = 1 : length(y_procibility)
                y_tworeceivetraining_temp(:,i_y_product) =y_procibility(i_y_product).* y_tworeceivetraining_temp(:,i_y_product); 
            end
            y_twoexception(:,i_y_twoexception) = sum(y_tworeceivetraining_temp,2);
        end
        
        vectorDetectionTwo = zeros(1,size(x_2,2));
        for i_vectorDetector = 1 : size(x_2,2)
            detectionData = norm(YR - y_twoexception(:,i_vectorDetector),2);
            vectorDetectionTwo(i_vectorDetector) = detectionData;
        end
        [minValueTwo,minIndexTwo] = min(vectorDetectionTwo);
         
        X_detectionTwo = x_2(:,minIndexTwo);
        X_finalDetection = zeros(NT,1);
        
        for i = 1 :NT_1
           X_finalDetection(I1(i)) = X_detectionOne(i);
        end
        
        for j = 1 :NT_2
           X_finalDetection(I2(j)) = X_detectionTwo(j);
        end
%       

   
      
       if(ii == 1)
        E_4 = [E_4,X_finalDetection];  
      end
      
      if(ii == 2)
        E_6 = [E_6,X_finalDetection];  
      end
      
      if(ii == 3)
        E_8 = [E_8,X_finalDetection];  
      end
      
      if(ii == 4)
        E_10 = [E_10,X_finalDetection];  
      end
      
      if(ii == 5)
        E_12 = [E_12,X_finalDetection];  
      end
       
%       if(ii == 1)
%         E_14 = [E_14,X_finalDetection];  
%       end
%       
%       if(ii == 2)
%         E_16 = [E_16,X_finalDetection];  
%       end
 end
    
        
        l4 =NT-length(find(E_4(:,jj)==X));
        errorNumber_4 = errorNumber_4 + l4;


       
        l6 =NT-length(find(E_6(:,jj)==X));
        errorNumber_6 = errorNumber_6 + l6;
            
        
                                                                                                                                                                                                                                                                                 
        l8 = NT-length(find(E_8(:,jj)==X));
        errorNumber_8 = errorNumber_8 + l8;

        
        
        l10 = NT-length(find(E_10(:,jj)==X));
        errorNumber_10 = errorNumber_10 + l10;
% 
%         
%         
        l12 = NT-length(find(E_12(:,jj)==X));
        errorNumber_12 = errorNumber_12 + l12;
        
%         l14 = NT-length(find(E_14(:,jj)==X));
%         errorNumber_14 = errorNumber_14 + l14;
%         
%         l16 = NT-length(find(E_16(:,jj)==X));
%         errorNumber_16 = errorNumber_16 + l16;
    
end

SER_4 = errorNumber_4/(NT*jj);
SER_6 = errorNumber_6/(NT*jj);
SER_12 = errorNumber_12/(NT*jj);
SER_8 = errorNumber_8/(NT*jj);
SER_10 = errorNumber_10/(NT*jj);
% SER_12 = errorNumber_12/(NT*jj);
% SER_16 = errorNumber_16/(NT*jj);
% 
SER = [SER_4,SER_6,SER_8,SER_10,SER_12];
x = -10:5:10;
semilogy(x,SER,'-o');
grid on