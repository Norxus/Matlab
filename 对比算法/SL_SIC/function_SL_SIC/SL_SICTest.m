function [X] = SL_SICTest(Y,H,M,NT_2,L,Eb)
% M为调制方式（进制），NT为发射端个数，NT_2为分割后第二部分的个数，NR为接收端个数，L为每个向量的训练次数
% IterationTimes为迭代次数，Eb为设定的SNR范围（目前只支持5个不同的信噪比进行比较）
NR = size(H,1);
NT = size(H,2);
NT_1 = NT - NT_2;

if M == 2
      X_S = [1,-1];
elseif M == 4
     X_S = [1+1i,1-1i,-1+1i,-1-1i]/sqrt(2);
elseif M ==16
    X_S = [3+1i,3-1i,3-3i,3+3i,-3-3i,-3+3i,-3-1i,-3+1i,-1+1i,-1-1i,-1+3i,-1-3i,1+3i,1-3i,1+1i,1-1i]/sqrt(10);
end

    
  %向量分割阶段
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
    
    if M == 2
        H_2R = [real(H2);imag(H2)];
    else
        H_2R = [real(H2),-imag(H2);imag(H2),real(H2)];
    end
    
    W = null(H_2R')';
    
%     X1 = X(I1);
%     X2 = X(I2);

        YR = [sign(real(Y)); sign(imag(Y))]; 
        YR_1 = W*YR;
        
          %数据训练阶段
        x_2 = zeros(NT_2,(M)^NT_2);
        for i = 1 : (M)^NT_2
            dectofour = Dec2(i,M,NT_2) + ones(NT_2,1);
            for j = 1 : NT_2
                x_2(j,i) = X_S(dectofour(j));
            end
        end
        
        x_1 = zeros(NT_1,(M)^NT_1);
        for i = 1 : (M)^NT_1
            dectofour = Dec2(i,M,NT_1) + ones(NT_1,1);
            for j = 1 : NT_1
                x_1(j,i) = X_S(dectofour(j));
            end
        end
        xk_temp = [];
        xk_temp_false = [];
        xk_2 = [];
        x_all = zeros(2*(NR-NT_2),(M)^NT*L);
        x_all_false = zeros(2*(NR-NT_2),(M)^NT*L);
        x_j_factor = 1;
        
        %生成所有可能的接收信号
%         for l = 1 : L
%             x_k_factor = 1;%条件索引
%             Z_train = sqrt(NT)*(10^(-Eb(ii)/20))*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);
% % Z_train = 0;
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
       %以上代码应该全部正确
       
       %生成所有可能的接收信号
%        x_all = zeros(2*(NR-NT_2),(M)^NT*L);
       x_all_sign = zeros(2*NR,(M)^NT*L);
       x_k_factor = 1;%条件索引
        for k_x1 = 1 : size(x_1,2)
            for j_x2 = 1 : size(x_2,2)
                for l = 1 : L
                Z_train = sqrt(NT)*(10^(-Eb/20))*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);
                y_temp = (H1*x_1(:,k_x1)+H2*x_2(:,j_x2)) + Z_train;
                y_1_kjl = [sign(real(y_temp)); sign(imag(y_temp))];
                x_all_sign(:,x_k_factor) = y_1_kjl;
%                 y_1_kjl = W*y_1_kjl;
%                 x_all(:,x_k_factor) = y_1_kjl;
                x_k_factor = x_k_factor + 1;
                end
            end
        end
        
%         x_all_sum = sum(x_all_false,1);
%         interval = M^(NT_2)*L;
% %         计算出有效向量的概率和索引
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
        
        %计算期望y_tk
        if M == 2
            y_exception = zeros(2*NR-NT_2,(M)^NT_1);
        else
            y_exception = zeros(2*(NR-NT_2),(M)^NT_1);
        end
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
        
        %对第一部分向量的检测
        
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
                Z_train = sqrt(NT)*10^(-Eb/20)*[randn(NR,1)+1i*randn(NR,1)]/sqrt(2);
                    y_temp = H1*X_detectionOne+H2*x_2(:,j_two) + Z_train;
                    y_2_jl = [sign(real(y_temp)); sign(imag(y_temp))];
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
        y_twoeffective = zeros(2*NR+1,(M)^NT*L);
        for i_calculate = 1 : size(x_2,2)
            Linterval = (i_calculate - 1)*L + 1;
            Uinterval =  i_calculate*L;
          [y_twoeffective_temp] = SearchSame( X_two_tranining(:,Linterval:Uinterval));
               
            y_twoeffective(:,Linterval:Linterval+size(y_twoeffective_temp,2)-1) = y_twoeffective_temp;
        end
        
         %计算向量2的期望
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
       X = X_finalDetection;
end