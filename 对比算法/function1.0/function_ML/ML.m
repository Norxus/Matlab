function [X] = ML(Y,H,M)
NR = size(H);
NT = size(H,2);
if M == 2
      X_S = [1,-1];
elseif M == 4
     X_S = [1+1i,1-1i,-1+1i,-1-1i]/sqrt(2);
elseif M ==16
    X_S = [3+1i,3-1i,3-3i,3+3i,-3-3i,-3+3i,-3-1i,-3+1i,-1+1i,-1-1i,-1+3i,-1-3i,1+3i,1-3i,1+1i,1-1i]/sqrt(10);
end
x = zeros(NT,(M)^NT);
for i = 1 : (M)^NT
    dec2 = Dec2(i,M,NT) + ones(NT,1);
    for j = 1 : NT
        x(j,i) = X_S(dec2(j));
    end
end
if M == 2
    X_ML = x;
else
    X_ML = [real(x);imag(x)];
end        
        if M == 2
            H_hat = [real(H);imag(H)];
        else
            H_hat = [real(H),-imag(H);imag(H),real(H)];
        end
        
        if M == 2
            G = zeros(NT,2*NR);
        else
            G = zeros(2*NT,2*NR);
        end
        
        x_ML_detection = zeros(1,size(X_ML,2));
        
        for i = 1 : 2*NR
           G(:,i) =  H_hat(i,:)*Y_hat(i);
        end
        
        for i = 1 : size(X_ML,2)
           x_ML_detection(i) = prod(normcdf(sqrt(2*p)*(G')*X_ML(:,i),0,1));
        end
        
        Y_mid_OneBit  = X_ML(:,find(x_ML_detection == max(x_ML_detection)));
        if M == 2
            X = Y_mid_OneBit;
        else
            X = [ Y_mid_OneBit(1:NT)+1i* Y_mid_OneBit(NT+1:2*NT)];
        end
end