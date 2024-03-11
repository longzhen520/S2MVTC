function [x,objV] = wshrinkObj_weight_up(x,rho,sX, isWeight,mode,rank_low)
%rho is the weighted vector
if isWeight == 1
    %     C = 2*sqrt(2)*sqrt(sX(3)*sX(2));
    C = sqrt(sX(3)*sX(2));
end
if ~exist('mode','var')
    % mode = 1 --> lateral slice
    % mode = 2 --> front slice
    % mode = 3 --> top slice
    mode = 1;
end

X=reshape(x,sX);
if mode == 1
    Y=X2Yi(X,3);
elseif mode == 3
    Y=shiftdim(X, 2);
else
    Y = X;
end

% [n1,n2,n3]=size(Y);
% Yhat=zeros(n1,n2,n3);
% for c=1:n2
%     Yhat(:,c,:)=fft(Y(:,c,:));
% end

Yhat = fft(Y,[],3);

objV = 0;
if mode == 1
    n3 = sX(3);
elseif mode == 3
    n3 = sX(1);
else
    n3 = min(sX(3),rank_low);
end







if isinteger(n3/2)
    endValue = int16(n3/2+1);
    
    for i = 1:endValue
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');%SVD 分解
        
        if isWeight         %权重计算
            weight = C./(diag(shat) + eps);
            tau = rho*weight;
            shat = soft(shat,diag(tau));
        else
            tau = rho;
            shat = max(shat - tau,0);
        end
        
        objV = objV + sum(shat(:));%目标值累加
        Yhat(:,:,i) = uhat*shat*vhat';
        
        if i > 1%对称处理
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';
            objV = objV + sum(shat(:));
        end
    
    end


    [uhat,shat,vhat] = svd(full(Yhat(:,:,endValue+1)),'econ');%SVD分解 economic SVD

    if isWeight
        weight = C./(diag(shat) + eps);
        tau = rho*weight;
        shat = soft(shat,diag(tau));
    else
        tau = rho;
        shat = max(shat - tau,0);

    end
    

    objV = objV + sum(shat(:));
    Yhat(:,:,endValue+1) = uhat*shat*vhat';
else
    endValue = int16(n3/2+1);
    for i = 1:endValue
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');%SVD分解
        
        if isWeight
            weight = C./(diag(shat) + eps);
            tau = rho*weight;
            shat = soft(shat,diag(tau));
        else
            tau = rho;
            shat = max(shat - tau,0);

        end
        objV = objV + sum(shat(:));
        Yhat(:,:,i) = uhat*shat*vhat';
        if i > 1 
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';%%之前是n3
            objV = objV + sum(shat(:));
        end
    end
end

%% Raw IFFT
Y = ifft(Yhat,[],3);
%% New IFFT
% for c=1:n2
%     Y(:,c,:)=ifft(Yhat(:,c,:));
% end
if mode == 1
    X = Yi2X(Y,3);
elseif mode == 3
    X = shiftdim(Y, 1);
else
    X = Y;
end

x =X(:);
% x(x==0)=1;
end