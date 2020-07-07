function [err,jac,w ] =mfrvobj_adl(atotal,x,y,polyconstr)
%Errors for Exponential MIDAS polynomial
% xmidas=x.xmidas;
% xols=x.ylagf
[dayLag]=size(x.xmidas,2);
%a=atotal(1:end-size(x.xols,2)); % Compressing [ B_0, B_1, theta_1, theta_2, b_2, b_3];
alpha=atotal(1); 
n_hvar = size(x.xmidas,3);
beta=atotal([2,2+2*n_hvar+1:2+3*n_hvar-1]);              % Extracting factor coefficients
%k1=atotal(3);                                          % Hyperparameter: Theta (1)
aols=atotal(end-size(x.xols,2)+1:end);                  % Coefficient of lag y(p)
%k2=atotal(4);                                          % Hyperparameter: Theta(2) 
%{
if length(a)==4
    k2=a(4);
end
%}
params=atotal(3:2+2*n_hvar);                                % Converting into row vector

% w=exp_weights(dayLag,k1,k2);
[xx w]=midas_X(x,'exp',params,polyconstr);
if isa(y, 'logical')
    yfit = glmval(beta,xx,'probit');                      % Generate predicted value using generalized linear model.
else
    if ~isempty(x.xols)
        yfit = alpha+xx*beta+x.xols*aols(:);
    else
        yfit = alpha + xx*beta;
    end
end

err=(y-yfit);


% Be aware that the dimension of the Jacobian matrix varies according the
% the dimension of the input data (X).
if size(x.xmidas,3)<2
    jac=-[ones(size(xx)) xx beta*jacob(@midas_X,3,x,'exp',params,polyconstr) x.xols];
else
    beta = reshape(repelem(beta,1,2),2*length(beta),1);                     % Aligning the dimension of beta
    jacob_mat = jacob(@midas_X,3,x,'exp',params,polyconstr).*beta';           % Computing the first order Jacobian matrix with respect to changing thetas.
    jac=-[ones(size(xx,1),1) xx jacob_mat x.xols];                          
end
%RSS=err'*err;
%err'*err
%a
