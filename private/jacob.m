%Helper function for MIDAS toolbox. Computes Jacobian of function handler
%fname. num is a place of parameters in fname.
%Arthur Sinko July 2010

function out=jacob(fname,num,varargin)          %num represents the hyperparameters in the poly function.
params0=varargin;
out=[];
eps=1e-6;
%warning('Current stepsize = %.9f \n ', eps)
dim_x = size(params0{1}(1).xmidas,3);

if dim_x < 2
    for i=1:length(params0{num})
        params1=params0;
        params2=params0;
        params1{num}(i)=params0{num}(i)+eps/2;
        params2{num}(i)=params0{num}(i)-eps/2;
        temp=(feval(fname,params1{:})-feval(fname,params2{:}))/eps;     % This computes the weight and 
        out=[out temp];
    end
else
    theta_vec = reshape(params0{num},2, dim_x);                         % Reshape for manipulation
    for i = 1:2                                                         % Here: assuming only two hyperparameters only ( Further improvement needed ).
        params1=params0;
        params2=params0;
        theta_up = theta_vec;
        theta_down = theta_vec;
        theta_up(i,:) = theta_vec(i,:) + eps/2;                                  % Increasing theta
        params1{num}=reshape(theta_up, dim_x*2,1);
        theta_down(i,:) = theta_vec(i,:) - eps/2;
        params2{num}=reshape(theta_down, dim_x*2,1);
        x_up = feval(fname,params1{:});
        x_down = feval(fname,params2{:});
        temp=(x_up - x_down)/eps;
        out=[out temp];
    end
end