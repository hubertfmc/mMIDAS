function [OutputForecast,OutputEstimate,MixedFreqData,ExtendedForecast] = m1MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,varargin)

% Parse inputs and set defaults
callerName = 'MIDAS_ADL';
parseObj = inputParser;
parseObj.addParameter('Xlag',9,@(x)validateattributes(x,{'numeric','char'},{},callerName));
parseObj.addParameter('Ylag',1,@(x)validateattributes(x,{'numeric','char','cell'},{},callerName));
parseObj.addParameter('Horizon',1,@(x)validateattributes(x,{'numeric','char'},{},callerName));
parseObj.addParameter('EstStart',[],@(x)validateattributes(x,{'numeric','char'},{},callerName));
parseObj.addParameter('EstEnd',[],@(x)validateattributes(x,{'numeric','char'},{},callerName));
parseObj.addParameter('ExoReg',[],@(x)validateattributes(x,{'numeric'},{'2d'},callerName));
parseObj.addParameter('ExoRegDate',{},@(x)validateattributes(x,{'cell'},{},callerName));
parseObj.addParameter('Method','FixedWindow',@(x)validateattributes(x,{'char'},{},callerName));
parseObj.addParameter('Polynomial','Beta',@(x)validateattributes(x,{'char'},{},callerName));
parseObj.addParameter('PolyStepFun',[],@(x)validateattributes(x,{'numeric'},{'2D'},callerName));
parseObj.addParameter('AlmonDegree',[],@(x)validateattributes(x,{'numeric'},{'scalar','integer','nonnegative'},callerName));
parseObj.addParameter('Discount',0.9,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<=',1},callerName));
parseObj.addParameter('DiscountIncrease',true,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
parseObj.addParameter('Display','full',@(x)validateattributes(x,{'char'},{},callerName));
parseObj.addParameter('PlotWeights',true,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
parseObj.parse(varargin{:});
xlag = parseObj.Results.Xlag;
ylag = parseObj.Results.Ylag;
horizon = parseObj.Results.Horizon;
estStart = parseObj.Results.EstStart;
estEnd = parseObj.Results.EstEnd;
ExoReg = parseObj.Results.ExoReg;
ExoRegDate = parseObj.Results.ExoRegDate;
estMethod = parseObj.Results.Method;
polynomial = parseObj.Results.Polynomial;
stepfun = parseObj.Results.PolyStepFun;
almonDegree = parseObj.Results.AlmonDegree;
discount = parseObj.Results.Discount;
discountIncrease = parseObj.Results.DiscountIncrease;
dispFlag = parseObj.Results.Display;
plotWeightsFlag = parseObj.Results.PlotWeights;
if isempty(almonDegree)
    almonDegree = 3;
end
if strcmpi(dispFlag,'full') || strcmpi(dispFlag,'time')
    dispTime = 1;
else
    dispTime = 0;
end
if strcmpi(dispFlag,'full') || strcmpi(dispFlag,'estimate')
    dispNLS = 1;
else
    dispNLS = 0;
end

% Handle lagged Y
iscellYlag = iscell(ylag);
if iscellYlag    
    ylagVec = cell2mat(ylag);
    ylagVec = ylagVec(:)';
    ylag = max(ylagVec);    
else
    ylagVec = 1:ylag;
end

% Prepare mixed frequency regression data
MixedFreqData = mMixFreqData(DataY,DataYdate,DataX,DataXdate,xlag,ylag,horizon,estStart,estEnd,dispTime);
EstY = MixedFreqData.EstY;
EstX = MixedFreqData.EstX;
EstYdate = MixedFreqData.EstYdate;
EstXdate = MixedFreqData.EstXdate;
EstLagY = MixedFreqData.EstLagY;
OutY = MixedFreqData.OutY;
OutX = MixedFreqData.OutX;
OutYdate = MixedFreqData.OutYdate;
OutXdate = MixedFreqData.OutXdate;
OutLagY = MixedFreqData.OutLagY;
xlag = MixedFreqData.Xlag;
nobs = size(EstY,1);
nforecast = size(OutY,1);

% MixFreqData works on 1:ylag, but we extract the required lagged data
if iscellYlag
    EstLagY = EstLagY(:,ylagVec);
    OutLagY = OutLagY(:,ylagVec);
    if dispTime && length(ylagVec) < ylag
        disp(['Note: Ylag = ',num2str(ylagVec),'. Regression time frame shows more lags, but it does not affect estimation.'])
    end
end

% Handle exogenous regressors
nExoReg = size(ExoReg,2);
if nExoReg > 0
    
    if size(ExoReg,1) ~= size(ExoRegDate,1)
        error('Exogenous data and dates must match.')
    end
    
    if size(ExoReg,1) < size(DataY,1)
        error('Exogenous time series data cannot be shorter than DataY.')
    end
    
    EstExoReg = zeros(nobs,nExoReg);
    OutExoReg = zeros(nforecast,nExoReg);
    
    % Exogenous regressors mimic x by setting Horizon = 0
    date_ck = sum(datevec(ExoRegDate,'dd/mm/yyyy')== datevec(DataYdate,'dd/mm/yyyy'),1)==size(DataYdate,1);
    if sum(sum(date_ck,1),2)==6
        if length(ExoReg)~= length(EstY)
            EstExoReg(:,:) = ExoReg(2:end,:);
        else
            EstExoReg(:,:) = ExoReg(:,:);
        end
        %{
        for t = 1:nExoReg
            Mimic = mMixFreqData(DataY,DataYdate,ExoReg(:,t),ExoRegDate,1,ylag,0,estStart,estEnd,false);
            if any(EstY ~= Mimic.EstY)
                error('Unable to process exogenous regressors. Make sure the time series is long enough.')
            end
            if nforecast < length(Mimic.OutY)
                OutExoReg(:,t) = Mimic.OutX(1:nforecast,:);
            else
                OutExoReg(:,t) = Mimic.OutX;
            end        
            EstExoReg(:,t) = Mimic.EstX;
        end
        %}
    else
        error('Exogenous variables and explanatory variables do not have the same length of time.');
    end
    % Append exogenous regressors to the lagged Y
    EstLagY = [EstLagY,EstExoReg];
    OutLagY = [OutLagY,OutExoReg];
end

switch lower(estMethod)
    case 'fixedwindow'        
        OutputEstimate = MIDAS_Estimate(EstY,EstX,EstLagY,EstXdate,polynomial,stepfun,almonDegree,nExoReg,ylagVec);
        %OutputForecast = MIDAS_Forecast(OutY,OutX,OutLagY,OutputEstimate.estParams,OutputEstimate.estWeights,discount,OutYdate,discountIncrease);               
    %{
    case {'rollingwindow','recursive'}        
        nobs = size(EstY,1);
        nroll = size(OutY,1);
        if nroll == 0
            error('Rolling window does not apply because there are no rolling periods. Decrease "EstEnd".')
        end
        Ybig = [EstY;OutY];
        Xbig = [EstX;OutX];
        LagYbig = [EstLagY;OutLagY];
        Xdatebig = [EstXdate;OutXdate];
        Ydatebig = [EstYdate;OutYdate];
        Yf = zeros(nroll,1);
        YfDate = zeros(nroll,1);
        for t = 1:nroll
            if strcmpi(estMethod,'rollingwindow')
                EstYroll = Ybig(t:nobs-1+t,:);
                EstXroll = Xbig(t:nobs-1+t,:);
                EstLagYroll = LagYbig(t:nobs-1+t,:);
                EstXdateroll = Xdatebig(t:nobs-1+t,:);                
            else
                EstYroll = Ybig(1:nobs-1+t,:);
                EstXroll = Xbig(1:nobs-1+t,:);
                EstLagYroll = LagYbig(1:nobs-1+t,:);
                EstXdateroll = Xdatebig(1:nobs-1+t,:);
            end            
            OutYroll = Ybig(nobs+t,:);
            OutXroll = Xbig(nobs+t,:);
            OutLagYroll = LagYbig(nobs+t,:);
            OutYdateroll = Ydatebig(nobs+t,:);
            OutputEstimate = MIDAS_Estimate(EstYroll,EstXroll,EstLagYroll,EstXdateroll,polynomial,stepfun,almonDegree,nExoReg,ylagVec);
            %OutputForecast = MIDAS_Forecast(OutYroll,OutXroll,OutLagYroll,OutputEstimate.estParams,OutputEstimate.estWeights,discount,OutYdateroll,discountIncrease);
            Yf(t) = OutputForecast.Yf;
            YfDate(t) = datenum(OutputForecast.YfDate);
        end
                
        gap = OutY - Yf;
        seq = (0:length(gap)-1)';
        if ~discountIncrease
            seq = flipud(seq);
        end
        gapDiscount = gap .* discount.^seq;
        MSFE = mean(gap.*gap);
        RMSE = sqrt(MSFE);
        DMSFE = mean(gapDiscount.*gapDiscount);
        %OutputForecast = struct('Yf',Yf,'YfDate',datestr(YfDate),'RMSE',RMSE,'MSFE',MSFE,'DMSFE',DMSFE);
      %}          
    otherwise
        error('Estimation method must be fixedWindow, rollingWindow or recursive')
        
end

% A copy of AIC, BIC for OutputForecast
%OutputForecast.aic = OutputEstimate.aic;
%OutputForecast.bic = OutputEstimate.bic;

if dispNLS
    if any(isnan(OutputEstimate.se))
        warning('The estimated parameter covariance matrix is not positive definite.')        
    end
        
    nreg = length(OutputEstimate.estParams);   
    result = cell(nreg+1,3);
    result(1,:)={' Estimator','SE','t-stat'};
    for t = 1:nreg
        %result(t+1,1) = OutputEstimate.paramName(t);
        result(t+1,:) = {OutputEstimate.estParams(t),OutputEstimate.se(t),OutputEstimate.tstat(t)};
    end
    disp(' ')
    disp(OutputEstimate.model)
    disp(result)
    
    disp(['Goodness of fit:   ',num2str(OutputEstimate.r2)])
    disp(['Noise variance:    ',num2str(OutputEstimate.sigma2)])
    disp(['Log likelihood:    ',num2str(OutputEstimate.logL)])
    disp(['Akaike criteria:   ',num2str(OutputEstimate.aic)])
    disp(['Bayesian criteria: ',num2str(OutputEstimate.bic)])
    disp(' ')
end

var_name={};
for nb_fact = 1: size(DataX,2)
    var_name =[var_name, strcat('HiFreq Var',num2str(nb_fact))];
end

if plotWeightsFlag
    w_plot = figure(10);
    set(w_plot,'name','PolyWeight')
    weights = OutputEstimate.estWeights;
    plot(1:xlag,weights);
    title(OutputEstimate.model);
    xlabel('Lag');
    ylabel('Weights');
    xlim([1,xlag]);
    legend(var_name);
end


%-----------------
% Extended out-of-sample forecast
% [ Hubert: Temporary suspend forecasting function.]
%-----------------
ExtendedForecast ='NA';
OutputForecast = 'NA';
%{
ExtendedForecast = MIDAS_ExtendForecast(DataY,DataYdate,DataX,DataXdate,ExoReg,ExoRegDate,...
    OutputEstimate.estParams,OutputEstimate.estWeights,xlag,ylag,horizon,estStart,ylagVec);

%TEMP SUSPENSE of forecasting function to examine the regression result

Reformat serial dates to string dates
if ~isempty(MixedFreqData.EstYdate)
    MixedFreqData.EstYdate = reshape(cellstr(datestr(MixedFreqData.EstYdate)),size(MixedFreqData.EstYdate));
end
if ~isempty(MixedFreqData.EstXdate)
    MixedFreqData.EstXdate = reshape(cellstr(datestr(MixedFreqData.EstXdate)),size(MixedFreqData.EstXdate));
end
if ~isempty(MixedFreqData.EstLagYdate)
    MixedFreqData.EstLagYdate = reshape(cellstr(datestr(MixedFreqData.EstLagYdate)),size(MixedFreqData.EstLagYdate));
end
if ~isempty(MixedFreqData.OutYdate)
    MixedFreqData.OutYdate = reshape(cellstr(datestr(MixedFreqData.OutYdate)),size(MixedFreqData.OutYdate));
end
if ~isempty(MixedFreqData.OutXdate)
    MixedFreqData.OutXdate = reshape(cellstr(datestr(MixedFreqData.OutXdate)),size(MixedFreqData.OutXdate));
end
if ~isempty(MixedFreqData.OutLagYdate)
    MixedFreqData.OutLagYdate = reshape(cellstr(datestr(MixedFreqData.OutLagYdate)),size(MixedFreqData.OutLagYdate));
end
if ~isempty(OutputForecast.YfDate)
    OutputForecast.YfDate = cellstr(OutputForecast.YfDate);
end
if ~isempty(ExtendedForecast.YfDate)
    ExtendedForecast.YfDate = cellstr(ExtendedForecast.YfDate);
end
%}
end


%-------------------------------------------------------------------------
function OutputEstimate = MIDAS_Estimate(EstY,EstX,EstLagY,EstXdate,polynomial,stepfun,almonDegree,nExoReg,ylagVec)

% Parameter estimation of ADL-MIDAS model
% Support weight polynomials of beta, betaNN,expAlmon,umidas,step,Almon

% Convert to original MIDAS package data format
y = EstY;
x.xmidas = EstX;
x.xmidasd = EstXdate;
x.xols = EstLagY;
output_type = 'struct';
polyconstr = 'es';

% Regressor dimension
xlag = size(EstX,2);
ylag = size(EstLagY,2);

% Estimate ADL-MIDAS model using functions in the original MIDAS package
switch lower(polynomial)
    %{
    case 'beta'        
        MdlName = 'MIDAS: Normalized beta density with a zero last lag';
        ParamName = {'Const','HighFreqSlope','Beta1','Beta2'}; 
        numofparams = 2;
        Results = bnls_adl_new(y,x,numofparams,output_type,polyconstr);
        weightsNN = Results.weightsNN(1:xlag);
    case 'betann'        
        MdlName = 'MIDAS: Normalized beta density with a non-zero last lag';
        ParamName = {'Const','HighFreqSlope','Beta1','Beta2','Beta3'};
        numofparams = 2;
        Results = bnlsNN_adl_new(y,x,numofparams,output_type,polyconstr);
        weightsNN = Results.weightsNN(1:xlag);
    %}
    case 'expalmon'                                                         % Focus application (1)
        MdlName = 'MIDAS: Normalized exponential Almon lag polynomial';
        ParamName = {'Const','HighFreqSlope','ExpAlmon1','ExpAlmon2'};
        numofparams = 2; % Not sure if we need to change this number for more hyperparameter
        exp0=[-1,0]';                                                        % Initial value of theta
        if size(x.xmidas,3)>1
            exp0=repmat(exp0,size(x.xmidas,3),1);                            % Create extra hyperparameters for other variables
        end
        if isa(y,'logical')
            fprintf('Initiating probit regression \n')
            a00 =  glmfit([midas_X(x,'exp',exp0,polyconstr) x.xols],y,'binomial','link','probit'); % When dealing with binary explanatory inputs.
        elseif isa(y,'numeric')
            a00=ols(y,[midas_X(x,'exp',exp0,polyconstr) x.xols]);               % Creating variables matrix with regressors and lag_Y.
        end
        fprintf("Preliminary OLS result: %.5f \n", a00)
        LB0=[-10000 -10000 -100 -100]';
        UB0=[10000 10000 10 0.5]';
        if numofparams==1
            exp0=exp0(1);
            LB0=LB0(1:3);
            UB0=UB0(1:3);
        end
        a0 =[a00(1:2); exp0; a00(3:end)];
        % LB = [LB0(1:2); -Inf; LB0(3:end)];
        % UB = [UB0(1:2);  Inf; UB0(3:end)];
        %{
        nEndReg = size(x.xmidas,3)+ylag+1;                                  %Indicating the number of endorgenous variables { y_lag, X, Beta_0}
        if nEndReg ~= length(a00)
            nExoReg = length(a00(nEndReg+1:end));
            fprintf('Found %.1f exogenous variables \n',nExoReg);
        else
            nExoReg = 0;
            fprintf('Found NO exogenous variables \n');
        end
        %}
        LB = [LB0;-Inf(nExoReg,1)];
        UB = [UB0;Inf(nExoReg,1)];
        options =optimset('Display','on','TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'MaxFunEvals',1000,'Jacobian','on');
        Results = enls1_adl_new(y,x,numofparams,output_type,polyconstr,a0,options,LB,UB);
        weightsNN = Results.weightsNN(1:xlag,:);                            % Not sure why we only take out the first xlag coefficients?
    %{
    case 'umidas'        
        MdlName = 'MIDAS: Unrestricted coefficients (U-MIDAS)';
        ParamName = cell(1,1+xlag);
        ParamName{1} = 'Const';
        for m = 1:xlag
            ParamName{m+1} = ['Umidas',num2str(m)];
        end
        stepfun = 1:xlag;
        Results = midas_sf_adl_new(y,x,stepfun,output_type,polyconstr);
        estParams = Results.aopt;
        weightsNN = zeros(xlag,1);
        weightsNN(1:stepfun(1)) = estParams(2);
        for m = 2:length(stepfun)
            weightsNN(stepfun(m-1)+1:stepfun(m)) = estParams(m+1);
        end
    case 'step'        
        if isempty(stepfun)
            stepfun = 3:3:xlag;
        end
        stepfun(stepfun > xlag) = xlag;
        if stepfun(end) < xlag
            stepfun = [stepfun,xlag];
        end
        MdlName = ['MIDAS: Step function at lags ',num2str(stepfun)];
        nstep = length(stepfun);
        ParamName = cell(1,1+nstep);
        ParamName{1} = 'Const';
        for m = 1:nstep
            ParamName{m+1} = ['Step',num2str(m)];
        end
        Results = midas_sf_adl_new(y,x,stepfun,output_type,polyconstr);
        estParams = Results.aopt;
        weightsNN = zeros(xlag,1);
        weightsNN(1:stepfun(1)) = estParams(2);
        for m = 2:length(stepfun)
            weightsNN(stepfun(m-1)+1:stepfun(m)) = estParams(m+1);
        end        
    case 'almon'        
        MdlName = ['MIDAS: Almon lag polynomial of order ',num2str(almonDegree)];
        ParamName = cell(1,2+almonDegree);
        ParamName{1} = 'Const';
        ParamName{2} = 'AlmonDegree0';
        for m = 1:almonDegree
            ParamName{m+2} = ['AlmonDegree',num2str(m)];
        end
        Results = almon_adl_new(y,x,almonDegree,output_type,polyconstr);
        estParams = Results.aopt;
        AlmonMat = zeros(xlag,almonDegree+1);
        for m = 1:xlag
            AlmonMat(m,:) = m .^ (0:almonDegree);
        end
        weightsNN = AlmonMat * estParams(2:end-ylag);
        %}
end

% Collect the results
estParams = Results.aopt;
EstParamCov = Results.estvar;
E_x = 0;
for nb_fact = 1:size(EstX,3)
    E_x = E_x + EstX(:,:,nb_fact)*weightsNN(:,nb_fact);     % Compute expected value of regressors
end
yfit = estParams(1) + E_x + EstLagY * estParams(end-ylag+1:end);
resid = EstY - yfit;
rss = resid' * resid;
nobs = length(resid);
nvar = length(estParams);
sigma2 = rss / (nobs-nvar-1);
se2 = diag(EstParamCov);
if any(se2<0)
    se2(se2<0) = NaN;
end
se = sqrt(se2);
tstat = estParams(:) ./ se;
logMSE = log(rss/nobs);
logL = -nobs/2 * (1 + log(2*pi) + logMSE);
aic = nobs*logMSE + 2*nvar;
bic = nobs*logMSE + nvar*log(nobs);
r2 = 1 - var(resid)/var(EstY);

% Complete parameter names
ParaNamesPart = cell(1,ylag);
count = 0;
for m = ylagVec
    count = count + 1;
    ParaNamesPart{count} = ['Ylag',num2str(m)];
end
for m = 1:nExoReg
    count = count + 1;
    ParaNamesPart{count} = ['ExoReg',num2str(m)];
end    
ParamName = [ParamName,ParaNamesPart];

% Create an output struct
OutputEstimate = struct('model',MdlName,'paramName',{ParamName},...
    'estParams',estParams,'EstParamCov',EstParamCov,'se',se,'tstat',tstat,...
    'sigma2',sigma2,'yfit',yfit,'resid',resid,'estWeights',weightsNN,...
    'logL',logL,'r2',r2,'aic',aic,'bic',bic);

end


%-------------------------------------------------------------------------
function OutputForecast = MIDAS_Forecast(OutY,OutX,OutLagY,estParams,estWeights,discount,OutYdate,discountIncrease)

% Out of sample forecast of the ADL-MIDAS model

ylag = size(OutLagY,2);
if ~isempty(OutY)
    Yf = estParams(1) +  OutX * estWeights + OutLagY * estParams(end-ylag+1:end);
    gap = OutY - Yf;
    seq = (0:length(gap)-1)';
    if ~discountIncrease
        seq = flipud(seq);
    end
    gapDiscount = gap .* discount.^seq;
    MSFE = mean(gap.*gap);
    RMSE = sqrt(MSFE);
    DMSFE = mean(gapDiscount.*gapDiscount);
else
    Yf = [];
    MSFE = [];
    RMSE = [];
    DMSFE = [];
end

OutputForecast = struct('Yf',Yf,'YfDate',datestr(OutYdate),'RMSE',RMSE,'MSFE',MSFE,'DMSFE',DMSFE);

end

%-------------------------------------------------------------------------
function ExtendedForecast = MIDAS_ExtendForecast(DataY,DataYdate,DataX,DataXdate,ExoReg,ExoRegDate,estParams,estWeights,xlag,ylag,horizon,estStart,ylagVec)

warning('off','MixFreqData:EstStartOutOfBound');
warning('off','MixFreqData:EstEndOutOfBound');

DataYBig = DataY;
DataYdateBig = datenum(DataYdate);
step = mean(diff(DataYdateBig));

Yf = zeros(1000,1);
YfDate = zeros(1000,1);
for t = 1:1000
    
    % Pad one-period ahead observation with zero 
    currentDate = DataYdateBig(end);
    newDate = currentDate + step;
    DataYBig = [DataYBig;0];     %#ok<AGROW>
    DataYdateBig = [DataYdateBig;newDate]; %#ok<AGROW>
    
    % Extract lagged and high-freq regressors for one-period forecast
    MixedFreqData = mMixFreqData(DataYBig,DataYdateBig,DataX,DataXdate,xlag,ylag,horizon,estStart,currentDate,false);    
    OutX = MixedFreqData.OutX;    
    OutLagY = MixedFreqData.OutLagY;
    OutYdate = MixedFreqData.OutYdate;
    if size(OutX,1) > 1
        error('Bug!')
    end
    if size(OutX,1) == 0
        Yf = Yf(1:t-1);
        YfDate = YfDate(1:t-1);
        break
    end
    if length(ylagVec) ~= max(ylagVec)         
        OutLagY = OutLagY(:,ylagVec);          
    end
    
    % Extract exogenous regressors for one-period forecast, mimic x by setting Horizon = 0 
    nExoReg = size(ExoReg,2);
    OutExoReg = zeros(1,nExoReg);
    if nExoReg > 0
        Mimic = mMixFreqData(DataYBig,DataYdateBig,ExoReg(:,1),ExoRegDate,1,ylag,0,estStart,currentDate,false);
        if size(Mimic.OutX,1) == 0
            Yf = Yf(1:t-1);
            YfDate = YfDate(1:t-1);
            break
        end        
        OutExoReg(:,1) = Mimic.OutX;
    end
    for m = 2:nExoReg
        Mimic = mMixFreqData(DataYBig,DataYdateBig,ExoReg(:,m),ExoRegDate,1,ylag,0,estStart,currentDate,false);                  
        OutExoReg(:,m) = Mimic.OutX;
    end   
    OutLagY = [OutLagY,OutExoReg]; %#ok<AGROW>
    
    % One-period forecast
    Yf(t) = estParams(1) +  OutX * estWeights + OutLagY * estParams(end-length(ylagVec)-nExoReg+1:end);
    YfDate(t) = OutYdate;
    
    % For AR type recursive forecast
    DataYBig(end) = Yf(t);    
end

ExtendedForecast = struct('Yf',Yf,'YfDate',datestr(YfDate));

warning('on','MixFreqData:EstStartOutOfBound');
warning('on','MixFreqData:EstEndOutOfBound');

end


