function Output = mMixFreqData(DataY,DataYdate,DataX,DataXdate,xlag,ylag,horizon,estStart,estEnd,dispFlag)


% Complete data
mask = ~isnan(DataY);
DataY = DataY(mask);
DataYdate = DataYdate(mask);
%mask = ~isnan(DataX);

e_r = sum(~isnan(DataX),2)==size(DataX,2);
DataX = DataX(e_r,:);
DataXdate(e_r,:);
%{
if isempty(e_r) ~= 1
    DataX(e_r,:)=[];
    DataXdate(e_r,:)=[];            %Deleting any rows with 'NAN'
else
    fprintf('Found no "NAN" entries.\n')
end
%}
DataY = DataY(:);
DataYdate = DataYdate(:);
DataX = DataX(:,:);                             % H: Changed
DataXdate = DataXdate(:);

% Convert dates to numerical vectors
formatIn = 'yyyy/MM/dd';
DataYdateVec = datevec(DataYdate,formatIn);
DataXdateVec = datevec(DataXdate,formatIn);
DataYdateNum = datenum(DataYdateVec);
DataXdateNum = datenum(DataXdateVec);
estStart = datenum(estStart);
estEnd = datenum(estEnd);

% Frequency of Y and X
dateFormat = {'year(s)','month(s)','day(s)','hour(s)','minute(s)','second(s)'};
[periodY,unitY] = DateFreq(DataYdateVec);
[periodX,unitX] = DateFreq(DataXdateVec);

% Convert descriptive lags to numeric lags
ylag = lagNum(ylag,periodY,unitY);
xlag = lagNum(xlag,periodX,unitX);
horizon = lagNum(horizon,periodX,unitX);
if ylag < 0
    error('ylag cannot be negative.')
end
if xlag < 0
    error('xlag cannot be negative')
end

% Minimum and maximum dates that data support
%(Making sure the alignment of date in different regressors.)
minDateY = DataYdateNum(ylag+1,:);
minDateX = DataXdateNum(max(1,xlag+horizon),:);
if minDateY > minDateX
    minDate = minDateY;
else
    minDate = minDateX;
end

maxDateY = DataYdateNum(end,:);
maxDateX = DataXdateNum(end,:);
if horizon < 0
    maxDateX = DataXdateVec(end,:);
    maxDateX(unitX) = maxDateX(unitX) + periodX * horizon;
    maxDateX = datenum(maxDateX);
    
    % Test: not out of bound
    %{
    for m = 1:50
        locEnd = find(DataYdateNum >= maxDateX-1e-12, 1);
        loc = find(DataXdateNum >= DataYdateNum(locEnd)-1e-12, 1);
        if loc-horizon > size(DataX,1)
            maxDateX = DataXdateVec(end,:);
            maxDateX(unitX) = maxDateX(unitX) + periodX * horizon - m;
            maxDateX = datenum(maxDateX);
        else
            break
        end
    end
    %}
end
if maxDateY > maxDateX
    maxDate = maxDateX;
else
    maxDate = maxDateY;
end

% Check and set default sample period
if isempty(estStart)
    estStart = minDate;
elseif estStart < minDate
    warning('MixFreqData:EstStartOutOfBound',...
    'Start date cannot be earlier than %s due to lagged regressors. Reset start date to %s',...
        datestr(minDate), datestr(minDate))
    estStart = minDate;
end

if isempty(estEnd)
   estEnd = maxDate;
elseif estEnd > maxDate
   warning('MixFreqData:EstEndOutOfBound',...
       'Terminal date cannot be later than %s. Reset terminal date to %s',...
       datestr(maxDate),datestr(maxDate))
   estEnd = maxDate;
end

% Construct Y data
tol = 1e-10;                                        % H: Abitrary?
locStart = find(DataYdateNum >= estStart-tol, 1);
locEnd = find(DataYdateNum >= estEnd-tol, 1);
EstY = DataY(locStart:locEnd);
EstYdate = DataYdateNum(locStart:locEnd);           % Focus on certain time period.

locForecastEnd = find(DataYdateNum >= maxDate-tol, 1);
OutY = DataY(locEnd+1:locForecastEnd);
OutYdate = DataYdateNum(locEnd+1:locForecastEnd);
nobs = locEnd - locStart + 1;
nforecast = size(OutY,1);

% Construct lagged Y data
EstLagY = zeros(nobs,ylag);
EstLagYdate = zeros(nobs,ylag);
for m = 1:ylag
    EstLagY(:,m) = DataY(locStart-m:locEnd-m);              % Each Y_t represents a column vector.
    EstLagYdate(:,m) = DataYdateNum(locStart-m:locEnd-m);
end
OutLagY = zeros(nforecast,ylag);
OutLagYdate = zeros(nforecast,ylag);
for m = 1:ylag
    OutLagY(:,m) = DataY(locEnd+1-m:locForecastEnd-m);
    OutLagYdate(:,m) = DataYdateNum(locEnd+1-m:locForecastEnd-m);
end

% Construct lagged X data
if size(DataX,2)<2
    EstX = zeros(nobs,xlag);
else
    EstX = zeros(nobs, xlag, size(DataX,2)); %Forming a high dimension matrix
end

EstXdate = zeros(nobs,xlag);
for t = 1:nobs
    loc = find(DataXdateNum >= EstYdate(t)-tol, 1);
    %loc
    if isempty(loc)
        loc = length(DataXdateNum);
    end
    
    if loc-horizon > size(DataX,1)        
        nobs = t - 1;
        EstY = EstY(1:nobs,:);
        EstYdate = EstYdate(1:nobs,:);
        EstLagY = EstLagY(1:nobs,:);
        EstLagYdate = EstLagYdate(1:nobs,:);
        if size(EstX,3)>1
            EstX = EstX(1:nobs,:,:);
        else
            EstX = EstX(1:nobs,:);
        end
        EstXdate = EstXdate(1:nobs,:);
        maxDate = EstYdate(end);
        warning('MixFreqData:EstEndOutOfBound',...
            'Horizon is a large negative number. Observations are further truncated to %s',datestr(maxDate))
        break
    else
        if size(DataX,2)<2
            
            EstX(t,:) = DataX(loc-horizon:-1:loc-horizon-xlag+1);
        else
            for nb_col = 1:size(DataX,2)
                EstX(t,:,nb_col)= DataX(loc-horizon:-1:loc-horizon-xlag+1, nb_col);
            end
        end
        %fprintf("Embedding high dimension matrix: \n");
        EstX;                                                                % Printing function
        EstXdate(t,:) = DataXdateNum(loc-horizon:-1:loc-horizon-xlag+1);
    end    
end

if size(DataX, 2)<2 %Specifying the condition for high dimension matrix
    OutX = zeros(nforecast,xlag);
else
    OutX = zeros(nforecast, xlag, size(DataX,2));
end

OutXdate = zeros(nforecast,xlag);
for t = 1:nforecast
    loc = find(DataXdateNum >= OutYdate(t)-tol, 1);
    if isempty(loc)
        loc = length(DataXdateNum);
    end
    if loc-horizon > size(DataX,1)        
        nforecast = t - 1;
        OutY = OutY(1:nforecast,:);
        OutYdate = OutYdate(1:nforecast,:);
        OutLagY = OutLagY(1:nforecast,:);
        OutLagYdate = OutLagYdate(1:nforecast,:);
        OutX = OutX(1:nforecast,:);
        OutXdate = OutXdate(1:nforecast,:);
        break
    else
        if size(DataX,2)<2
            OutX(t,:) = DataX(loc-horizon:-1:loc-horizon-xlag+1);
        else
            for nb_col = 1:size(DataX,2)
                OutX(t,:,nb_col)= DataX(loc-horizon:-1:loc-horizon-xlag+1,nb_col); % Construct high-dimension forecasting matrix.
            end
        end
        OutXdate(t,:) = DataXdateNum(loc-horizon:-1:loc-horizon-xlag+1);
    end   
end

if dispFlag
    % Display mixed frequency data
    fprintf('Frequency of Data Y: %d %s\n',periodY,dateFormat{unitY});
    fprintf('Frequency of Data X: %d %s\n',periodX,dateFormat{unitX});
    fprintf('Start    Date: %s\n',datestr(estStart));
    fprintf('Terminal Date: %s\n',datestr(estEnd));
    fprintf('\n');
    
    % Display timeframe of mixed frequency regression
    fprintf('Mixed frequency regression time frame:\n')
    for m = [1,2,nobs]
        fprintf('Reg Y(%s) on ',datestr(EstYdate(m),2));
        if ylag == 1
            fprintf('Y(%s),',datestr(EstLagYdate(m,1),2));
        end
        if ylag == 2
            fprintf('Y(%s),Y(%s),',datestr(EstLagYdate(m,1),2), datestr(EstLagYdate(m,end),2));
        end
        if ylag > 2
            fprintf('Y(%s),...,Y(%s),',datestr(EstLagYdate(m,1),2), datestr(EstLagYdate(m,end),2));
        end
        if xlag == 1
            fprintf('X(%s)',datestr(EstXdate(m,1),2));
        end
        if xlag == 2
            fprintf('X(%s),X(%s)',datestr(EstXdate(m,1),2), datestr(EstXdate(m,end),2));
        end
        if xlag == 3
            fprintf('X(%s),X(%s),X(%s)',datestr(EstXdate(m,1),2),datestr(EstXdate(m,2),2),datestr(EstXdate(m,end),2));
        end
        if xlag > 3
            fprintf('X(%s),X(%s),...,X(%s)',datestr(EstXdate(m,1),2),datestr(EstXdate(m,2),2),datestr(EstXdate(m,end),2));
        end
        fprintf('\n');
        if m == 2
            fprintf('...\n');
        end
    end
end

Output = struct('EstY',EstY,'EstYdate',EstYdate,'EstX',EstX,'EstXdate',EstXdate,...
    'EstLagY',EstLagY,'EstLagYdate',EstLagYdate,...
    'OutY',OutY,'OutYdate',OutYdate,'OutX',OutX,'OutXdate',OutXdate,...
    'OutLagY',OutLagY,'OutLagYdate',OutLagYdate,'Xlag',xlag,'Ylag',ylag,...
    'MinDate',minDate,'MaxDate',maxDate);

end


%-------------------------------------------------------------------------
function [period,unit] = DateFreq(DateVec)

% DateFreq: Identify data frequency
%
% Input Arguments:
%
% DateVec: T-by-6 MATLAB vector format data: [year,month,day,hour,min,sec]
%
% Output Arguments:
%
% period: length of two consecutive dates
%
% unit: unit of length measure
%       o 1 = year
%       o 2 = month
%       o 3 = day
%       o 4 = hour
%       o 5 = minutes
%       o 6 = seconds
%
% Notes:
%
% Frequency   period   unit
% yearly         1      1  
% semiannual     6      2 
% quarterly      3      2 
% monthly        1      2 
% biweekly       14     3
% weekly         7      3
% daily          1      3
% hourly         1      4
% minutely       1      5
% secondly       1      6

DateDiff = diff(DateVec); % H: This finds out the time difference between two rows

% Check annual or lower frequency
modeUse = mode(DateDiff(:,1));
if modeUse >= 1 %If the first column (Year) distance is mostly larger than 1, 
    period = modeUse; % Then it should be a yearly data
    unit = 1;
    return
end

% Check monthly frequency, quarter = 3 months, semiannual = 6 months
modeUse = mode(DateDiff(:,2));
mask = modeUse < 0;
modeUse(mask) = modeUse(mask) + 12;
if modeUse >= 1
    period = modeUse;
    unit = 2;
    return
end

% Check daily frequency, week = 7 days, biweekly = 14 days
modeUse = mode(DateDiff(:,3));
mask = modeUse < 0;
modeUse(mask) = modeUse(mask) + 30;
if modeUse >= 1
    period = modeUse;
    unit = 3;
    return
end

% Check hourly frequency
modeUse = mode(DateDiff(:,4));
mask = modeUse < 0;
modeUse(mask) = modeUse(mask) + 24;
if modeUse >= 1
    period = modeUse;
    unit = 4;
    return
end

% Check minutely frequency
modeUse = mode(DateDiff(:,5));
mask = modeUse < 0;
modeUse(mask) = modeUse(mask) + 60;
if modeUse >= 1
    period = modeUse;
    unit = 5;
    return
end

% Check secondly frequency
elapse = etime(DateVec(2:end,:),DateVec(1:end-1,:));
period = mean(elapse);
unit = 6;

end

%-------------------------------------------------------------------------
function [dataMode,countMax] = mode(data)

% mode: mode (most frequent value) of a vector
%
% Input Arguments:
%
% data: a vector of data
%
% Output Arguments:
%
% dataMode: mode of the vector
%
% countMax: frequency count at the mode
%

nobs = numel(data);
data = sort(data);
count = 1;
countMax = 1;
dataMode = data(1);
for t = 2:nobs
    if data(t) == data(t-1)
        count = count + 1;        
    else        
        if count > countMax
            countMax = count;
            dataMode = data(t-1);
        end
        count = 1;
    end
end

if count > countMax
    countMax = count;
    dataMode = data(nobs);
end

end

%------------------------------------------------------------------------
function xlag = lagNum(xlag,period,unit)
% Not entirely sure about the usage of this function?
if isnumeric(xlag)
    validateattributes(xlag,{'numeric'},{'scalar','integer'},'','lag');
    return
end

multiplier = str2double(xlag(1:end-1));
if isnan(multiplier)
    error('The description "%s" cannot be recognized. The format should be 3m, 1q, etc.',xlag)
end
    
% Convert multiplier to daily frequency (business days)
ndaysPerYear = 260;
ndaysPerQuarter = 65;
ndaysPerMonth = 22;
nhoursPerDay = 8;
switch xlag(end)
    case 'y'
        multiplier = multiplier * ndaysPerYear;
    case 'q'
        multiplier = multiplier * ndaysPerQuarter;
    case 'm'
        multiplier = multiplier * ndaysPerMonth;
    case 'd'
        multiplier = multiplier * 1;
    case 'h'
        multiplier = multiplier / nhoursPerDay;
    case 's'
        multiplier = multiplier / (nhoursPerDay*60*60);
end

% Compute number of lags
switch unit
    case 1
        xlag = round(multiplier / (ndaysPerYear * period));
    case 2
        xlag = round(multiplier / (ndaysPerMonth * period));
    case 3
        xlag = round(multiplier / period);
    case 4
        xlag = round(multiplier / (period / nhoursPerDay));
    case 5
        xlag = round(multiplier / (period / nhoursPerDay / 60));
    case 6
        xlag = round(multiplier / (period / nhoursPerDay / 60 / 60));
end
        
end
