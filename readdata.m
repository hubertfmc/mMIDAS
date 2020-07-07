function [dt, sp,ip, ft, id, nl_var] = readdata(c_id,freq)
%==========================================================================
% Function: Parse excel worksheets in corresponding directory and return
% corresponding date and date matrix.
% Input: Country name
% Output: pt: price vector (presumably weekly / daily observations)
%         dt: dividend vector (Quarterly/Semi-annual/annual data)
%         ft: financial characteristics (Quarterly/semi-annual/annual data)
%         id: Name of sample
% By FMC, 20 Aug 2019
%==========================================================================

%
nan_tol = 0.95 ;
addpath (pwd + "\Data\"+freq);
wb= c_id+".xlsx";
if exist(wb,'file')~=2
    warning("Cannot find data file in the search directory! \n");
    return
else
    fprintf('Parsing from %s \n',wb);
end
% Name each Variables
[~,n_ws]=xlsfinfo(wb);
nl_var = [];
for i = 1:length(n_ws)
    if n_ws{i}~=string('REQUEST_TABLE')
        d_var = strcat('d_',string(n_ws{i}));
        t_var = strcat('t_', string(n_ws{i}));
        if n_ws{i}== string('IP')|n_ws{i}== string('IDY')
            [index.(d_var), index.(t_var)]=xlsread(wb,n_ws{i});
            index.(t_var) = index.(t_var)(2:end,1); 
        elseif n_ws{i} == string('ID')
            [fvar.(d_var), fvar.(t_var)]=xlsread(wb,n_ws{i});
            fvar.(t_var) = fvar.(t_var)(2:end,1:2);
            nl_var = [nl_var; [d_var, t_var]];
        else
            [fvar.(d_var), fvar.(t_var)]=xlsread(wb,n_ws{i});
            fvar.(t_var) = fvar.(t_var)(2:end,1);
            nl_var = [nl_var; [d_var, t_var]];
        end
    end
end
%% Verify whether data have the same frequency
%{
size_list = [];
for i = 1:length(nl_var)
    size_list = [size_list;[[size(variable.(nl_var(i,1)),1),size(variable.(nl_var(i,1)),2)], ...
                            [size(variable.(nl_var(i,2)),1),size(variable.(nl_var(i,2)),2)]]];
end
%}
%% Clearing NAN data

for i = 1:length(nl_var)
    if nl_var(i,1)~= string('d_ID')
        var_d = fvar.(nl_var(i,1));
            if nl_var(i,1)== string('d_D_q')
                nnan_l = sum(~isnan(var_d))== size(var_d,1);
            else
                nnan_l = sum(~isnan(var_d)) > nan_tol*size(var_d,1);
            end
            for j = 1:length(nl_var)
                if nl_var(j,1)~= string('d_ID')
                    fvar.(nl_var(j,1))= fvar.(nl_var(j,1))(:,nnan_l);
                else
                    fvar.(nl_var(j,1))= fvar.(nl_var(j,1))(nnan_l',:);
                    fvar.(nl_var(j,2))= fvar.(nl_var(j,2))(nnan_l',:);
                end
            end
    end        
end
%%
%{
[d_d,t_dt] = xlsread(wb, 'D');
t_dt = t_dt(2:end,1);


[p_t,t_pt] = xlsread(wb, 'SP');
t_pt = t_pt(6:end,1);

[id,i_id] = xlsread(wb, 'ID');
id = id(2:end,1:2);

[d_re,t_re] = xlsread(wb, 'RE');
t_re = t_re(2:end,1);

[d_te, t_te] = xlsread(wb, 'TE');
t_te = t_te(2:end,1);

re_r = d_re./d_te;

[d_roa, t_roa] = xlsread(wb, 'ROA');
t_roa = t_roa(2:end,1);

[d_pbr, t_pbr] = xlsread(wb, 'PBr');
t_pbr = t_pbr(2:end,1);

[d_eps, t_eps] = xlsread(wb, 'EPS');
t_eps = t_eps(2:end,1);

[d_eps, t_eps] = xlsread(wb, 'EPS');
t_eps = t_eps(2:end,1);





d_ft = [];
if size(d_roa,1)==size(re_r,1)&& ...
        size(d_pbr,1)==size(re_r,1)&& size(d_eps,1)==size(re_r,1)
        d_ft(:,:,1) = re_r;
        d_ft(:,:,2) = d_roa;
        d_ft(:,:,3) = d_pbr;
        d_ft(:,:,4) = d_eps;
else
    warning('Data have different length!');
    return
end

if size(t_roa,1)== size(d_re,1) && size(t_pbr,1)==size(d_re,1) ...
        size(t_eps,1)== size(d_re,1) && size(t_pbr,1)==size(d_re,1)
    t_ft = t_re;
else
    warning('Size of time vectors are different')
    return
end

%%
%Check if date vector matches among factors.

% Check if factors with the same number of samples
if size(p_t,2)~=size(id,1)|size(d_ft,2)~=size(id,1) 
    warning('Factors with different samples!');
    return
end
%%
% Data cleasing
% Excluding companies with NAN data.
col_vec = sum(~isnan(d_d))== size(d_d,1);                               % Extract column with NAN entries
fprintf('Excluding...')
id(logical(1-col_vec)',:)
d_d = d_d(:, col_vec);                                                  % Retain useful data.
p_t = p_t(:, col_vec);
d_ft = d_ft(:, col_vec,:);
id = id(col_vec',:);
%}
%%
for i = 1:length(nl_var)
    if nl_var(i,1)~= string('d_ID')
        var_d = fvar.(nl_var(i,1));
            nz_list = [];
            for nb_c = 1:size(var_d,2)
                if nnz(var_d(:,nb_c))> nan_tol*size(var_d,1)                                  %   Retain sample with more than half observations being non-zero dividend, otherwise 
                    nz_list=[nz_list,nb_c];
                end
            end
            for j = 1:length(nl_var)
                if nl_var(j,1)~= string('d_ID')
                    fvar.(nl_var(j,1))= fvar.(nl_var(j,1))(:,nz_list);
                else
                    fvar.(nl_var(j,1))= fvar.(nl_var(j,1))(nz_list',:);
                    fvar.(nl_var(j,2))= fvar.(nl_var(j,2))(nz_list',:);
                end
            end
    end        
end
%% Erasing NAN data from price

ip_nan = logical(sum(~isnan(index.d_IP),2)==size(index.d_IP,2));
index.d_IP = index.d_IP (ip_nan,:);
index.t_IP = index.t_IP (ip_nan,:);
fvar.d_SP = fvar.d_SP(ip_nan,:);
fvar.t_SP = fvar.t_SP(ip_nan,:);

sp_nan = logical(sum(~isnan(fvar.d_SP),2)== size(fvar.d_SP,2));
index.d_IP = index.d_IP (sp_nan,:);
index.t_IP = index.t_IP (sp_nan,:);
fvar.d_SP = fvar.d_SP(sp_nan,:);
fvar.t_SP = fvar.t_SP(sp_nan,:);

%%

% Compute log return of the price data
fvar.d_SP = (fvar.d_SP(2:end,:)- fvar.d_SP(1:end-1,:))./fvar.d_SP(1:end-1,:);                               % Compute the log return.
fvar.t_SP = fvar.t_SP(2:end);

index.d_IP = (index.d_IP(2:end,:)-index.d_IP(1:end-1,:))./index.d_IP(1:end-1,:);
index.t_IP = index.t_IP(2:end);

%%
% Plot return graph and dividend payment
%{
figure;
plot(datetime(fvar.t_SP,'InputFormat','dd/MM/yyyy'),[fvar.d_SP, index.d_IP]);
ts = "Time series plot of log return -"+ string(c_id); 
title(ts);
inst_list = [fvar.t_ID(:,1);string(c_id)];
lgd = legend(inst_list);
lgd.Location = 'south';
lgd.Orientation='horizontal';
lgd.NumColumns = 4;
ylabel('Log return (%)');


figure
plot(datetime(fvar.t_D_q,'InputFormat','dd/MM/yyyy'), fvar.d_D_q);
ts2 = "Time series plot of dividend payout -"+c_id;
title(ts2);
legend(fvar.t_ID(:,1));
ylabel('Dividend per share');
%}
%%
% Combine result
ip.data=index.d_IP;
ip.time=index.t_IP;

for i = 1:length(nl_var)
    if nl_var(i,1) == string('d_ID')
        id.sym = fvar.t_ID(:,1);
        id.fname = fvar.t_ID(:,2);
        id.dfreq = fvar.d_ID(:,1);
        id.acctfreq = fvar.d_ID(:,2);
        id.industry = fvar.d_ID(:,3);
    elseif nl_var(i,1) == string('d_SP')
        sp.time = fvar.t_SP;
        sp.data = fvar.d_SP;
    elseif nl_var(i,1) == string('d_D_q')|nl_var(i,1) == string('d_D_yr')
        dt.annual = fvar.t_D_yr;
        dt.d_yr = fvar.d_D_yr;
        dt.quarter = fvar.t_D_q;
        dt.d_q = fvar.d_D_q;
    else
        ft.time=fvar.t_RE;
        ft.re= fvar.d_RE;
        ft.te=fvar.d_TE;
        ft.roa = fvar.d_ROA;
        ft.d_roa = fvar.d_ROA(2:end,:) - fvar.d_ROA(1:end-1,:);
        ft.cash = fvar.d_Cash;
        ft.c_ratio = fvar.d_Cash./fvar.d_TA;
        ft.pbr=fvar.d_PBr;
        ft.rer=fvar.d_RE./fvar.d_TE;
        ft.ta = fvar.d_TA;
        ft.log_ta = log(fvar.d_TA);
    end
end
end
