clear all
nruns = 1;
alpha_idx = 0;
alpha_bracket = 2.0:0.01:3.0;
for alpha_temp=alpha_bracket
    tic
    alpha_idx = alpha_idx + 1;
    for run = 1:nruns
        N  = 100;
        x0 = randn(1,length(0:2:N));
        LB = -Inf(1,length(0:2:N));
        UB = Inf(1,length(0:2:N));
        f = @(params)loss(T_g_handle(g_handle(params,N),alpha_temp),g_handle(params,N),alpha_temp);
        [paramsEst(alpha_idx,run,:),fvalue(alpha_idx,run)] = fminsearchbnd(f, x0, LB, UB);
    end
    [fvalue_min(alpha_idx),idx] = min(fvalue(alpha_idx,:));
    params_Est_min(alpha_idx,:) = squeeze(paramsEst(alpha_idx,idx,:));
    toc
end

[function_value,f_idx] = min(fvalue_min);
parameter = params_Est_min(f_idx,:);
alpha_est = alpha_bracket(f_idx);


function alpha_true = alpha_true(g)
    f = @(x) g(x-0.5)-x-g(0.5);
    if f(0.5)*f(1)<0
        temp = fzero(f,[0.5,1]);
        alpha_true = 1/(2*temp-1);
    else
        alpha_true = NaN;
    end
end
function MSE = loss(T_g,g,alpha_temp)
        T_g_handle = @(y) -alpha_temp.*g(g(-y./alpha_temp));

    SE = 0;
    ybracket = -0.5:0.01:0.5;
    for y = ybracket
        SE = (T_g(y)-g(y)).^2+SE;
    end
    MSE = SE./length(ybracket);
    alpha = alpha_true(g);
    if ~isnan(alpha)
        MSE = Inf;
    elseif abs(alpha -alpha_temp)>1e-2
        MSE = Inf;
    end
end
function T_g_handle = T_g_handle(g,alpha_temp)
end
function g_handle = g_handle(params,N)
    order = 0:2:N;
    g_handle = @(y) sum(params.*y.^order);
end