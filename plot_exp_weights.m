
function w_mat =plot_exp_weights(dayLag,lB, uB, delta)


rng = lB:delta:uB;
t1_vec = rng;
t2_vec = rng;
w_mat=[];
    for c_t1 = 1: length(rng)
        for c_t2 = 1: length(rng)
                   w_vec = cal_exp_weight(dayLag,t1_vec(c_t1),t2_vec(c_t2));
                   w_mat = [w_mat; w_vec];
        end
    end

[a,b]= find(w_mat(:,3:end)~=0);
w_mat = w_mat(a,:);
w_mat= w_mat(sum(~isinf(w_mat(:,3:end)),2)==size(w_mat(:,3:end),2),:);

plot(1:size(w_vec(1,3:end),2),w_vec(:,3:end));
%label(string(:,1:2));

end

function w_vec = cal_exp_weight(L,t1,t2)
        iii=1:L;
        iii=repelem(iii,2,1);
        iii(2,:)= iii(2,:).^2;
        theta=[t1,t2];
        w_vec(1,1:2)= theta;
        result = theta*iii;
        w=exp(result)./exp(sum(result));
        w_vec (1,3:3+size(w,2)-1) = w;
end