function [R_opt, t_opt, best_set] = robustSolver_KS(adj_, pts_3d, pts_3d_, noise)

n_ele=size(pts_3d, 1);

best_set=zeros(1,1);co=0;
for i=1:n_ele
    if sum(adj_(i,:))>0
        co=co+1;
        best_set(co)=i;
    end
end

if co<3
    R_opt=eye(3);
    t_opt=zeros(3,1);
    return
end

%% RANSAC procedures

max_itr=1e+3;
opt_size=0;
opt_set=[];
itr=0;
while 1
    
    itr=itr+1;
    
    this_set=best_set(randperm(co,3));
    
    [R_this,t_this]=minimal_PCR_KS(pts_3d(this_set,:),pts_3d_(this_set,:));

    re=sqrt(sum((R_this*pts_3d(best_set,:)'+t_this-pts_3d_(best_set,:)').^2));

    this_set=best_set(re<=1.5*3.5*noise);

    p_=mean(pts_3d(this_set,:),1)';
    q_=mean(pts_3d_(this_set,:),1)';
    
    H=(pts_3d(this_set,:)'-p_)*(pts_3d_(this_set,:)'-q_)';
    
    [U,~,V]=svd(H);
    
    R_opt=V*U';
    
    t_opt=q_-R_opt*p_;
    
    re=sqrt(sum((R_opt*pts_3d(best_set,:)'+t_opt-pts_3d_(best_set,:)').^2));
    
    this_set_=best_set(re<=3.5*noise);

    count=numel(this_set);
    
    if count>opt_size
        opt_set=this_set_;
        opt_size=count;
        max_itr=log(0.01)/log(1-(opt_size/co)^3);
    end
    
    if itr>=max_itr || count==co
        break
    end
    
end

p_=mean(pts_3d(opt_set,:),1)';
q_=mean(pts_3d_(opt_set,:),1)';

H=(pts_3d(opt_set,:)'-p_)*(pts_3d_(opt_set,:)'-q_)';

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-R_opt*p_;

end
