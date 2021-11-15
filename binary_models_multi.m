function symbols= binary_models_multi(nbird,p_leader,TMAX,type)
switch(type)
    case 'A'
        symbols = binary_model_A_multi(nbird,p_leader,TMAX);
    case 'B'
        symbols = binary_model_B_multi(nbird,p_leader,TMAX);
end
end

function symbols = binary_model_A_multi(nbird,p_leader,TMAX)

x=zeros(1,TMAX);
followers=zeros(nbird-1,TMAX);
status=rand(1,TMAX);
status_follower=rand(nbird-1,TMAX);
for ind=2:TMAX
    if status(ind)<p_leader
        x(ind)=x(ind-1);
    else
        x(ind)=randi(2)-1;
    end
    
    for f_ind=1:nbird-1
        
        curr_status=status_follower(f_ind,ind);
        if curr_status<p_leader
            followers(f_ind,ind)=x(ind-1);
            
        else
            follower_num=randi(nbird-1);
            
            followers(f_ind,ind)=followers(follower_num,ind-1);
        end
    end
end

symbols(:,1)=uint8(x);
symbols(:,2)=uint8(followers(1,:));


end

function symbols = binary_model_B_multi(nbird,p_leader,TMAX)


x=zeros(1,TMAX);
followers=zeros(nbird-1,TMAX);
status=rand(1,TMAX);
status_follower=rand(nbird-1,TMAX);
for ind=2:TMAX
    if status(ind)<p_leader
        x(ind)=x(ind-1);
    else
        x(ind)=randi(2)-1;
    end
    
    
    for f_ind=1:nbird-1
        
        curr_status=status_follower(f_ind,ind);
        if curr_status<p_leader
            followers(f_ind,ind)=x(ind-1);
            
        else
            
            follower_num=f_ind;
            followers(f_ind,ind)=followers(follower_num,ind-1);
        end
    end
end


symbols(:,1)=uint8(x);
symbols(:,2)=uint8(followers(1,:));


end
