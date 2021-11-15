
function [theta_out] = modified_VM(lbox,vs,r,eta,TMAX,leader_weight,randomize_F,randomize_L,randomize_L_every_two)

%This function computes the modified Vicsek Models of types A, B, C, D, C' and D'
%in [1] Sattari et. al, Sci. Adv (2021)

%This function relies on the external MEX function mex_calc_nearby. 
%To compile the mex, run 
%mex mex_calc_nearby.cpp
%in the MATLAB command window in the same folder as this .m file and 
%the mex_calc_nearby.cpp file. A gcc compiler is required.

%Inputs:
%lbox - the square box side length in arb. units. lbox=10.
%vs - the speed in arb. units of length over time steps. vs= 0.3 in [1].
%r - the interaction radius in arb. units. r=3 in [1].
%eta - the noise value  between 0 and 2pi. eta varies in [1].
%TMAX - the time step length of the simulation. TMAX=2e6 in [1].
%leader_weight - the weight of leaders. varies from 0 to 10 in [1].

%The inputs randomize_F, randomize_L, and randomize_L_every_two determine
%the types as described in [1]
%(randomize_F,randomize_L,randomize_L_every_two) is set as
%Type A:  (1,1,0)
%Type B:  (0,1,0)
%Type C:  (1,0,0)
%Type D:  (0,0,0)
%Type C': (1,0,1)
%Type D': (0,0,1)

%Outputs:
%theta_out - the theta values of the leader (first row) and follower
%(second row)

%If you have any question or comment about this code please 
%email sulimon.sattari@gmail.com

%Initialize wmatrix: Leader influences follower with weight leader_weight,
%follower does not influence leader in this model.
wmat=[1 0; leader_weight 1];

%Initialize number of agents
ncell=2;

%Initialize arrays
vx=zeros(ncell,TMAX);
vy=zeros(ncell,TMAX);
xb=rand(ncell,1).*lbox;
yb=rand(ncell,1).*lbox;
ang=pi.*2.*rand(ncell,1);

%current step number is 1
nsteps=1;
%initial conditions
vxb=vs.*cos(ang);
vyb=vs.*sin(ang);
vx(:,nsteps)=vxb;
vy(:,nsteps)=vyb;
theta_out(:,nsteps)=ang;
theta=ang;

%cell 1 is the leader, cell 2 is the follower
is_follower(1)=0;
is_follower(2)=1;


for nsteps=2:TMAX
%We use the previous positions calculate the next velocity
    x_prev=xb;
    y_prev=yb;
%step the positions forward
    xb=xb+vxb;
    yb=yb+vyb;
    
    for cell1=ncell:-1:1
        %Map cells back onto the domain using periodic boundary coditions
        [xb,yb]=boundary_conditions(xb,yb,cell1,lbox);
        
        %Calculate nearby cells
        [nearang]=calc_nearby(x_prev,y_prev,cell1,lbox,r);
        
        
        %For each respective type, compute interaction
        %calc_nearby_average_A randomizes that agent before computing the
        %interaction, calc_nearby_average does not
        
        if randomize_L && randomize_F
            
            [vxtemp,vytemp]=calc_nearby_average_A(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1),vs);
        end
        
        
        if randomize_F && ~randomize_L
            if is_follower(cell1)
                
                [vxtemp,vytemp]=calc_nearby_average_A(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1),vs);
            else
                if randomize_L_every_two
                    if mod(nsteps,2)
                        [vxtemp,vytemp]=calc_nearby_average(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1));
                    else
                        [vxtemp,vytemp]=calc_nearby_average_A(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1),vs);
                    end
                else
                    [vxtemp,vytemp]=calc_nearby_average(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1));
                end
            end
        end
        
        
        if randomize_L && ~randomize_F
            if is_follower(cell1)
                
                [vxtemp,vytemp]=calc_nearby_average(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1));
                
            else
                
                [vxtemp,vytemp]=calc_nearby_average_A(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1),vs);
            end
            
        end
        
        
        if ~randomize_L && ~randomize_F
            if is_follower(cell1)
                [vxtemp,vytemp]=calc_nearby_average(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1));
            else
                
                if randomize_L_every_two
                    if mod(nsteps,2)
                        [vxtemp,vytemp]=calc_nearby_average(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1));
                    else
                        [vxtemp,vytemp]=calc_nearby_average_A(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1),vs);
                    end
                else
                    [vxtemp,vytemp]=calc_nearby_average(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1));
                end
                
            end
            
        end
        
        
        vtemp=[vxtemp,vytemp];
        %Compute theta from vx vy
        curr_theta=calc_curr_theta(vxtemp,vytemp,vtemp);
        
        theta(cell1)=curr_theta;
        
        
        
    end
    
    
    %add noise
    
    ang=theta+eta.*(rand(ncell,1)-0.5);
    
    
    
    
    
    %compute the angle modulo 2pi
    ang(ang>2*pi)=ang(ang>2*pi)-2*pi;
    ang(ang<0) = ang(ang<0) + 2*pi;
    
    %compute the speed vx vy from the angle
    
    vxb=vs.*cos(ang);
    vyb=vs.*sin(ang);
    
    %save the speed and angle
    vx(:,nsteps)=vxb;
    vy(:,nsteps)=vyb;
    theta_out(:,nsteps)=ang;
end




end

function [xb,yb]=boundary_conditions(xb,yb,cell1,lbox)

if(xb(cell1)<0);xb(cell1)=xb(cell1)+lbox; end
if (yb(cell1)<0);yb(cell1)=yb(cell1)+lbox;end
if (xb(cell1)>lbox);xb(cell1)=xb(cell1)-lbox;end
if (yb(cell1)>lbox);yb(cell1)=yb(cell1)-lbox;end

end

function [nearang]=calc_nearby(all_x,all_y,cell1,lbox,r)
[nearang]=mex_calc_nearby(all_x,all_y,cell1,lbox,r);

nearang=find(nearang);

end

function [vxtemp,vytemp]=calc_nearby_average_A(wmat,nearang,cell1,curr_vx_nearby,curr_vy_nearby,vs)
curr_weights=wmat(cell1,:);
curr_cell=find(nearang==cell1);
weights0=curr_weights(nearang)';
weights=weights0;
[random_vx,random_vy]=erase_memory(vs);
curr_vx_nearby(curr_cell)=random_vx;
curr_vy_nearby(curr_cell)=random_vy;



    vxtemp=(sum(weights.*curr_vx_nearby,1))./(sum(weights,1));
    vytemp=(sum(weights.*curr_vy_nearby,1))./(sum(weights,1));
    


end

function [vxtemp,vytemp]=calc_nearby_average(wmat,nearang,cell1,curr_vx_nearby,curr_vy_nearby)
curr_weights=wmat(cell1,:);
weights0=curr_weights(nearang)';
weights=weights0;


    vxtemp=(sum(weights.*curr_vx_nearby,1))./(sum(weights,1));
    vytemp=(sum(weights.*curr_vy_nearby,1))./(sum(weights,1));
    

end

function curr_theta=calc_curr_theta(vxtemp,vytemp,vtemp)
if vytemp>=0
    curr_theta=acos(vxtemp/norm(vtemp));
else
    if vytemp<0
        curr_theta=2*pi-acos(vxtemp/norm(vtemp));
        
    end
end
end

function [random_vx,random_vy]=erase_memory(vs)
rand_theta=2*pi*rand();
random_vx=vs*cos(rand_theta);
random_vy=vs*sin(rand_theta);
end
