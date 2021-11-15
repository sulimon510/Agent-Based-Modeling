
function [theta_out] = modified_VM_multi(lbox,vs,r,eta,ncell,TMAX,wmat)

%This function computes the modified Vicsek Models of Models A and B in
%[1] Sattari et. al, Sci. Adv (2021)

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

%The input wmat determines the interaction strengths and model type.
%wmat is an ncell x ncell matrix whose ith row and jth column determines
%the interaction between particle i and j. For example, if ncell=3, 
%wmat = [1 0 0; 0 1 0; 0 0 1] means that each particle only depends
%on its own past and does not interact with others. In model B where the
%leaders and followers mutually interact but followers do not interact with
%each other, and the leader to follower weight is set to 4, the wmat for
%ncell=4 (1 leader, 3 followers) is [1 1 1 1; 4 1 0 0; 4 0 1 0; 4 0 0 1]. 


%Outputs:
%theta_out - the theta values of the leader (first row) and follower
%(second row)

%If you have any question or comment about this code please 
%email sulimon.sattari@gmail.com


vx=zeros(ncell,TMAX);
vy=zeros(ncell,TMAX);

xb=rand(ncell,1).*lbox;
yb=rand(ncell,1).*lbox;

ang=pi.*2.*rand(ncell,1);
nsteps=1;
vxb=vs.*cos(ang);
vyb=vs.*sin(ang);
vx(:,nsteps)=vxb;
vy(:,nsteps)=vyb;


theta_out(:,1)=ang;

theta=ang;



for nsteps=2:TMAX
    x_prev=xb;
    y_prev=yb;
    xb=xb+vxb;
    yb=yb+vyb;
    
    for cell1=ncell:-1:1
        
        
        [xb,yb]=boundary_conditions(xb,yb,cell1,lbox);
        
        
        
        [nearang]=calc_nearby(x_prev,y_prev,cell1,lbox,r);
        
        
        
        
        [vxtemp,vytemp]=calc_nearby_average(wmat,nearang,cell1,vx(nearang,nsteps-1),vy(nearang,nsteps-1));
        
        
        vtemp=[vxtemp,vytemp];
        curr_theta=calc_curr_theta(vxtemp,vytemp,vtemp);
        
        
        theta(cell1)=curr_theta;
        
        
        
    end
    
    ang=theta+eta.*(rand(ncell,1)-0.5);
    
    
    ang(ang>2*pi)=ang(ang>2*pi)-2*pi;
    ang(ang<0) = ang(ang<0) + 2*pi;
    
    vxb=vs.*cos(ang);
    vyb=vs.*sin(ang);
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

