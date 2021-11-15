function symbols=binary_models(TMAX,c,type)

switch(type)
    case('A')
        symbols= binary_A(TMAX,c);
    case('B')
        
        symbols= binary_B(TMAX,c)
    case('C')
        
        symbols= binary_C(TMAX,c)
    case('D')
        
        symbols= binary_D(TMAX,c)
    case('CPrime')
        
        symbols= binary_CPrime(TMAX,c)
    case('DPrime')
        
        symbols= binary_DPrime(TMAX,c)
end

end
%Type A
function symbols = binary_A(TMAX,c)
status=rand(1,TMAX);
for ind=1:TMAX
    if status(ind)>c
        x(ind)=0;
    else
        x(ind)=1;
    end
end

y(1)=0;y(2:TMAX)=x(1:end-1);
symbols(:,1)=uint8(x);
symbols(:,2)=uint8(y);

end


%Type B
function symbols = binary_B(TMAX,c)

x=randi(2,[1 TMAX])-1;
status=rand(1,TMAX);
y(1)=0;
for ind=2:TMAX
    if status(ind)>c
        y(ind)=x(ind-1);
    else
        y(ind)=y(ind-1);
    end
end
symbols(:,1)=uint8(x);
symbols(:,2)=uint8(y);

end


%Type C
function symbols = binary_C(TMAX,c)

x(1)=0;
status=rand(1,TMAX);
for ind=2:TMAX
    if status(ind)>c
        x(ind)=x(ind-1);
    else
        x(ind)=randi(2)-1;
    end
end

y(1)=0;
y(2:TMAX)=x(1:end-1);
symbols(:,1)=uint8(x);
symbols(:,2)=uint8(y);

end

%Type D
function symbols = binary_D(TMAX,c)

x(1)=0;
status=rand(1,TMAX);
for ind=2:TMAX
    if status(ind)>c
        x(ind)=x(ind-1);
    else
        x(ind)=randi(2)-1;
    end
end
status=rand(1,TMAX);
y(1)=0;
for ind=2:TMAX
    if status(ind)>c
        y(ind)=x(ind-1);
    else
        y(ind)=y(ind-1);
    end
end
symbols(:,1)=uint8(x);
symbols(:,2)=uint8(y);

end
%Type C'
function symbols = binary_CPrime(TMAX,c)

x(1)=0;
status=rand(1,TMAX);
for ind=2:TMAX
    if mod(ind,2)
        if status(ind)>c
            x(ind)=x(ind-1);
        else
            x(ind)=randi(2)-1;
        end
    else
        x(ind)=randi(2)-1;
    end
end

y(1)=0;
y(2:TMAX)=x(1:end-1);
symbols(:,1)=uint8(x);
symbols(:,2)=uint8(y);



end

%Type D'
function symbols = binary_DPrime(TMAX,c)

x(1)=0;
status=rand(1,TMAX);

for ind=2:TMAX
    
    if mod(ind,2)
        if status(ind)>c
            x(ind)=x(ind-1);
        else
            x(ind)=randi(2)-1;
        end
    else
        x(ind)=randi(2)-1;
    end
end
status=rand(1,TMAX);
y(1)=0;
for ind=2:TMAX
    if status(ind)>c
        y(ind)=x(ind-1);
    else
        y(ind)=y(ind-1);
    end
end
symbols(:,1)=uint8(x);
symbols(:,2)=uint8(y);

end

