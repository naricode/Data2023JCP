function z=MOTP01(x)
    NFE=NFE+1;
    
    n=numel(x);

    z1=x(1);
    
    g=1+9/(n-1)*sum(x(2:end));
    
    h=1-sqrt(z1/g);
    
    z2=g*h;
    
    z=[z1 z2]';

end