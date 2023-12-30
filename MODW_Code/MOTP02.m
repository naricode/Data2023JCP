function z=MOTP02(x)
    global NFE
    NFE=NFE+1;
    
    n=numel(x);
    
    z1=1-exp(-sum((x(1:end)-1/sqrt(n)).^2));
    
    z2=1-exp(-sum((x(2:end)+1/sqrt(n)).^2));
   
    z=[z1 z2]';

end