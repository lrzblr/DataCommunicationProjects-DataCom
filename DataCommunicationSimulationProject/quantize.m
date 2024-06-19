function [y, x2, errorcuantizacion] = quantize(x, nivel)

    % Levels of uniform quantization 
    dif=(max(x)-min(x))/(nivel-1);
    val=[min(x):dif:max(x)];
    % y = vector to plot levels of quantization
    y=kron(val',ones(1,size(x,1)));
    %No amplification
    xp=x;
    % Value decision: it decides where the value corresponds to one level
    % according minimum distance
    ref=repmat(xp',nivel,1); 
    x1=abs(y-ref);   
    [~, x2]=min(x1);

    %------------ QUANTIZATION ERROR
    %Signal normalization 
    m=2/(nivel-1);
    b=(-nivel-1)/(nivel-1);
    xerror=m.*x2+b;

    for i=1:size(xp,1)
        if xp(i)==0
            xp(i)=0.001;
        end    
        er(i)=abs((xp(i)-xerror(i))/xp(i));
    end

    errorcuantizacion=mean(er);
end