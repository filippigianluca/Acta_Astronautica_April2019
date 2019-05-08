function [upperFunctionValue upperEqualityConstrVals upperInequalityConstrVals]=upperLevelSMD(upperLevelMember, lowerLevelMember, testProblemName)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function call here
    fhandle = str2func(testProblemName);

    equalityConstrVals = [];
    inequalityConstrVals = [];
    
    [upperFunctionValue upperEqualityConstrVals upperInequalityConstrVals] = fhandle(upperLevelMember, lowerLevelMember);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd1(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);
    
    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2).^2) + sum((xu2 - tan(xl2)).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd2(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);   
    
    functionValue = sum((xu1).^2) ...
                    - sum((xl1).^2) ...
                    + sum((xu2).^2) - sum((xu2 - log(xl2)).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd3(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    c = ones(size(xu2));    
    
    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - c).^2) + sum((xu2.^2 - tan(xl2)).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd4(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);
    
    functionValue = sum((xu1).^2) ...
                - sum((xl1).^2) ...
                + sum((xu2).^2) - sum((abs(xu2) - log(1+xl2)).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd5(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r); 
    
    term2 = 0;
    for i=1:q-1
        term2 = term2 + (xl1(i+1) - xl1(i).^2).^2 + (xl1(i) - 1).^2;
    end
    
    %Same as smd5
    functionValue = sum((xu1).^2) ...
                - term2 ...
                + sum((xu2).^2) - sum((abs(xu2) - xl2.^2).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd6(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = floor((length(xl) - r)/2 - eps);
    s = ceil((length(xl) - r)/2 + eps);

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q+s);
    xl2 = xl(q+s+1:q+s+r);   
    
    functionValue = sum((xu1).^2) ...
                    - sum(xl1(1:q).^2) + sum(xl1(q+1:q+s).^2) ...
                    + sum((xu2).^2) - sum((xu2 - xl2).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd7(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);  
    
    m = [1:p];
    functionValue = 1+1/400*sum((xu1).^2) - prod(cos(xu1./sqrt(m))) ...
                    - sum((xl1).^2) ...
                    + sum((xu2).^2) - sum((xu2 - log(xl2)).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd8(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);  
    
    term2 = 0;
    for i=1:q-1
        term2 = term2 + (xl1(i+1) - xl1(i).^2).^2 + (xl1(i) - 1).^2;
    end
    
    functionValue = 20+exp(1)-20*exp(-0.2*sqrt(1/p*sum((xu1).^2))) - exp(1/p*sum(cos(2*pi*xu1)))  ...
                - term2 ...
                + sum((xu2).^2) - sum((xu2 - xl2.^3).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd9(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);  
    
    functionValue = sum((xu1).^2) ...
                    - sum((xl1).^2) ...
                    + sum((xu2).^2) - sum((xu2 - log(1+xl2)).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals(1) = sum(xu1.^2)+sum(xu2.^2) - floor(sum(xu1.^2)+sum(xu2.^2)+0.5);
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd10(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    a = 2*ones(size(xu1));
    c = 2*ones(size(xu2));    

    
    functionValue = sum((xu1 - a).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - c).^2) - sum((xu2 - tan(xl2)).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    for i=1:p
        inequalityConstrVals(i) = xu1(i) + xu1(i).^3 - sum(xu1.^3) - sum(xu2.^3);
    end

    for i=1:r
        inequalityConstrVals(p+i) = xu2(i) + xu2(i).^3 - sum(xu2.^3) - sum(xu1.^3);
    end

    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd11(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);
    
    functionValue = sum((xu1).^2) ...
                    - sum((xl1).^2) ...
                    + sum((xu2).^2) - sum((xu2 - log(xl2)).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = xu2 - 1/sqrt(r) - log(xl2);
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd12(xu, xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;

    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    a = 2*ones(size(xu1));
    c = 2*ones(size(xu2));    

    functionValue = sum((xu1 - a).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - c).^2) + sum(tan(abs(xl2))) - sum((xu2 - tan(xl2)).^2);

    functionValue = -functionValue;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    for i=1:p
        inequalityConstrVals(i) = xu1(i) + xu1(i).^3 - sum(xu1.^3) - sum(xu2.^3);
    end

    for i=1:r
        inequalityConstrVals(p+i) = xu2(i) + xu2(i).^3 - sum(xu2.^3) - sum(xu1.^3);
    end

    inequalityConstrVals(p+r+1:p+2*r) = xu2 - tan(xl2);
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
