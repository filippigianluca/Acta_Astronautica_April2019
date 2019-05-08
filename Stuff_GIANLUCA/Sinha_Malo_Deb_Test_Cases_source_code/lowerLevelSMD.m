function [lowerFunctionValue lowerEqualityConstrVals lowerInequalityConstrVals]=lowerLevelSMD(upperLevelMember, lowerLevelMember, testProblemName)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function call here
    fhandle = str2func(testProblemName);

    equalityConstrVals = [];
    inequalityConstrVals = [];
    
    [lowerFunctionValue lowerEqualityConstrVals lowerInequalityConstrVals] = fhandle(upperLevelMember, lowerLevelMember);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd1(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - tan(xl2)).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd2(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - log(xl2)).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd3(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + q + sum(xl1.^2 - cos(2*pi*xl1)) ...
                    + sum((xu2.^2 - tan(xl2)).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd4(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                        + q + sum(xl1.^2 - cos(2*pi*xl1)) ...
                        + sum((abs(xu2) - log (1+xl2)).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd5(xu,xl)

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
    
    functionValue = sum((xu1).^2) ...
                        + term2 ...
                        + sum((abs(xu2) - xl2.^2).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [functionValue equalityConstrVals inequalityConstrVals] = smd6(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = floor((length(xl) - r)/2 - eps);
    s = ceil((length(xl) - r)/2 + eps);
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q+s);
    xl2 = xl(q+s+1:q+s+r);

    term2 = sum(xl1(1:q).^2);
    for i=q+1:2:q+s-1
        term2 = term2 + (xl1(i+1) - xl1(i)).^2;
    end
    
    functionValue = sum((xu1).^2) ...
                    + term2 ...
                    + sum((xu2 - xl2).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd7(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^3) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - log(xl2)).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd8(xu,xl)

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
    
    functionValue = sum(abs(xu1)) ...
                        + term2 ...
                        + sum((xu2 - xl2.^3).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals = [];
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd9(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - log(1+xl2)).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals(1) = sum(xl1.^2)+sum(xl2.^2) - floor(sum(xl1.^2)+sum(xl2.^2)+0.5);
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd10(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    b = 2*ones(size(xl1));

    functionValue = sum((xu1).^2) ...
                    + sum((xl1 - b).^2) ...
                    + sum((xu2 - tan(xl2)).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    for i=1:q
        inequalityConstrVals(i) = xl1(i) + xl1(i).^3 - sum(xl1.^3);
    end
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd11(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - log(xl2)).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals(1) = sum((xu2 - log(xl2)).^2) - 1;
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [functionValue equalityConstrVals inequalityConstrVals] = smd12(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    b = 2*ones(size(xl1));

    functionValue = sum((xu1).^2) ...
                    + sum((xl1 - b).^2) ...
                    + sum((xu2 - tan(xl2)).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    for i=1:q
        inequalityConstrVals(i) = xl1(i) + xl1(i).^3 - sum(xl1.^3);
    end
    inequalityConstrVals(q+1) = sum((xu2 - tan(xl2)).^2) - 1;
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
