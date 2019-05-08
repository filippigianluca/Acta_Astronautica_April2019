Readme file for SMD test problems.



SMD test-suite
------------------------------------------------------------------------
SMD test-suite contains 12 bilevel test problems that are scalable in terms of number of variables at upper and lower levels. All the problems contain single-objective to be minimized at both levels. A description of the SMD test problems can be found in the attached PDF file. You can refer to the following paper for further details about the SMD test-suite.

Authors: Ankur Sinha, Pekka Malo, and Kalyanmoy Deb.
Paper Title: Test problem construction for single-objective bilevel optimization.
Journal: Evolutionary Computation Journal.
Year: 2013.
------------------------------------------------------------------------



Files in the package
------------------------------------------------------------------------
There are the following Matlab files in this package:
upperLevelSMD.m
lowerLevelSMD.m
getOptimalSolutionSMD.m
------------------------------------------------------------------------



Executing the code
------------------------------------------------------------------------
Files upperLevelSMD.m and lowerLevelSMD.m contain the source code for upper and lower levels of SMD test problem. They can be called by using the following command on Matlab:
[upperFunctionValue upperEqualityConstrVals upperInequalityConstrVals]=upperLevelSMD(upperLevelMember, lowerLevelMember, testProblemName)
[lowerFunctionValue lowerEqualityConstrVals lowerInequalityConstrVals]=lowerLevelSMD(upperLevelMember, lowerLevelMember, testProblemName)

Both commands take three inputs, i.e. upper level decision vector, lower level decision vector and name of the test problem. The name of the test problem is a string. The output again consists of three variables that are function value, equality constraint values and inequality constraint values. In case there are no equality constraints or inequality constraints in the problem then the code returns an empty set.

For example, if the upper and lower level decision vectors are given as:
upperLevelMember=[6.425 5.245 8.656 2.943 4.142];
lowerLevelMember=[2.246 1.873 7.456 4.734 8.335];
then SMD1 can be evaluated as follows:
[upperFunctionValue upperEqualityConstrVals upperInequalityConstrVals]=upperLevelSMD(upperLevelMember, lowerLevelMember, 'smd1')
[lowerFunctionValue lowerEqualityConstrVals lowerInequalityConstrVals]=lowerLevelSMD(upperLevelMember, lowerLevelMember, 'smd1')

This returns:
upperFunctionValue = -2.6919e+03
upperEqualityConstrVals = []
upperInequalityConstrVals = []
lowerFunctionValue = 2.6660e+03
lowerEqualityConstrVals = []
lowerInequalityConstrVals = []

The file getOptimalSolutionSMD.m can be used to find the exact optimal solution to the SMD test problems. It can be called by using the following command:
[upperMember lowerMember upperFunctionValue lowerFunctionValue]=getOptimalSolutionSMD(upperDimensions,lowerDimensions,testProblemName)

It takes as input the number of upper level dimensions, the number of lower level dimensions and test problem name. It returns the optimal upper level member, the optimal lower level member, the optimal upper level function value and the optimal lower level function value as output.

For example, if one wants to find out the optimal solution of SMD1 test problem with 5 dimensions at the upper level and 5 dimensions at the lower level, then one needs to call the following command:
[upperMember lowerMember upperFunctionValue lowerFunctionValue]=getOptimalSolutionSMD(5,5,'smd1')

This returns:
upperMember = [0     0     0     0     0]
lowerMember = [0     0     0     0     0]
upperFunctionValue = 0
lowerFunctionValue = 0

Please note that the file getOPtimalSolutionSMD.m makes a call to upperLevelSMD.m and lowerLevelSMD.m in order to compute the function values at the optimum.
------------------------------------------------------------------------


Contact
------------------------------------------------------------------------
In case you have any questions, comments, suggestions, or you want to report any bugs, you can send an email to Ankur Sinha (Ankur.Sinha@aalto.fi)

Ankur Sinha, PhD
Aalto University School of Business
Helsinki Finland
------------------------------------------------------------------------