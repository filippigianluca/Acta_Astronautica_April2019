function angle=qck(angle)

% %***********************************************************************
% *     This function takes any angle and reduces it, if necessary,
% *     so that it lies in the range from 0 to 2 PI radians.
% *
% *     INPUTS TO THE FUNCTION
% *          ANGLE   = The ange to be reduced (in radians)
% *
% *     OUTPUTS FROM THE FUNCTION
% *          QCK     = The angle reduced, if necessary, to the range
% *                    from 0 to 2 PI radians (in radians)
% *
% *     MISSION PLANNING SUBROUTINES AND FUNCTIONS CALLED
% *          PI
% *
% *     PROGRAMMER:    W.T. Fowler
% *
% *     DATE:          July, 1978
% *
% *     VERIFIED BY:   Darrel Monroe, 8/20/90
% *                    Matteo Ceriotti - 10-01-2007
% *
% ************************************************************************

twopi = 2*pi;
 
diff = twopi * (fix(angle/twopi) + min([0,sign(angle)]));

angle = angle -diff;
