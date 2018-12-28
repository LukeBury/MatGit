%-Abstract
%
%   CSPICE_ROTMAT calculates the rotation matrix generated by
%   a rotation of a specified angle about a specified axis applied
%   to a matrix. This rotation is thought of as rotating the
%   coordinate system.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY,
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING,
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      m1      the matrix to apply the rotation. In matrix algebra, the
%              components of the matrix are relative to one particular
%              coordinate system. Applying cspice_rotmat changes the
%              components of 'm1' so that they are relative to a rotated
%              coordinate system.
%
%              [3,3]   = size(m1); double = class(m1)
%
%      angle   the angle in radians through which to rotate the original
%              coordinate system.
%
%              [1,1]   = size(angle); double = class(angle)
%
%      iaxis   the index for the axis of the original coordinate system
%              about which to perform the rotation by 'angle'.
%              'iaxis' = 1,2 or 3 respectively designates the x-, y-,
%              or z-axis.
%
%              [1,1]   = size(iaxis); int32 = class(iaxis)
%
%   the call:
%
%      mout = cspice_rotmat( m1, angle, iaxis)
%
%   returns:
%
%      mout    matrix resulting from the application of 'angle' to the
%              input matrix 'm1'.
%
%              [3,3]   = size(mout); double = class(mout)
%
%              If
%
%                 [angle]
%                       iaxis
%
%              denotes the rotation matrix by 'angle' radians about 'iaxis',
%              (see the Rotations Required Reading document) then 'mout' is
%              given by the following matrix equation:
%
%                 mout = [angle]      * m1
%                               iaxis
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Rotate the 3x3 identity matrix by 90 degrees about
%      % the y axis.
%      %
%
%      %
%      % Create the 3x3 identity matrix.
%      %
%      ident = eye(3);
%
%      %
%      % Rotate 'ident' by Pi/2 about the Y axis.
%      %
%      r = cspice_rotmat( ident, cspice_halfpi, 2 )
%
%   MATLAB outputs:
%
%      r =
%
%          0.0000         0   -1.0000
%               0    1.0000         0
%          1.0000         0    0.0000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine rotmat_c.
%
%   MICE.REQ
%   ROTATION.REQ
%
%-Version
%
%   -Mice Version 1.1.1, 10-MAR-2015, EDW (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.1.0, 24-JAN-2009, EDW (JPL)
%
%      Corrected the function definition name. This wrapper had a
%      the function name "cspice_rotate" instead of "cspice_rotmat."
%
%   -Mice Version 1.0.0, 17-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   rotate a matrix
%
%-&

function [mout] = cspice_rotmat( m1, angle, iaxis )

   switch nargin
      case 3

         m1    = zzmice_dp(m1);
         angle = zzmice_dp(angle);
         iaxis = zzmice_int(iaxis);

      otherwise

         error ( [ 'Usage: [mout(3,3)] = ' ...
                   'cspice_rotmat( m1(3,3), angle, iaxis )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [mout] = mice('rotmat_c', m1, angle, iaxis );
   catch
      rethrow(lasterror)
   end

