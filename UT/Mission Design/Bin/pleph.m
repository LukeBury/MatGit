function [r,v]=pleph(et,ntarg,ncent,kmflag)

% PLEPH  Read JPL planetary ephemeris and give body states
%    [R,V]=PLEPH(ET,NTARG,NCENT,KMFLAG) returns the position R
%    and velocity V of body NTARG relative to body NCENT at Julian
%    date ET.  The numbering convention for NTARG and NCENT is:
%
%     1 -> Mercury          8 -> Neptune
%     2 -> Venus            9 -> Pluto
%     3 -> Earth           10 -> Moon
%     4 -> Mars            11 -> Sun
%     5 -> Jupiter         12 -> Solar system barycenter
%     6 -> Saturn          13 -> Earth-Moon barycenter
%     7 -> Uranus          14 -> Nutations (longitude and obliquity)
%               15 -> Librations, if on file
%
%    If nutations or librations are desired, set NCENT = 0.
%
%    Units for position and velocity are AU and AU/day or, if KMFLAG
%    is nonzero, km and km/s.  Units for nutations and librations
%    are radians and radians/day.

% Declare variables

fname='binephem';
FID=fopen(fname,'r');

% EPH_HEADER_SIZE Set sizes for JPL Ephemeris header file parameters

% DE200 values
ttlsize  = [3 80];
cnamsize = 8;
sssize   = 3;
iptsize  = [3 13];

TTL=char(fread(FID,ttlsize,[num2str(ttlsize(2)) '*uchar=>uchar']));
SST=fread(FID,sssize,'float64=>float64');
NCON=fread(FID,1,'int32=>float64');
CNAM=char(fread(FID,[NCON cnamsize],[num2str(cnamsize) '*uchar=>uchar']));
CVAL=fread(FID,NCON,'float64=>float64');
AU=fread(FID,1,'float64=>float64');
EMRAT=fread(FID,1,'float64=>float64');
IPT=fread(FID,iptsize,[num2str(iptsize(2)),'*int32=>float64']);
NUMDE=fread(FID,1,'int32=>float64');
NCOEF=fread(FID,1,'int32=>float64');
HEADSIZE=ftell(FID);

% Handle NTARG=NCENT

%if ntarg==ncent   % Degenerate case
%    r = zeros(3,1);
%    v = zeros(3,1);
%    return
%end

% Set up LIST for call to STATE_EPH

list=zeros(1,12);

% Check for nutation call

if ntarg==14 || ncent==14
    list(11)=2;
    [r,v]=state_eph(et,list,1,kmflag,FID,SST,AU,IPT,NCOEF,HEADSIZE);
    r=r(1:2,12);
    v=v(1:2,12);
    return
end

%Check for libration call

if ntarg==15 || ncent==15
    list(12)=2;
    [r,v]=state_eph(et,list,1,kmflag,FID,SST,AU,IPT,NCOEF,HEADSIZE);
    r=r(:,13);
    v=v(:,13);
    return
end

% Check for Earth call

for i=1:2
    k = ntarg;
    if i ==2, k = ncent; end
    if k <= 10, list(k)  = 2; end % Planets and Moon
    if k == 10, list(3)  = 2; end % If Moon, get EM-Bary
    if k ==  3, list(10) = 2; end % If Earth, get Moon
    if k == 13, list(3)  = 2; end % If EM-bary, get EM-Bary
end

% Make call to STATE_EPH

[pv,vv]=state_eph(et,list,1,kmflag,FID,SST,AU,IPT,NCOEF,HEADSIZE);

if ntarg == 13 || ncent == 13 % EM-Bary
    pv(:,13) = pv(:,3);
    vv(:,13) = vv(:,3);
end

if ntarg*ncent == 30 && ntarg+ncent == 13 % Earth and Moon
    pv(:,3)=zeros(3,1);
    vv(:,3)=zeros(3,1);
else
    if list(3) == 2 % Get Earth from EM-Bary and Moon
        pv(:,3)=pv(:,3)-pv(:,10)/(1+EMRAT);
        vv(:,3)=vv(:,3)-vv(:,10)/(1+EMRAT);
    end
    if list(10) ==2 % Get Moon from EM-Bary and Geocentric Moon
        pv(:,10)=pv(:,3)+pv(:,10);
        vv(:,10)=vv(:,3)+vv(:,10);
    end
end

% Get NTARG relative to NCENT
r=pv(:,ntarg)-pv(:,ncent);
v=vv(:,ntarg)-vv(:,ncent);
fclose(FID);
end

function [pv,vv] = state_eph(jd,list,baryflag,kmflag,FID,SST,AU,IPT,NCOEF,HEADSIZE)

% STATE_EPH Read and interpolate JPL Planetary Ephemeris File
%     STATE_EPH(JD,LIST,BARYFLAG,KMFLAG,EPHFILE) reads and
%     interpolates a JPL planetary ephemeris file
%
%     INPUTS: JD      Julian Date at desired interpolation epoch
%
%             LIST    1-by-12 array indicating which bodies and what
%                     type of interpolation is desired:
%                     LIST(i) = 0 -> No interpolation for body i
%                             = 1 -> Position only
%                             = 2 -> Position and velocity
%
%                     Designation of astronomical bodies
%                      i = 1 -> Mercury
%                        = 2 -> Venus
%                        = 3 -> Earth-Moon barycenter
%                        = 4 -> Mars
%                        = 5 -> Jupiter
%                        = 6 -> Saturn
%                        = 7 -> Uranus
%                        = 8 -> Neptune
%                        = 9 -> Pluto
%                        = 10-> Geocentric Moon
%                        = 11-> Nutations in Earth longitude and 
%                               obliquity
%                        = 12-> Lunar librations (if on file)
%
%            BARYFLAG  = 0 -> Planetary states are heliocentric
%                     ~= 0 -> Planetary states are solar-system 
%                             barycentric
%
%            KMFLAG    = 0 -> Units are AU and AU/day
%                     ~= 0 -> Units are km and km/s
%
%            EPHFILE  Name of ephemeris file
%
%     OUTPUTS: PV       3-by-13 vector of positions
%              
%              VV       3-by-13 vector of velocities
%
%     States are relative to Earth mean equator and equinox of J2000
%     if the DE number is >= 200; of B1950 if the DE number is < 200.

% Declare global variables
% eph_global
% Check to see if JD is in the ephemeris file interval
if jd < SST(1) || jd > SST(2), error('Specified date not in ephemeris interval.'), end

% Determine which interval jd is within and read in the coefficients for that
% interval

iint=max([1 ceil((jd-SST(1))/SST(3))]);
fseek(FID,HEADSIZE+(iint-1)*NCOEF*8,-1);
db=fread(FID,NCOEF,'*float64'); % buf is ncoef-by-1 coefficient vector


% Normalize time to be on [0 1]

ts=SST(3);
t=(jd-db(1))/ts;

% Pre-allocate pv and vv matrices

pv=zeros(3,12);
vv=zeros(3,12);
au2km=AU; % define au2km conversion factor


% Get barycentric state of Sun

buf=getcoef(t,db,3,11,IPT);

[pvsun,vvsun]=int_eph(buf,t,ts,IPT(2,11),IPT(3,11));

% Get barycentric state of selected planets

for i=1:10
    if list(i) ~= 0
       buf=getcoef(t,db,3,i,IPT);
       [pv(:,i),vv(:,i)]=int_eph(buf,t,ts,IPT(2,i),IPT(3,i));
   end
end

% Put state in desired units

if ~kmflag
    pv=pv/au2km;
    vv=vv/au2km;
    pvsun=pvsun/au2km;
    vvsun=vvsun/au2km;
else
    vv=vv/86400;
    vvsun=vvsun/86400;
end

% If heliocentric state is desired (baryflag is false), subtract solar state from
% planetary state

if ~baryflag
    pv=pv-pvsun*[(list(1:10)>0) 0 0];
    vv=vv-vvsun*[(list(1:10)>0) 0 0];
end

% Get nutations and librations if requested and available

for i=11:12
    if (list(i) > 0 && IPT(2,i+1) > 0)
        buf=getcoef(t,db,i-9,i+1,IPT);
        [pnut,vnut]=int_eph(buf,t,ts,IPT(2,i+1),IPT(3,i+1));
        pv(1:size(pnut,1),i)=pnut;
        vv(1:size(vnut,1),i)=vnut;
    end
end

% Append barycentric solar state to pv and vv

pv=[pv(:,1:10) pvsun pv(:,11:12)];
vv=[vv(:,1:10) vvsun vv(:,11:12)];



end

function [pv,vv]=int_eph(buf,t,ts,ncf,na)% Subfunction INT_EPH

% INT_EPH Interpolate a set of chebyshev coefficients to give position and velocity
%

% Transform T from [0 1] to [-1 1]

tc=2*(mod(na*t,1)+fix(t))-1;

%tc=2*t-1;

% Evaluate Chebyshev polynomials at TC

pc=chebyval(tc,ncf-1);

% Multiply and add to get position vector

pv=fliplr(buf')*flipud(pc); % Flip matrices to sum from smalles to largest

% Differentiate to get velocity

pcp=chebyder(tc,pc);
vv=fliplr(buf')*flipud(pcp)*2/ts*na;

end

function buf=getcoef(t,db,ndim,ibod,IPT) % Subfunction getcoef

% eph_global;

l=fix(IPT(3,ibod)*t-fix(t));
ist=IPT(1,ibod)+l*ndim*IPT(2,ibod);
ien=IPT(1,ibod)+(l+1)*ndim*IPT(2,ibod)-1;

buf=reshape(db(ist:ien),[IPT(2,ibod),ndim]);
end

function t=chebyval(x,n)

% CHEBYVAL Evaluate Chebyshev polynomial
%     CHEBYVAL(X,N) returns the zeroth through the Nth Chebyshev
%     polynomial evaluated at X.  If X is scalar, CHEBYVAL returns
%     an N+1-by-1 vector whose elements are the zeroth through Nth
%     Chebyshev polynomials evaluated at X. If X is a p-by-1 or 
%     1-by-p vector, CHEBYVAL returns a p-by-N+1 or N+1-by-p matrix.
%     If X is a matrix, CHEBYVAL returns a matrix of size(X)-by-N+1.

c=size(x);
x=x(:); % Unwind x for easier evaluation
t=[ones(size(x)) x zeros(size(x,1),n-1)]; % pre-allocate t

for i=3:n+1
    t(:,i)=2*x.*t(:,i-1)-t(:,i-2); % Recurrence relation for Cheby. pol.
end

% Reshape t for output

if length(c)==2 & c(1)==1
    t=t';
elseif length(c)==2 & c(2)==1
else
    t=reshape(t,[c,n+1]);
end
end

function tp=chebyder(x,t)

% CHEBYVAL Evaluate Chebyshev polynomial derivatives
%     CHEBYVAL(X,T) returns the Chebyshev polynomial derivatives
%     corresponding to the Chebyshev polynomials evaluated at X.
%     T is the matrix of Chebyshev polynomials evaluated at X
%     (see CHEBYVAL).

% Determine dimensions

c=size(t);
cc=size(x);
n=c(length(c))-1;
if length(c)==2
    n=max(c)-1;
end

% Unwind pc and x for easier evaluation

x=x(:);
t=t(:)';

tp=[zeros(size(x)) ones(size(x)) zeros(size(x,1),n-1)]; % pre-allocate tp

% Differentiate recurrence relation for Chebyshev polynomials

for i=3:n+1
    tp(:,i)=2*(t(:,i-1)+x.*tp(:,i-1))-tp(:,i-2);
end

% Reshape tp to be the same dimensions as t

if length(cc)==2 & cc(1)==1
    tp=tp';
elseif length(cc)==2 & cc(2)==1
else
    tp=reshape(tp,[cc,n+1]);
end
end
