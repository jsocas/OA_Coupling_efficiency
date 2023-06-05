function atm = ogsAtmosphere(ID,corrF,varargin)

p = inputParser;
p.addRequired('ID', @isnumeric);
p.addRequired('corrF', @isnumeric); % Correction factor to r0 and fr0 depending on the telescope elevation
p.addParamValue('L0',30,@isnumeric);
p.addParamValue('vs',[5 10 15 20],@ismatrix);
p.addParamValue('vd',[0 pi/4 pi pi],@ismatrix);
p.addParamValue('z',[0 2000 4000 10000],@ismatrix);
p.addParamValue('semilla', 0,@isnumeric);
p.parse(ID, corrF, varargin{:});
ID = p.Results.ID;
corrF = p.Results.corrF;
L0 = p.Results.L0;
vs = p.Results.vs;
vd = p.Results.vd;
z = p.Results.z;
semilla = p.Results.semilla;

s = RandStream('mt19937ar', 'seed', semilla);

% z  = [0 2000 4000 10000];
% vs = [5 10 15 20];
% vd = [0 pi/4 pi pi];
fprintf(' *** OGS ATMOSPHERE \n')
switch ID
    case 10
        fprintf('Case 10: PERFECT ATM\n ***\n')
        r0 = corrF * 100;
        fr0 = corrF * [1 1 1 1];
        tag = 'PERFECT';
    case 1
        fprintf('Case 1: NIGHT Oct-Nov-Dec-Feb\n ***\n')
        r0 = corrF * 0.1894;
        fr0 = corrF * [0.5600 0.1291 0.0533 0.2482];
        tag = 'OGS-ATM1';
    case 2
        fprintf('Case 2: NIGHT Mar-Apr-May\n ***\n')
        r0 = corrF * 0.1744;
        fr0 = corrF * [0.5927 0.1583 0.0532 0.1856];
        tag = 'OGS-ATM2';
    case 3
        fprintf('Case 3: NIGHT Jun-Jul-Aug-Sep\n ***\n')
        r0 = corrF * 0.1914;
        fr0 = corrF * [0.5650 0.1297 0.0974 0.1964];
        tag = 'OGS-ATM3';
    case 4
        fprintf('Case 4: DAY Jun-Jul-Aug Sun-elevation:60�\n ***\n')
        r0 = corrF * 0.0373;
        fr0 = corrF * [0.9284 0.0170 0.0113 0.0427];
        tag = 'OGS-ATM4';
    case 5
        fprintf('Case 5: DAY Jun-Jul-Aug Sun-elevation:35�\n ***\n')
        r0 = corrF * 0.0280;
        fr0 = corrF * [0.6038 0.0204 0.0773 0.2912];
        tag = 'OGS-ATM5';
    case 6
        fprintf('Case 6: DAY Sep-Oct-Nov Sun-elevation:72�\n ***\n')
        r0 = corrF * 0.0917;
        fr0 = corrF * [0.9946 0.0021 0.0009 0.0036];
        tag = 'OGS-ATM6';
    case 7
        fprintf('Case 7: DAY Sep-Oct-Nov Sun-elevation:40�\n ***\n')
        r0 = corrF * 0.0293;
        fr0 = corrF * [0.9943 0.0019 0.0014 0.0052];
        tag = 'OGS-ATM7';
    case 8
        fprintf('Case 8: DAY STRONG ATM\n')
        r0 = corrF * 0.0200;
        fr0 = corrF * [0.9943 0.0019 0.0014 0.0052];
        tag = 'OGS-ATM8';
    otherwise
        error('Oups! only 7 cases, so far!')
end

atm = atmosphere(photometry.V,r0,L0,...
    'altitude',z,...
    'fractionnalR0',fr0,...
    'windSpeed',vs,...
    'windDirection',vd,...
    'randStream',s);
atm.tag = tag;