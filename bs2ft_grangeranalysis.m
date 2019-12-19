function [data] =  bs2ft_grangeranalysis(Fxy, Fyx, varargin)

% BS2FT_GRANGERANALYSIS convert bsmart output data of granger causality
% analysis into fieldTrip structure.
% 
% Usage:
%   data =  bs2ft_grangeranalysis(Fxy, Fyx, time,freq,label,param)
%   data =  bs2ft_grangeranalysis(Fxy, Fyx, cfg)
% 
% Output:
%   data     - converted granger causality data in fieldTrip format
% 
% Inputs: 
%   Fxy      - matrix, the causality measure from x to y (combined_chs x frequencies x times)
%   Fyx      - matrix, the causality measure from y to x (combined_chs x frequencies x times)
%              The order of Fx2y/Fy2x is 1 to 2:N; 2 to 3:N; ...; N-1 to N,
%               where N is the number of channels.
%              That is,
%               1st row=1&2; 2nd=1&3; ...; (N-1)th=1&N; ...; (N(N-1)/2)th=(N-1)&N.
% 
%   cfg      - structure, you can input all the following parametes as 'cfg.***'
% 
%   time     - vector, physical time of each time point [s] (1 x times)
%   freq     - vector, frequencies of each frequency point [Hz] (1 x frequencies)
%   label    - cellstr, labels of each channel (label_name x 1)
%   param    - (optional)string, desired field name of granger causality data [default: 'grangerspctrm']

% 20170329 Yuasa

narginchk(3,6);
if ~isstruct(varargin{1})
    assert(nargin >= 5, message('MATLAB:narginchk:notEnoughInputs'));
    time    = varargin{1};
    freq    = varargin{2};
    label   = varargin{3};
    try param   = varargin{4}; catch, param = 'grangerspctrm'; end
else
    cfg     = varargin{1};
    time    = cfg.time;
    freq    = cfg.freq;
    label   = cfg.label;
    try param   = cfg.param; catch, param = 'grangerspctrm'; end
end

nlabel      = length(label);
cmbdata     = zeros(nlabel,nlabel,length(freq),length(time));
ixy = 1; 
for ilp = 1:nlabel
  for jlp = (ilp+1):nlabel
      cmbdata(ilp,jlp,:,:)  = Fxy(ixy,:,:);
      cmbdata(jlp,ilp,:,:)  = Fyx(ixy,:,:);
      ixy = ixy + 1;
  end
end

data   = [];
data.label           = label;
data.dimord          = 'chan_chan_freq_time';
data.(param)         = cmbdata;
data.time            = time;
data.freq            = freq;
