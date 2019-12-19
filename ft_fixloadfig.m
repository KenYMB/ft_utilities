function    ft_fixloadfig(handlenum)

% It fix appdata name to enable interactive mode in fieldTrip.
% You need to run this function after you load a saved figure.
% Both of figure handle and axis handle is ok.
% Multi-plot figure do not need this fix.
% 
% Usage:
%   uiopen(FIGUREPATH,1)
%   ft_fixloadfig(gcf)

% using: fieldtrip, ft_private(fixname), cellstrfind

% 20170203 Yuasa: create

%-- set figure handle
figh    = handlenum;
while(~strcmp(get(figh,'Type'),'figure'))
    figh     = get(figh,'Parent');
end
%-- axis figure handle
axish   = handlenum;
isbreak = false;
if ~strcmp(get(axish,'Type'),'axes')
    axish    = get(axish,'Children');
    axesType = get(axish,'Type');
    axesidx  = cellstrfind(axesType,'axes');
    if isempty(axesidx)
        warning('No figure axis is found.');
        isbreak = true;                     % skip ft_fixloadfig
    elseif length(axesidx) > 1
        warning('Too many figure axes are found. Cannot specify the corresponding axis.');
        isbreak = true;                     % skip ft_fixloadfig
    else
        axish   = axish(axesidx);
    end
end

if ~isbreak
  %-- get all appdata
  appdatas        = getappdata(figh);
  %-- get names of appdata
  appdatalist     = fieldnames(appdatas);
  %-- find the index of appdata for fieldTrip
  appdataidx      = cellstrfind(appdatalist,'x*_*');
  
  if isempty(appdataidx)                    % skip ft_fixloadfig
      warning('No data is saved for fieldTrip interactive mode.');
  elseif length(appdataidx) ~= 1            % skip ft_fixloadfig
      warning('Too many datas are saved. Cannot specify the corresponding data.');
  else    % length(appdataidx) == 1
      ft_private('set','silent');
      
      olddataname = appdatalist{appdataidx};
      dataname    = fixname(num2str(double(axish)));
      if ~strcmp(dataname, olddataname)
          data        = appdatas.(olddataname);
          setappdata(figh,dataname,data)
          rmappdata(figh,olddataname)
      end
      
      ft_private('remove','silent');
  end
end

    