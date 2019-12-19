function bs_plot_grangernet(Fxy,Fyx,location,thre,time,foilim,chan,label)
% BS_PLOT_GRANGERNET visualize the Granger causality network.
% The line of granger network is judged by the max peak causality along
% with frequencies in each time point.
% 
% Syntax:
%   bs_plot_grangernet(Fxy,Fyx,location,thre,time,foilim,chan,label)
% 
% Input(s):
%   Fxy,Fyx     - Granger causality 
%   location    - Location of the sites
%   thre        - Threshold
%   time        - Specify the window number
%   foilim(1)   - Starting frequency
%   foilim(2)   - Ending frequency
%   chan        - Channels of interest
%   label       - channel name
% 
% See also: ganetwork, conetwork, mov_bi_ga, one_bi_ga, granger_causality_network.

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.2$ $Date: 14-Sep-2007 10:38:30$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang

% Using: draw_layout (, make_layout, arrow, textoval, ellipse)

% 20170331 Yuasa: modified "ganetwork" to use gui

%-- parameter settings
% time  = time;
% foilim = foilim;
% Fxy      = Fxy;
% Fyx      = Fyx;
% location = location;

%-- error checking
if ( foilim(1)>size(Fxy,2) )
    errordlg('please input correct start frequency','parameter lost');
    return
end
if ( foilim(2)>size(Fxy,2) )
    errordlg('please input correct end frequency','parameter lost');
    return
end
if(ndims(Fxy)>2)
    if ( time>size(Fxy,3) )
        errordlg('please input correct window','parameter lost');
        return
    end
end
if(ndims(Fxy)==2)
    if ( time~=1 )
        errordlg('please input correct window','parameter lost');
        return
    end
end

%-- find channel numbers
Ncmb = size(Fxy,1);
Nchans = (1+sqrt(1+8*Ncmb))/2;    % Nchans = (1+sqrt(1+8Ncmb))/2, Ncmb = number of pairs of coherence = (Nchans-1)Nchans/2

if nargin < 7 || isempty(chan)
    chan = 1:Nchans;
end
if nargin < 8 || isempty(label)
    label = num2cellstr(1:Nchans);
end
if isempty(location)
    location = prepare_layout(Nchans);
end

wholetime = size(Fxy,3);
wholefreq = size(Fxy,2);
gcpowmax  = round(nanmax([Fxy(:);Fyx(:)]),3,'significant');

%-- main
gch = figure('Name','Granger Causality Network','NumberTitle','off','Menubar','none','ToolBar','none');
gpa = plot_gc_layout(Fxy,Fyx,thre,time,foilim,chan,label,location);
gpa.Tag = 'Main';

setappdata(gch,'Fxy',Fxy);
setappdata(gch,'Fyx',Fyx);
setappdata(gch,'thre',thre);
setappdata(gch,'time',time);
setappdata(gch,'foilim',foilim);
setappdata(gch,'chan',chan);
setappdata(gch,'location',location);
setappdata(gch,'label',label);

axpos = gpa.Position;
axrsp = 1-(axpos(1)+axpos(3));
axusp = 1-(axpos(2)+axpos(4));
htl = uicontrol('Style','text','Tag','Tlabel',...
                'Units','normalized','Position',[axpos(1)-0.09 axpos(2)/3 0.08 axpos(2)/3],...
                'String','Time','HorizontalAlignment','right');
hsl = uicontrol('Style','slider','Tag','Tslider',...
                'Units','normalized','Position',[axpos(1) axpos(2)/3 axpos(3) axpos(2)/3],...
                'Min',1,'Max',wholetime,'SliderStep',[1/(wholetime-1) 2/(wholetime-1)],'Value',time,...
                'Callback',@timesld);
hfl = uicontrol('Style','text','Tag','Flabel',...
                  'Units','normalized','Position',[(axpos(1)+axpos(3))+axrsp/8 (axpos(2)+axpos(4))+axusp/20 axrsp*3/4 axusp/2],...
                'String','Freq');
hfd = uicontrol('Style','slider','Tag','FMINslider',...
                'Units','normalized','Position',[(axpos(1)+axpos(3))+axrsp/6 axpos(2) axrsp/3 axpos(4)],...
                'Min',1,'Max',wholefreq,'SliderStep',[1/(wholefreq-1) 5/(wholefreq-1)],'Value',foilim(1),...
                'Callback',@freqsld);
hfu = uicontrol('Style','slider','Tag','FMAXslider',...
                'Units','normalized','Position',[(axpos(1)+axpos(3))+axrsp/2 axpos(2) axrsp/3 axpos(4)],...
                'Min',1,'Max',wholefreq,'SliderStep',[1/(wholefreq-1) 5/(wholefreq-1)],'Value',foilim(2),...
                'Callback',@freqsld);
hthr = uicontrol('Style','pushbutton','Tag','THREb',...
                'Units','normalized','Position',[axpos(1)/5 (axpos(2)+axpos(4))+axusp/10 axpos(1)*4/5 axusp*3/5],...
                'String','Threshold','Callback',@threshedit);
hthsl = uicontrol('Style','slider','Tag','THREslider',...
                'Units','normalized','Position',[axpos(1)/2 axpos(2) axrsp/3 axpos(4)],...
                'Min',0,'Max',gcpowmax,'SliderStep',[round(gcpowmax/20,1,'significant') round(gcpowmax/5,1,'significant')],'Value',thre,...
                'Callback',@threshsld);

end%bs_plot_grangernet

function hax = plot_gc_layout(Fxy,Fyx,thre,time,foilim,chan,label,location)
%-- find channel numbers
Ncmb = size(Fxy,1);
Nchans = (1+sqrt(1+8*Ncmb))/2;    % Nchans = (1+sqrt(1+8Ncmb))/2, Ncmb = number of pairs of coherence = (Nchans-1)Nchans/2

%-- main
foiidx = foilim(1):foilim(2);
[dag, ispeak] = gc_data_processing(Fxy,Fyx,thre,time,foiidx,chan);

%-- plot figure
tstr = sprintf('Window: %d, Frequency: %.0f - %.0f Bins, Threshold: %.2f',time,foilim(1),foilim(2),thre);
if ~ispeak
    tstr = sprintf('No peak found. %s',tstr);
end
%-- set parameters for draw
circle = zeros(Nchans,1);
circle(chan) = 1; 
left   = location(1,:);
right  = location(2,:);
%-- plot main
cla;
draw_layout(dag,label,circle,left,right,'gc');
hax = gca;
title(tstr,'Tag','Title');

set(hax.Children.findobj('Type','text'),'ButtonDownFcn', @chansel);
set(hax.Children.findobj('Type','patch','Tag',''),'ButtonDownFcn', @chansel);

end%plot_gc_layout

function [dag, ispeak] = gc_data_processing(Fxy,Fyx,thre,time,foiidx,chan)

Ncmb = size(Fxy,1);
Nchans = (1+sqrt(1+8*Ncmb))/2;    % Nchans = (1+sqrt(1+8Ncmb))/2, Ncmb = number of pairs of coherence = (Nchans-1)Nchans/2

ispeak  = true;
dag     = zeros(Nchans,Nchans);
dag2    = dag;
dataF   = squeeze(Fxy(:,:,time))';
dataB   = squeeze(Fyx(:,:,time))';

%-- get upper triangular part index along row direction
 triidx = [];
 [triidx(:,1), triidx(:,2)] = find(tril(ones(Nchans),-1));
 triidx = num2cell(fliplr(triidx));
 
for icmbch = 1:Ncmb
    %-- x2y
    peak=findpeaks(dataF(foiidx,icmbch));       % need SignalProcessingToolbox
    if ~isempty(peak) && max(peak)>thre
        dag(triidx{icmbch,[1 2]})=1;
    end
    %-- y2x
    peak=findpeaks(dataB(foiidx,icmbch));       % need SignalProcessingToolbox
    if ~isempty(peak) && max(peak)>thre
        dag(triidx{icmbch,[2 1]})=1;
    end
end

%-- set channels of interests
chanary = false(1,Nchans);
chanary(chan) = true;
dag(~chanary,:) = 0;
dag(:,~chanary) = 0;

if max(dag(:)==dag2(:))
%-- if no peak is found
    ispeak  = false;
    for icmbch=1:Ncmb
        %-- x2y
        if max(dataF(foiidx,icmbch))>thre
            dag(triidx{icmbch,[1 2]})=1;
        end
        %-- y2x
        if max(dataB(foiidx,icmbch))>thre
            dag(triidx{icmbch,[2 1]})=1;
        end
    end
    %-- remove channels of no-interests
    dag(~chanary,:) = 0;
    dag(:,~chanary) = 0;
end

end%gc_data_processing

function layout_data = prepare_layout(num)
% make default layout map

%-- set square or rectangle
nx = ceil(sqrt(num));
lstclmn  = mod(num,nx);
lstclmn2 = mod(num,nx+1);
if lstclmn && (~lstclmn2 || lstclmn2 >= lstclmn) 
  nx = nx + 1;
end
ny = ceil(num/nx);

posx = 0.9/(nx+1)*(1:nx)+0.05;
posy = 0.9/(ny+1)*sort([1:2:ny (2:2:ny)+0.1])+0.05;

[meshx, meshy] = meshgrid(posx, posy);
if nx >= 3
  meshy(:,1:2:nx) = meshy(:,1:2:nx) - (meshy(1,1)-0.05)*0.5/2;
  meshy(:,2:2:nx) = meshy(:,2:2:nx) + (meshy(1,1)-0.05)*0.5/2;
end
if ny >= 3
  meshx(1:2:ny,:) = meshx(1:2:ny,:) - (meshx(1,1)-0.05)*0.2/2;
  meshx(2:2:ny,:) = meshx(2:2:ny,:) + (meshx(1,1)-0.05)*0.2/2;
end

layout_data   = [reshape(meshx,1,[]); reshape(meshy,1,[])];

end%prepare_layout

function c = num2cellstr(a)
    c = cell(size(a));
    for i=1:numel(a)
        c{i} = num2str(a(i));
    end 
end%num2cellstr


%----------------- GUI sub function -----------------------%
function timesld(hsl, cld)
    gch = hsl.Parent;
    
    Fxy      = getappdata(gch,'Fxy');
    Fyx      = getappdata(gch,'Fyx');
    thre     = getappdata(gch,'thre');
    time     = getappdata(gch,'time');
    foilim   = getappdata(gch,'foilim');
    chan     = getappdata(gch,'chan');
    location = getappdata(gch,'location');
    label    = getappdata(gch,'label');
        
    hsl.Value   = round(hsl.Value);
    if time ~= hsl.Value
        time = hsl.Value;
        setappdata(gch,'time',time);
        
        gpa = plot_gc_layout(Fxy,Fyx,thre,time,foilim,chan,label,location);
    end

end%timesld

function freqsld(hfreq, cld)
    gch = hfreq.Parent;
    hfu = gch.Children.findobj('Tag','FMAXslider');
    hfd = gch.Children.findobj('Tag','FMINslider');
    
    Fxy      = getappdata(gch,'Fxy');
    Fyx      = getappdata(gch,'Fyx');
    thre     = getappdata(gch,'thre');
    time     = getappdata(gch,'time');
    foilim   = getappdata(gch,'foilim');
    chan     = getappdata(gch,'chan');
    location = getappdata(gch,'location');
    label    = getappdata(gch,'label');
        
    switch hfreq.Tag
        case 'FMAXslider'
            hfreq.Value = max(hfd.Value+2,round(hfreq.Value));
        case 'FMINslider'
            hfreq.Value = min(hfu.Value-2,round(hfreq.Value));
    end
    if max(foilim ~= [hfd.Value hfu.Value])
        foilim = [hfd.Value hfu.Value];
        setappdata(gch,'foilim',foilim);
        
        plot_gc_layout(Fxy,Fyx,thre,time,foilim,chan,label,location);
    end

end%freqsld

function threshedit(hthr, cld)
    gch = hthr.Parent;
    hthsl = gch.Children.findobj('Tag','THREslider');
    
    Fxy      = getappdata(gch,'Fxy');
    Fyx      = getappdata(gch,'Fyx');
    thre     = getappdata(gch,'thre');
    time     = getappdata(gch,'time');
    foilim   = getappdata(gch,'foilim');
    chan     = getappdata(gch,'chan');
    location = getappdata(gch,'location');
    label    = getappdata(gch,'label');    
    
    hthr.Value = str2double(cell2mat(inputdlg('Threshold?','',1,num2cellstr(thre))));
    if isnan(hthr.Value), hthr.Value = thre;  end
    if hthr.Value < hthsl.Min, hthr.Value = hthsl.Min;  end
    if hthr.Value > hthsl.Max, hthr.Value = hthsl.Max;  end

    if thre ~= hthr.Value
        thre = hthr.Value;
        hthsl.Value = thre;
        setappdata(gch,'thre',thre);
        
        plot_gc_layout(Fxy,Fyx,thre,time,foilim,chan,label,location);
    end

end%threshedit

function threshsld(hthsl, cld)
    gch = hthsl.Parent;
    
    Fxy      = getappdata(gch,'Fxy');
    Fyx      = getappdata(gch,'Fyx');
    thre     = getappdata(gch,'thre');
    time     = getappdata(gch,'time');
    foilim   = getappdata(gch,'foilim');
    chan     = getappdata(gch,'chan');
    location = getappdata(gch,'location');
    label    = getappdata(gch,'label');
    
    hthsl.Value = round(hthsl.Value,3,'significant');

    if thre ~= hthsl.Value
        thre = hthsl.Value;
        setappdata(gch,'thre',thre);
        
        plot_gc_layout(Fxy,Fyx,thre,time,foilim,chan,label,location);
    end

end%threshsld

function chansel(hobj, cld)
    gpa = hobj.Parent;
    gch = gpa.Parent;
    Mpos = gpa.CurrentPoint;
    
    Fxy      = getappdata(gch,'Fxy');
    Fyx      = getappdata(gch,'Fyx');
    thre     = getappdata(gch,'thre');
    time     = getappdata(gch,'time');
    foilim   = getappdata(gch,'foilim');
    chan     = getappdata(gch,'chan');
    location = getappdata(gch,'location');
    label    = getappdata(gch,'label');
    
    Mdist    = mean((location - repmat(Mpos(1,1:2)',1,size(location,2))).^2,1);
    selchidx = find(Mdist == min(Mdist),1);
    
    if max(chan == selchidx)
        chan(chan == selchidx) = [];
    else
        chan = sort([chan selchidx]);
    end
    
    setappdata(gch,'chan',chan);
    plot_gc_layout(Fxy,Fyx,thre,time,foilim,chan,label,location);

end%chansel
% [EOF]