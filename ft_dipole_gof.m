function gof = ft_dipole_gof(dip)
% 
% GOF = FT_DIPOLE_GOF(dip)
% 
% 'dip' is output variable of ft_dipolefit
% GOF is calculated by dip.Vdata and dip.Vmodel.
% 
% GOF based on residual variance is calculated by following equation, 
%   GOF = 1 - abs(Vdata - Vmodel)^2 / sum(Vdata^2)

% R2 GOF is calculated by following equation, 
%   GOF = 1 - norm(Vdata - Vmodel) / norm(Vdata - mean(Vdata))

% 20160907: Yuasa

%-- R2 gof 
% residual = mean((dip.Vdata - dip.Vmodel).^2,1);
% gof      = 1- (residual ./ var(dip.Vdata, 1, 1));

%-- Xi2 gof
% Xi2 = sum((dip.Vdata - dip.Vmodel).^2 ./ repmat(var(dip.Vdata, 0, 1),size(dip.Vdata,1),1),1);
% gof      = 1- Xi2 / (size(dip.Vdata,1) -1);

%-- rv gof
rv = sum((dip.Vdata - dip.Vmodel).^2,1) ./ sum(dip.Vdata.^2,1);
gof      = 1- rv;