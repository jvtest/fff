function ThermalRadiationSELFTEST()
%
% SELFTEST for ThermalRadiation.

interface();

integrateAll();

displacementLaw();

end


function interface()
%
% Tests that class can be constructed and gives sensible outputs.
%

temp = 400;
wl = linspace(1, 20, 1000);
tr = simulation.ThermalRadiation(wl, temp);

assert(tr.temperature == temp, 'JV:SelftestFailed', ...
    'Temperature property not correct');
assert(all(tr.lambda == wl), 'JV:SelftestFailed', ...
    'Wavelength property not correct');
assert(numel(tr.l) == numel(tr.lambda),'JV:SelftestFailed', ...
    'Output radiation intensity vector has wrong size');


tr2 = simulation.ThermalRadiation();
il = tr2.planckLaw(wl, temp);
assert(all(il == tr.l), 'JV:SelftestFailed', ...
    'Expected the two class usages to produce the same result');

end

function integrateAll()
%
% Tests that we reproduce the Stefan-Boltzmann law.
%

sb = 5.670373*1e-8; % Stefan-Boltzmann constant
tol = 3.7e-6;

temp = 50:50:1000;
ir = temp;
for q = 1:length(temp)
    ir(q) = simulation.ThermalRadiation.integrate(0, inf, temp(q));
end

warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');

pp = polyfit(temp, ir, 4);
assert(abs(1 - sb/(pi*pp(1))) <= tol, 'JV:SelftestFailed', ...
    'Integration did not reproduce Stefan-Boltzmann law with adequate precision');
assert(all(abs(pp(2:end) ./ pp(1)) <= 1e-4), 'JV:SelftestFailed', ...
    'Integration did not reproduce Stefan-Boltzmann law with adequate precision');


warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');


end

function displacementLaw()
%
% Tests that the displacement law is accurately reproduced
%

b = Constants.b(1);

t = 50:5:2000;
nt = length(t);

l = 0.1:0.01:70;

lPeak = zeros(size(t));

tr = simulation.ThermalRadiation();

for q = 1:nt
    y = tr.planckLaw(l, t(q));
    [~, k] = max(y);
    lPeak(q) = l(k);
end


pp = polyfit(1./t, lPeak .* 1e-6, 1);
% First test: require the constant term to be less than the error we would
% have seen if we had two measurements of lPeak and t, and solved the
% linear equations.
TOL = (l(2) - l(1))^2 * (1 - 2*t(end) + 2*t(end)^2) / (t(end) - t(1))^2;
TOL = sqrt(TOL) * 1e-6;
assert(abs(pp(2)/pp(1)) < TOL, 'JV:SelftestFailed', ...
    'Wien displacement law is not reproduced with satisfactory precision');
% Second test: require the linear constant term to be less than the error
% we would have seen if we had two measurements of lPeak and t, and solved
% the linear equations.
TOL = 2*((l(2)-l(1))^2 * 1e-12*t(1)*t(end)/(t(end) - t(1)));
TOL = sqrt(TOL);
assert(abs(pp(1) - b) < TOL, 'JV:SelftestFailed', ...
    'Wien displacement law is not reproduced with satisfactory precision');


end




