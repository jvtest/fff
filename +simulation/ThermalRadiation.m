classdef ThermalRadiation < handle
    %
    % Numerically evaluates Planck's law
    %
    % Static members provide integration and sensitivity calculations. 
    %
    
    % Author: Åge Andreas Falnes Olsen, Justervesenet
    % Date: 25-oct-2012
    
    
    properties (SetAccess = protected)
        
        % Spectral radiance in [W/m^2/sr] per wavelength in micrometers.
        % The choice of units is a bit weird but tries to accomodate the
        % fact that units for wavelength are most practical in micrometers.
        % 
        l; 
        lambda; % wavelength in um
        temperature; % the temperature in Kelvin
                
        c1; % first radiation constant.
        c2; % Second radiation constant
        
    end
    
    methods
        
        function this = ThermalRadiation(wavelengths, temperature)
            %
            % this = ThermalRadiation()
            % this = ThermalRadiation(wavelengths, temperature)
            %
            % wavelengths    - wavelengths to compute radiation intensity
            %                  for. Micrometer. Vector or scalar.
            % temperature    - Temperature. Vector or scalar. In Kelvin.
            %
            % wavelengths and temperature can not both be a vector at the
            % same time; one of them must be a scalar. If wavelengths is a
            % vector, the values in class property "l" are intensity as a
            % function of wavelength for a fixed temperature. Conversely,
            % if temperature was input as a vector the resulting "l" values
            % are intensity at a fixed wavelength versus temperature. 
            %
            
            
            % Convert the constants so that wavelengths are given in
            % micrometers. The radiation intensity is W per wavelength and
            % hence the unit conversion for the first radiation constant
            % seems a bit unmotivated at first.
            
            this.c1 = Constants.c1(1) .* 1e24;
            this.c2 = Constants.c2(1) .* 1e6;
            
            
            if nargin == 2
                this.temperature = temperature;
                this.lambda = wavelengths;
                this.l = this.planckLaw(wavelengths, temperature);
            end
            
        end
        
        function irr = planckLaw(this, w, t)
            %
            % Computes the radiation intensity at (wavelength, temperature)
            %
            % irr = planckLaw(w, t)
            %
            % w     - wavelengths (vector, um)
            % t     - temperature (scalar, Kelvin)
            %
            % The output is in [W/m^2] per micrometer wavelength per solid
            % angle (sr).
            %
            
            irr = this.c1 ./ w .^ 5;
            irr = irr ./ (exp(this.c2 ./ (w .* t)) - 1);
            
        end
        
        function irr = wienLaw(this, w, t)
            %
            % Computes the Wien approximation to the Planck law.
            %
            
            irr = this.c1 ./ w .^ 5;
            irr = irr ./ (exp(this.c2 ./ (w .* t)));
            
        end
        
    end
    
    methods (Static = true)
        
        function res = sensitivity(l0, l1, t)
            %
            % Computes sensitivity of imaginary sensor to temperature
            %
            % l0    - minimum wavelength (scalar, um)
            % l1    - maximum wavelength (scalar, um)
            % t     - temperatures (vector)
            %
            % res   - output sensitivity, size(res) = size(t)
            %
            % The sensitivity is RELATIVE. The computation is based on
            % manual differentiation of the Planck law. Numeric
            % differentation leads to numeric noise.
            %
            
            obj = simulation.ThermalRadiation();

            res = obj.integrate(l0, l1, t);
            n = length(res);
            for q = 1:n
                fn = @(x)obj.planckLaw(x, t(q)) ./ (x .* (1 - exp(-obj.c2 ./ (x.*t(q)))));
                res(q) = integral(fn, l0, l1) / res(q);
            end
            res = res .* obj.c2 ./ t;
            
        end
        
        function res = integrate(l0, l1, t, rl)
            %
            % Integrates the Planck law in a specified band
            %
            % res = integrate(l0, l1, t)
            % res = integrate(l0, l1, t, rl)
            %
            % l0    - lower wavelength range (um)
            % l1    - upper wavelength range (um)
            % t     - temperatures (Kelvin)
            % rl    - optional function of wavelength to be multiplied with
            %         the Planck law
            %
            % res   - output. Size(res) = size(t).
            %
            
            obj = simulation.ThermalRadiation();
            if nargin == 3
                rl = @(x)1;
            end
            
            res = zeros(size(t));
            n = length(res);
            for q = 1:n
                fn = @(x)obj.planckLaw(x, t(q)) .* rl(x);
                res(q) = integral(fn, l0, l1);
            end
            
        end
        
        
        
    end
    
end