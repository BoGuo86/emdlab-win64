% physical constants

classdef emdlab_phy_constants
    
    properties (Constant = true)
        
        % permeability of vaccume [H/m]
        mu0 = 4*pi*1e-7;  
        
        % reluctivity of vaccume [m/H]
        nu0 = 1/(4*pi*1e-7);  
        
        % permittivity of vaccume [H/m]
        e0 = 8.85e-12;  
        
    end
    
end