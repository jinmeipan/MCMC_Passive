%snowpit structure
classdef prior

properties
    %basic
    swe_mean
    swe_std
    density_mean
    density_std
    T_mean
    T_std
    dmax_mean
    dmax_std
    pex_mean
    pex_std
    %directly given or calculated if given swe
    sd_mean
    sd_std
    %lognormal fit
    swe_mu
    swe_sigma
    density_mu
    density_sigma
    T_mu
    T_sigma
    dmax_mu
    dmax_sigma
    pex_mu
    pex_sigma
    sd_mu
    sd_sigma
    %note
    note
end

methods
    function prior=calc_sdprior(prior)
        SWE_mean=prior.swe_mean; %mm
        SWE_std=prior.swe_std; %to be revised
        density_mean=prior.density_mean; %kg/m^3
        density_std=prior.density_std;   %kg/m^3
        
        sd_mean=SWE_mean/density_mean; %tundra density
        sd_cov=SWE_mean^2./density_mean^2 * (SWE_std^2/SWE_mean^2 + density_std^2/...
            density_mean^2);
        prior.sd_mean=sd_mean;
        prior.sd_std=sd_cov^0.5;
    end
    
    function prior=calc_logparam(prior)
        m=prior.swe_mean;
        v=prior.swe_std.^2;
        prior.swe_mu = log((m^2)/sqrt(v+m^2));
        prior.swe_sigma = sqrt(log(v/(m^2)+1));

        m=prior.density_mean;
        v=prior.density_std.^2;
        prior.density_mu = log((m^2)/sqrt(v+m^2));
        prior.density_sigma = sqrt(log(v/(m^2)+1));
        
        m=prior.T_mean;
        v=prior.T_std.^2;
        prior.T_mu = log((m^2)/sqrt(v+m^2));
        prior.T_sigma = sqrt(log(v/(m^2)+1));
        
        m=prior.dmax_mean;
        v=prior.dmax_std.^2;
        prior.dmax_mu = log((m^2)/sqrt(v+m^2));
        prior.dmax_sigma = sqrt(log(v/(m^2)+1));
        
        m=prior.pex_mean;
        v=prior.pex_std.^2;
        prior.pex_mu = log((m^2)/sqrt(v+m^2));
        prior.pex_sigma = sqrt(log(v/(m^2)+1));
        
        m=prior.sd_mean;
        v=prior.sd_std.^2;
        prior.sd_mu = log((m^2)/sqrt(v+m^2));
        prior.sd_sigma = sqrt(log(v/(m^2)+1));
    end
end


end