function default_val = get_default(name)
    % This funtion returns default settings of some parameters.

    % name: the name of the parameter.
    %
    % default_val: the default setting value of the parameter "name".


    switch name
        case {"pcm_eps"}
            default_val = 1e-8;
        case {"pcm_nu"}
            default_val = 0.9;
        case {"pcm_mu"}
            default_val = 0.4;
        case {"pcm_vrho"}
            default_val = 1.5;
        case {"ifprint"}
            default_val = "no";
        otherwise
            default_val = "Unknown parameter name";
    end
end