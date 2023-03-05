% Load or unload spice kernel files
function [] = ephem(flg)
    switch flg
        case 'load'
            
            path2data = "/Users/purnanandelango/Documents/trajopt-util/mutil/data/";
            % path2data = "C:\Users\elango\Documents\FreeMotion_v2\EPHEM\ephemeris_data\";

            cspice_furnsh( cellfun(@char,{ path2data + "naif0011.tls.pc",...
                                           path2data + "de421.bsp",...
                                           path2data + "pck00010.tpc" },'UniformOutput',false) );            
        case 'unload'
            cspice_kclear;
        otherwise
            error("Invalid flag.")
    end
end
