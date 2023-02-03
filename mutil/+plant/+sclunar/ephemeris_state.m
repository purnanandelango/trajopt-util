% Function to find ephemeris state of a target wrt an observer
% Output:
%   ephm_state = position and velocity of target body wrt observer 
% Input:
%   t               = time in days since start_JD 
%   start_JD        = Julian Date at the start of the trajectory
%   target_body     = description of the target body; eg. 301 or 'MOON'
%   observer_body   = description of the observer body; eg. 399 or 'EARTH'
function ephm_state = ephemeris_state(t,start_JD,target_body,observer_body)
    
    % Calling this at each call to ephemeris_state will lead to considerable slowdown
    % cspice_furnsh( { 'naif0011.tls.pc',...
    %                  'de421.bsp',...
    %                  'pck00010.tpc' } );

    J2000           = 2451545.0;                    % [day]
    t_past_J2000    = (t+start_JD-J2000)*24*3600;   % time past Julian Date [s]
    ptarg           = mice_spkezr( num2str(target_body), t_past_J2000, 'J2000', 'none', num2str(observer_body));
    ephm_state      = [ptarg.state]';

    % This is important for unloading ephemeris files from memory
    % cspice_kclear
end