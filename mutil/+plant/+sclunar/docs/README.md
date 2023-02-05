## 3-DoF Spacecraft Dynamics in Cislunar Space

The equations of motion are represented in the Earth-centered inertial frame (ECI) defined by Earth's Mean Equator and Mean Equinox (MEME) at 12:00 Terrestrial Time on 1 January 2000, with the origin at the instantaneous Moon 
center. This frame is labelled as `J2000` in SPICE.

### Accessing Ephemeris via SPICE

[NAIF SPICE](https://naif.jpl.nasa.gov/naif/toolkit.html) toolkit for MATLAB, called `mice`, is necessary for querying ephemeris state of celestial bodies in specific coordinate frames.
Download `mice` from [here](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html) and include `mice\lib` and `mice\src\mice` to MATLAB path.

### Data

 - `pck00010.tpc` &ndash; Orientation and size/shape data for natural bodies.
 - `naif0011.tls.pc` &ndash; Leap second kernel.
 - `de421.bsp` &ndash; Position of planets and Moon between 1900 and 2050. See [this document](https://ipnpr.jpl.nasa.gov/progress_report/42-178/178C.pdf) for more details.

### Solar Radiation Pressure

SRP acceleration from canonball model:

$$
    a_{\scriptscriptstyle\text{SRP}} = -\frac{1}{2}\beta C_R \frac{GM_{S}}{\|r_{\scriptscriptstyle sc-S}\|_2^3}r_{\scriptscriptstyle sc-S}, 
$$

where 

$$
\begin{aligned}
    \beta ={} & \frac{\sigma^\star_{\scriptscriptstyle m/A}}{\sigma_{\scriptscriptstyle m/A}}\\
    \sigma^\star_{\scriptscriptstyle m/A} ={} & \frac{2\Phi L^2}{cGM_S}
\end{aligned}
$$

 - $r_{\scriptscriptstyle sc-S}$ is the vector pointing from spacecraft to Sun, 
 - $C_R = 1.15$ is the radiation pressure coefficient, 
 - $G$ is the universal gravitational constant,
 - $M_S$ is the mass of Sun,
 - $\Phi = 1360~\text{kg/s}^3$ is solar flux,
 - $L = 1.495978707 \times 10^8~\text{km}$ is the astronomical unit,
 - $c = 299792458~\text{m/s}$ is speed of light. 

For a spacecraft with perfectly reflective surface, the mass to area ratio, $\sigma_{\scriptscriptstyle  m/A}^\star = 1.53 \times 10^{-3}~\text{kg/m}^2$, is the critical value at which the solar radiation pressure acting on the spacecraft is equal to the gravitational force due to the Sun.

More details and assumptions about this model are provided in Section 2.3.2 of [this thesis](https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2021_Muralidharan.pdf). 