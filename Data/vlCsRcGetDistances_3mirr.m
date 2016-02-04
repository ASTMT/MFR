function dist = vlCsRcGetDistances_3mirr

global IM

% Values taken from Ritchey-Chretion Baseline Design (TMT.SEN.SPE.06.001.DRF02)

% Distance from M1 to M2 is 27.09375 meters
dist.M1M2 = 27.09375 * IM.UnitsPerMeter;

% M3 is 3.5 meters ABOVE the primary so we use a negative value.
dist.M1M3 = -3.50 * IM.UnitsPerMeter;

% Distance from M3 to focal plane is 20.0 meters
dist.M3Focus = 20.0 * IM.UnitsPerMeter;

% Distance from azimuth plane to elevation axis is 20 meters.
dist.AzEl = 20.0 * IM.UnitsPerMeter;
