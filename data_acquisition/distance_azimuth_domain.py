#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Apr 05 2022
@author: Sylvain Brisson sylvain.brisson@ens.fr

Domain obspy object allowing to select the stations based on distances bounds and azimuth bounds.
"""

# Plotting routines
import matplotlib.pyplot as plt
from obspy.imaging.beachball import beach
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Maths
import numpy as np
from numpy import meshgrid, pi,sin,cos,arccos,arctan

# Obspy inherited class
from obspy.clients.fdsn.mass_downloader import Domain

# Obspy distance and azimuth computation routines
from obspy.geodetics.base import gps2dist_azimuth,locations2degrees,degrees2kilometers
import geopy
from geopy.distance import geodesic

# to have better defined great circles
# https://stackoverflow.com/questions/40270990/cartopy-higher-resolution-for-great-circle-distance-line
class LowerThresholdOthographic(ccrs.Orthographic):
    @property
    def threshold(self):
        return 1e3

class DistanceAzimuthDomain(Domain):
    """Define a surface on a sphere bounded by distances to a center point (an event) and azimuth (angles from the center with reference to a pivot).
    """
    
    def __init__(self, lat_center, lon_center, lat_pivot, lon_pivot, dmin, dmax, azmin, azmax):
        Domain.__init__(self)
        self.lat_center = lat_center
        self.lon_center = lon_center
        self.lat_pivot = lat_pivot
        self.lon_pivot = lon_pivot
        self.dmin = dmin 
        self.dmax = dmax 
        self.azmin = azmin 
        self.azmax = azmax

    def get_query_parameters(self):
        return {
        "latitude" : self.lat_center,
        "longitude": self.lon_center,
        "minradius": self.dmin,
        "maxradius": self.dmax,
        }

    def is_in_domain(self, lat, lon):
        
        # compute distances between event - pivot - station
        Δes = locations2degrees(self.lat_center, self.lon_center, lat, lon)
        Δep = locations2degrees(self.lat_center, self.lon_center, self.lat_pivot, self.lon_pivot)
        Δps = locations2degrees(self.lat_pivot, self.lon_pivot, lat, lon)
        
        # compute azimuth
        az = arccos((cos(Δps*pi/180.) - cos(Δep*pi/180.)*cos(Δes*pi/180.)) / (sin(Δep*pi/180.)*sin(Δes*pi/180.)))*180./pi
                
        return (self.azmin <= az and az <= self.azmax and self.dmin <= Δes and Δes <= self.dmax)
    
    def __repr__(self):
        out = ">> distance-azimuth bounded domain\n"
        out += f"- centered upon ({self.lat_center:.1f}°lat,{self.lon_center:.1f}°lon)\n"
        out += f"- distance in [{self.dmin:.1f},{self.dmax:.1f}]°\n"
        out += f"- azimuth in  [{self.azmin:.1f},{self.azmax:.1f}]°, with respect to the pivot point ({self.lat_pivot:.1f}°lat,{self.lon_pivot:.1f}°lon)"
        return out
    
    def plot(self):
            
        _,az_ep,_ = gps2dist_azimuth(self.lat_center,self.lon_center,self.lat_pivot, self.lon_pivot)
        origin = geopy.Point(self.lat_center,self.lon_center)
        
        fig = plt.figure()
            
        proj = LowerThresholdOthographic(
            central_latitude    = self.lat_pivot,
            central_longitude   = self.lon_pivot
        )   

        ax = fig.add_subplot(1, 1, 1, projection=proj)
        ax.stock_img()
        
        ax.plot(self.lon_center, self.lat_center, ".r",  transform = ccrs.PlateCarree())
        ax.plot(self.lon_pivot,  self.lat_pivot,  ".g",  transform = ccrs.PlateCarree())

        lats = []
        lons = []

        for dist in (self.dmin, self.dmax):
            for az in (self.azmin, self.azmax):
            
                az_np = az_ep + az
                dist_km = degrees2kilometers(dist)
                point = geodesic(kilometers=dist_km).destination(origin, bearing=az_np)
                lats.append(point.latitude)
                lons.append(point.longitude)
        
        lons[2],lons[3] = lons[3],lons[2]
        lats[2],lats[3] = lats[3],lats[2]
        lons.append(lons[0])
        lats.append(lats[0])

        ax.plot(lons, lats, "b", transform=ccrs.Geodetic(), lw=1)
        ax.plot(lons, lats, ".b", transform=ccrs.PlateCarree())


        # and plot great circles between source and domain boundaries
        for lat,lon in zip(lats[:2],lons[:2]):
            ax.plot([lon,self.lon_center], [lat,self.lat_center], "k",ls="--", transform=ccrs.Geodetic(), lw=1, alpha=0.5)
        ax.plot([self.lon_pivot,self.lon_center], [self.lat_pivot,self.lat_center], "k", transform=ccrs.Geodetic(), lw=1, alpha=0.5)

        return ax
            
    
class DistanceNorthAzimuthDomain(Domain):
    """Define a surface on a sphere bounded by distances to a center point (an event) and azimuth (angles from the center with reference to the North pole).
    """
    
    def __init__(self, lat, lon, dmin, dmax, azmin, azmax):
        Domain.__init__(self)
        self.lat = lat
        self.lon = lon
        self.dmin = dmin 
        self.dmax = dmax 
        self.azmin = azmin 
        self.azmax = azmax

    def get_query_parameters(self):
        return {
        "latitude" : self.lat,
        "longitude": self.lon,
        "minradius": self.dmin,
        "maxradius": self.dmax,
        }

    def is_in_domain(self, lat, lon):

        _,az,_ = gps2dist_azimuth(self.lat, self.lon, lat, lon)
        dist = locations2degrees(self.lat, self.lon, lat, lon)   
        
        return (self.azmin <= az and az <= self.azmax and self.dmin <= dist and dist <= self.dmax)
    
    def __repr__(self):
        
        return f"Domain centered upon {self.lat:.1f}°lat {self.lon:.1f}°lon, distance in [{self.dmin:.1f},{self.dmax:.1f}]°, azimuth (with respect to the North, in clock wise order) [{self.azmin:.1f},{self.azmax:.1f}]°"
    


    
    
    
# Routines to verify that the domain is working, not to been used
    
def latlon2xyz(lat,lon):
    """
    Convert (lat,lon) position into (x,y,z) on the unit sphere
    """
    theta = (90 - lat)*np.pi/180
    phi = (90 + lon)*np.pi/180
        
    x = np.sin(theta)*np.sin(phi)
    y = np.cos(theta)
    z = np.sin(theta)*np.cos(phi)
    
    return x,y,z
    
def set_axes_equal(ax: plt.Axes):
    """Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)
    
def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])


if __name__ == "__main__":
    
    # Testing
    
    from mpl_toolkits.mplot3d import Axes3D
    
    lat_center, lon_center = 10.0,60.0
    dmin,dmax = 20.0, 70.0
    
    # test 1
    # lat_pivot, lon_pivot = 90.0, 0.0
    # azmin = 20.0
    # azmax = 40.0
    # domain = DistanceNorthAzimuthDomain(lat_center, lon_center, dmin, dmax, azmin, azmax)
    
    # test 2
    lat_pivot, lon_pivot = 65.0, -20.0
    azmin = -20.0
    azmax = 20.0
    domain = DistanceAzimuthDomain(
        lat_center  = lat_center,
        lon_center  = lon_center,
        lat_pivot   = lat_pivot,
        lon_pivot   = lon_pivot,
        dmin        = dmin,
        dmax        = dmax,
        azmin       = azmin,
        azmax       = azmax
    )
    
    
    domain.plot()
    plt.show()
    exit()
    
    lon_1d,lat_1d = np.linspace(-180,180.0,100), np.linspace(-90.0, 90.0, 100)
    
    lat_mesh, lon_mesh = np.meshgrid(lat_1d, lon_1d)
        
    # plot 3D
    fig, ax = plt.subplots(1, 1, figsize = (8,8),subplot_kw={'projection':'3d'})
    plt.axis('off')
    
    # compute X,Y,Z position of points and color
    
    X_in,Y_in,Z_in = [],[],[]
    X_out,Y_out,Z_out = [],[],[]
    for lat in lat_1d:
        for lon in lon_1d:
            x,y,z = latlon2xyz(lat,lon)
            if domain.is_in_domain(lat,lon):
                X_in.append(x)
                Y_in.append(y)
                Z_in.append(z)
            else:
                X_out.append(x)
                Y_out.append(y)
                Z_out.append(z)
 
    X_in,Y_in,Z_in = map(np.array, [X_in,Y_in,Z_in])
    X_out,Y_out,Z_out = map(np.array, [X_out,Y_out,Z_out])
            
    points_in, = ax.plot([],[],[],'r.', markersize=5, alpha=0.9)
    points_out, = ax.plot([],[],[],'k.', markersize=5, alpha=0.9)
        
    azimuth, elev = -90, 85
    ax.view_init(elev, azimuth )
    
    # routine to plot only visible points
    def plot_visible(azimuth, elev):
        #transform viewing angle to normal vector in data coordinates
        a = azimuth*np.pi/180. -np.pi
        e = elev*np.pi/180. - np.pi/2.
        A = [ np.sin(e) * np.cos(a),np.sin(e) * np.sin(a),np.cos(e)]  
        # concatenate coordinates
        B = np.c_[X_in,Y_in,Z_in]
        # calculate dot product 
        # the points where this is positive are to be shown
        cond = (np.dot(B,A) >= 0)
        # filter points by the above condition
        x_c = X_in[cond]
        y_c = Y_in[cond]
        z_c = Z_in[cond]
        # set the new data points
        points_in.set_data(x_c, y_c)
        points_in.set_3d_properties(z_c, zdir="z")
        # idem for points out
        B = np.c_[X_out,Y_out,Z_out]
        # calculate dot product 
        # the points where this is positive are to be shown
        cond = (np.dot(B,A) >= 0)
        # filter points by the above condition
        x_c = X_out[cond]
        y_c = Y_out[cond]
        z_c = Z_out[cond]
        # set the new data points
        points_out.set_data(x_c, y_c)
        points_out.set_3d_properties(z_c, zdir="z")
        fig.canvas.draw_idle()
        
    plot_visible(azimuth, elev)

    def rotate(event):
        if event.inaxes == ax:
            plot_visible(ax.azim, ax.elev)
    fig.canvas.mpl_connect('motion_notify_event', rotate)
    
    # plotting source position
    x,y,z = latlon2xyz(lat_center,lon_center)
    ax.plot(x.tolist(),y.tolist(),z.tolist(),"*b", ms=10, label = "Center")
    
    x,y,z = latlon2xyz(lat_pivot, lon_pivot)
    ax.plot(x.tolist(),y.tolist(),z.tolist(),"ob", ms=10, label = "Pivot")
    
    ax.set_box_aspect([1,1,1])
    set_axes_equal(ax)
    ax.legend()
    
    plt.show()
    
    
