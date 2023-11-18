import rasterio
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from skimage import measure
import scipy.ndimage as ndimage
from osgeo import gdal, osr
from math import radians, sin, cos, sqrt, atan2


def degrees_to_direction(degrees):
    
    # Normalize degrees to be within 0 to 360 range
    degrees %= 360
    
    # range of angles for a finite eight neighbor grid 
    direction_ranges = {
        (337.5, 22.5): 'N',
        (22.5, 67.5): 'NE',
        (67.5, 112.5): 'E',
        (112.5, 157.5): 'SE',
        (157.5, 202.5): 'S',
        (202.5, 247.5): 'SW',
        (247.5, 292.5): 'W',
        (292.5, 337.5): 'NW'
    }
    
    # Check which range the degrees fall into
    for angle_range, direction in direction_ranges.items():
        if angle_range[0] <= degrees < angle_range[1]:
            return direction    
        else
            print('degree not valid')
            return None

def direction_to_step(direction):
    
    if direction == 'N':
    elif direction == 'NE':
    elif direction == 'E':
    elif direction == 'SE':
    elif direction == 'S':
    elif direction == 'SW':
    elif direction == 'W':
    elif direction == 'NW':
    

    
def extract_transformation(file_name):
    # FUNCTION TO EXTRACT 

    dataset = gdal.Open(file_name)
    
    if dataset is None:
        print("Unable to open the file")
        geotransform = None

    geotransform = dataset.GetGeoTransform()

    if geotransform is None:
        print("No geotransform found")

    return geotransform

def pixel_to_geo(x, y, geotransform):
    # FUNCTION TO CONVERT PIXEL COORDINATE TO GEOGRAPHIC COORDINATE
    lon = geotransform[0] + x * geotransform[1] + y * geotransform[2]
    lat = geotransform[3] + x * geotransform[4] + y * geotransform[5]
    return lat, lon

def geo_to_pixel(lat, lon, geotransform):
    # FUNCTION TO CONVERT GEOGRAPHIC COORDINATE TO PIXEL COORDINATE

    inv_geotransform = gdal.InvGeoTransform(geotransform)
    x, y = gdal.ApplyGeoTransform(inv_geotransform, lon, lat)
    x_int = int(x + 0.5)  # Round to the nearest integer
    y_int = int(y + 0.5)  # Round to the nearest integer
    return x_int, y_int

def calculate_distance(lat1, lon1, lat2, lon2):
    # Convert latitude and longitude from degrees to radians
    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    radius_earth = 6371  # Radius of the Earth in kilometers

    distance = radius_earth * c  # Distance in kilometers
    distance = distance*1000 # DIstance in meters
    
    return distance 

def extract_tiff(tiff_file):
    # FUNCTION WILL EXTRACT THE DATA, METADATA, AND TRANSFORM INFORMATION FROM TIFF FILE 
    
    # Open the TIFF file
    with rasterio.open(tiff_file) as src:
        
        # Read the image and its metadata
        img = src.read()
        metadata = src.meta
        transform = src.transform
        return img,metadata,transform

def extract_all_shoreline(depth,trs):
    # FUNCTION WILL EXTRACT ALL THE BOUNDARIES (AKA SHORELINES) IN THE IMAGE
    
    ploton = False

    shore = measure.find_contours(depth,0.5)
    

    if ploton:
        for c in shore:
            lon, lat = trs * (c[:, 1], c[:, 0])
            plt.plot(lon,lat,linewidth=2,color='black')

    return shore
    
def extract_main_shoreline(depth,trs):
    # FUNCTION WILL EXTRACT JUST THE LONGEST BOUNDARY (AKA SHORELINE) IN THE IMAGE
    
    ploton = False
    
    shore = extract_all_shoreline(depth,trs)
    
    main_shore = max(shore,key=len) # longest continuos shorelin
    lon1,lat1 = trs*(main_shore[:,1],main_shore[:,0])
    #main_shore2 = sorted(shore,key=len)[-2]
    #lon2,lat2 = trs*(main_shore2[:,1],main_shore2[:,0])
    if ploton:
        plt.plot(lon1,lat1,linewidth=2,color='red')
     #   plt.plot(lon2,lat2,linewidth=2,color='green')
        plt.title('shorelines')
        plt.show()
    
    return lon1,lat1,main_shore


def zero_outside_basin(depth,shore,trs): 
    # FUNCTION WILL ZERO OUT ALL THE DEPTH INFORMATION OUTSIDE THE MAIN SHORELINE 
    
    ploton = False
    

    xx = shore[:,1]
    yy = shore[:,0]
    
    
    mask = np.zeros_like(depth,dtype=bool)
    mask[np.round(shore[:, 0]).astype('int'), np.round(shore[:, 1]).astype('int')] = 1
    mask = ndimage.binary_fill_holes(mask)


    depth_inside = depth
    depth_inside[~mask] = 0

    if ploton:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        
        axes[0].imshow(depth,cmap='binary')
        axes[0].plot(xx,yy,linewidth=2,color='red')
        axes[0].set_title('All Depth')
        
        
        axes[1].imshow(depth_inside,cmap='binary')
        axes[1].set_title('Depth Inside Main Basin')
        
        plt.tight_layout()
        plt.show()

    return depth_inside

def find_islands(depth,trs):
    # FUNCTION WILL FIND ALL THE BOUNDARIES (AKA SHORELINES FOR ISLANDS) WITHIN THE MAIN BASIN
    
    ploton = False
    
    mshore = extract_main_shoreline(depth,trs)
    mshore = mshore[2]
    mlon,mlat = trs * (mshore[:, 1], mshore[:, 0])

    islands = measure.find_contours(depth)
    islands = islands[1:] 
    if ploton:
        plt.plot(mlon,mlat,linewidth=2,color='black')
        for i in islands:
            lon, lat = trs * (i[:, 1], i[:, 0])
            plt.plot(lon,lat,linewidth=2,color='red')
        plt.title('islands within lake')
        plt.show()
    
    return islands

def calculate_fetch(buoy_lonlat,depth_binary,geotransform,direction): # < ------------------------------------------------------------ PROBLEM HERE
    # FUNCTION WILL FIND THE FETCH FROM THE BUOY TO THE NEAREST SHORELINE IN SPECIFIED DIRECTION
        
    xb,yb = geo_to_pixel(buoy_lonlat[1],buoy_lonlat[0],geotransform)
    dy,dx = direction
    #dy,dx = degree_to_vector(direction)

    

    distance = 0
    
    maxWhile = 1e9
    whileCount = 0
    
    x = xb
    y = yb
    
    print(dy,dx)
    
    while 0 <= y < depth_binary.shape[0] and 0 <= x < depth_binary.shape[1] and whileCount <= maxWhile:
        
        x0 = x
        y0 = y
        y = y0 + dy
        x = x0 + dx
        lat1,lon1 = pixel_to_geo(x0,y0,geotransform)
        lat2,lon2 = pixel_to_geo(x,y,geotransform)
        
        step_size = calculate_distance(lat1, lon1, lat2, lon2)
        distance += step_size

        if depth_binary[y, x] == 0:
            return distance,lat2,lon2

        if whileCount + 1 == maxWhile:
            print('max number of while loop reached')
            return None,None,None
        else:
            whileCount = whileCount + 1

def plot_lake(img,metadata,transform,buoy_loc,shoreline,islands,fetch_loc):
    # FUNCTION WILL PLOT THE LAKE, THE MAIN SHORELINE, AND BUOY LOCATION

    # Plot the first band of the multiband image
    plt.figure(figsize=(8, 8))
    img_plot = plt.imshow(img[0],
            extent=[transform[2], transform[2] + transform[0] * img.shape[2],transform[5] + transform[4] * img.shape[1], transform[5]],
            cmap='viridis') 
    
    plt.scatter(buoy_loc[0],buoy_loc[1],color='red', marker='o', s=100) 
    plt.scatter(fetch_loc[0],fetch_loc[1],color='blue', marker='o', s=100)
    plt.annotate('', xy=buoy_loc, xytext=fetch_loc, arrowprops=dict(arrowstyle='->', color='black'))   
     
    plt.plot(shoreline[0],shoreline[1],linewidth=2,color='black')
    for i in islands:
        lon, lat = transform * (i[:, 1], i[:, 0])
        plt.plot(lon,lat,linewidth=2,color='blue')

    # plot details
    cbar = plt.colorbar(img_plot,orientation='horizontal')
    cbar.set_label('Depth [m]')  # Set your colorbar label here
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Lake Superior')
    plt.grid(True)
    plt.show()
    

def main():
    
    
    file_name = 'LS.tiff'  # lake superior bathy data from noaa

    LS,LSm,LSt = extract_tiff(file_name)
    geo_trs = extract_transformation(file_name)
    
    # Station 45005 (East Superior) https://www.ndbc.noaa.gov/station_history.php?station=45004
    lat = 47.585
    lon = -86.585
    buoy_location = [lon,lat]

    depth = LS[0] # multiband for some reason? only look at first band
    depth = (depth < -10).astype(int) # zero out shallow depths to make shoreline cleaner
    
    main_shore = extract_main_shoreline(depth,LSt)
        
    dd = zero_outside_basin(depth,main_shore[2],LSt)
    ii = find_islands(dd,LSt)
    

    #wind_dir = 45 # from North [deg]
    wind_dir = (0,1)
    f_dist,flat,flon = calculate_fetch(buoy_location,dd,geo_trs,wind_dir)
    print(f"{f_dist/1000} kilometers")
    
    fetch_pt = [flon,flat]
    
    plot_lake(LS,LSm,LSt,buoy_location,main_shore,ii,fetch_pt)
    
    
    

if __name__ == '__main__':
    main()
