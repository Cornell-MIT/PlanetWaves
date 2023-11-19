#!/usr/bin/env python3

import sys
import rasterio
import pandas as pd
import numpy as np
import math
import logging
import colorlog
import scipy.ndimage as ndimage

from skimage import measure
from matplotlib import pyplot as plt
from osgeo import gdal, osr
from math import radians, sin, cos, sqrt, atan2

'''
find_fetch.py: Find the fetch to a buoy station given the bathymetry of a lake (.tiff) for any arbitrary direction

Author: Una Schneck (schneck.una@gmail.com)

'''

# Define log_level as a global variable
log_level = logging.DEBUG

def setup_custom_logger():
    
    global log_level
    
    logger = logging.getLogger("LOGGING")

    logger.setLevel(log_level)

    if not logger.handlers:
        stdout = colorlog.StreamHandler(stream=sys.stdout)
        fmt = colorlog.ColoredFormatter(
            "%(name)s: %(white)s%(asctime)s%(reset)s | %(log_color)s%(levelname)s%(reset)s | %(blue)s%(filename)s:%(lineno)s%(reset)s | %(process)d >>> %(log_color)s%(message)s%(reset)s"
        )
        stdout.setFormatter(fmt)
        logger.addHandler(stdout)

    return logger

def make_log():
    
    global log_level
    mylog = setup_custom_logger()
    
    return mylog
    
def extract_transformation(file_name):
    # FUNCTION TO EXTRACT THE TRANFORMATION MATRIX FROM TIFF FILE
    
    mylog = make_log()
     
    try:
        dataset = gdal.Open(file_name)
        if dataset is None:
            mylog.critical(f"Unable to open the file --> {file_name} <-- . Double-check spelling.")
            return None
        else:
            mylog.info(f"File {file_name} imported successfully")
            
        try:
            geotransform = dataset.GetGeoTransform()
            mylog.info("extract_transformation() sucessful ")
            return geotransform
        except Exception as e:
            mylog.critical(f"Error while getting geotransform: {e}")
            return None
        
    except Exception as e:
        mylog.critical(f"Exception occurred: {e}")
        return None


def pixel_to_geo(x, y, geotransform):
    # FUNCTION TO CONVERT PIXEL COORDINATE TO GEOGRAPHIC COORDINATE
    
    mylog = make_log()

    try:
        lon = geotransform[0] + x * geotransform[1] + y * geotransform[2]
        lat = geotransform[3] + x * geotransform[4] + y * geotransform[5]
        mylog.info("pixel_to_geo() sucessful ")
        return lat, lon
    except Exception as e:
        mylog.critical(f"Error in pixel_to_geo(): {e}")
        return None,None

def geo_to_pixel(lat, lon, geotransform):
    # FUNCTION TO CONVERT GEOGRAPHIC COORDINATE TO PIXEL COORDINATE
    
    mylog = make_log()
    
    try:
        inv_geotransform = gdal.InvGeoTransform(geotransform)
        x, y = gdal.ApplyGeoTransform(inv_geotransform, lon, lat)
        x_int = int(x + 0.5)  # Round to the nearest integer
        y_int = int(y + 0.5)  # Round to the nearest integer
        mylog.info("geo_to_pixel successful")
        return x_int, y_int
    except Exception as e:
        mylog.critical(f"Error in geo_to_pixel(): {e}")
        return None, None

def calculate_distance(lat1, lon1, lat2, lon2):
    # Convert latitude and longitude from degrees to radians
    
    mylog = make_log()
    
    try:
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
        
        mylog.info("calculate_distance() sucessful")
        return distance 
    except Exception as e:
        mylog.critical(f"Error in calculate_distance(): {e}")
        return None

def extract_tiff(tiff_file):
    # FUNCTION WILL EXTRACT THE DATA, METADATA, AND TRANSFORM INFORMATION FROM TIFF FILE 
    
    mylog = make_log()
    
    try:
        # Open the TIFF file
        with rasterio.open(tiff_file) as src:
            
            # Read the image and its metadata
            img = src.read()
            metadata = src.meta
            transform = src.transform
            
            mylog.info("extract_tiff() succesful")
            
            return img,metadata,transform
    except Exception as e:
        mylog.critical(f"Error in calculate_distance(): {e}")
        return None,None,None

def clean_depth(depth_tiff, threshold):
    # EXTRACT DEPTH THRESHOLDED TO SOME VALUE TO MAKE THE ARRAY AND SUBSEQUENT SHORELINE EXTRACTION CLEANER
    
    mylog = make_log()
    
    try:
        depth = depth_tiff[0] # multiband for some reason? only look at first band
        depth = (depth < threshold).astype(int) # zero out shallow depths to make shoreline cleaner        
        mylog.info("clean_depth() sucessful")
        return depth
    except Exception as e:
        mylog.critical(f"Error in clean_depth(): {e}")
        return None
    
def extract_all_shoreline(depth,trs):
    # FUNCTION WILL EXTRACT ALL THE BOUNDARIES (AKA SHORELINES) IN THE IMAGE
    
    mylog = make_log()
    
    ploton = False
    
    try:
        shore = measure.find_contours(depth,0.5)
        
        if ploton:
            for c in shore:
                lon, lat = trs * (c[:, 1], c[:, 0])
                plt.plot(lon,lat,linewidth=2,color='black')

        return shore
    except Exception as e:
        mylog.critical(f"Error in extract_all_shoreline(): {e}")
        return None
        
def extract_main_shoreline(depth,trs):
    # FUNCTION WILL EXTRACT JUST THE LONGEST BOUNDARY (AKA SHORELINE) IN THE IMAGE
    
    mylog = make_log()
    
    ploton = False
    
    try:
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
        mylog.info("extract_main_shoreline() sucessful")
        return lon1,lat1,main_shore
    except Exception as e:
        mylog.critical(f"Error in extract_main_shoreline(): {e}")
        return None,None,None
        
def zero_outside_basin(depth,shore,trs): 
    # FUNCTION WILL ZERO OUT ALL THE DEPTH INFORMATION OUTSIDE THE MAIN SHORELINE 
    
    mylog = make_log()
    
    ploton = False

    try:
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
        mylog.info("zero_outside_basin() sucessful")
        return depth_inside
    except Exception as e:
        mylog.critical(f"Error in zero_outside_basin(): {e}")
        return None

def find_islands(depth,trs):
    # FUNCTION WILL FIND ALL THE BOUNDARIES (AKA SHORELINES FOR ISLANDS) WITHIN THE MAIN BASIN
    
    mylog = make_log()
    
    ploton = False
    
    try:
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
        mylog.info("find_island() successful")
        return islands
    except Exception as e:
        mylog.critical(f"Error in find_island(): {e}")
        return None
    

def degrees_to_vector(degrees):
    # Convert degrees to radians
    
    radians = math.radians(degrees)
    
    # Calculate the x and y components of the vector
    x = math.cos(radians)
    y = math.sin(radians)
    
    # Normalize the vector
    vector_length = math.sqrt(x**2 + y**2)
    normalized_x = x / vector_length
    normalized_y = y / vector_length
    
    return normalized_x, normalized_y
    
def grid_calculate_fetch(buoy_lonlat,depth_binary,geotransform,direction): # 
    # FUNCTION WILL FIND THE FETCH FROM THE BUOY TO THE NEAREST SHORELINE IN SPECIFIED DIRECTION
    # only works for cardinal directions
    
    my_log = make_log()
    
    mylog.debug("grid_calculate_fetch() only works for cardinal direction. Need true ray casting for fetch calculation")
    
    try:
        xb,yb = geo_to_pixel(buoy_lonlat[1],buoy_lonlat[0],geotransform)
        dy,dx = direction
        #dy,dx = degree_to_vector(direction)

        distance = 0
        
        maxWhile = 1e9
        whileCount = 0
        
        x = xb
        y = yb
            
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
                mylog.error('max number of while loop reached')
                return None,None,None
            else:
                whileCount = whileCount + 1
    except Exception as e:
        mylog.critical(f"Error in grid_calculate_fetch(): {e}")
        return None

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
        plt.plot(lon,lat,linewidth=2,color='black')

    # plot details
    cbar = plt.colorbar(img_plot,orientation='horizontal')
    cbar.set_label('Depth [m]')  # Set your colorbar label here
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Lake Superior')
    plt.grid(True)
    plt.show()
    
def find_fetch(buoy,wind_dir):
    # CALCULATE THE FETCH AND INTERSECTION POINT WITH THE SHORELINE FROM THE BUOY FOR A GIVEN WIND DIRECTION
    
    mylog = make_log()
    
    ploton = True
    
    blat = buoy[0]#47.585
    blon = buoy[1]#-86.585
        
    depth_file = 'LS.tiff'  # lake superior bathy data from noaa

    # RUN CALCULATIONS TO FIND FETCH FOR BUOY GIVEN WIND DIRECTION AND MAKE PLOT
    LS,LSm,LSt = extract_tiff(depth_file)
    geo_trs = extract_transformation(depth_file)
    depth = clean_depth(LS,-10)
    main_shore = extract_main_shoreline(depth,LSt)
    dd = zero_outside_basin(depth,main_shore[2],LSt)
    ii = find_islands(dd,LSt)
    #f_dist,flat,flon = grid_calculate_fetch([blon,blat],dd,geo_trs,wind_dir)
    #print(f"{f_dist/1000} kilometers")
    mylog.error("Cannot calculate fetch yet. Ray casting is not working")
    #f_dist,flat,flon = ray_casting() # <---- not workiung
    
    if ploton:
        plot_lake(LS,LSm,LSt,[blon,blat],main_shore,ii,[flon,flat])
        

    return f_dist

def ray_casting(vRayStart,direction_deg,depth):
    
    # DDA Algorithm based on 
    # https://github.com/OneLoneCoder/Javidx9/blob/master/PixelGameEngine/SmallerProjects/OneLoneCoder_PGE_RayCastDDA.cpp
    
    mylog = make_log()
    
    mylog.debug("DDA ray casting is not working yet. Need to debug here.")
    
    vMapCheck = vRayStart
    vecMap = np.where(depth == 1, 0, 1) # invert so that shores are 1s and liquid is 0
    vMapSize = depth.shape
    
    vRayDir = degrees_to_vector(direction_deg)
    
    if vRayDir[0] != 0 and vRayDir[1] != 0:
        vRayUnitStepSize = [ sqrt(1 + (vRayDir[1] / vRayDir[0]) * (vRayDir[1] / vRayDir[0])), sqrt(1 + (vRayDir[0] / vRayDir[1]) * (vRayDir[0] / vRayDir[1])) ];
    else:
        vRayUnitStepSize = vRayDir
            
    vMapCheck = vRayStart
    vRayLength1D = [0, 0]
    vStep = [0, 0]

    # Establish Starting Conditions
    if vRayDir[0] < 0:
        vStep[0] = -1
        vRayLength1D[0] = (vRayStart[0] - float(vMapCheck[0])) * vRayUnitStepSize[0]
    else:
        vStep[0] = 1
        vRayLength1D[0] = (float(vMapCheck[0] + 1) - vRayStart[0]) * vRayUnitStepSize[0]

    if vRayDir[1] < 0:
        vStep[1] = -1
        vRayLength1D[1] = (vRayStart[1] - float(vMapCheck[1])) * vRayUnitStepSize[1]
    else:
        vStep[1] = 1
        vRayLength1D[1] = (float(vMapCheck[1] + 1) - vRayStart[1]) * vRayUnitStepSize[1]

    # Perform "Walk" until collision or range check
    bTileFound = False
    fMaxDistance = max(vMapSize)
    fDistance = 0.
    
    mylog.debug(vMapCheck[1])
    mylog.debug(vMapSize[0])
    mylog.debug(vMapCheck[0])
    mylog.debug(vecMap[vMapCheck[1] * vMapSize[0] + vMapCheck[0]])
    
    return 
    
    while not bTileFound and fDistance < fMaxDistance:
        # Walk along shortest path
        if vRayLength1D[0] < vRayLength1D[1]:
            vMapCheck[0] += vStep[0]
            fDistance = vRayLength1D[0]
            vRayLength1D[0] += vRayUnitStepSize[0]
        else:
            vMapCheck[1] += vStep[1]
            fDistance = vRayLength1D[1]
            vRayLength1D[1] += vRayUnitStepSize[1]

        # Test tile at new test point
        if (
            0 <= vMapCheck[0] < vMapSize[0]
            and 0 <= vMapCheck[1] < vMapSize[1]
            and vecMap[vMapCheck[1] * vMapSize[0] + vMapCheck[0]] == 1
        ):
            bTileFound = True

    # Calculate intersection location
    vIntersection = None
    if bTileFound:
        vIntersection = [vRayStart[0] + vRayDir[0] * fDistance, vRayStart[1] + vRayDir[1] * fDistance]


    return vIntersection, fDistance
    
def main():

    station = 45004 # Station 45004 (East Superior) https://www.ndbc.noaa.gov/station_history.php?station=45004

    if station == 45004:
        buoy_lat = 47.585
        buoy_lon = -86.585

    #wind_direction = (1,0) # < -------------------------------------------- problem here! This is quantize so won't work in most directions
    #fetch_m = find_fetch([buoy_lat,buoy_lon],wind_direction)

    depth_file = 'LS.tiff' 
    LS,LSm,LSt = extract_tiff(depth_file)
    geo_trs = extract_transformation(depth_file)
    depth = clean_depth(LS,-10)
    main_shore = extract_main_shoreline(depth,LSt)
    dd = zero_outside_basin(depth,main_shore[2],LSt)
    
    
    wind_dir = 190
    bx,by = geo_to_pixel(buoy_lat, buoy_lon, geo_trs)
    ray_casting([bx,by],wind_dir,dd)


if __name__ == '__main__':
    main() 
    
    

    
    
    
