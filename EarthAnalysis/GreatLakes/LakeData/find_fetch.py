import rasterio
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from skimage import measure
import scipy.ndimage as ndimage
from osgeo import gdal, osr
from math import radians, sin, cos, sqrt, atan2
import math
   
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

def clean_depth(depth_tiff, threshold):
    
    depth = depth_tiff[0] # multiband for some reason? only look at first band
    depth = (depth < threshold).astype(int) # zero out shallow depths to make shoreline cleaner
    
    return depth
    
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
    
def ray_casting(vRayStart,direction_deg,depth):
    
    # DDA Algorithm based on https://github.com/OneLoneCoder/Javidx9/blob/master/PixelGameEngine/SmallerProjects/OneLoneCoder_PGE_RayCastDDA.cpp
    
    vMapCheck = vRayStart
    vecMap = depth
    vMapSize = np.array(depth).shape
    vRayDir = degrees_to_vector(direction_deg)
    
    vRayUnitStepSize = (
    math.sqrt(1 + (vRayDir[1] / vRayDir[0]) * (vRayDir[1] / vRayDir[0])),
    math.sqrt(1 + (vRayDir[0] / vRayDir[1]) * (vRayDir[0] / vRayDir[1]))
    )
    
    vStep = [0, 0]
    vRayLength1D = [0, 0]

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

    # Walk until a collision
    bTileFound = False
    fMaxDistance = 100000.0
    fDistance = 0.0
    
    while not bTileFound and fDistance < fMaxDistance:
        if vRayLength1D[0] < vRayLength1D[1]:
            vMapCheck[0] += vStep[0]
            fDistance = vRayLength1D[0]
            vRayLength1D[0] += vRayUnitStepSize[0]
        else:
            vMapCheck[1] += vStep[1]
            fDistance = vRayLength1D[1]
            vRayLength1D[1] += vRayUnitStepSize[1]

        if (
            0 <= vMapCheck[0] < vMapSize[0]
            and 0 <= vMapCheck[1] < vMapSize[1]
            and vecMap[vMapCheck[1] * vMapSize[0] + vMapCheck[0]] == 1
        ):
            bTileFound = True
    

    if bTileFound:
        vIntersection = vRayStart + vRayDir * fDistance
    else:
        vIntersection = None
        fDistance = None
        print('no intersection found')
        
    return vIntersection, fDistance
    
    
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
    f_dist,flat,flon = calculate_fetch([blon,blat],dd,geo_trs,wind_dir)
    print(f"{f_dist/1000} kilometers")
    
    if ploton:
        plot_lake(LS,LSm,LSt,[blon,blat],main_shore,ii,[flon,flat])
    
    return f_dist

if __name__ == '__main__':
    
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
    
    
    wind_dir = 0
    print(degrees_to_vector(wind_dir))
    bx,by = geo_to_pixel(buoy_lat, buoy_lon, geo_trs)
    ray_casting([bx,by],wind_dir,dd)
    
