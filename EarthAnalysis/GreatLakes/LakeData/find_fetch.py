import rasterio
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from skimage import measure
import scipy.ndimage as ndimage    


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
    #for x,y in zip(xx,yy):
    #    if 0 <= x < depth.shape[1] and 0 <= y < depth.shape[0]:
    #        print(x)
    #        print(y)
    #        mask[y,x] = True


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


def plot_lake(img,metadata,transform,buoy_loc,shoreline,islands):
    # FUNCTION WILL PLOT THE LAKE, THE MAIN SHORELINE, AND BUOY LOCATION

    # Plot the first band of the multiband image
    plt.figure(figsize=(8, 8))
    img_plot = plt.imshow(img[0],
            extent=[transform[2], transform[2] + transform[0] * img.shape[2],transform[5] + transform[4] * img.shape[1], transform[5]],
            cmap='viridis') 
    plt.scatter(buoy_loc[1],buoy_loc[0],color='red', marker='o', s=100) 
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
    
def calculate_fetch(buoy,depth_binary,direction): # < ------------------------------------------------------------ PROBLEM HERE
    # FUNCTION WILL FIND THE FETCH FROM THE BUOY TO THE NEAREST SHORELINE IN SPECIFIED DIRECTION
    
    y,x = buoy[1],buoy[0]
    dy,dx = direction
    distance = 0
    
    maxWhile = 1e7
    whileCount = 0
    
    while 0 <= y < depth_binary.shape[0] and 0 <= x < depth_binary.shape[1] and whileCount <= maxWhile:
        distance += 1
        y += dy
        x += dx
        if depth_binary[y, x] == 0:
            return distance
        if whileCount + 1 == maxWhile:
            print('max number of while loop reached')
            distance = None
            return distance
        else:
            whileCount = whileCount + 1

def main():
    
    
    file_name = 'LS.tiff'  # lake superior bathy data from noaa
    
    # Station 45005 (East Superior) https://www.ndbc.noaa.gov/station_history.php?station=45004
    lat = 47.585
    lon = -86.585
    buoy_location = [lat,lon]

    LS,LSm,LSt = extract_tiff(file_name)

    depth = LS[0] # multiband for some reason? only look at first band
    depth = (depth < -10).astype(int) # zero out shallow depths to make shoreline cleaner
    
    main_shore = extract_main_shoreline(depth,LSt)
    
    dd = zero_outside_basin(depth,main_shore[2],LSt)
    ii = find_islands(dd,LSt)
    
    #plot_lake(LS,LSm,LSt,buoy_location,main_shore,ii)
    
    wind_dir = (1, 0)
    x_dist = calculate_fetch(buoy_location,dd,wind_dir)
    print(x_dist)

if __name__ == '__main__':
    main()
