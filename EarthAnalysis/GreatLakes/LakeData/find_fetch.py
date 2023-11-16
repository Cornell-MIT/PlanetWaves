import rasterio
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from skimage import measure

def extract_tiff(tiff_file):
    # Open the TIFF file
    with rasterio.open(tiff_file) as src:
        
        # Read the image and its metadata
        img = src.read()
        metadata = src.meta
        transform = src.transform
        return img,metadata,transform

def extract_all_shoreline(depth,trs):
    
    ploton = False

    shore = measure.find_contours(depth,0.5)

    if ploton:
        for c in shore:
            lon, lat = trs * (c[:, 1], c[:, 0])
            plt.plot(lon,lat,linewidth=2,color='black')

    return shore
    
def extract_main_shoreline(depth,trs):
    
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


def zero_outside_basin(depth,trs):
    
    ploton = True
    
    shore = extract_main_shoreline(depth,trs)
    shore = shore[2]
    xx = shore[:,1]
    yy = shore[:,0]
    
    mask = np.zeros_like(depth,dtype=bool)
    
    for x,y in zip(xx,yy):
        mask[y,x] = True


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

def find_islands(depth,main_shoreline,trs):
    
    #islands = measure.find_contours(main_shoreline)
    #plt.plot(main_shoreline[:,1],main_shoreline[:,0],linewidth=2,color='black')
    #for i in islands:
    #    plt.plot(i[:,1],i[:,0],linewidth=2,color='red')
    #plt.title('islands within lake')
    #plt.show()
    pass
    

def plot_lake(img,metadata,transform,buoy_loc,shoreline):
    # Plot the first band of the multiband image
    plt.figure(figsize=(8, 8))
    img_plot = plt.imshow(img[0],
            extent=[transform[2], transform[2] + transform[0] * img.shape[2],transform[5] + transform[4] * img.shape[1], transform[5]],
            cmap='viridis') 
    plt.scatter(buoy_loc[1],buoy_loc[0],color='red', marker='o', s=100) 
    plt.plot(shoreline[0],shoreline[1],linewidth=2,color='black')
    #plt.plot(shoreline[2],shoreline[3],linewidth=2,color='black')

    # plot details
    cbar = plt.colorbar(img_plot,orientation='horizontal')
    cbar.set_label('Depth [m]')  # Set your colorbar label here
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Lake Superior')
    plt.grid(True)
    plt.show()


def main():
    file_name = 'LS.tiff'  # lake superior data from noaa
    
    # Station 45005 (East Superior) https://www.ndbc.noaa.gov/station_history.php?station=45004
    lat = 47.585
    lon = -86.585
    buoy_location = [lat,lon]

    LS,LSm,LSt = extract_tiff(file_name)
    
    
    depth = LS[0]
    depth = (depth < -10).astype(int)
    
        #shoreline = extract_main_shoreline(depth,LSt)
    #plot_lake(LS,LSm,LSt,buoy_location,shoreline)
    zero_outside_basin(depth,LSt)
    #find_islands(shoreline[0:1],LSt)


if __name__ == '__main__':
    main()
