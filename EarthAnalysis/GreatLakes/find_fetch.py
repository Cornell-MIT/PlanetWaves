#!/usr/bin/env python3

import csv
import logging
import os
import sys
from math import atan2, cos, radians, sin, sqrt


import colorlog
import numpy as np
import rasterio
import scipy.ndimage as ndimage
from alive_progress import alive_bar
from matplotlib import pyplot as plt
from osgeo import (  # pip install --global-option=build_ext --global-option="-I/usr/include/gdal" GDAL==`gdal-config --version`
    gdal, osr)
from scipy import io
from skimage import measure

'''
find_fetch.py:

    Find the fetch to a buoy station given the bathymetry
    of a lake (.tiff) for any arbitrary wind direction.

Author: Una Schneck (schneck.una@gmail.com)

Last modified: 04/18/2024
'''

#####################################################################################################
#####################################################################################################
# LOGGING FUNCTIONS

# Define log_level as a global variable
log_level = logging.INFO


def make_log() -> logging.Logger:
    """
    MAKE LOGGER
    """

    global log_level

    logger = logging.getLogger(__name__)

    logger.setLevel(log_level)

    if not logger.handlers:
        stdout = colorlog.StreamHandler(stream=sys.stdout)
        fmt = colorlog.ColoredFormatter(
            "%(name)s: %(white)s%(asctime)s%(reset)s | %(log_color)s%(levelname)s%(reset)s | %(blue)s%(filename)s:%(lineno)s%(reset)s  >>> %(log_color)s%(message)s%(reset)s"
        )
        stdout.setFormatter(fmt)
        logger.addHandler(stdout)

    return logger


def not_used(func):
    """
    wrapper for function to not use for whatever reason
    """
    def wrapper():
        raise RuntimeError(f"Function {func.__name__} should not be used.")
    return wrapper

#####################################################################################################
#####################################################################################################
# CLASSES


class Buoy:
    """
    Class containing information on Buoy of interest:
    (1) Latitude [deg]
    (2) Longitude [deg]
    (3) Name of lake [string]
    """
    def __init__(self, lat, lon, lakename):
        self.lat = lat
        self.lon = lon
        self.lakename = lakename

#####################################################################################################
#####################################################################################################
# HELPER FUNCTIONS


def degrees_to_vector(degrees):
    """
    FUNCTION TO CONVERT DEGREES INTO A NORMALIZED VECTOR
    """
    myradians = radians(degrees)

    # Calculate the x and y components of the vector
    x = cos(myradians)
    y = sin(myradians)

    # Normalize the vector
    vector_length = sqrt(x**2 + y**2)
    normalized_x = x / vector_length
    normalized_y = y / vector_length

    return normalized_x, normalized_y


def calculate_distance(lat1, lon1, lat2, lon2):
    """
    FUNCTION TO CALCULATE DISTANCE BETWEEN TWO POINTS ON A GLOBE USING HAVERSINE
    """
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
        radius_earth = 6371.009  # Radius of the Earth in kilometers

        distance = radius_earth * c  # Distance in kilometers
        distance = distance*1000  # DIstance in meters

        mylog.info(f"{calculate_distance.__name__} ran successfully")
        return distance
    except Exception as e:
        mylog.critical(f"Error in calculate_distance(): {e}")
        return None


#####################################################################################################
#####################################################################################################
# TIFF FILE FUNCTIONS

def extract_tiff(tiff_file):
    """
    FUNCTION TO EXTRACT THE DATA, METADATA, AND TRANSFORM INFORMATION FROM TIFF FILE
    """
    mylog = make_log()

    try:
        # Open the TIFF file
        with rasterio.open(tiff_file) as src:

            # Read the image and its metadata
            img = src.read()
            metadata = src.meta
            transform = src.transform

            mylog.info(f"{extract_tiff.__name__} ran successfully")

            return img, metadata, transform
    except Exception as e:
        mylog.critical(f"Error in extract_tiff(): {e}")
        return None, None, None


def extract_transformation(file_name):
    """
    FUNCTION TO EXTRACT THE TRANFORMATION MATRIX FROM TIFF FILE
    """
    mylog = make_log()

    try:
        dataset = gdal.Open(file_name)
        if dataset is None:
            mylog.critical(f"Unable to open the file --> {file_name} <-- . Double-check spelling.")
            return None
        else:
            mylog.info(f"File <{file_name}> imported successfully")

        try:
            geotransform = dataset.GetGeoTransform()
            mylog.info(f"{extract_transformation.__name__} ran successfully")
            return geotransform
        except Exception as e:
            mylog.critical(f"Error while getting geotransform: {e}")
            return None

    except Exception as e:
        mylog.critical(f"Exception occurred: {e}")
        return None


def pixel_to_geo(x, y, geotransform):
    """
    FUNCTION TO CONVERT PIXEL COORDINATE TO GEOGRAPHIC COORDINATE
    """
    mylog = make_log()

    try:
        lon = geotransform[0] + x * geotransform[1] + y * geotransform[2]
        lat = geotransform[3] + x * geotransform[4] + y * geotransform[5]
        mylog.info(f"{pixel_to_geo.__name__} ran successfully")
        return lat, lon
    except Exception as e:
        mylog.critical(f"Error in pixel_to_geo(): {e}")
        return None, None


def geo_to_pixel(lat, lon, geotransform):
    """
    FUNCTION TO CONVERT GEOGRAPHIC COORDINATE TO PIXEL COORDINATE
    """
    mylog = make_log()

    try:
        inv_geotransform = gdal.InvGeoTransform(geotransform)
        x, y = gdal.ApplyGeoTransform(inv_geotransform, lon, lat)
        x_int = int(x + 0.5)  # Round to the nearest integer
        y_int = int(y + 0.5)  # Round to the nearest integer
        mylog.info(f"{geo_to_pixel.__name__} ran successfully")
        return x_int, y_int
    except Exception as e:
        mylog.critical(f"Error in geo_to_pixel(): {e}")
        return None, None


def return_grid_cell_loc_of_buoy(MyBuoy, map_meta):
    """
    Convert latitude and longitude to array indices
    """
    col = int((MyBuoy.lon - map_meta['transform'].c) / map_meta['transform'].a)
    row = int((MyBuoy.lat - map_meta['transform'].f) / map_meta['transform'].e)

    return col, row


def return_grid_cell_size(map_meta):
    """
    find the actual size of the grid cell of the array given the transformation
    """
    scale_x = map_meta['transform'].a
    scale_y = map_meta['transform'].e

    grid_cell_width = 1 / scale_x  # has units of the coordinate system (e.g. degrees)
    grid_cell_height = 1 / scale_y  # has units of the coordinate system (e.g. degrees)

    return grid_cell_width, grid_cell_height

#####################################################################################################
#####################################################################################################
# SHORELINE AND DEPTH FUNCTIONS


def clean_depth(depth_tiff, threshold):
    """
     FUNCTION TO EXTRACT DEPTH THRESHOLDED TO SOME VALUE TO MAKE
     THE ARRAY AND SUBSEQUENT SHORELINE EXTRACTION CLEANER
    """
    mylog = make_log()

    try:
        depth = depth_tiff[0]  # multiband for some reason? only look at first band
        depth = (depth < threshold).astype(int)  # zero out shallow depths to make shoreline cleaner
        mylog.info(f"{clean_depth.__name__} ran successfully")
        return depth
    except Exception as e:
        mylog.critical(f"Error in clean_depth(): {e}")
        return None


def extract_all_shoreline(depth, trs):
    """
    FUNCTION TO EXTRACT ALL THE BOUNDARIES (AKA SHORELINES) IN THE IMAGE
    """
    mylog = make_log()

    ploton = False

    try:
        shore = measure.find_contours(depth, 0.5)

        if ploton:
            for c in shore:
                lon, lat = trs * (c[:, 1], c[:, 0])
                plt.plot(lon, lat, linewidth=2, color='black')
        mylog.info(f"{extract_all_shoreline.__name__} ran successfully")
        return shore
    except Exception as e:
        mylog.critical(f"Error in extract_all_shoreline(): {e}")
        return None


def extract_main_shoreline(depth, trs):
    """
    FUNCTION TO EXTRACT JUST THE LONGEST BOUNDARY
    (AKA SHORELINE) IN THE IMAGE
    """
    mylog = make_log()

    ploton = False

    try:
        shore = extract_all_shoreline(depth, trs)

        main_shore = max(shore, key=len)  # longest continuos shorelin
        lon1, lat1 = trs*(main_shore[:, 1], main_shore[:, 0])
        # main_shore2 = sorted(shore,key=len)[-2]
        # lon2,lat2 = trs*(main_shore2[:,1],main_shore2[:,0])
        if ploton:
            plt.plot(lon1, lat1, linewidth=2, color='red')
            # plt.plot(lon2,lat2,linewidth=2,color='green')
            plt.title('shorelines')
            plt.show()
        mylog.info(f"{extract_main_shoreline.__name__} ran successfully")
        return lon1, lat1, main_shore
    except Exception as e:
        mylog.critical(f"Error in extract_main_shoreline(): {e}")
        return None, None, None


def zero_outside_basin(depth, shore):
    """
    FUNCTION TO ZERO OUT ALL THE DEPTH INFORMATION OUTSIDE THE MAIN SHORELINE
    """
    mylog = make_log()

    ploton = False

    try:
        xx = shore[:, 1]
        yy = shore[:, 0]

        mask = np.zeros_like(depth, dtype=bool)
        mask[np.round(shore[:, 0]).astype('int'), np.round(shore[:, 1]).astype('int')] = 1
        mask = ndimage.binary_fill_holes(mask)

        depth_inside = depth
        depth_inside[~mask] = 0

        if ploton:
            fig, axes = plt.subplots(1, 2, figsize=(10, 5))

            axes[0].imshow(depth, cmap='binary')
            axes[0].plot(xx, yy, linewidth=2, color='red')
            axes[0].set_title('All Depth')

            axes[1].imshow(depth_inside, cmap='binary')
            axes[1].set_title('Depth Inside Main Basin')

            plt.tight_layout()
            plt.show()
        mylog.info(f"{zero_outside_basin.__name__} ran successfully")
        return depth_inside
    except Exception as e:
        mylog.critical(f"Error in zero_outside_basin(): {e}")
        return None


def find_islands(depth, trs):
    """
    FUNCTION TO FIND ALL THE BOUNDARIES (AKA SHORELINES FOR ISLANDS) WITHIN THE MAIN BASIN
    """
    mylog = make_log()

    ploton = False

    try:
        mshore = extract_main_shoreline(depth, trs)
        mshore = mshore[2]
        mlon, mlat = trs * (mshore[:, 1], mshore[:, 0])

        islands = measure.find_contours(depth)
        islands = islands[1:]
        if ploton:
            plt.plot(mlon, mlat, linewidth=2, color='black')
            for i in islands:
                lon, lat = trs * (i[:, 1], i[:, 0])
                plt.plot(lon, lat, linewidth=2, color='red')
            plt.title('islands within lake')
            plt.show()
        mylog.info(f"{find_islands.__name__} ran successfully")
        return islands
    except Exception as e:
        mylog.critical(f"Error in {find_islands.__name__}: {e}")
        return None

#####################################################################################################
#####################################################################################################
# FETCH FUNCTIONS


@not_used
def grid_calculate_fetch(MyBuoy, depth_binary, geotransform, direction):
    """
    FUNCTION WILL FIND THE FETCH FROM THE BUOY TO THE NEAREST SHORELINE IN SPECIFIED DIRECTION
    NOTE: not used because only works for cardinal directions
    """
    mylog = make_log()

    mylog.warning(f"{grid_calculate_fetch.__name__} only works for cardinal direction. Need true ray casting for fetch calculation")

    try:
        xb, yb = geo_to_pixel(MyBuoy.lat, MyBuoy.lon, geotransform)
        dy, dx = direction

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
            lat1, lon1 = pixel_to_geo(x0, y0, geotransform)
            lat2, lon2 = pixel_to_geo(x, y, geotransform)

            step_size = calculate_distance(lat1, lon1, lat2, lon2)
            distance += step_size

            if depth_binary[y, x] == 0:
                return distance, lat2, lon2

            if whileCount + 1 == maxWhile:
                mylog.error('max number of while loop reached')
                return None, None, None
            else:
                whileCount = whileCount + 1
    except Exception as e:
        mylog.critical(f"Error in {grid_calculate_fetch.__name__}: {e}")
        return None


def dda_ray_casting(start_pt, direction_deg, depth):
    """
    FUNCTION TO FIND FETCH USING RAY CASTING
    DDA Algorithm for ray casting
    https://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm)
    """

    mylog = make_log()

    direction_deg = direction_deg % 360

    # invert so that shores are 1s and liquid is 0
    depth = np.where(depth == 1, 0, 1)

    x_start = start_pt[0]
    y_start = start_pt[1]

    max_fetch = sqrt(depth.shape[0]**2 + depth.shape[1]**2)

    x_end = x_start + max_fetch * cos(radians(direction_deg))
    y_end = y_start + max_fetch * sin(radians(direction_deg))

    dx = x_end - x_start
    dy = y_end - y_start

    steps = int(abs(dx)) if abs(dx) > abs(dy) else int(abs(dy))

    x_increment = dx / steps
    y_increment = dy / steps

    x = x_start
    y = y_start

    grid_collision_point = None

    for _ in range(steps):
        x += x_increment
        y += y_increment

        x_rounded = round(x)
        y_rounded = round(y)

        # Check if the line has hit a boundary in the array
        if not (0 <= x_rounded < depth.shape[1]) or not (0 <= y_rounded < depth.shape[0]):
            mylog.warning(f"{dda_ray_casting.__name__}:  no collision point found")
            return None

        if depth[y_rounded][x_rounded]:
            grid_collision_point = (x_rounded, y_rounded)
            mylog.info(f"{dda_ray_casting.__name__}:  collision point found")
            mylog.info(f"{dda_ray_casting.__name__} ran successfully")
            break

    return grid_collision_point


def loop_thru_direction(depth_file, MyBuoy):
    """
    PERFORM CALCULATIONS OF FETCH FOR ALL DIRECTIONS
    SPECIFIED
    """

    LS, LSm, LSt = extract_tiff(depth_file)
    geo_trs = extract_transformation(depth_file)
    depth = clean_depth(LS, -10)
    main_shore = extract_main_shoreline(depth, LSt)
    dd = zero_outside_basin(depth, main_shore[2])
    ii = find_islands(dd, LSt)
    bx, by = geo_to_pixel(MyBuoy.lat, MyBuoy.lon, geo_trs)

    return bx, by, dd, geo_trs, LS, LSm, LSt, main_shore, ii


def make_table(bx, by, winds, dd, geo_trs, MyBuoy):
    """
    MAKE A DICTIONARY OF FETCHES AND THE INTERSECTION WITH SHORELINE
    """
    wind_fetch = {}

    with alive_bar(len(winds), bar='smooth') as bar:
        for wind_dir in winds:
            dda_pt = dda_ray_casting([bx, by], wind_dir, dd)
            flat, flon = pixel_to_geo(dda_pt[0], dda_pt[1], geo_trs)
            fetch_dist = calculate_distance(MyBuoy.lat, MyBuoy.lon, flat, flon)
            wind_fetch[wind_dir] = fetch_dist
            bar()

    return wind_fetch, flat, flon, fetch_dist


def write_to_csv(my_csv_filename, wind_fetch) -> None:
    """
    WRITE DICTIONARY OF FETCHES TO CSV
    """
    mylog = make_log()

    with open(my_csv_filename, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Wind [deg]', 'Fetch [m]'])  # Header row
        for key, value in wind_fetch.items():
            writer.writerow([key, value])

    mylog.info(f'Fetch table saved to: {csv_filename}')


def save_bathymetry_as_mat(LS) -> None:
    """
    SAVE THE CLEANED DEPTHS TO A MAT FILE
    """
    mylog = make_log()
    mat_filename = 'LakeSuperior_cleaned.mat'
    io.savemat(mat_filename, {'LS': LS})
    mylog.info(f'Depth profile saved to {mat_filename}')

#####################################################################################################
#####################################################################################################
# PLOTTING RESULTS FUNCTIONS


def plot_lake(img, transform, MyBuoy, shoreline, islands, fetch_loc, fetch_dist, wind_dir):
    """
    FUNCTION WILL PLOT THE LAKE, THE MAIN SHORELINE, AND BUOY LOCATION
    """

    # Plot the first band of the multiband image
    plt.figure(figsize=(8, 8))
    img_plot = plt.imshow(img[0],
                          extent=[transform[2], transform[2] + transform[0] * img.shape[2],
                                  transform[5] + transform[4] * img.shape[1], transform[5]],
                          cmap='viridis')

    plt.scatter(MyBuoy.lon, MyBuoy.lat, color='red', marker='o', s=100)
    plt.scatter(fetch_loc[0], fetch_loc[1], color='magenta', marker='o', s=100)
    plt.annotate('', xy=[MyBuoy.lon, MyBuoy.lat], xytext=fetch_loc, arrowprops=dict(arrowstyle='->', color='black'))

    plt.plot(shoreline[0], shoreline[1], linewidth=2, color='black')
    for i in islands:
        lon, lat = transform * (i[:, 1], i[:, 0])
        plt.plot(lon, lat, linewidth=2, color='black')

    # plot details
    cbar = plt.colorbar(img_plot, orientation='horizontal')
    cbar.set_label('Depth [m]')  # Set your colorbar label here
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f'{MyBuoy.lakename}\nWind Direction: {wind_dir[-1]} deg, Fetch : {fetch_dist} km')
    plt.grid(True)
    # plt.savefig(f'Lake_Superior_{wind_dir}deg.jpg', format='jpg')
    plt.show()

#####################################################################################################
#####################################################################################################
# MAIN FUNCTION

# MAIN FUNCTION:
#   INPUTS:
#       (1) BUOY LOCATION
#       (2) BATHYMETRY OF GIVEN LAKE
#       (3) WIND DIRECTION
#   OUTPUTS:
#       (1) FETCH (aka DISTANCE)
#       (2) INTERSECTION POINT WITH NEAREST SHORELINE


def main(depth_file_name, csv_filename) -> None:

    # >> python find_fetch.py --info will print out the description of the overall code
    if "--info" in sys.argv:
        print("find_fetch.py:\n"
              "\n"
              "    Find the fetch to a buoy station given the bathymetry\n"
              "    of a lake (.tiff) for any arbitrary wind direction.\n"
              "\n"
              "Author: Una Schneck (schneck.una@gmail.com)\n"
              "\n"
              "Last modified: 04/18/2024\n")

    mylog = make_log()
##########################################################################
# MAIN INPUTS
    last_dir = 1
    winds = list(range(0, last_dir, 1))
    station = 45004
    # winds = list(range(0, 359, 1))

###########################################################################
# RUN CALCULATION

    if station == 45004:
        # Station 45004 (East Superior) https://www.ndbc.noaa.gov/station_history.php?station=45004
        mylog.warning('Note to self: only classified buoy is 45004')
        blat: float = 47.585
        blon: float = -86.585
        lakename: str = "Lake Superior"
        deepwater_buoy = Buoy(blat, blon, lakename)
    else:
        mylog.critical(f'Station is unknown. Only 45004 currently supported')
        return

    mylog.info(f"Calculating fetch for {len(winds)} directions at station {station} in {lakename}")

    here = os.getcwd()

    depth_file = os.path.join(here, "LakeData", depth_file_name)

    if os.path.exists(depth_file):
        mylog.info(f"{depth_file_name} successfully found")
    else:
        mylog.critical(f"{depth_file_name} not found. Check file location.")
        return

    bx, by, dd, geo_trs, LS, LSm, LSt, main_shore, ii = loop_thru_direction(depth_file, deepwater_buoy)
    wind_fetch, flat, flon, fetch_dist = make_table(bx, by, winds, dd, geo_trs, deepwater_buoy)
    # csv_filename = os.path.join(here, csv_filename)
    # write_to_csv(csv_filename,wind_fetch)

    plot_lake(LS, LSt, deepwater_buoy, main_shore, ii, [flon, flat], round(fetch_dist / 1000, 2), winds)


# TO RUN: >> python find_fetch.py
if __name__ == '__main__':

    # bathymetry profile tiff
    depth_file_name = 'LS.tiff'
    # name of csv to save the fetches list for each direction
    csv_filename = 'WindFetchLS.csv'

    main(depth_file_name, csv_filename)