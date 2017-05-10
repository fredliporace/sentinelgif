#!/usr/bin/python
# -*- coding: utf-8 -*-

import click

import os
import sys
import re
import urllib2
import numpy as np
import requests
import shutil
import uuid
import json
import math
import time
import pprint
from PIL import Image, ImageFont, ImageDraw
from skimage import exposure
from osgeo import gdal, ogr, osr
from shapely.wkt import loads

from matplotlib import cm

def query_builder(lat=None, lon=None, start_date=None, end_date=None,
                  cloud_min=None, cloud_max=None):
    """ Builds the proper search syntax (query) for SAT API """

    query = []
    and_string = ''
    search_string = ''

    if start_date and end_date:
        query.append(date_range_builder(start_date, end_date))
    elif start_date:
        query.append(date_range_builder(start_date, '2100-01-01'))
    elif end_date:
        query.append(date_range_builder('2009-01-01', end_date))

    if cloud_min and cloud_max:
        query.append(cloud_cover_prct_range_builder(cloud_min, cloud_max))
    elif cloud_min:
        query.append(cloud_cover_prct_range_builder(cloud_min, '100'))
    elif cloud_max:
        query.append(cloud_cover_prct_range_builder('-1', cloud_max))

    if lat and lon:
        query.append(lat_lon_builder(lat, lon))

    if query:
        and_string = '&'.join(map(str, query))

    search_string = and_string

    return search_string

def date_range_builder(start='2013-02-11', end=None):
    """
    Builds date range query
    Accepts start and end date in this format YYYY-MM-DD
    """
    if not end:
        end = time.strftime('%Y-%m-%d')

    return 'date_from=%s&date_to=%s' % (start, end)

def cloud_cover_prct_range_builder(min=0, max=100):
    """
    Builds cloud cover percentage range query
    Accepts bottom and top range in float, e.g. 1.00
    """
    return 'cloud_from=%s&cloud_to=%s' % (min, max)

def lat_lon_builder(lat=0, lon=0):
    """ Builds lat and lon query """

    return( 'intersects={"type":"Point","coordinates":[%f,%f]}' % (lon, lat) )

def search(quer, limit=200):
    """ Call sat api and return s2 scenes"""
    req = '%s?satellite_name=sentinel-2&%s&limit=%s' % \
    (sat_api_url, quer, limit)
    print req
    r = requests.get(req)

    r_dict = json.loads(r.text)
    result = {}

    if 'error' in r_dict:
        result['status'] = u'error'
        result['code'] = r_dict['error']['code']
        result['message'] = r_dict['error']['message']

    elif 'meta' in r_dict:
        result['status'] = u'SUCCESS'
        result['total'] = r_dict['meta']['found']
        result['results'] = [{'sceneID': i['scene_id'],
                              'path': '{:03d}'.format(ord(i['grid_square'][0])),
                              'row': '{:03d}'.format(ord(i['grid_square'][1])),
                              'date': i['date'],
                              # Assuming the points order in geometry
                              # is always the same, it would be better to
                              # check lower/upper left/right from coordinates
                              # values
                              'lowerLeftCornerLatitude':i['tile_geometry']['coordinates'][0][1][1],
                              'lowerLeftCornerLongitude':i['tile_geometry']['coordinates'][0][1][0],
                              'lowerRightCornerLatitude':i['tile_geometry']['coordinates'][0][2][1],
                              'lowerRightCornerLongitude':i['tile_geometry']['coordinates'][0][2][0],
                              'upperRightCornerLatitude':i['tile_geometry']['coordinates'][0][3][1],
                              'upperRightCornerLongitude':i['tile_geometry']['coordinates'][0][3][0],
                              'upperLeftCornerLatitude':i['tile_geometry']['coordinates'][0][4][1],
                              'upperLeftCornerLongitude':i['tile_geometry']['coordinates'][0][4][0],
                              '2':i['download_links']['aws_s3'][1],
                              '3':i['download_links']['aws_s3'][2],
                              '4':i['download_links']['aws_s3'][3],
                              '8':i['download_links']['aws_s3'][7],
                              'metadata':i['product_meta_link'],
                              'processing_level':i['processing_level'],
                              'cloud': i['cloud_coverage']}
                             for i in r_dict['results']]

    return result
##########################   

##############
sat_api_url = 'https://api.developmentseed.org/satellites'
##############

@click.group()
def cli():
    pass

@cli.command()
@click.option('--lat', type=float, default=None,
    help='Latitude of the query, between 90 and -90.')
@click.option('--lon', type=float, default=None,
    help='Longitude of the query, between 180 and -180.')
@click.option(
    '--cloud', type=float, default=20.,
    help='Maximum cloud percentage (%) allowed.')
@click.option(
    '--start_date', type=str, default='2015-01-01',
    help='Start date of the query in the format YYYY-MM-DD.')
@click.option(
    '--end_date', type=str, default=time.strftime('%Y-%m-%d'),
    help='End date of the query in the format YYYY-MM-DD.')
@click.option(
    '--buffer', type=int, default=10000,
    help='Buffer size around lat/lon point for image creation.')
@click.option(
    '--taskid', type=str, default=str(uuid.uuid1()),
    help='UUID of task.')
@click.option(
    '--ndvi', is_flag=True,
    help='Create NDVI animation instead of RGB')
@click.option(
    '--path', type=click.Path(exists=True), default='.',
    help='Set the path where the file will be saved.')

def worker(lat, lon, cloud, start_date, end_date, buffer, taskid, ndvi, path):
    """ Create animated GIF from landsat 8 data"""

    #Test 
    #lat lon has to be defined, path_row not supported
    if (not lat) | (not lon):
        print "No defined lat-lon"
        sys.exit(1)

    #Query Scenes
    print
    print "Building Landsat-API request"
    landsat_query = query_builder(lat=lat, lon=lon, start_date=start_date, end_date=end_date, cloud_max=cloud)
    print landsat_query
    print "Searching Landsat 8 images"
    candidate_scenes = search(landsat_query)
    #pprint.pprint(candidate_scenes)

    if not candidate_scenes.has_key('results'):
        print "Landsat-API Querry returned with 'Not Found message'"
        sys.exit(1)
        
    im2process = candidate_scenes['results']
    all_ids = [i['sceneID'] for i in im2process]
    
    print '{} Landsat scene found'.format(len(all_ids))
    print 'landsat ids: {}'.format(", ".join(all_ids))

    #Construct AOI  (square in WebMercator)
    wgs = osr.SpatialReference()  
    wgs.ImportFromEPSG(4326)
     
    wmerc = osr.SpatialReference()  
    wmerc.ImportFromEPSG(3857)
    wgsTowm = osr.CoordinateTransformation(wgs, wmerc)    
    wmTowgs = osr.CoordinateTransformation(wmerc, wgs)    
    
    #Create AOI - 10km buffer square (WebMercator) around point 
    pt = ogr.Geometry(ogr.wkbPoint)
    pt.AddPoint(lon, lat)
    
    pt.Transform(wgsTowm)
    shPt = loads(pt.ExportToWkt())
    polB = shPt.buffer(buffer, cap_style=3)
    
    aoi = ogr.CreateGeometryFromWkt(polB.wkt)
    aoi.Transform(wmTowgs) #Transform AOI in WGS84
    pt = pol = polB = None
    
    print "Excluding Landsat 8 scenes not covering the Entire AOI and not TOA if NDVI is requested"
    proc_images = []
    for ii in range(len(im2process)):
        imgMeta = im2process[ii]

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(imgMeta['lowerLeftCornerLongitude'], imgMeta['lowerLeftCornerLatitude'])
        ring.AddPoint(imgMeta['upperLeftCornerLongitude'], imgMeta['upperLeftCornerLatitude'])
        ring.AddPoint(imgMeta['upperRightCornerLongitude'], imgMeta['upperRightCornerLatitude'])
        ring.AddPoint(imgMeta['lowerRightCornerLongitude'], imgMeta['lowerRightCornerLatitude'])
        ring.AddPoint(imgMeta['lowerLeftCornerLongitude'], imgMeta['lowerLeftCornerLatitude'])
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        if aoi.Within(poly):

            if( not ndvi or imgMeta['processing_level'] == 'Level-1C'):
                proc_images.append(imgMeta)
        
        ring = poly = None

    if len(proc_images) == 0:
        print 'No Image found covering the AOI - change buffer size or change lat-lon'  
    else:

        #Check Only if ROW is the same (same date) 
        all_pr = ['{:03d},{:03d}'.format(int(i['path']),int(i['row'])) for i in proc_images]
        all_row = [i['row'] for i in proc_images]
        if len(list(set(all_row))) > 1:
            print '''AOI covering more than one Row : 
            Please choose one of the following: {}
            Using --path_row option'''.format(' | '.join(list(set(all_pr))))
            sys.exit(1)

        workdir = os.path.join(path, taskid)
        if not os.path.exists(workdir):
            os.makedirs(workdir, 0775)

        font = ImageFont.load_default().font
        
        l8_images = []
        date_array = []        
        for ii in range(len(proc_images)):
             
            im = proc_images[ii]
            print 'Processing Landsat image {}'.format(im['sceneID'])
            
            out_im = os.path.join(workdir, '{}.tif'.format(im['date']))
            print out_im

            try:
                WRSPath = im['path']
                WRSRow = im['row']

                # Get B03 scene geographic metadata
                bqa = '/vsicurl/{addr_name}'.format(addr_name=im['3'])
                src_ds = gdal.Open(bqa, gdal.GA_ReadOnly)
                geoT = src_ds.GetGeoTransform()
                proj = src_ds.GetProjection()
                src_ds = None
          
                imSpatialRef = osr.SpatialReference()
                imSpatialRef.ImportFromWkt(proj)
          
                aoiSpatialRef = osr.SpatialReference()
                aoiSpatialRef.ImportFromEPSG(4326)
                coordTransform = osr.CoordinateTransformation(aoiSpatialRef, imSpatialRef)
      
                aoi.Transform(coordTransform) # reproject the aoi in UTM
                aoi_bounds = aoi.GetEnvelope()
          
                x_off = int((aoi_bounds[0] - geoT[0]) / geoT[1])
                y_off = int((aoi_bounds[3] - geoT[3]) / geoT[5])             
                x_size = int(((aoi_bounds[0] - geoT[3]) / geoT[5]) - ((aoi_bounds[1] - geoT[3]) / geoT[5]))
                y_size = int(((aoi_bounds[2] - geoT[3]) / geoT[5]) - ((aoi_bounds[3] - geoT[3]) / geoT[5]))
                #print 'x_off: %f, y_off: %f, x_size: %f, y_size: %f' % (x_off, y_off, x_size, y_size)

                #Create RGB file
                ngeo = list(geoT)
                ngeo[0] = aoi_bounds[0]
                ngeo[3] = aoi_bounds[3]
            
                if ndvi:
                    # Here we use band5 and band4 as it would be for LS8,
                    # but the image is obtained from equivalent Sentinel-2
                    # bands (band8 and band4, respectively)
                    band5_address = '/vsicurl/{0}'.format(im['8'])
                    awsim5 = gdal.Open(band5_address, gdal.GA_ReadOnly)
                    arr5 = awsim5.GetRasterBand(1).ReadAsArray(x_off, y_off, x_size, y_size) 
                    # No need conversion since Level-1C is already TOA
                    # reflectance
                    
                    band4_address = '/vsicurl/{0}'.format(im['4'])
                    awsim4 = gdal.Open(band4_address, gdal.GA_ReadOnly)
                    arr4 = awsim4.GetRasterBand(1).ReadAsArray(x_off, y_off, x_size, y_size) 

                    ratio = np.where( arr5*arr4 > 0, np.nan_to_num((arr5 - arr4) / (arr5 + arr4)), 0)
                    awsim4 = awsim5 = arr4 = arr5 = None

                    #Use winter colormap (http://matplotlib.org/examples/color/colormaps_reference.html)
                    img = Image.fromarray(np.uint8(cm.winter((ratio + 1.) / 2.) * 255)).convert('RGB')
                    ratio = None
                    draw = ImageDraw.Draw(img)
                    xs,ys = draw.textsize(im['date'],  font=font)
                    draw.rectangle([ (5, 5), (xs+15, ys+15) ], fill=(255,255,255))
                    draw.text((10, 10), im['date'], (0,0,0), font=font)
                    out_jpg = out_im.replace('.tif','.jpg')
                    img.save(out_jpg)           

                else:
                    driver = gdal.GetDriverByName("GTiff")
                    dst_ds = driver.Create(out_im, x_size, y_size, 3, gdal.GDT_Byte)
                    dst_ds.SetGeoTransform(tuple(ngeo))
                    dst_ds.SetProjection(proj)
                    rgb = [4,3,2]
                    for b in range(len(rgb)):
                        band_address = '/vsicurl/{0}'.format(im[str(rgb[b])])
                        awsim = gdal.Open(band_address, gdal.GA_ReadOnly)
                        arr = awsim.GetRasterBand(1).ReadAsArray(x_off, y_off, x_size, y_size)                                       
                        p2, p98 = np.percentile(arr[arr > 0], (2, 98))
                        dst_ds.GetRasterBand(b+1).WriteArray(np.where(arr > 0, exposure.rescale_intensity(arr, in_range=(p2, p98), out_range=(1,255)), 0))
                        dst_ds.GetRasterBand(b+1).SetNoDataValue(0)
                        awsim = arr = None
                    dst_ds = None    # save, close
              
                    img = Image.open(out_im)
                    draw = ImageDraw.Draw(img)
                    xs,ys = draw.textsize(im['date'],  font=font)
                    draw.rectangle([ (5, 5), (xs+15, ys+15) ], fill=(255,255,255))
                    draw.text((10, 10), im['date'], (0,0,0), font=font)
                    out_jpg = out_im.replace('.tif','.jpg')
                    img.save(out_jpg)                
                    os.remove(out_im)
                    
                date_array.append(im['date'])
                l8_images.append(out_jpg)
            except:
                #print "Unexpected error", sys.exc_info()[0]
                #raise
                print 'Failed to process Landsat image {}'.format(im['sceneID'])

        if len(date_array) > 0:
            #Sort image by date and rename with number
            sorted_index = np.argsort(date_array)
            l8sort = [l8_images[i] for i in sorted_index]
            for i in range(len(l8sort)):
                os.rename(l8sort[i], os.path.join(workdir, '{:05d}.jpg'.format(i)))
        
            #This part can be replace in pure python
            #Create GIF
            gif_file = os.path.join(path, "%s.gif" % taskid)
            inJpg = os.path.join(workdir, "*.jpg")
            os.system('convert -delay 30 -depth 8 -layers optimize -quality 80 -loop 0 {0} {1}'.format(inJpg, gif_file))
        
        shutil.rmtree(workdir)

if __name__ == '__main__':
    worker()
