#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 17:03:40 2019

@author: leonidas

# Some Code from here: https://git.earthdata.nasa.gov/projects/LPDUR/repos/nasa_viirs_surfacereflectance/browse
# Πληροφορίες σχετικά με τα Black Marble για τα bits, scale values, hdf dataset: https://viirsland.gsfc.nasa.gov/PDF/VIIRS_BlackMarble_UserGuide.pdf
"""

import h5py, os
import numpy as np
from osgeo import gdal, gdal_array
from pyresample import image, geometry
import osgeo.osr
import glob
import os.path
from scipy import stats
import datetime
import numpy.ma as ma
import math
from lunar_irrad_DNB import mt2009
import calendar
import sys

def reproject(dataset, inAreaDefinition, outAreaDefinition):
   
    msg_con_nn = image.ImageContainerNearest(dataset, inAreaDefinition, radius_of_influence=2000)
    area_con_nn = msg_con_nn.resample(outAreaDefinition)
    return(area_con_nn)


def geotiff(outDir,outName,resolution, fill_value,ImageContainer, array):
    
    result_data_nn=array

    outputName = os.path.normpath(os.path.join(outDir, '{}_2100.tif'.format(outName))) # Generate output filename
    print(outputName)
    nRow, nCol = result_data_nn.shape                                           # Define rows/cols from array size
    dataType = gdal_array.NumericTypeCodeToGDALTypeCode(result_data_nn.dtype)   #result_data_nn.dtype # Define output data type or np.int16 for smaller files
    if dataType==3:
        fill_value=32767

    driver = gdal.GetDriverByName('GTiff')                              
    
    outFile = driver.Create(outputName, nCol, nRow, 1, dataType)                        # Specify parameters of the GeoTIFF
    band = outFile.GetRasterBand(1)                                                     # Get band 1
    band.WriteArray(result_data_nn)                                                     # Write data array to band 1
    band.FlushCache   
                                                         # Export data
    band.SetNoDataValue(fill_value)
    uul=ImageContainer.geo_def.upper_left_extent
                                              # Set fill value
    geoInfo = (uul[0], resolution, 0, uul[1], 0, -resolution)         # Define geotransform parameters
    outFile.SetGeoTransform(geoInfo)
    
    srs = osgeo.osr.SpatialReference()
    srs.ImportFromEPSG(2100)                                                    # Set Geotransform
    outFile.SetProjection(srs.ExportToWkt())                                                              # Set projection
    outFile = None          
    


def BRDF_data(BRDF_file,greek_def):
    # https://lpdaac.usgs.gov/products/vnp43ma4v001/
    with h5py.File(BRDF_file, 'r') as hdf:
    
        fileMetadata = hdf['HDFEOS INFORMATION']['StructMetadata.0'][()].split() # Read file metadata
        fileMetadata = [m.decode('utf-8') for m in fileMetadata]                 # Clean up file metadata
        
        ulc = [i for i in fileMetadata if 'UpperLeftPointMtrs' in i][0]    # Search file metadata for the upper left corner of the file
        ulcLon = float(ulc.split('=(')[-1].replace(')', '').split(',')[0]) # Parse metadata string for upper left corner lon value
        ulcLat = float(ulc.split('=(')[-1].replace(')', '').split(',')[1]) # Parse metadata string for upper left corner lat value
        
        lrc = [i for i in fileMetadata if 'LowerRightMtrs' in i][0]  
        lrcLon = float(lrc.split('=(')[-1].replace(')', '').split(',')[0]) # Parse metadata string for upper left corner lon value
        lrcLat = float(lrc.split('=(')[-1].replace(')', '').split(',')[1]) # Parse metadata string for upper left corner lat value
                       
        GROUP_DNB_SDR = 'HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/'
        M4            = hdf[GROUP_DNB_SDR]['Nadir_Reflectance_M4'][...] #Albedo_BSA_M4
        M5            = hdf[GROUP_DNB_SDR]['Nadir_Reflectance_M5'][...] 
        M7            = hdf[GROUP_DNB_SDR]['Nadir_Reflectance_M7'][...] 
        
        scale_factor=1.0E-4 #  0.0001
        M_avg = np.nanmean( [ M4, M5, M7 ], axis=0 )*scale_factor
        fill_value= 32767
        M_avg[M4==fill_value]=fill_value
        M_avg[M5==fill_value]=fill_value
        M_avg[M7==fill_value]=fill_value
        sinu_def = geometry.AreaDefinition('sinu', 'MODIS Sinusoidal', 'sinu',
                                       {'a': '6371007.181', 'b': '6371007.181',
                                        'lat_0': '0', 'lat_ts': '0',
                                        'lon_0': '0', 'proj': 'sinu', 'units':'m'},
                                        1200, 1200,
                                        [ulcLon, lrcLat,
                                         lrcLon,ulcLat])   
    
        M_avg_brdf=reproject(M_avg, sinu_def, greek_def)
        brdf= M_avg_brdf.image_data
        return(brdf)
        
        
        
def LunarCorrection(theta, brdf, Em):
  #from Millers formula
    SRF_INTEG = 0.32560294 #VIIRS spectral response function
    #print(Em*SRF_INTEG)
    Lm = ((Em / math.pi) * brdf * np.cos((theta * math.pi)/180)) #Angles should be in radians, not degrees.
    Lm = Lm *  100 * SRF_INTEG   # *  100 to convert to nW per sq centimeter per stradian
    return (Lm)


def modeUTC_time(timearray, ymd):
    utc_mode=stats.mode(np.round(timearray,4).ravel()).mode[0]
    hours = str(int(utc_mode)).zfill(2)
    minutes = (utc_mode*60) % 60
    #seconds = (utc_mode*3600) % 60
    modetime="%s%s%02d" % (ymd,hours, minutes)
    return(modetime)
    
def julian_days(julian_day):
    julian_day_str= julian_day[5:]
    julian_day_integer = int(julian_day[5:])
    year=int(julian_day[1:5])
    if julian_day_integer==1:
        if calendar.isleap(year-1):
            previous_julian_day_integer = 366
        else:
            previous_julian_day_integer = 365
        previous_day_year=year-1
    else:
        previous_julian_day_integer=julian_day_integer-1
        previous_day_year=year
    result={"year":str(year), 
            "julian_day":julian_day,
            "julian_day_str":julian_day_str,
            "julian_day_integer":julian_day_integer,
            "previous_julian_day_str":"%03d" % (previous_julian_day_integer), 
            "previous_day_year":str(previous_day_year)
            }
    return(result)
    
    
def correctDNB(mt_time_format,DNB, LZA, brdf ):
    LunarIrradiance=mt2009(mt_time_format)
    Lm=LunarCorrection(LZA,brdf,LunarIrradiance['lun_irrad_scl'])
    Lm[LZA>90]=0 #set LM=0 and do not substract if LZA >90, moon under horizon
    DNB_corrected = np.subtract(DNB, Lm)
    DNB_corrected[DNB_corrected < 0] = 0
    return(DNB_corrected)
    
    
def generateGeotiffs(inDir, outDir, BRDF_Dir, exportOriginalDNB = False):    
    
    # ---------------------------------------------------------------------------------------------------------- 
  
    #BRDF_Dir='../data/BRDF_vnp43ma4v001_h20v05'
    
       
    resolution= 500  

    # ---------------------------------------------------------------------------------------------------------- 
    xres = 527
    yres = 645


    
    xmin = 369503	
    ymin = 4019222	
    xmax = 633003	
    ymax = 4341722
    
    greek_def = geometry.AreaDefinition('ggrs', 'Greek Grid', 'ggrs',
                               {  'x_0':'500000', 'y_0':'0', 'lat_0': '0','k':'0.9996', 'lon_0': '24', 'proj': 'tmerc', 'ellps':'GRS80', 'units':'m'},
                                xres, yres,
                                [xmin,ymin,xmax,ymax])  
    
    
    
    
    
    
    #os.chdir(inDir)   
    if not os.path.exists(outDir): os.makedirs(outDir)
    #import re
    pattern='VNP46A1.A2018*.*h5'
    hdf_files = glob.glob(os.path.join(inDir,pattern))
    #hdf_files = [os.path.join(inDir,f) for f in os.listdir(inDir) if re.search(r'(^VNP46A1.A201[7-9](15[2-9]|16[0-9]|17[0-9]|18[0-9]|19[0-9]|20[0-9]|21[0-9]|22[0-9]|23[0-9]|24[0-3]).*\.h5$)', f)]
    hdf_files.sort()
    


    for hdf in hdf_files:
        inFile = os.path.basename(hdf)
        print(inFile)
    
        jd = julian_days(inFile.split(".")[1])
        julian_day=jd['julian_day']
        
        
        BRDF_files = glob.glob("{}/VNP43MA4.{}.*.h5".format(BRDF_Dir,julian_day))
        if len(BRDF_files)==0:#σε περίπτωση που δεν βρεθεί αντίστοιχο αρχείο BRDF, πήγαινε στο επόμενο VNP46A1
            continue
        BRDF_file = BRDF_files[0]
        outName = inFile.rsplit('.', 1)[0]                      # Parse out the file extension, keep file name
        
        tiff_files = glob.glob(os.path.join(outDir,"{}.*.tif".format(outName)))
        if len(tiff_files)>0:
            print ("Geotiff already exist: {}".format(outName))
            continue

        
        ymd = datetime.datetime.strptime(inFile.split(".")[1][1:8], '%Y%j').date().strftime('%Y%m%d')
        
        brdf=BRDF_data(BRDF_file,greek_def)
      
        
    
        with h5py.File(hdf, 'r') as hdf:
            try: 
            
                fileMetadata = hdf['HDFEOS']['GRIDS']['VNP_Grid_DNB'].attrs
                lrcLon=fileMetadata['EastBoundingCoord'][0]
                ulcLon=fileMetadata['WestBoundingCoord'][0]
                ulcLat=fileMetadata['NorthBoundingCoord'][0]
                lrcLat=fileMetadata['SouthBoundingCoord'][0]
                
                
                # ================== DNB ===============================
                #At-sensor DNB radiance, nW·cm-2·sr-1 (16-bit unsigned integer )
                scale_factor=0.1 #DNB scale factor
                #DNB_fillvalue=65535 δεν χρειάζεται να κάνω masking με το fillvalue γιατί κάτά το export με το geotiff σε band.SetNoDataValue(fill_value) τα pixels με 65535 γίνονται Nodata
                DNB = hdf['HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/']['DNB_At_Sensor_Radiance_500m'][...]*scale_factor
                
                # ================== LZA ===============================
                LZ_fillvalue=-32768
                LZ_scale_factor=0.01
                Lunar_Zenith= hdf['HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/']['Lunar_Zenith'][...]*LZ_scale_factor
                LZ_mask = np.isin(Lunar_Zenith, LZ_fillvalue)
            
                 # ================== QF_DNB ===============================
                #QF_DNB_fillvalue=65535
                QF_DNB  = hdf['HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/']['QF_DNB'][...]
                mask_QF_DNB =np.isin(QF_DNB, np.array([1,2,4,8,16,256,512,1024,2048]))
            
                # ================== QF_Cloud_Mask ===============================
                QF_Cloud_Mask= hdf['HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/']['QF_Cloud_Mask'][...]    
                #DNB_mask = np.isin(DNB, DNB_fillvalue)
                
                Cloud_Detection_Results_Confidence_Indicator=(QF_Cloud_Mask & 192)>>6
                mask_Cloud_Detection_Results_Confidence_Indicator=np.logical_not(np.isin(Cloud_Detection_Results_Confidence_Indicator, np.array([0])))
                
                Shadow=(QF_Cloud_Mask & 256)>>8
                mask_Shadow=np.isin(Shadow, np.array([1]))
                
                Snow_Ice_Surface=(QF_Cloud_Mask & 1024)>>10
                mask_Snow_Ice_Surface=np.isin(Snow_Ice_Surface, np.array([1]))
                
                
                Land_Water_Background=(QF_Cloud_Mask & 14)>>1
                mask_Land_Water_Background= np.isin(Land_Water_Background, np.array([2,3])) # Inland and Sea Water
                
                # ================== UTC_time ===============================
                UTC_time_fillvalue=-999.9
                UTC_time= hdf['HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/']['UTC_Time'][...]
                UTC_time_mask = np.isin(UTC_time, UTC_time_fillvalue)
        
                
                # ================== Combine all Masks ===============================
                Mask=np.logical_or.reduce((LZ_mask,
                                           UTC_time_mask,
                                           mask_Land_Water_Background,
                                           mask_Cloud_Detection_Results_Confidence_Indicator,
                                           mask_Shadow,
                                           mask_Snow_Ice_Surface,
                                           mask_QF_DNB
                                          ))    
                
                # ================== apply fill_value to masked pixels ================
                fill_value = 65535
                DNB[Mask]=fill_value
                UTC_time[Mask]=UTC_time_fillvalue
                Lunar_Zenith[Mask]=LZ_fillvalue
                
                DNB = ma.masked_array(DNB, mask=(Mask))
                UTC_time = ma.masked_array(UTC_time, mask=(Mask))
                Lunar_Zenith = ma.masked_array(Lunar_Zenith, mask=(Mask))
                
                # ================== Export to Geotiffs ===============================
                #export DNB as geotiff
         
                XDim, YDim, = 2400,2400
                #https://pyresample.readthedocs.io/en/latest/geo_def.html
                wgs84_def = geometry.AreaDefinition('WGS 84"', 'WGS_1984', 'WGS 84',
                                               "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                                                XDim, YDim,
                                                [ulcLon, lrcLat,
                                                 lrcLon,ulcLat])  
            
            
            
            
        
            
                DNB_2100 = reproject(DNB, wgs84_def, greek_def)
                DNB_2100_ma = ma.masked_array(DNB_2100.image_data, mask=(DNB_2100.image_data==65535))
                Lunar_Zenith_2100 = reproject(Lunar_Zenith, wgs84_def, greek_def)
                UTC_time_2100 = reproject(UTC_time, wgs84_def, greek_def).image_data
                UTC_time_2100_filtered = UTC_time_2100.compressed()#UTC_time_2100[UTC_time_2100>=0] #avoid fill value -999.9 #compressed: Return all the non-masked data as a 1-D array.
        
                print("Unique hours Values:{}".format(np.unique(UTC_time_2100_filtered.astype('int'))))
                
                
                if np.unique(UTC_time_2100_filtered.astype('int')).size==1: #αν τα δεδομένα δεν σπάνε σε δύο βραδιές
                    
                    mt_time_format = modeUTC_time(UTC_time_2100_filtered, ymd)
                    
                    
                    DNB_corrected = correctDNB(mt_time_format,DNB_2100_ma, Lunar_Zenith_2100.image_data, brdf )
                    #export DNB
                    mydate="{}-{}".format(jd['year'],jd['julian_day_str'])
                    geotiff(outDir,
                            "{}.DNB_{}_#{}#".format(outName,mt_time_format, mydate),
                            resolution,
                            fill_value,
                            DNB_2100,
                            np.ma.filled(DNB_corrected, fill_value=fill_value))
                    if exportOriginalDNB:
                        geotiff(outDir,
                            "{}.original_DNB_{}_#{}#".format(outName,mt_time_format, mydate),
                            resolution,
                            fill_value,
                            DNB_2100,
                            np.ma.filled(DNB_2100_ma, fill_value=fill_value))
                        
                else: # αν τα δεδομένα σπάνε σε δυο βραδίες
                    for i in np.unique(UTC_time_2100_filtered.astype('int')):
                        hour_limit = 20
                        if i>hour_limit: #current date
                            #apply filter only for >22
                            mydate="{}-{}".format(jd['year'],jd['julian_day_str'])
                            mt_time_format = modeUTC_time(UTC_time_2100_filtered[UTC_time_2100_filtered>hour_limit],ymd)
                            DNB_2100_ma_by_time = ma.masked_array(DNB_2100_ma, mask=(UTC_time_2100.astype('int')<hour_limit))   
        
                        else: #previous date
                            #apply filter only for >22
                            mt_time_format = modeUTC_time(UTC_time_2100_filtered[UTC_time_2100_filtered<hour_limit],ymd)
                            DNB_2100_ma_by_time = ma.masked_array(DNB_2100_ma, mask=(UTC_time_2100.astype('int')>hour_limit))
                            mydate="{}-{}".format(jd['previous_day_year'],jd['previous_julian_day_str'])
                            
                        DNB_corrected = correctDNB(mt_time_format,DNB_2100_ma_by_time, Lunar_Zenith_2100.image_data, brdf)
                        
                        geotiff(outDir,
                                "{}.DNB_{}_#{}#".format(outName,mt_time_format, mydate),
                                resolution,
                                fill_value,
                                DNB_2100,
                                np.ma.filled(DNB_corrected, fill_value=fill_value))
                        if exportOriginalDNB:
                             geotiff(outDir,
                                "{}.original_DNB_{}_#{}#".format(outName,mt_time_format, mydate),
                                resolution,
                                fill_value,
                                DNB_2100,
                                np.ma.filled(DNB_2100_ma_by_time, fill_value=fill_value))

            except Exception as e: 
                print(e)
                print("An exception occurred.File:{}".format(str(inFile)))            
                f = open('errors.txt', 'a')
                f.write("{}\n".format(str(inFile))) # python will convert \n to os.linesep
                f.close() 
            finally:
                print("Done")


if __name__ == '__main__':
    inDir  = sys.argv[1] # directory with VNP46A1
    outDir = sys.argv[2] # directory to export lunar corrected (and original) VNP46A1 as geotiff
    BRDF_Dir = sys.argv[3] # 
    exportOriginalDNB=sys.argv[4] # True/False
    generateGeotiffs(inDir,outDir,BRDF_Dir,exportOriginalDNB)
    
    
