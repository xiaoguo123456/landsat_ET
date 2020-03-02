from functions import *
import os
from read_MTL import landsat_metadata
from osgeo import gdal
import gdal
from ET_algorithm import *

def get_EA(meta_path,Ta,RH,Rn,WS,DT, NIR_name=None, RED_name=None,outdir=None):
    """

    :param meta_path:
    :param merra_path:
    :param NIR_name:
    :param RED_name:
    :param outdir:输出路径
    :return:
    """
    meta_path = os.path.abspath(meta_path)
    path=os.path.split(meta_path)[0]
    meta_info=landsat_metadata(meta_path)

    if NIR_name == None:
        band5_path = os.path.join(path+'/'+meta_info.FILE_NAME_BAND_5)
    else:
        band5_path = NIR_name
    if RED_name == None:
        band4_path = os.path.join(path+'/'+meta_info.FILE_NAME_BAND_4)
    else:
        band4_path = RED_name
    band3_path = os.path.join(path+'/'+meta_info.FILE_NAME_BAND_3)

    inds = gdal.Open(band4_path)
    in_band=inds.GetRasterBand(1)
    x_size=in_band.XSize
    y_size=in_band.YSize
    red_data = gdal.Open(band4_path).ReadAsArray().astype(np.float)
    nir_data = gdal.Open(band5_path).ReadAsArray().astype(np.float)
    green_data=gdal.Open(band3_path).ReadAsArray().astype(np.float)

    if outdir is not None:
        outdir = os.path.abspath(outdir)
        outname = create_outname(outdir, meta_path, "EA", "tif")
    else:
        folder = os.path.split(meta_path)[0]
        outname = create_outname(folder, meta_path, "EA", "tif")
    driver = inds.GetDriver()
    outimg = driver.Create(outname, x_size, y_size, 1, gdal.GDT_Float32)
    outimg.SetGeoTransform(inds.GetGeoTransform())
    outimg.SetProjection(inds.GetProjection())
    outBand = outimg.GetRasterBand(1)
    np.seterr(divide='ignore', invalid='ignore')
    ndvi = np.where(red_data + nir_data == 0, 0, (nir_data - red_data) / (nir_data + red_data))

    ndwi=np.where(green_data+nir_data==0,0,(green_data-nir_data)/(nir_data + red_data))
    # NDVI,LEcor,Ta,RH,Rn,WS,DT,Rs

    xBlockSize = 128
    yBlockSize = 128
    for i in range(0, y_size, yBlockSize):
        if i + yBlockSize < y_size:
            numRows = yBlockSize
        else:
            numRows = y_size - i
        for j in range(0, x_size, xBlockSize):
            if j + xBlockSize < x_size:
                numCols = xBlockSize
            else:
                numCols = x_size - j
            ndvi_data = ndvi[i:i + numRows, j:j + numCols]
            ndwi_data = ndwi[i:i + numRows, j:j + numCols]
            Ta1=Ta[i:i + numRows, j:j + numCols]
            RH1=RH[i:i + numRows, j:j + numCols]
            Rn1=Rn[i:i + numRows, j:j + numCols]
            WS1=WS[i:i + numRows, j:j + numCols]
            DT1=DT[i:i + numRows, j:j + numCols]
            e = np.where(ndvi_data==0,0,EA(NDVI=ndvi_data, NDWI=ndwi_data,Ta=Ta1, RH=RH1, Rn=Rn1, WS=WS1, DT=DT1 ))
            outBand.WriteArray(e, j, i)
    outBand.SetNoDataValue(0)
    outBand.FlushCache
    outimg = None

    print("Saved output at {0}".format(outname))