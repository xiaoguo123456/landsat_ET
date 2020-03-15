from functions import *
import os
from read_MTL import landsat_metadata
from osgeo import gdal
import gdal
from ET_algorithm import *
import multiprocessing
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
    muti_caculate(128,128,ndvi,ndwi,Ta,RH,Rn,WS,DT,outBand)
    outBand.SetNoDataValue(0)
    outBand.FlushCache
    outimg = None

    print("Saved output at {0}".format(outname))

def gen_tub(x_bloc,y_bloc,arr):
    tub=[]
    for i in range(0,arr.shape[0],y_bloc):
        if i+y_bloc<arr.shape[0]:
            num_row=y_bloc
        else:
            num_row=arr.shape[0]-i
        for j in range(0,arr.shape[1],x_bloc):
            if j+x_bloc<arr.shape[1]:
                num_col=x_bloc
            else:
                numcol=arr.shape[1]-j
            tub.append((arr[i:i+num_row,j:j+num_col],j,i))   
    return tub
def muti_caculate(x_block,y_block,ndvi,ndwi,ta,rh,rn,ws,dt,outBand):
    result=[]
    row_col=[]
    ndvi_list=gen_tub(x_block,y_block,ndvi)
    ndwi_list=gen_tub(x_block,y_block,ndwi)
    ta_list=gen_tub(x_block,y_block,ta)
    rh_list=gen_tub(x_block,y_block,rh)
    rn_list=gen_tub(x_block,y_block,rn)
    ws_list=gen_tub(x_block,y_block,ws)
    dt_list=gen_tub(x_block,y_block,dt)
    pool=multiprocessing.Pool(multiprocessing.cpu_count())
    for i in range(len(ndvi_list)):
        temp_result=pool.apply_async(EA,(ndvi_list[i][0],ndwi_list[i][0],ta_list[i][0],rh_list[i][0],rn_list[i][0],ws_list[i][0],dt_list[i][0],))
        result.append(temp_result)
        row_col.append((ndvi_list[i][1],ndvi_list[i][2]))
    for j,data in enumerate(result):
        outBand.WriteArray(data.get(),row_col[j][0],row_col[j][1])