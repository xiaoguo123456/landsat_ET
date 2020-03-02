from osgeo import osr
import os
import gdal
import numpy as np
# from scipy import interpolate
import cv2
#投影坐标转换为地理坐标，在转换为图像坐标
#参考http://www.voidcn.com/article/p-xtxfepsi-bgc.html
def getSRSPair(dataset):
    """ 获得给定数据的投影参考系和地理参考系
    :param dataset: GDAL地理数据 :return: 投影参考系和地理参考系 """
    prosrs = osr.SpatialReference()
    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()
    return prosrs, geosrs
def geo2lonlat(dataset, x, y):
    """将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
    :param dataset: GDAL地理数据
    :param x: 投影坐标x
    :param y: 投影坐标y
    :return: 投影坐标(x, y)对应的经纬度坐标(lon, lat)
    """
    prosrs, geosrs = getSRSPair(dataset)
    ct = osr.CoordinateTransformation(prosrs, geosrs)
    coords = ct.TransformPoint(x, y)
    return coords[:2]
def lonlat2geo(dataset, lon, lat):
    """将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定）
    :param dataset: GDAL地理数据
    :param lon: 地理坐标lon经度
    :param lat: 地理坐标lat纬度
    :return: 经纬度坐标(lon, lat)对应的投影坐标 """
    prosrs, geosrs = getSRSPair(dataset)
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    coords = ct.TransformPoint(lon, lat)
    return coords[:2]
def imagexy2geo(dataset, row, col):
    """
    :param dataset:
    :param row:像素的行号
    :param col:像素的列号
    :return:
    """
    trans = dataset.GetGeoTransform()
    px = trans[0] + row * trans[1] + col * trans[2]
    py = trans[3] + row * trans[4] + col * trans[5]
    return px, py
def geo2imagexy(dataset, x, y):
    """
    :param dataset:
    :param x:投影或地理坐标x
    :param y: 投影或地理坐标y
    :return:
    """
    trans = dataset.GetGeoTransform()
    a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
    b = np.array([x - trans[0], y - trans[3]])
    return np.linalg.solve(a, b)  # 使用numpy的linalg.solve进行二元一次方程的求解
def create_outname(outdir, inname, suffix, ext=False):
    """
    Quick way to create unique output file names within iterative functions
    This function is built to simplify the creation of output file names. Function allows
    ``outdir = False`` and will create an outname in the same directory as inname. Function will
    add a the user input suffix, separated by an underscore "_" to generate an output name.
    this is useful when performing a small modification to a file and saving new output with
    a new suffix. Function merely returns an output name, it does not save the file as that name.
    :param outdir:      either the directory of the desired outname or False to create an outname
                        in the same directory as the inname
    :param inname:      the input file from which to generate the output name "outname"
    :param suffix:      suffix to attach to the end of the filename to mark it as output
    :param ext:         specify the file extension of the output filename. Leave blank or False
                        and the outname will inherit the same extension as inname.
    :return outname:    the full filepath at which a new file can be created.
    """
    # isolate the filename from its directory and extension
    if os.path.isfile(inname):
        head, tail = os.path.split(inname)
        noext = tail.split('.')[:-1]
        noext = '.'.join(noext)
    else:
        head = ""
        tail = inname
        if "." in inname:
            noext = tail.split('.')[:-1]
            noext = '.'.join(noext)
        else:
            noext = inname
    # create the suffix
    if ext:
        suffix = "_{0}.{1}".format(suffix, ext)
    else:
        ext = tail.split('.')[-1:]
        suffix = "_{0}.{1}".format(suffix, ''.join(ext))
    if outdir:
        outname = os.path.join(outdir, noext + suffix)
        return outname
    else:
        outname = os.path.join(head, noext + suffix)
        return outname

def mask(merra_data,landsat):
    """
    用landsat数据截取merra数据
    :param lansat_data: landsat某个波段的数据
    :param merra:
    :return:
    """

    landsat_data=gdal.Open(landsat)
    xsize=landsat_data.GetRasterBand(1).XSize
    ysize = landsat_data.GetRasterBand(1).YSize
    data=merra_data

    land_geo = imagexy2geo(landsat_data, 0, 0)
    land_lonlan = geo2lonlat(landsat_data, land_geo[0], land_geo[1])

    merra_data = np.zeros((5,5))
    a=int((90-land_lonlan[0])/0.5)+1
    b=int((180+land_lonlan[1])/0.625)
    merra_data = data[a:a + 5, b:b+ 5]
    # y, x = np.mgrid[-1:1:5j, -1:1:5j]
    # data=interpolate.interp2d(y,x,merra_data, kind='linear')
    # xnew = np.linspace(-1, 1, xsize)  # x
    # ynew = np.linspace(-1, 1, ysize)
    # data=data(xnew,ynew)
    merra_data=np.array(merra_data).astype(np.float)
    data=cv2.resize(merra_data,dsize=(xsize,ysize))
    return data


def MTL_file_name(file_dir):
    L = []
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            if (os.path.splitext(file)[1] == '.txt')&( 'MTL' in file):
                L.append(os.path.join(root, file))
    return L

