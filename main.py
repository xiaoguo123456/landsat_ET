import os
from GET_ET import *
from functions import *
from pyhdf.SD import SD
from read_MTL import landsat_metadata
import time
from multiprocessing import Pool, cpu_count
#guoxiaozheng123456@163.com




if __name__=="__main__":
    """
    """
    pool=Pool(cpu_count())
    # 读取对应hdf文件中的全部landsat数据的头文件，即MTL_txt文件
    mtl_file = MTL_file_name('/mnt/lustre/users/boyyunjun/GXZ/Landsat_data/')
    start_time=time.time()
    print("总共{0}幅图像".format(len(mtl_file)))
    for i in range(0, len(mtl_file)):
 
        if not os.path.isfile(mtl_file[i]):
            print(mtl_file[i] + '不存在')
            continue
        path = os.path.split(mtl_file[i])[0]
        mtl_name = os.path.split(mtl_file[i])[1]

        meta_info = landsat_metadata(mtl_file[i])
        # time
        year = meta_info.DATE_ACQUIRED[0:4]
        month = meta_info.DATE_ACQUIRED[5:7]
        day = meta_info.DATE_ACQUIRED[8:10]
        # row col
        row = meta_info.WRS_ROW
        col = meta_info.WRS_PATH

        # 根据影像读取信息
        band5_path = os.path.join(path , meta_info.FILE_NAME_BAND_5)
        
        print('正在处理第{0}幅影像,行号为{1}，列号为{2}，时间{3}...'.format((i + 1),
                                                            meta_info.WRS_ROW,
                                                            meta_info.WRS_PATH,
                                                            meta_info.DATE_ACQUIRED))
        # 读取merra数据
        merra_path = os.path.join('/mnt/lustre/GLASS02/MERRA2_LT_daily' + '/{0}'.format(year) + '/{0}{1}'.format(year, month) + \
                                    '/MERRA2.LT.GEO.DAILY.{0}{1}{2}.hdf'.format(year, month, day))
        if not os.path.isfile(merra_path):
            print(mtl_name + '对应的merra数据不存在')
            continue
        merra_data = SD(merra_path)
        rn_data = merra_data.select('NET_RADIATION_LT_DAILY_MEAN')[:]
        ws_data = merra_data.select('WIND2M_LT_DAILY_MEAN')[:]
        rh_data = merra_data.select('RH_LT_DAILY_MEAN')[:]
        ta_data = merra_data.select('T2M_LT_DAILY_MEAN')[:] - 273.15
        dt_data = merra_data.select('T2M_LT_DAILY_MAX')[:] - merra_data.select('T2M_LT_DAILY_MIN')[:]

        ws1 = mask(merra_data=ws_data, landsat=band5_path)
        rh1 = mask(merra_data=rh_data, landsat=band5_path)
        ta1 = mask(merra_data=ta_data, landsat=band5_path)
        dt1 = mask(merra_data=dt_data, landsat=band5_path)
        rn1 = mask(merra_data=rn_data, landsat=band5_path)
        temp_name=mtl_file[i]
        get_EA(meta_path=mtl_file[i], Ta=ta1, RH=rh1, Rn=rn1, WS=ws1, DT=dt1, outdir='/mnt/lustre/users/boyyunjun/GXZ/landsat_ET/')

#            pool.apply_async(get_EA,(temp_name,ta1,rh1,rn1,ws1,dt1,'/mnt/lustre/users/boyyunjun/GXZ/landsat_ET/',))
       # pool.apply_async(get_EA(meta_path=mtl_file[i], Ta=ta1, RH=rh1, Rn=rn1, WS=ws1, DT=dt1, outdir='/mnt/lustre/users/boyyunjun/GXZ/landsat_ET/'))

    pool.close()
    pool.join()
    end_time=time.time()
    print('总共用时{}'.format(end_time-start_time))













