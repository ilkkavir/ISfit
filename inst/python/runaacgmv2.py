def runaacgmv2(inlat,inlon,height,year,month,day,hour,minute,second,method_code):
    import datetime
    import aacgmv2
    dt = datetime.datetime(year,month,day,hour,minute,second)
    coord = aacgmv2.convert_latlon_arr(inlat,inlon,height,dt,method_code)
    return coord
