def startflipchem(year,month,day,hour,minute,second):
    import datetime
    import flipchem
    dt = datetime.datetime(year,month,day,hour,minute,second)
    fc = flipchem.Flipchem(dt)
    return fc
