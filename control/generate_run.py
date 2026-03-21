import os
import string

import CalUtil

def createSeriesMonthly(beginYear, beginMonth, endYear, endMonth,
    dirBase, strBase):

    listMonths = CalUtil.iterateMonth(
        beginYear, beginMonth, endYear, endMonth)

    strFilesMonthly = ""

    for year, month in listMonths:
        nDay = CalUtil.lengthMonth(*(year, month))
        strYear, strMonth, strDay \
            = CalUtil.strDate(year, month, nDay)
        dirYearly = os.path.join(dirBase, strYear)
        dirMonthly = os.path.join(dirYearly, strMonth)
        if os.path.exists(dirMonthly):
            strFileBase = os.path.join(dirMonthly,
                strBase + "." + strYear + strMonth)
            strFileDailyWild = strFileBase + "??.nc"
            strFileBase = os.path.join(dirMonthly,
                strBase.replace("daily", "monthly") + "." + strYear + strMonth)
            strFileMonthly = strFileBase + ".nc"
            command = "ncra -O " + strFileDailyWild + " " + strFileMonthly
            strFilesMonthly += strFileMonthly + " "
            print(command)
            os.system(command)

    strBeginYear, strBeginMonth, strDay \
        = CalUtil.strDate(beginYear, beginMonth, 1)
    strEndYear, strEndMonth, strDay \
        = CalUtil.strDate(endYear, endMonth, 1)
    strCatFileMonthly = os.path.join(dirBase,
        strBase.replace("daily", "monthly")
        + "." + strBeginYear + strBeginMonth
        + "-" + strEndYear + strEndMonth + ".nc")

    command = "ncrcat -O " + strFilesMonthly + strCatFileMonthly
    os.system(command)

def doByMonth(beginYear, beginMonth, endYear, endMonth,
    dirControl, nameTemplate):

    listMonths = CalUtil.iterateMonth(
        beginYear, beginMonth, endYear, endMonth)

    for year, month in listMonths:
        nDay = CalUtil.lengthMonth(*(year, month))
        strYear, strMonth, strDay \
            = CalUtil.strDate(year, month, nDay)
        dirRun = os.path.join(dirControl, strYear)
        nameFile = os.path.join(dirRun,
            nameTemplate.replace("yyyy", strYear).replace("mm", strMonth))
        print(nameFile)

        # decrement month
        decYear, decMonth = CalUtil.addMonth(year, month, -1)
        nDecDay = CalUtil.lengthMonth(decYear, decMonth)
        strDecYear, strDecMonth, strDecDay \
            = CalUtil.strDate(decYear, decMonth, 1)

        handleTemplate = open(nameTemplate, "r")
        strTemplate = string.Template(handleTemplate.read())
        handleTemplate.close()
        print(strYear, strMonth, nDecDay + nDay)
        strFile = strTemplate.safe_substitute(
            yyyy=strYear, mm=strMonth, nn=("%2.2i" % (nDecDay + nDay)),
            yyyymmdd_begin=strDecYear+strDecMonth+strDecDay)
        # print(strFile)
        if (not os.path.exists(dirRun)):
            print("creating " + dirRun)
            os.mkdir(dirRun)
        handleFile = open(nameFile, "w")
        handleFile.write(strFile)
        handleFile.close()

def setHalfYear(beginYear, endYear,
    dirControl, nameTemplate):

    for year in range(beginYear, endYear + 1):
        for month in range(1, 7):
            strYear, strMonth, strDay \
                = CalUtil.strDate(year, month, 1)
            dirRun = os.path.join(dirControl, strYear)
            nameFile = os.path.join(dirRun,
                nameTemplate + strYear + "01-06.csh")

        for month in range(8, 13):
            strYear, strMonth, strDay \
                = CalUtil.strDate(year, month, 1)
            dirRun = os.path.join(dirControl, strYear)
            nameFile = os.path.join(dirRun,
                nameTemplate + strYear + "07-12.csh")

def stageByMonth(beginYear, beginMonth, endYear, endMonth,
    dirStage, strRun, listLinks):

    listMonths = CalUtil.iterateMonth(
        beginYear, beginMonth, endYear, endMonth)

    for year, month in listMonths:
        nDay = CalUtil.lengthMonth(*(year, month))
        strYear, strMonth, strDay \
            = CalUtil.strDate(year, month, nDay)
        dirRun = os.path.join(dirStage, strRun + strYear + strMonth)
        if (not os.path.exists(dirRun)):
            print("creating " + dirRun)
            os.mkdir(dirRun)
        else:
            print("exists " + dirRun)
        # create data links
        if (not os.listdir(dirRun)):
            for link in listLinks:
                command = "ln -s " \
                        + link + " " + dirRun 
                print(command)
                os.system(command)

if __name__ == "__main__":

    dirStage = os.path.join(os.getenv("DATA_DIR"), "MATCH")
    dirModel = os.path.join(os.getenv("HOME"), "MATCH")
    dirControl = os.path.join(dirModel, "runs")

    listLinks = [os.path.join(dirModel, "data", "*"),
                 os.path.join(dirModel, "grib", "*")]

    # doByMonth(2015, 1, 2015, 12, dirControl, "run_sarb4_hourly_yyyymm.csh")
    stageByMonth(2013, 1, 2013, 12, dirStage, "run_sarb4_viirs_hourly_", listLinks)

