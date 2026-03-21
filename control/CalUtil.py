"""
module: CalUtil
"""

import datetime
import calendar
import unittest

__all__ = ["IllegalMonthError", "IllegalYearRangeError",
           "strDate", "strDayNumber", "extractDate",
           "leapYear", "lengthYear", "lengthMonth",
           "dateToDayNumber", "midMonth"]

class IllegalMonthError(ValueError):
    def __init__(self, month):
        self.month = month
    def __str__(self):
        return "bad month number %r; must be 1-12" % self.month

class IllegalYearRangeError(ValueError):
    def __init__(self, beginYear, endYear):
        self.beginYear = beginYear
        self.endYear = endYear
    def __str__(self):
        return "bad year range %r - %r" % (self.beginYear, self.endYear)

def strDate(year, month, day):
    """
    year
    month
    day
    return - tuple of strings (yyyy, mm, dd)
    """

    strYear  = "%4.4i" % year
    strMonth = "%2.2i" % month
    strDay   = "%2.2i" % day

    return strYear, strMonth, strDay

def strDayNumber(year, dayNumber):
    """
    year
    dayNumber
    return - tuple of strings (yyyy, ddd)
    """

    strYear   = "%4.4i" % year
    strDayNum = "%3.3i" % dayNumber

    return strYear, strDayNum

def extractDate(date):
    """
    date
    return - tuple (year, month, day)
    """

    year = date // 10000
    month = (date - 10000 * year) // 100
    day = date - 10000 * year - 100 * month

    return year, month, day

def leapYear(year):
    """
    return - True if leap year
    """

    return (year % 4 == 0 and year % 100 != 0) \
        or (year % 400 == 0)

def lengthYear(year):
    """
    return - length of year
    """

    if (leapYear(year)):
        return 366
    else:
        return 365

def lengthMonth(year, month):
    """
    return - length of month
    """

    if not 1 <= month <= 12:
        raise IllegalMonthError(month)

    length = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if (leapYear(year)):
        length[1] = 29

    return length[month - 1]

def addMonth(year, month, nMonth):
    """
    year
    month
    nMonth - number of months to add, positive or negative
    return - tuple (incremented year, incremented month)
    """

    if not 1 <= month <= 12:
        raise IllegalMonthError(month)

    monthsSince2000 = (year - 2000) * 12 + month - 1
    incMonthsSince2000 = monthsSince2000 + nMonth
    incYear = incMonthsSince2000 // 12 + 2000
    # incMonth = (month - 1 + nMonth) % 12 + 1
    incMonth = incMonthsSince2000 % 12 + 1

    return (incYear, incMonth)

def dateToDayNumber(year, month, day):

    if not 1 <= month <= 12:
        raise IllegalMonthError(month)

    dayNumber = 0

    for iMonth in range(1, month):
        dayNumber += lengthMonth(year, iMonth)

    dayNumber += day

    return dayNumber

def midMonth(year, month):

    day = int(round(lengthMonth(year, month) / 2.0))

    return dateToDayNumber(year, month, day)

def iterateMonth(beginYear, beginMonth, endYear, endMonth,
    splitDate=False):

    if not beginYear <= endYear:
        raise IllegalYearRangeError(beginYear, endYear)
    if not 1 <= beginMonth <= 12:
        raise IllegalMonthError(beginMonth)
    if not 1 <= endMonth <= 12:
        raise IllegalMonthError(endMonth)

    listYearMonths = []
    listYears = []
    listMonths = []

    for year in range(beginYear, endYear + 1):
        bmonth, emonth = 1, 12
        if (year == beginYear):
            bmonth = beginMonth
        if (year == endYear):
            emonth = endMonth
        for month in range(bmonth, emonth + 1):
            listYearMonths.append((year, month))
            listYears.append(year)
            listMonths.append(month)

    if (splitDate):
        return listYears, listMonths
    else:
        return listYearMonths

# test classes
class TestCalendar(unittest.TestCase):

    def setUp(self):
        self.objDateToday = datetime.date.today()

    def testExtractDate(self):
        date = 10000 * self.objDateToday.year \
             + 100 * self.objDateToday.month \
             + self.objDateToday.day
        year, month, day = extractDate(date)
        self.assertEqual(year, self.objDateToday.year)
        self.assertEqual(month, self.objDateToday.month)
        self.assertEqual(day, self.objDateToday.day)

    def testLeapYear(self):
        for year in range(1900, 2101):
            self.assertEqual(leapYear(year), calendar.isleap(year))

    def testLengthYear(self):
        for year in range(1900, 2101):
            if (calendar.isleap(year)):
                length = 366
            else:
                length = 365
            self.assertEqual(lengthYear(year), length)

    def testLengthMonth(self):
        print("\n")
        for year in range(2000, self.objDateToday.year + 1):
            print("testLengthMonth: %4.4i" % year)
            for month in range(1, 13):
                self.assertEqual(lengthMonth(year, month), \
                    calendar.monthrange(year, month)[1])

    def testAddMonth(self):
        print("\n")
        for year in range(2000, self.objDateToday.year + 1):
            for month in range(1, 13):
                nMonth = -3
                incYear, incMonth = addMonth(year, month, nMonth)
                print("testAddMonth: %4.4i %2.2i %4.4i %2.2i" \
                      % (year, month, incYear, incMonth))

    def testDateToDayNumber(self):
        print("\n")
        for year in range(2000, self.objDateToday.year + 1):
            print("testDateToDayNumber: %4.4i" % year)
            for month in range(1, 13):
                for day in range(lengthMonth(year, month) + 1):
                    dayNumber = dateToDayNumber(year, month, day)

    def tearDown(self):
        del self.objDateToday

class TestBuiltInCalendar(unittest.TestCase):

    def setUp(self):
        self.year = datetime.date.today().year
        self.calendar = calendar.Calendar()

    def testYearDatesCalendar(self):
        yearCalendar = self.calendar.yeardatescalendar(self.year)

    def tearDown(self):
        del self.year
        del self.calendar

if __name__ == "__main__":
    suiteCalendar \
        = unittest.TestLoader().loadTestsFromTestCase(TestCalendar)
    suiteBuiltInCalendar \
        = unittest.TestLoader().loadTestsFromTestCase(TestBuiltInCalendar)
    tests = unittest.TestSuite([suiteCalendar, suiteBuiltInCalendar])
    unittest.TextTestRunner(verbosity=2).run(tests)

