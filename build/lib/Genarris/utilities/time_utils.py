
"""
Created on Tue Feb  6 19:26:58 2018

@author: timcrose
"""

from datetime import datetime
from datetime import date
from time import sleep
from time import mktime
from time import time as gtime
from time import ctime
from time import localtime
from time import strftime

def get_greg_time_from_time_str(time_str='00:00:00', ds=None, DATE_FMT='%Y-%m-%d', TIME_FMT='%H:%M:%S'):

    if ds is None:
        #default to today
        ds = datetime.today().strftime(DATE_FMT)
        
    dtt = datetime.strptime(ds + ' ' + time_str, DATE_FMT + ' ' + TIME_FMT)
    greg_time = float(mktime(dtt.timetuple()))

    return greg_time

def get_greg_from_mdYHMS(mon, day, yr, hr, minute, sec):
    fmt = '%m-%d-%Y %H:%M:%S'
    s = str(mon).zfill(2) + '-'  + str(day).zfill(2) + '-' + str(yr).zfill(4) + ' ' + str(hr).zfill(2) + ':' + str(minute).zfill(2) + ':' + str(int(sec)).zfill(2)
    dt = datetime.strptime(s, fmt)
    tt = dt.timetuple()
    greg_time = mktime(tt)
    return greg_time

def delay_start(time_of_day_to_start):
    start_time_greg = get_greg_time_from_time_str(time_str=time_of_day_to_start)

    if gtime() > start_time_greg:
        #trying to start tomorrow so need to get amount of delay until start time
        delay_time = get_greg_time_from_time_str(time_str='23:59:59') +1 - gtime() + start_time_greg - get_greg_time_from_time_str(time_str='00:00:00')
    else:
        delay_time = start_time_greg - gtime()

    sleep(delay_time)
    
def get_now_MM_DD_YYYY():
    return str(datetime.now().month) + '_' + str(datetime.now().day) + '_' + str(datetime.now().year)


def wait_til_weekday(time_of_day_to_start='09:00:00'):
    weekno = datetime.today().weekday()

    #Monday is 0 and Sunday is 6
    while not ((weekno in range(5)) and gtime() < get_greg_time_from_time_str(time_str=time_of_day_to_start)):
        sleep(3666)
        weekno = datetime.today().weekday()


def wait_til_weekday_old(time_of_day_to_start='09:15:00'):
    weekno = datetime.today().weekday()
    weekno_original = datetime.today().weekday()

    if weekno_original == 4:
        start_time_greg = get_greg_time_from_time_str(time_str=time_of_day_to_start)
        if gtime() < start_time_greg:
            #trying to start today
            return

    # check if will be a weekday tmw
    #Monday is 0 and Sunday is 6
    #Wait until today and tmw is a weekday
    #today is a weekday if weekno is 0, 1, 2, 3, 4
    #tmw is a weekday if weekno is 0, 1, 2, 3, 6
    #So today is a weekday and tmw is a weekday if weekno is 0, 1, 2, 3
    #But if we started the program on Friday morning (before the requested start time), then we should also allow it to exit the loop
    while weekno not in [0, 1, 2, 3]:
        sleep(3666)
        weekno = datetime.today().weekday()


def get_date_time_str_from_greg(greg):
    return ctime(int(greg))
def get_time_str_from_greg(greg):
    return strftime("%H:%M:%S", localtime(greg))

def day_of_week_from_date_str(date_str, delim='_'):
    '''
    date_str: str
        must be of format month_day_year. e.g. 12_20_2018 (or other delimiter)
    delim: delimiter of passed input string. See date_str above
    
    return:
    day-of-the-week number 0 - 6 where 0 is Monday, 6 is Sunday
    '''
    month, day, year = [int(i) for i in date_str.split(delim)]
    born = date(year, month, day)
    #return born.strftime('%A') Day name like "Sunday"
    return born.weekday()
