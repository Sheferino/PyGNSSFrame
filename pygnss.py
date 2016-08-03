# -*- coding: utf-8 -*-
"""
Created on Fri May 20 15:20:33 2016

@author: Dmitry
@version: 1.0

Python GNSS framework

GNSS data is converted to separate files and moved to folder with same name as data file

Attributes:
    __pathext__ - folder name extension
    __obsext__ - observation files extension
    path - data-folder path
    time - time array, datetime type
    continious_time - time array, float type
    SV_list - filtered satellite names list
    SV_list_full - full satellite names list
    G01,G02...R07...E15...C12...(RINEX names) - containers with GNSS data(dict):
        .C1C, .L2P... (RINEX names) - measurements array, numpy.array type
        .x, .y, .z - satellite coordinates, numpy.array type
        .dt - satellite clock, numpy.array type
        .num - satellite number, int
        .prn - satellite PRN (frequency number/liter), int
        .name - RINEX name
        .fullname - full RINEX name (with liter)
        
Methods:
    __init__ - creating.
    
    __iter__, __next__ - SV  attributes iteration from .SV_list
    
    
    
--->    getobs(timestamp) - return a dict with observations for specified time point
        input: timestamp, datetime type
        output: dict with keys:
            .C1C, .L2P... (RINEX names) - array of measurements for all SV
            .x, .y, .z - array of coordinate for all SV
            coordinates - matrice of satellite coordinates
            name, num, prn - arrays of SV attributes
        
    
    filt - filtering by parameters:
        gnss='G'/'R'/'E'/'C'
        prn=-7:6
        num=1:32
======>     time

    reset - reset .SV_list to .SV_list_full
    
    __parseRINEX__ - parsing RINEX files
    
    
        

"""

from datetime import datetime
import os
import re
import numpy as np
import shutil
from scipy import interpolate


class PyGNSS:
    __pathext__ = '.pygnss'  # расширение папки
    __obsext__ = '.obs.csv'  # расширение файлов с измерениями
    __ephext__ = '.xyz.csv'

    def __init__(self, file_name, time_interval=30):
        '''Creating PyGNSS object'''
        self.path = file_name + self.__pathext__

        # проверяем существование папки с парсеными измереними
        if not (os.path.exists(self.path)):
            # если нет, парсим файлы данных
            if (file_name[-1] == 'o') or (file_name[-1] == 'O'):
                # читаем RINEX 3
                self.__parseRINEX__(file_name, time_interval)
                # читаем эфемериды sp3
                if os.path.exists(os.getcwd() + '/' + file_name[:-3] + 'sp3'):
                    self.__parseEphemeris__(file_name[:-3] + 'sp3')
            elif file_name[-3:-1] == 'csv':
                # читаем измерения GSS
                pass
            elif file_name[-3:-1] == 'mrk':
                # читаем измерения МРК
                pass
        ## Открытие существующего
        # Формирование списка НКА
        self.SV_list_full = []
        for filename in os.listdir(self.path):
            if re.search(self.__obsext__, filename):
                self.SV_list_full.append(filename[:3])
        self.SV_list = self.SV_list_full
        # Формирование вектора времени
        self.time = []
        self.continious_time = []
        file_time = open('{}/time.csv'.format(self.path))
        for line in file_time:
            self.time.append(datetime.fromtimestamp(float(line)))
            self.continious_time.append(float(line))
        self.continious_time = np.array(self.continious_time)
        file_time.close()

        # Создание атрибутов, соответствующих НКА и их заполнение
        # Чтение информации о наблюдениях
        for SV_name in self.SV_list_full:
            file_SVobs = open('{}/{}{}'.format(self.path, SV_name, self.__obsext__))
            flag_firstline = True
            num_obs = 0
            for line in file_SVobs:
                if flag_firstline:
                    flag_firstline = False
                    # создание атрибутов у каждого SV
                    obs_types = line.strip().split(';')
                    SV = dict.fromkeys(obs_types)
                    for key in SV:
                        SV[key] = np.nan * np.arange(len(self.time))
                    SV['name'] = SV_name  # имя НКА
                    SV['num'] = int(SV_name[1:])  # номер НКА
                    # заполнение литер, ГНСС и полных имён
                    if SV_name[0] == 'G':
                        SV['gnss'] = 'G'
                        SV['prn'] = SV['num']
                        SV['fullname'] = SV['name']
                    elif SV_name[0] == 'R':
                        SV['gnss'] = 'R'
                        # подтягивание номера литеры
                        file_prntable = open(self.path + '/litertable.csv')
                        for count in range(SV['num'] - 1):
                            file_prntable.readline()
                        SV['prn'] = float(file_prntable.readline())
                        if not(np.isnan(SV['prn'])):
                            SV['prn']=int(SV['prn'])
                        file_prntable.close()
                        SV['fullname'] = '{}({:+2})'.format(SV['name'], SV['prn'])
                    elif SV_name[0] == 'E':
                        SV['gnss'] = 'E'
                        SV['prn'] = SV['num']
                        SV['fullname'] = SV['name']

                else:
                    dat = line.strip().split(';')
                    for k in range(len(dat) - 1):
                        SV[obs_types[k]][num_obs] = np.nan if dat[k] == '' else float(dat[k])
                    num_obs += 1
            file_SVobs.close()

            # Чтение информации об эфемеридах
            if os.path.exists('{}/{}{}'.format(self.path, SV_name, self.__ephext__)):
                file_SVeph = open('{}/{}{}'.format(self.path, SV_name, self.__ephext__))
                flag_firstline = True
                num_obs = 0
                for line in file_SVeph:
                    if flag_firstline:
                        flag_firstline = False
                        # создание атрибутов у каждого SV
                        eph_types = line.strip().split(';')
                        SV_eph = dict.fromkeys(eph_types)
                        for key in SV_eph:
                            SV_eph[key] = np.nan * np.arange(len(self.time))
                    else:
                        dat = line.strip().split(';')
                        for k in range(len(dat)):
                            SV_eph[eph_types[k]][num_obs] = np.nan if dat[k] == '' else float(dat[k])
                        num_obs += 1
                file_SVeph.close()
                # интерполяция пропусков
                for eph_type in eph_types:
                    ind = ~(np.isnan(SV_eph[eph_type]))
                    f_interp = interpolate.interp1d(self.continious_time[ind], SV_eph[eph_type][ind],
                                                    bounds_error=False, kind='cubic')
                    SV_eph[eph_type] = f_interp(self.continious_time)

                SV.update(SV_eph)  # добавляем словарь с эфемеридами в основной словарь SV

            setattr(self, SV_name, SV)

    def __iter__(self):
        self.__pointer__ = 0
        return (self)

    def __next__(self):
        if self.__pointer__ == len(self.SV_list):
            raise StopIteration
        else:
            self.__pointer__ += 1
            return (getattr(self, self.SV_list[self.__pointer__ - 1]))

    def remove(self):
        shutil.rmtree(self.path)

    def getobs(self, SV_name):
        # потерял значение, задача метода изменена, необходимо переписать
        file_SVobs = open('{}/{}{}'.format(self.path, SV_name, self.__obsext__))
        flag_firstline = True
        num_obs = 0
        for line in file_SVobs:
            if flag_firstline:
                obs_types = line.split(';')
                SV = dict.fromkeys(obs_types)
                for key in SV:
                    SV[key] = np.nan * np.arange(len(self.time))
                flag_firstline = False
            else:
                dat = line.strip().split(';')
                for k in range(len(dat)):
                    SV[obs_types[k]][num_obs] = np.nan if dat[k] == '' else float(dat[k])
                num_obs += 1
        file_SVobs.close()
        return (SV)

    def filt(self, **options):
        ''' Фильтрация по параметрам. Изменяет список НКА для выдачи итератором'''
        regexp = ''
        self.SV_list = self.SV_list_full.copy()

        # собираем регулярное выражение для фильтрации имен файлов(НКА) в списке
        regexp = '[{}]'.format(options.get('gnss', 'GREC'))
        for SV in self.SV_list_full:
            if not (re.search(regexp, SV)):
                self.SV_list.remove(SV)  # удаляем элементы, которые не проходят фильтрацию

        # фильтрация по PRN
        prnlist = options.get('prn', None)
        if prnlist:
            for SV in self.SV_list_full:
                if not (getattr(self, SV)['prn'] in prnlist) and (getattr(self, SV)['name'] in self.SV_list):
                    self.SV_list.remove(SV)  # удаляем элементы, которые не проходят фильтрацию

        # фильтрация по NUM
        numlist = options.get('num', None)
        if numlist:
            for SV in self.SV_list_full:
                if not (getattr(self, SV)['num'] in numlist) and (getattr(self, SV)['name'] in self.SV_list):
                    self.SV_list.remove(SV)  # удаляем элементы, которые не проходят фильтрацию

    def reset(self):
        ''' Сброс параметров фильтра и возвращение исходного списка НКА'''
        self.SV_list = self.SV_list_full.copy()

    def cuttime(self, start_time, end_time):
        """"
        :param start_time: datetime type
        :param end_time: datetime type
        :return: PyGNSS обрезанный по времени
        """
        self.reset()
        ind_not_removing = (self.continious_time >= start_time.timestamp()) & (
        self.continious_time <= end_time.timestamp())
        for SV in self:
            for key in SV.keys():
                if not (type(SV[key]) == np.ndarray):
                    # защита от невекторных свойств
                    continue
                SV[key] = SV[key][ind_not_removing]
        self.continious_time = self.continious_time[ind_not_removing]
        self.time = [datetime.fromtimestamp(temp) for temp in self.continious_time]

    def __parseRINEX__(self, rinex_file_name, time_interval):
        file_rinexobs = open(rinex_file_name)
        os.mkdir(self.path)
        ## Чтение шапки и временных меток
        ObsTypes = [[] for k in range(3)]
        list_time = []
        for line in file_rinexobs:
            if 'RINEX VERSION' in line:  # проверка формата
                if line[5] != '3':
                    print('Incompatible RINEX version')
                    break
            elif 'OBS TYPES' in line:  # виды измерений
                ObsTypes_temp = line[7:60].split()
                if line[0] == 'G':
                    gnss = 1
                elif line[0] == 'R':
                    gnss = 2
                elif line[0] == 'E':
                    gnss = 3
                elif line[0] == ' ':
                    pass
                else:
                    gnss = None

                if gnss == None:
                    print('Unknown GNSS type')
                    continue
                else:
                    ObsTypes[gnss - 1].extend(ObsTypes_temp)
            elif '>' in line:  # метки времени
                list_time.append(datetime.strptime(line[2:21], '%Y %m %d %H %M %S'))
        file_rinexobs.close()

        ## Подготовка вектора времени
        time_interval = np.maximum((list_time[2].timestamp() - list_time[1].timestamp()), time_interval)
        vct_timestamp = np.arange(time_interval * np.ceil(list_time[0].timestamp() / time_interval),
                                  list_time[-1].timestamp() + 0.001,
                                  time_interval)  # поправка 0,001 чтобы последний элемент тоже был включён
        list_time = []
        for timestmp in vct_timestamp:
            list_time.append(datetime.fromtimestamp(timestmp))
        # Запись в файл     
        file_time = open('{}/time.csv'.format(self.path), 'w')
        for timestamp in list_time:
            file_time.write(str(datetime.timestamp(timestamp)) + '\n')
        file_time.close()

        ## Создание файла инициализации        
        pass
        ## Чтение измерений и запись в отдельные файлы
        file_rinexobs = open(rinex_file_name)
        num_obs = 0
        dict_SepObsFiles = {}
        list_ObsFilesForRemoving = []
        flag_Header = False
        for line in file_rinexobs:
            # метка времени
            if (line[0] == '>'):
                # если нет такого элемента в векторе времени
                if not (datetime.strptime(line[2:21], '%Y %m %d %H %M %S') in list_time):
                    flag_Header = False
                    continue
                kol_missed_obs = list_time.index(datetime.strptime(line[2:21],
                                                                   '%Y %m %d %H %M %S')) - num_obs - 1  # для отработки пропущенных моментов времени
                # поиск соответствующего номера наблюдения
                num_obs = list_time.index(datetime.strptime(line[2:21], '%Y %m %d %H %M %S'))
                # Заполнение пустыми строками файлов, в которые не писали на прошлом такте
                for SepObsFile in list_ObsFilesForRemoving:
                    SepObsFile.write('\n')
                # Заполнение пустыми строками пропусков в измерениях
                for SepObsFile in dict_SepObsFiles.values():
                    for k in range(kol_missed_obs):
                        SepObsFile.write('\n')
                # новая инициализация списка имён файлов для механизма записи пустых строк
                list_ObsFilesForRemoving = list(dict_SepObsFiles.values())
                flag_Header = True
                continue
                # измерения GPS
            elif (line[0] == 'G') and (flag_Header):
                line = line.replace('G ', 'G0')  # защита от файлов с пробелом вместо нуля
                gnss = 0
            # измерения GLO
            elif (line[0] == 'R') and (flag_Header):
                line = line.replace('R ', 'R0')  # защита от файлов с пробелом вместо нуля
                gnss = 1
            # измерения GAL
            elif (line[0] == 'E') and (flag_Header):
                line = line.replace('E ', 'E0')  # защита от файлов с пробелом вместо нуля
                gnss = 2
            # остальные строки
            else:
                continue

                # запись наблюдений в отдельный файл
            # если соответствующего файла нет, то создаём
            SV = line[:3]
            if not (dict_SepObsFiles.get(SV)):
                dict_SepObsFiles[SV] = open('{}/{}{}'.format(self.path, SV, self.__obsext__), 'w')
                # заполнение шапки
                dict_SepObsFiles[SV].write(';'.join(ObsTypes[gnss]) + '\n')
                # заполнение пустых моментов времени
                dict_SepObsFiles[SV].write((';' * len(ObsTypes[gnss]) + '\n') * num_obs)
                # добавление в перечень для механизма записи пустых строк
                list_ObsFilesForRemoving.append(dict_SepObsFiles[SV])
            # записываем
            SepObsFile = dict_SepObsFiles.get(SV)
            list_obs = list(line[k:k + 13].strip() for k in range(4, len(line), 16))  # observations list
            SepObsFile.write(';'.join(list_obs) + '\n')
            list_ObsFilesForRemoving.remove(SepObsFile)  # удаляем имена файлов, в которых писали на этом такте

        # closing
        file_rinexobs.close()
        for SepObsFile in dict_SepObsFiles.values():
            SepObsFile.close()

        # prn reading      
        prntable = np.arange(30) * np.nan
        if os.path.exists(os.getcwd() + '/' + rinex_file_name[:-1] + 'g'):
            file_navinfo = open(rinex_file_name[:-1] + 'g')
        elif os.path.exists(os.getcwd() + '/' + rinex_file_name[:-1] + 'G'):
            file_navinfo = open(rinex_file_name[:-1] + 'G')
        else:
            file_navinfo = None
        if file_navinfo != None:
            flag_3str = 0
            flag_header = True
            for line in file_navinfo:
                if flag_header:  # проматываем заголовок
                    if 'END OF HEADER' in line:
                        flag_header = False
                    continue
                if line[0] == 'R':  # начало блока с информацией о НКА ГЛОНАСС
                    line = line.replace('R ', 'R0')  # защита от файлов с пробелом вместо нуля
                    dat = line.split()
                    num_SV = int(dat[0][1:])
                    flag_3str = 3
                flag_3str -= 1
                if flag_3str == 0:  # строка с номером литеры
                    try:
                        prntable[num_SV - 1] = int(float(line[61:80].replace('D','E')))
                    except ValueError:
                        prntable[num_SV - 1] = int(float(line[61:63]))
                    else:
                        prntable[num_SV - 1]=np.nan

            file_navinfo.close()  # закрываем файл с нав. информацией
        # creation file with liter table
        file_prntable = open(self.path + '/litertable.csv', 'w')
        for prn in prntable:
            file_prntable.write(str(prn) + '\n')
        file_prntable.close()

    def __parseEphemeris__(self, file_name):
        ''' Ephemeris reading from sp3 file'''

        # list time reading
        list_time = []
        file_time = open('{}/time.csv'.format(self.path))
        for line in file_time:
            list_time.append(datetime.fromtimestamp(float(line)))

        # ephemeris reading
        dict_SepObsFiles = {}
        list_ObsFilesForRemoving = []
        num_obs = 0
        file_ephemeris = open(file_name)
        for line in file_ephemeris:
            if line[0] == '*':
                # timestamp
                if not (datetime.strptime(line[3:22], '%Y %m %d %H %M %S') in list_time):
                    continue
                # поиск соответствующего номера наблюдения
                kol_missed_obs = list_time.index(datetime.strptime(line[3:22],
                                                                   '%Y %m %d %H %M %S')) - num_obs - 1  # для отработки пропущенных моментов времени
                num_obs = list_time.index(datetime.strptime(line[3:22], '%Y %m %d %H %M %S'))
                # Заполнение пустыми строками файлов, в которые не писали на прошлом такте
                for SepObsFile in list_ObsFilesForRemoving:
                    SepObsFile.write('\n')
                # Заполнение пустыми строками пропусков в измерениях
                for SepObsFile in dict_SepObsFiles.values():
                    for k in range(kol_missed_obs):
                        SepObsFile.write('\n')
                # новая инициализация списка имён файлов для механизма записи пустых строк
                list_ObsFilesForRemoving = list(dict_SepObsFiles.values())
                continue
            # остальные строки
            elif line[0] == 'P':
                # запись наблюдений в отдельный файл
                SV = line[1:4]
                if SV[1] == ' ':
                    SV = SV[0] + '0' + SV[2:]
                # если соответствующего файла нет, то создаём
                if not (dict_SepObsFiles.get(SV)):
                    dict_SepObsFiles[SV] = open('{}/{}{}'.format(self.path, SV, self.__ephext__), 'w')
                    # заполнение шапки
                    dict_SepObsFiles[SV].write(';'.join(['x', 'y', 'z', 'dt']) + '\n')
                    # заполнение пустых моментов времени
                    dict_SepObsFiles[SV].write((';' * 4 + '\n') * num_obs)
                    # добавление в перечень для механизма записи пустых строк
                    list_ObsFilesForRemoving.append(dict_SepObsFiles[SV])
                # записываем
                SepObsFile = dict_SepObsFiles.get(SV)
                list_eph = list(line[k:k + 13].strip() for k in range(5, len(line), 14))  # ephemeris list
                SepObsFile.write(';'.join(list_eph) + '\n')
                list_ObsFilesForRemoving.remove(SepObsFile)  # удаляем имена файлов, в которых писали на этом такте

        # closing
        file_ephemeris.close()
        for SepObsFile in dict_SepObsFiles.values():
            SepObsFile.close()

# ---------------- конец класса
