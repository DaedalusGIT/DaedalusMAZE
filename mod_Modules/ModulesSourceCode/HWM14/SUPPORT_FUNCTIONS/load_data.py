from datetime import datetime, timedelta



def load_data_func_past(input_file, dataframe, size):
    print('import past started...')
    print(input_file.head())

    for i in range(size):
        # dataframe[i].date = input_file['Time (UTCG)'].iloc[i]
        # dataframe[i].date_Julian = input_file['DATE_JULIAN'].iloc[i]
        # dataframe[i].g_lat = input_file['LAT'].iloc[i]
        # dataframe[i].g_long = input_file['LON'].iloc[i]
        # dataframe[i].alt = input_file['ALT(km)'].iloc[i]
        # dataframe[i].datetime = datetime.strptime(dataframe[i].date,
        #                                           '%d %b %Y %H:%M:%S.%f')  # convert to datetime from text
        try:
            dataframe[i].date = input_file['Time (UTCG)'].iloc[i]
            dataframe[i].g_lat = input_file['Lat (deg)'].iloc[i]
            dataframe[i].g_long = input_file['Lon (deg)'].iloc[i]
            dataframe[i].alt = input_file['Alt (km)'].iloc[i]
        except:
            dataframe[i].date = input_file['Epoch(UTCG)'].iloc[i]
            dataframe[i].g_lat = input_file['Lat_GEOD(deg)'].iloc[i]
            dataframe[i].g_long = input_file['Lon_GEOD(deg)'].iloc[i]
            dataframe[i].alt = input_file['Height_WGS84 (km)'].iloc[i]
        
        try:
            dataframe[i].datetime = datetime.strptime(dataframe[i].date,'%d %b %Y %H:%M:%S.%f')  # convert to datetime from text day, month, year, hour,min,sec,msec
        except:
            #round nanoseconds to 6 digits and convert to dataframe
            time_string=dataframe[i].date
            date_time,microseconds=time_string.split('.')
            rounding = len(microseconds[:3])
            divisor=10**rounding
            new_micros=int(round(int(microseconds)/divisor,0))
            time_string=date_time+'.'+str(new_micros)
            timestring = datetime.strptime(time_string,'%b %d %Y %H:%M:%S.%f')
            dataframe[i].datetime=timestring
            #dataframe[i].datetime = datetime.strptime(dataframe[i].date,'%b %d %Y %H:%M:%S.%f')  # convert to datetime from text month, day, year, hour,min,sec,msec

    print('import completed...')

def load_data_func_future(input_file, dataframe, size):
    print('import past started...')
    print(input_file.head())



    for i in range(size):
        dataframe[i].date = input_file['Time (UTCG)'].iloc[i]
        # dataframe[i].date_Julian = input_file['DATE_JULIAN'].iloc[i]
        dataframe[i].g_lat = input_file['Lat (deg)'].iloc[i]
        dataframe[i].g_long = input_file['Lon (deg)'].iloc[i]
        dataframe[i].alt = input_file['Alt (km)'].iloc[i]
        dataframe[i].datetime = datetime.strptime(dataframe[i].date,
                                                  '%d %b %Y %H:%M:%S.%f')  # convert to datetime from text
        dataframe[i].f107 = input_file['F107'].iloc[i]
        dataframe[i].f107a = input_file['F107A'].iloc[i]
        dataframe[i].f107p = input_file['F107P'].iloc[i]
        dataframe[i].ap_daily = input_file['AP_DAILY'].iloc[i]
        dataframe[i].ap = input_file['AP'].iloc[i]
        dataframe[i].kp_daily = input_file['KP_DAILY'].iloc[i]
        dataframe[i].kp = input_file['KP'].iloc[i]
        dataframe[i].dst = input_file['DST'].iloc[i]
        dataframe[i].ae = input_file['AE'].iloc[i]

    print('import completed...')

def load_error_parameters(input_file,parameter):
    parameter.irine_error = input_file["IRI_ne_error"][0]/100
    parameter.iriTe_error = input_file["IRI_Te_error"][0]/100
    parameter.iriTi_error = input_file["IRI_Ti_error"][0]/100
    parameter.iriOplus_error = input_file["IRI_Oplus_error"][0]/100
    parameter.iriO2plus_error = input_file["IRI_O2plus_error"][0]/100
    parameter.iriNOplus_error = input_file["IRI_NOplus_error"][0]/100
    parameter.msiseO_error = input_file["MSISE_O_error"][0]/100
    parameter.msiseO2_error = input_file["MSISE_O2_error"][0]/100
    parameter.msiseN2_error = input_file["MSISE_N2_error"][0]/100
    parameter.B_error = input_file["B_error"][0] / 100
