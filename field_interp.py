import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.set_printoptions(threshold=np.inf)
from scipy.interpolate import griddata

depth = 10
depth_list = np.arange(0, 5000, 10)
CTD = pd.read_csv('data/CTD.txt',
                  sep='\t',
                  header=0)

CTD = CTD[['Station',
           'Longitude [degrees_east]',
           'Latitude [degrees_north]',
           'Bot. Depth [m]',
           'P',
           'T',
           'S',
           'Depth from Pressure (EOS80) [m]']]

CTD.columns = [['St',
                'Lon',
                'Lat',
                'BotDepth',
                'P',
                'T',
                'S',
                'Depth']]


def horizon_interpolate(df, st, step):
    depth_new = np.arange(0, 5000, step)
    t_out = np.interp(depth_new,
                      df[df['St'] == st]['Depth'].values,
                      df[df['St'] == st]['T'].values)
    s_out = np.interp(depth_new,
                      df[df['St'] == st]['Depth'].values,
                      df[df['St'] == st]['S'].values)

    lat = [df[df['St'] == st]['Lat'].values[0]] * len(s_out)
    lon = [df[df['St'] == st]['Lon'].values[0]] * len(s_out)
    st = [st] * len(s_out)
    a = np.column_stack((st, lat, lon, depth_new, t_out, s_out))
    return a


df_out = pd.DataFrame(columns=['St',
                               'Lat',
                               'Lon',
                               'Depth',
                               'T',
                               'S'], data=[])

for item in pd.unique(CTD['St']):
    a = horizon_interpolate(CTD, item, 10)
    df_temp = pd.DataFrame(data=a,
                           columns=['St', 'Lat', 'Lon', 'Depth', 'T', 'S'])
    df_out = df_out.append(df_temp,
                           ignore_index=True)

x = pd.unique(df_out['Lon'])
y = pd.unique(df_out['Lat'])

def grid(x_min, x_max, y_min, y_max, x_step, y_step):
    xi = np.arange(x_min, x_max, x_step)
    yi = np.arange(y_min, y_max, y_step)
    xx, yy = np.meshgrid(xi, yi)
    return xi, yi, xx, yy


xi, yi, xx, yy = grid(-42.5, -39.5, 10, 11.5, 0.002, 0.002)

def get_layer(ga, depth, xx, yy, dfCTD):

    z_t = dfCTD[dfCTD['Depth'] == depth]['T'].values
    #z_s = dfCTD['S'].values - 35
    LADCP_data = np.column_stack((x, y, z_t))

    setlev = 'set lev {}'.format(depth)
    ga(setlev, Quiet=True)
    vel = ga.exp('tt')
    lat = vel.grid.lat
    lon = vel.grid.lon
    T = ga.exp('tt').data.flatten()

    lat_ = []
    lon_ = []
    for j in range(len(lat)):
        for k in range(len(lon)):
            lat_.append(lat[j])
            lon_.append(lon[k])

    grads_data = np.column_stack((lon_, lat_, T.ravel()))



    def cutting(array, xmin, xmax, ymin, ymax):
        array[:, 2][(array[:, 1] > xmin)
                    & (array[:, 1] < xmax)
                    & (array[:, 0] >= ymin)
                    & (array[:, 0] <= ymax)] = np.nan
        data = grads_data[~np.isnan(grads_data[:, 2])]
        return data

    grads_data = cutting(grads_data, -41.5, -40.5, 10.4, 11.1)

    #cut_lon = (-41.5, -40.5)
    #cut_lat = (10.4, 11.1)

    full_data = np.concatenate((grads_data, LADCP_data))

    zi = griddata((full_data[:, 0], full_data[:, 1]),
                  full_data[:, 2],
                  (xx, yy),
                  method='linear')
    return zi


from grads.ganum import GaNum

ga = GaNum(Bin='grads')
ga.open('data/woa_vemafz_tt_10m.ctl', Quiet=True)

arr = get_layer(ga, depth_list[0], xx, yy, df_out)
arr = arr[np.newaxis,...]

for i in range(1, len(depth_list)-1):

    layer1 = get_layer(ga, depth_list[i], xx, yy, df_out)
    arr = np.vstack([arr, layer1[np.newaxis,...]])
    print np.shape(arr)
    print i

np.save('t', arr)