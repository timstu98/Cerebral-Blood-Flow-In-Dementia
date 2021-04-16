import pandas as pd
import numpy as np 
# pd.set_option("display.max_rows", None)
# pd.set_option("display.max_columns", None)

db = pd.read_csv('db.csv')
print(db)

def search(range_1,other_ranges=None,sort=None):
    
    filtered = db[ ( (db[range_1[0]] > range_1[1]) & (db[range_1[0]]< range_1[2]) ) ]

    if other_ranges:
        for range in other_ranges:
            filtered = filtered[ ( (filtered[range[0]] > range[1]) & (filtered[range[0]]< range[2]) ) ]

    if not sort:

        print(filtered)
        return filtered
    else:
        print(filtered.sort_values(by=sort))
        return filtered.sort_values(by=sort)


def filter(column_string):
    out = db.sort_values(by=column_string)
    return out


alpha_range = ['alpha',0.32,0.42]
Q_range = ['final Q',-0.05,1]
delay_range = ['delay',-2,5]
time_range = ['time',150,10000]

# other_ranges = [Q_range,delay_range]

# print(db)


filtered = search(time_range,sort=['final Q','final c'])
# filtered = search(alpha_range,other_ranges = other_ranges,sort=['delay','c_50'])
filtered.to_csv('filtered.csv',index=False)

filter_by = ['delay','kp_constant_s']
# out = filter(filter_by)
# out.to_csv('temp.csv',index=False)

