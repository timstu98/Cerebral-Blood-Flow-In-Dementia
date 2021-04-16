def interpolate(array,i):
    if i % 1 > 0.01:
        value = ( array[int(i+0.5)] + array[int(i-0.5)] ) / 2
    else:
        value = array[i]
    return value

array = [0,1,2,3,4,5]

print(interpolate(array,1))

print(interpolate(array,1.5))

print(interpolate(array,2.5))

print(interpolate(array,4))