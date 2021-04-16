import pandas as pd
import numpy  as np



# so the pt stuff, done in a baseline fraction form!
def test_p(frac,p_50,p_5):
    print('pt ',frac,':')
    p_p = (1-np.tanh( (frac-p_50) / p_5))
    print('p_p: ',p_p)
    p_n = (1+np.tanh( (frac-p_50) / p_5 ))
    print('p_n: ',p_n)

# c stuff for neurons, done in absolute value (easier to read off graphs)
def test_n(val,c_50,c_5):
    print('c ',val,':')
    n_p = (1+np.tanh( (val-c_50*0.909) / (c_5*0.909)))
    print('n_p: ',n_p)
    n_n = (1-np.tanh( (val-c_50*0.909) / (c_5*0.909 )))
    print('n_n: ',n_n)

# test_p(1,0.83,0.075)
test_p(0.95,1,2)
test_p(1.1,1,2)


# test_n(0.909,1.4,-0.8)
# test_n(1.1,1.4,-0.8)


#done in fraction form. you put in what you want the 1+/- tan to be
# p_or_n =#positive or negative
def return_vals_p(p_or_n,frac_1,sum_val_1,frac_2,sum_val_2):
    if p_or_n == 'p':
        p_5 = (frac_2-frac_1)/(  np.arctanh(1-sum_val_2) - np.arctanh(1-sum_val_1)    )
        p_50 = frac_1 - p_5*(np.arctanh(1-sum_val_1))
    elif p_or_n == 'n':
        p_5 = (frac_2-frac_1)/(  np.arctanh(-1+sum_val_2) - np.arctanh(-1+sum_val_1)    )
        p_50 = frac_1 - p_5*(np.arctanh(-1+sum_val_1))

    print('pericyte p/n: ',p_or_n)
    print('at frac ',frac_1,' output ',sum_val_1,' and at frac ',frac_2,' output ',sum_val_2)
    print('p_50: ',p_50)
    print('p_5: ',p_5)

# return_vals_p('p', 0.95,0.07833144559352,1.1,0.00149205)
# return_vals_p('n', 0.95,1.92166855,1.1,1.9985)

#avsolute values as easier with c.
def return_vals_n(p_or_n,val_1,sum_val_1,val_2,sum_val_2):
    if p_or_n == 'n':
        c_5 = (val_2-val_1)/( 0.909*(  np.arctanh(1-sum_val_2) - np.arctanh(1-sum_val_1)    ))
        c_50 = (1/0.909)*(val_1 - c_5*0.909*(np.arctanh(1-sum_val_1)))
    elif p_or_n == 'p':
        c_5 = (val_2-val_1)/( 0.909*( np.arctanh(-1+sum_val_2) - np.arctanh(-1+sum_val_1)    ))
        c_50 = (1/0.909)*(val_1 - c_5*0.909*(np.arctanh(-1+sum_val_1)))

    print('neuron p/n: ',p_or_n)
    print('at c value ',val_1,' output ',sum_val_1,' and at c value ',val_2,' output ',sum_val_2)
    print('c_50: ',c_50)
    print('c_5: ',c_5)

# return_vals_n('p', 0.909,1.09400,1.1,0.043259)
# return_vals_n('n', 0.909 ,0.9050488 ,1.1 , 1.9567409)