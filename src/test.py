import emission


def j_c(ts, f, ef, w):
    return emission.j_c(ts, f, ef, w)[4]
    
print(j_c(3000,1,7,4.5))