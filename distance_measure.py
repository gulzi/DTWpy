

def distancefunc(name="manhattan"):
    if name=="manhattan":
        return lambda x, y:abs(x - y)
    elif name == "euclidean":
        return lambda x,y:pow(x-y, 2)
    elif name == "canberra":
        return lambda x,y:(abs(x - y) / (x + y))

    return lambda x, y:abs(x - y)